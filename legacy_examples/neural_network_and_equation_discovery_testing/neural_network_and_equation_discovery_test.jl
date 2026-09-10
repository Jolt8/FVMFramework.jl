using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using NonlinearSolve
using Sparspak
using Dates
using DataInterpolations
using Random
using Lux

using XLSX
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using ForwardDiff
using DataFrames
using CSV
using JLD2
using OptimizationOptimisers

using FVMFramework

mesh_path = joinpath(@__DIR__, "meshes", "dialysis_tubing_cone_output.msh")

grid = togrid(mesh_path)

grid.cellsets
grid.facetsets

n_cells = length(grid.cells)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

#properties
Revise.includet(joinpath(@__DIR__, "properties", "steel_properties.jl"))
heated_properties = get_steel_properties()
surrounding_fluid_properties = merge_properties(heated_properties, ComponentVector(k = 1.0u"W/(m*K)"))
surrounding_fluid_properties.k

function update_solid_properties!(du, u, p, t, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.solid_rho[cell_id]
end

function solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol) #this will hopefully be handled by the neural network
end

u_proto = ComponentVector(
    temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto);

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])

rng = Random.default_rng()
nn = Lux.Chain(
    Lux.Dense(6, 10, tanh),
    Lux.Dense(10, 1)
)

# Initialize neural network parameters and state
p_nn_named_tuple, st_nn = Lux.setup(rng, nn)

p_nn = ComponentVector(p_nn_named_tuple)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        k = u"W/(m*K)",
        cp = u"J/(kg*K)",
        rho = u"kg/m^3",
        mass = u"kg",
        mass_face = u"kg",
    ),
    special_caches = ComponentArray(
        mass_face = zeros(n_cells, n_faces)u"kg",
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector()
)

add_region!(
    config, "surrounding_fluid";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 21.0u"°C",
    ),
    properties = surrounding_fluid_properties,
    region_function =
    function inlet!(du, u, p, t, cell_id, vol)
        #update_fluid_properties!(du, u, cell_id, vol, system)

        #du.heat[cell_id] += 10.0

        solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "dialysis_tubing_interior";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 21.0u"°C",
    ),
    properties = heated_properties,
    region_function =
    function heated_region!(du, u, p, t, cell_id, vol)
        #update_fluid_properties!(du, u, cell_id, vol, system)

        du.heat[cell_id] += 0.1

        solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

#Connection functions
function normal_solid_solid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes, st_nn
)
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function nn_solid_solid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes, st_nn
)
    T_ref = 300.0
    k_ref = 50.0
    d_ref = 0.2
    A_ref = 1.0

    du_heat_nn, _ = nn([
        u.temp[idx_a] / T_ref,
        u.temp[idx_b] / T_ref,
        u.k[idx_a] / k_ref,
        u.k[idx_b] / k_ref,
        cell_neighbor_distances[idx_a][face_idx] / d_ref,
        cell_neighbor_areas[idx_a][face_idx] / A_ref,
    ], p, st_nn)

    du.heat[idx_a] += du_heat_nn[1]
end

solid_solid_flux_closure!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes,
) = nn_solid_solid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes,
    st_nn
)

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Solid && typeof(phys_b) <: Solid && return solid_solid_flux_closure!
end

function solve_system!(du, u, p, t, geo, system)
    for cell_id in eachindex(grid.cells)
        update_solid_properties!(du, u, p, t, cell_id, geo.cell_volumes[cell_id], system)
    end

    solve_connection_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

function build_prob(f_closure, du0_vec, u0_vec, tMax, p_guess)
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
    )

    ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

    t0 = 0.0
    tspan = (t0, tMax)

    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

    return implicit_prob
end

experimental_tMax = 1000.0
saveat = experimental_tMax / 100

implicit_prob = build_prob(f_closure, du0_vec, u0_vec, experimental_tMax, p_nn)

@time sol = solve(implicit_prob, FBDF(linsolve = SparspakFactorization()), saveat = saveat, callback = approximate_time_to_finish_cb)
#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)

u_named = [ComponentVector(sol.u[i], system.state_axes) for i in 1:length(sol.u)];

# Create a DataFrame where each cell's properties are stored in separate columns
df_fake_data = DataFrame(t = sol.t)
for propertyname in propertynames(u_named[1])
    val = getproperty(u_named[1], propertyname)
    if val isa AbstractVector && length(val) == n_cells
        for cell_id in 1:n_cells
            df_fake_data[!, Symbol("$(propertyname)_$cell_id")] = [getproperty(u_named[i], propertyname)[cell_id] for i in eachindex(u_named)]
        end
    end
end

csv_path = joinpath(@__DIR__, "fake_experimental_data", "fake_experimental_data.csv")
#CSV.write(csv_path, df_fake_data)
dataframe_fake_experimental_data = CSV.read(csv_path, DataFrame)

# Create a dictionary of arrays of interpolations, one array per property name
fake_experimental_data_interps_dict = Dict{Symbol, Any}()
for propertyname in propertynames(u_named[1])
    val = getproperty(u_named[1], propertyname)
    fake_experimental_data_interps_dict[propertyname] = [
        LinearInterpolation(dataframe_fake_experimental_data[!, Symbol("$(propertyname)_$cell_id")], dataframe_fake_experimental_data.t)
        for cell_id in 1:n_cells
    ]
end

property_interps = NamedTuple(pairs(fake_experimental_data_interps_dict))

property_interps.temp[1](3.0)

observed_cells = [2]

function loss(θ)
    prob, system_copy, geo_copy = get!(task_local_storage(), :implicit_prob_setup) do
        system_copy = deepcopy(system)
        geo_copy = deepcopy(geo)
        
        f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo_copy, system_copy)
        
        prob = build_prob(f_closure, du0_vec, u0_vec, experimental_tMax, θ)
        return (prob, system_copy, geo_copy)
    end

    loss_prob = remake(prob, p = θ)

    sol = solve(
        loss_prob,
        FBDF(linsolve = SparspakFactorization()), 
        sensealg = ForwardSensitivity(),
        saveat = saveat,
        #callback = approximate_time_to_finish_cb
    )

    n_saves = Int(experimental_tMax / saveat)

    if length(sol.t) < n_saves + 1
        return nothing, nothing, 1e10
    end

    u_named = [ComponentVector(sol.u[i], system.state_axes) for i in eachindex(sol.u)]

    mean_squared_error = 0.0

    for (i, t) in enumerate(sol.t)
        for cell_id in observed_cells
            #@show u_named[i].temp[cell_id]
            #@show property_interps.temp[cell_id](t)
            mean_squared_error += abs2(u_named[i].temp[cell_id] - property_interps.temp[cell_id](t))
            #@show mean_squared_error
        end
    end

    return mean_squared_error / length(sol.t)
end

loss(p_nn)

root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files"

#du_complete, u_complete = regenerate_fvm_state(sol, system, solve_system!, geo, p_nn, u_additional_information = ComponentVector());

#sol_to_vtk(sol, du_complete, u_complete, grid, geo, @__FILE__, root_dir; include_zeros_fields = false) 

callback = function (state, loss; doplot = false)
    println("Current Loss: ", loss)
    return false  # Return true if you want to stop optimization early
end

adtype = Optimization.AutoForwardDiff()
opt_func = OptimizationFunction((p, l) -> loss(p), adtype)

opt_prob = OptimizationProblem(opt_func, p_nn)

println("--- Starting Adam Optimization ---")
opt_sol = Optimization.solve(opt_prob, OptimizationOptimisers.Adam(0.01), callback = callback, maxiters = 10)

println("--- Starting L-BFGS Optimization ---")
opt_prob_bfgs = OptimizationProblem(opt_func, opt_sol.u)
opt_sol_bfgs = Optimization.solve(opt_prob_bfgs, LBFGS(), callback = callback, maxiters = 2) 
#this takes forever, we really need a reverse differentiation library to be avaliable again (looking at YOU Enzyme.jl)

p_trained = opt_sol_bfgs.u

model_path = joinpath(@__DIR__, "learned_model.jld2")
jldsave(model_path; nn, p_trained, st_nn)

#recover later
model_data = load(model_path)
loaded_nn = model_data["nn"]
loaded_p = model_data["p_trained"]
loaded_st = model_data["st_nn"]

x = ComponentVector(
    temp_a = collect(290:1:310),
    temp_b = collect(290:1:310),
    k_a = 50.0,
    k_b = 50.0,
    distance = 0.2,
    area = 1.0
)

T_ref = 300.0
k_ref = 50.0
d_ref = 0.2
A_ref = 1.0

test_data_1 = [loaded_nn([temp_a/T_ref, x.temp_b[1]/T_ref, x.k_a/k_ref, x.k_b/k_ref, x.distance/d_ref, x.area/A_ref], loaded_p, loaded_st)[1][1] for temp_a in x.temp_a]
test_data_2 = [loaded_nn([x.temp_a[1]/T_ref, temp_b/T_ref, x.k_a/k_ref, x.k_b/k_ref, x.distance/d_ref, x.area/A_ref], loaded_p, loaded_st)[1][1] for temp_b in x.temp_b]

using Plots

plot(x.temp_a, test_data_1, label = "varied_temp_a")
plot(x.temp_b, test_data_2, label = "varied_temp_b")