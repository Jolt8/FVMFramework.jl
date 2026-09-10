using Revise

using FVMFramework

using Ferrite
using FerriteGmsh
using OrdinaryDiffEq
using SparseArrays
using ComponentArrays   
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero
using StaticArrays
using PreallocationTools
using ForwardDiff
using Polyester
using DataInterpolations

using XLSX

using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using ForwardDiff

using BenchmarkTools
using Plots

using Unitful

mesh_path = joinpath(@__DIR__, "cone_end/dialysis_tubing_cone_output.msh")

grid = togrid(mesh_path)

grid.cellsets
grid.facetsets

n_cells = length(grid.cells)
u_proto = (
    mass_fractions = (methylene_blue = zeros(n_cells), water = zeros(n_cells)),
)

config = create_fvm_config(grid, u_proto)

#since each mass fraction is modified, they have to be vectors
function normalize_mass_fractions(mass_fractions)
    total_mass_fractions = 0.0

    for (species_name, mass_fraction) in pairs(mass_fractions)
        total_mass_fractions += mass_fractions[species_name][1]
    end

    for species_name in keys(mass_fractions)
        mass_fractions[species_name][1] /= total_mass_fractions
    end

    return NamedTuple{keys(mass_fractions)}(first.(values(mass_fractions)))
end

#these are just for classifying regions to make sure they do the right connection functions
struct Fluid <: AbstractPhysics end
struct WellMixed <: AbstractPhysics end

dialysis_tubing_initial_mass_fractions = (
    methylene_blue = [0.0004],
    water = [0.9996]
)

dialysis_tubing_initial_mass_fractions = normalize_mass_fractions(dialysis_tubing_initial_mass_fractions)

add_region!(
    config, "dialysis_tubing_interior";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = dialysis_tubing_initial_mass_fractions,
    ),
    properties = (
        temp = ustrip(21.13u"°C" |> u"K"),
        k = 0.6, # k (W/(m*K))
        cp = 4184, # cp (J/(kg*K))
        rho = 1000, # rho (kg/m^3)
        mu = 1e-3, # mu (Pa*s)
        species_ids = (methylene_blue = 1, water = 2), #we could use mass_fractions for species loops, but this is just more consistent
        #=diffusion_coefficients = (
            methylene_blue = 1e-3,
            water = 1e-3
        ), #diffusion coefficients (m^2/s)=# 
        #not needed due to the assumption that the dialysis tubing interior and surrounding fluid are perfectly mixed
        molecular_weights = (
            methylene_blue = 0.31985, 
            water = 0.01802
        ), #species_molecular_weights [kg/mol]
    ), 
    optimized_syms = [],
    cache_syms = [:heat, :molar_concentrations, :mw_avg, :rho], 
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        #cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

surrounding_fluid_initial_mass_fractions = (
    methylene_blue = [0.0],
    water = [1.0]
)

surrounding_fluid_initial_mass_fractions = normalize_mass_fractions(surrounding_fluid_initial_mass_fractions)

surrounding_fluid_properties = (
    temp = ustrip(21.13u"°C" |> u"K"),
    k = 0.6, # k (W/(m*K))
    cp = 4184, # cp (J/(kg*K))
    rho = 1000, # rho (kg/m^3)
    mu = 1e-3, # mu (Pa*s)
    species_ids = (methylene_blue = 1, water = 2), #we could use mass_fractions for species loops, but this is just more consistent
    molecular_weights = (
        methylene_blue = 0.31985, 
        water = 0.01802
    ), #species_molecular_weights [kg/mol]
)

add_region!(
    config, "surrounding_fluid";
    type = WellMixed(),
    initial_conditions = (
        mass_fractions = surrounding_fluid_initial_mass_fractions,
    ),
    properties = surrounding_fluid_properties,
    optimized_syms = [],
    cache_syms = [:heat, :molar_concentrations, :mw_avg, :rho], 
    region_function =
    function surrounding_fluid!(du, u, cell_id, vol)
        #property updating/retrieval

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        #cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

#=
add_region!(
    config, "sampled_volume";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = surrounding_fluid_initial_mass_fractions,
    ),
    properties = surrounding_fluid_properties,
    optimized_syms = [],
    cache_syms = [:heat, :molar_concentrations, :mw_avg, :rho], 
    region_function =
    function surrounding_fluid!(du, u, cell_id, vol)
        #property updating/retrieval

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        #cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)
    =#

include("specific_physics/arrenhius_mass_fraction_diffusion_meth_blue_and_water.jl") 
#for arrenhius_mass_fraction_diffusion_meth_blue_and_water!

#this does nothing for now, but it will probably be used to determine the length of the p_vector in the future

#we'll come back to this later, but I'm thinking that instead I'll just add optimized_ in front of 
#the property name to denote that it's an optimized parameter
#=
macro optimize(key_value_pair) 
    param_name = key_value_pair.args[1] #this seem to always return the correct name
    param_value = key_value_pair.args[2]

    #=
    if typeof(key_value_pair.args[2]) <: Expr
        param_value = key_value_pair.args[2].args
    else
        param_value = key_value_pair.args[2]
    end
    =#

    #println(param_name)
    println(param_value)
    
    push!(p_tracker, OptimizedParameterTracker(param_name, param_value))
    
    return :($param_name => $param_value[1])
end

properties = (
    @optimize :test = 1,
    @optimize :diffusion_pre_exponential_factor = 1e-5,
    @optimize :diffusion_activation_energy = 1000.0
)
=#

add_patch!(
    config, "dialysis_tubing_surface";
    properties = (
        diffusion_pre_exponential_factor = 0.0006951928076296322,
        diffusion_activation_energy = 28718.5536309222
    ), 
    optimized_syms = [
        :diffusion_pre_exponential_factor,
        :diffusion_activation_energy
    ],
    patch_function =
    function membrane_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        
        #zeroing out regular diffusion that would happen across this face (now not necessary due to different connection groups)
        #du.mass_fractions.methylene_blue[idx_a] *= 0.0
        #du.mass_fractions.water[idx_a] *= 0.0
        #du.mass_fractions.methylene_blue[idx_b] *= 0.0
        #du.mass_fractions.water[idx_b] *= 0.0

        arrenhius_mass_fraction_diffusion_meth_blue_and_water!(
            du, u,
            idx_a, idx_b, face_idx,
            cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
            cell_volumes[idx_a], cell_volumes[idx_b]
        )
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
    )
end

function fluid_well_mixed_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    #nothing happens because the add_patch! is already taking care of this
end

function well_mixed_well_mixed_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    #nothing happens here, but we have a custom method inside parameter_fitting_operator for this
end


#this is the smallest I could make this function
#I like using <: here because it makes it look nice with syntax highlighting
function connection_map_function(type_a, type_b)
    typeof(type_a) <: Fluid && typeof(type_b) <: Fluid && return fluid_fluid_flux!
    typeof(type_a) <: WellMixed && typeof(type_b) <: Fluid && return fluid_well_mixed_flux!
    typeof(type_a) <: Fluid && typeof(type_b) <: WellMixed && return fluid_well_mixed_flux!
    typeof(type_a) <: WellMixed && typeof(type_b) <: WellMixed && return well_mixed_well_mixed_flux!
end

n_faces = length(config.geo.cell_neighbor_areas[1])
n_cells = length(config.geo.cell_volumes)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
)

#START OF OPTIMIZATION SETUP (just to help differentiate between what is and isn't optimization related)
#after doing some work on this, it might be best to make it so that custom methods are used every single time
#I think making something that is generalizable to any problem would be almost impossible to code
#So, we'll just have some basic methods that help with setup but then everything else requires editing the operator itself 

#=stuff we still have to add
    - methods for assigning experimental data only to certain cellsets 
        - I wonder if the memory cost from passing in experiemental data of n_fields * n_timestamps * n_cells is too much
        - For example, if we had 1 field, 4000 timesteps, and 10000 cells (pretty typical), that's 40,000,000 data points
        - However if we assume the field is constant for each field, then we only need n_fields * n_timestamps data points
        - We could also reduce the resolution of the timesteps, but that's not ideal
    - 
=#

experimental_data_path = joinpath(@__DIR__, "dialysis_tubing_diffusion_data.xlsx")

xf = XLSX.readxlsx(experimental_data_path)

#start of logic that doesn't need to be seen by the user

struct TrialData
    state_data::NamedTuple
    state_time::NamedTuple
    compared_data::NamedTuple
    compared_time::NamedTuple
end

function initialize_trials()
    return Dict{String, TrialData}()
end

function input_experimental_data(trials, xf, name::String; 
    data::NamedTuple, 
    trial_processor_function::Function
)  
    trial_processor_function(trials, data, name, xf)
end

function dv(data) #column_to_vector/row_to_vector, called it dv because it's data_to_vector
    data_vec = float.(vec(data))

    filter!(x -> !isempty(x), data_vec)

    return data_vec
end
#end of logic that doesn't need to be seen by the user

function moles_to_mass_fraction(species_moles, species_molecular_weight, rho, mixture_volume)
    return ((species_moles / mixture_volume) * species_molecular_weight) / rho
end

function process_experimental_data!(trials, data::NamedTuple, name, xf)
    #state data (affects the conditions of the simulation over time)
    temperatures = data.temp_data.temp .* u"°C"
    state_data = (
        temp = CubicSpline(ustrip.((temperatures .|> u"K")), data.temp_data.timestamps),
    )

    state_time = (
        temp = data.temp_data.timestamps, #this is only needed for defining tMax
    )

    #compared data (what the simulated data is compared against to get loss)
    mixture_volume = 500u"ml"
    mixture_rho = 1000.0u"kg/m^3"

    calibration_volume_fractions = 0.0
    corresonding_absorbances = 0.0

    if name == "T1" || name == "T2"
        calibration_volume_fractions = dv(xf["calibration"]["A2:A6"]) .* u"ml/ml"
        corresonding_absorbances = dv(xf["calibration"]["B2:B6"])
    elseif name == "T3" #For some reason the colorimeter became miscalibrated when I did T3
        calibration_volume_fractions = dv(xf["calibration"]["A2:A6"]) .* u"ml/ml"
        corresonding_absorbances = dv(xf["calibration"]["E2:E6"])
    end

    calibration_concentration = 0.016u"g/L"
    
    stock_solution_density = 1000.0u"kg/m^3"

    calibration_mass_fractions = (calibration_volume_fractions .* calibration_concentration) / stock_solution_density

    calibration_absorbances_to_mass_fractions_interp = CubicSpline(ustrip.(calibration_mass_fractions .|> u"kg/kg"), corresonding_absorbances)

    methylene_blue_mass_fractions = calibration_absorbances_to_mass_fractions_interp.(data.absorbance_data.methylene_blue_absorbance)
    water_mass_fractions = 1.0 .- methylene_blue_mass_fractions

    compared_data = (
        mass_fractions = (
            methylene_blue = ustrip.((methylene_blue_mass_fractions .|> u"kg/kg")),
            water = ustrip.((water_mass_fractions .|> u"kg/kg")),
        ),
    )

    compared_time = (
        mass_fractions = data.absorbance_data.timestamps,
    )

    trial = TrialData(
        state_data,
        state_time,
        compared_data,
        compared_time
    )

    trials[name] = trial
end

trials = initialize_trials()

input_experimental_data(trials, xf, "T1"; 
    data = (
        temp_data = (
            timestamps = dv(xf["T1 Temps"]["A2:A2997"]),
            temp = dv(xf["T1 Temps"]["B2:B2997"]),
        ),
        absorbance_data = ( #I think it would be safe to assume that values that are measured in the same field have the same timestamps
            timestamps = dv(xf["T1"]["B2:B10"]),
            methylene_blue_absorbance = dv(xf["T1"]["C2:C10"]),
        )
    ),
    trial_processor_function = process_experimental_data!
)


input_experimental_data(trials, xf, "T2"; 
    data = (
        temp_data = (
            timestamps = dv(xf["T2 Temps"]["A2:A1693"]),
            temp = dv(xf["T2 Temps"]["B2:B1693"]),
        ),
        absorbance_data = ( #I think it would be safe to assume that values that are measured in the same field have the same timestamps
            timestamps = dv(xf["T2"]["B2:B11"]),
            methylene_blue_absorbance = dv(xf["T2"]["C2:C11"]),
        )
    ),
    trial_processor_function = process_experimental_data!
)


input_experimental_data(trials, xf, "T3"; 
    data = (
        temp_data = (
            timestamps = dv(xf["T3 Temps"]["A2:A1727"]),
            temp = dv(xf["T3 Temps"]["B2:B1727"]),
        ),
        absorbance_data = ( #I think it would be safe to assume that values that are measured in the same field have the same timestamps
            timestamps = dv(xf["T3 Dropped Data"]["B2:B8"]),
            methylene_blue_absorbance = dv(xf["T3 Dropped Data"]["C2:C8"]),  
        )
    ),
    trial_processor_function = process_experimental_data!
)

#END OF OPTIMIZATION SETUP

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, special_caches)

u_test = (; create_views_inline(u0_vec, system.u_proto_axes)..., create_views_inline(get_tmp(system.u_diff_cache_vec, 0.0), system.u_cache_axes)...)

f_closure_implicit_pre = (du, u, p, t, state_data, state_time) -> 
methylene_blue_diffuion_parameter_fitting_f!(
    du, u, p, t, 

    state_data, state_time,
    system.p_axes,
    
    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals, 

    system.connection_groups, system.patch_groups, system.region_groups,

    system.merged_properties, system.du_diff_cache_vec, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system
f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trials["T1"].state_data, trials["T1"].state_time)

tMax = trials["T1"].state_time.temp[end]

prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, tMax), system.p_vec)
#@time sol = solve(prob, Tsit5(), saveat = tMax / 100)#, callback = approximate_time_to_finish_cb)

#VSCodeServer.@profview sol = solve(prob, Tsit5())

#sol_u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)

#sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

#sol_u_named_0.mass_fractions.methylene_blue == sol_u_named_end.mass_fractions.methylene_blue

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, system.p_vec, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0 = 0.0
#tMax is already defined above based on experimental data
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, system.p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

@time sol = solve(implicit_prob, FBDF(linsolve = KLUFactorization()))

sol_u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)

sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

sol_u_named_0.mass_fractions.methylene_blue == sol_u_named_end.mass_fractions.methylene_blue

sim_file = @__FILE__

u_proto_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

sol_to_vtk(sol, u_proto_named, grid, sim_file)


#VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve = KLUFactorization()))
#VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback=approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

#observed_cell_id = grid.cellsets["molar_concentrations_observation_point"][1] #this is what we should do
#observed_cell_id = grid.cellsets["sampled_volume"][1] #1676
observed_cell_id = 734

trials["T3"].compared_data.mass_fractions.methylene_blue

observed = Dict(
    "T1" => zeros(length(trials["T1"].compared_time.mass_fractions)), 
    "T2" => zeros(length(trials["T2"].compared_time.mass_fractions)), 
    "T3" => zeros(length(trials["T3"].compared_time.mass_fractions))
)
predicted = Dict(
    "T1" => zeros(length(trials["T1"].compared_time.mass_fractions)), 
    "T2" => zeros(length(trials["T2"].compared_time.mass_fractions)), 
    "T3" => zeros(length(trials["T3"].compared_time.mass_fractions))
)

trial_probs = Dict()

for trial_name in keys(trials)
    trial = trials[trial_name]

    tspan = (trial.state_time.temp[1], trial.state_time.temp[end])

    p_adjusted = [0.0, 0.0]

    f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trial.state_data, trial.state_time)

    ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
    #since this is just a diffusion problem, I doubt implicit solving is necessary

    trial_probs[trial_name] = ode_func
end


#START OF OPTIMIZATION SOLVING
function check_predicted_against_observed(θ)
    total_loss = 0.0
    solutions = Dict()
    for trial_name in keys(trials)
        trial = trials[trial_name]

        tspan = (trial.state_time.temp[1], trial.state_time.temp[end])
        p_adjusted = [exp(θ[1]), θ[2] * 1000]

        #f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trial.state_data, trial.state_time)

        #ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
        #since this is just a diffusion problem, I doubt implicit solving is necessary

        t_interval = trial.compared_time.mass_fractions[end] / 200

        prob_trial = ODEProblem(trial_probs[trial_name], u0_vec, tspan, p_adjusted)

        sol = solve(
            prob_trial, Tsit5(), 
            saveat = t_interval,
            #reltol = 1e-10, abstol = 1e-10
        )

        solutions[trial_name] = sol

        #I FOUND THE PROBLEM: using an implicit solver like FBDF here does not seem to agree well with ForwardDiff, it makes using gradient on this loss function take 10 seconds
        #instead of 0.6-1 second

        u_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

        #println(trial_name)

        trial_error = 0.0
        for time_idx in eachindex(trial.compared_time.mass_fractions)
            pred = u_named[time_idx].mass_fractions.methylene_blue[observed_cell_id]
            obs = trial.compared_data.mass_fractions.methylene_blue[time_idx]
            #println("pred: ", pred)
            #println("obs: ", obs)
            #println("")
            predicted[trial_name][time_idx] = pred
            observed[trial_name][time_idx] = obs
            trial_error += sum(abs2, pred .- obs)
        end

        total_loss += trial_error
    end
    return total_loss, solutions
end

function loss(θ)
    total_loss = 0.0
    for trial_name in keys(trials)
        trial = trials[trial_name]

        tspan = (trial.state_time.temp[1], trial.state_time.temp[end])
        p_adjusted = [exp(θ[1]), θ[2] * 1000]

        #f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trial.state_data, trial.state_time)

        #ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
        #since this is just a diffusion problem, I doubt implicit solving is necessary

        prob_trial = ODEProblem(trial_probs[trial_name], u0_vec, tspan, p_adjusted)

        sol = solve(
            prob_trial, Tsit5(), 
            saveat = trial.compared_time.mass_fractions,
            #reltol = 1e-14, abstol = 1e-14
        )

        #I FOUND THE PROBLEM: using an implicit solver like FBDF here does not seem to agree well with ForwardDiff, it makes using gradient on this loss function take 10 seconds
        #instead of 0.6-1 second

        u_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

        #println(trial_name)

        trial_error = 0.0
        for time_idx in eachindex(sol.u)
            pred = u_named[time_idx].mass_fractions.methylene_blue[observed_cell_id]
            obs = trial.compared_data.mass_fractions.methylene_blue[time_idx]
            #println("pred: ", pred)
            #println("obs: ", obs)
            #println("")
            #predicted[trial_name][time_idx] = pred
            #observed[trial_name][time_idx] = obs
            trial_error += sum(abs2, pred .- obs)
        end

        total_loss += trial_error
    end
    return total_loss
end

#= This is for recording once you've done your optimization
total_loss, solutions = check_predicted_against_observed([log(0.0006951928076296322), 28.7185536309222])

T1_u_named = [create_views_inline(solutions["T1"].u[i], system.u_proto_axes) for i in eachindex(solutions["T1"].u)]
T2_u_named = [create_views_inline(solutions["T2"].u[i], system.u_proto_axes) for i in eachindex(solutions["T2"].u)]
T3_u_named = [create_views_inline(solutions["T3"].u[i], system.u_proto_axes) for i in eachindex(solutions["T3"].u)]

sim_file = @__FILE__

t_interval = solutions["T1"].t[end] 
test = trials["T1"].compared_time.mass_fractions[end] 

sol_to_vtk(solutions["T1"], T1_u_named, grid, sim_file)
sol_to_vtk(solutions["T2"], T2_u_named, grid, sim_file)
sol_to_vtk(solutions["T3"], T3_u_named, grid, sim_file)
=#

grad = ForwardDiff.gradient(loss, [log(1.355e-5), 25.0]) #while using an implicit solver, this takes forever
#VSCodeServer.@profview ForwardDiff.gradient(loss, [log(1.355e-5), 19.0])

@time loss([log(1.355e-5), 19.0])

#VSCodeServer.@profview ForwardDiff.gradient(loss, [log(0.00010695283811760264), 25.0])
A_range = exp.(range(log(1e-4), log(5e-4), length = 20))
Ea_range = exp.(range(log(26.0), log(30.0), length = 20))

time_to_finish = (@timed loss([log(1e-7), 10.0])).time
approximate_time_to_finish = time_to_finish * length(A_range) * length(Ea_range)

losses = zeros(length(A_range), length(Ea_range))

#=
for (A_idx, A) in enumerate(A_range)
    for (Ea_idx, Ea) in enumerate(Ea_range)
        curr_loss = loss([log(A), Ea])
        losses[A_idx, Ea_idx] = curr_loss
        #=
        if Ea_idx > 1
            if losses[A_idx, Ea_idx] > losses[A_idx, Ea_idx-1]
                losses[A_idx, Ea_idx:end] .= losses[A_idx, Ea_idx-1]
                break
            end
        end
        if A_idx > 1
            if losses[A_idx, Ea_idx] > losses[A_idx-1, Ea_idx]
                losses[A_idx:end, Ea_idx] .= losses[A_idx-1, Ea_idx]
                break
            end
        end
        =#
        println("loss: ", curr_loss, "    A: ", A, "    Ea: ", Ea)
    end
end
=#

losses

idx = argmin(losses)

losses[idx]
A_range[idx[1]] #1.0e-7
Ea_range[idx[2]] #4638.456927175568

surface(A_range, Ea_range, losses)

total_loss, solutions = check_predicted_against_observed([log(A_range[idx[1]]), Ea_range[idx[2]]])

predicted
observed

plot(trials["T1"].compared_time.mass_fractions, observed["T1"])
plot!(trials["T1"].compared_time.mass_fractions, predicted["T1"])

plot(trials["T2"].compared_time.mass_fractions, observed["T2"])
plot!(trials["T2"].compared_time.mass_fractions, predicted["T2"])

plot(trials["T3"].compared_time.mass_fractions, observed["T3"])
plot!(trials["T3"].compared_time.mass_fractions, predicted["T3"])

loss_history = [] #loss accumulator
parameter_history = [] #parameter accumulator

optimization_callback_tracker = function (state, l)
    println("loss: ", l)
    println("optimized_parameters, ", state.u)
    println("grad, ", state.grad)
    #println("parameters: ", state.p)
    push!(loss_history, l)
    push!(parameter_history, state.p)
    return false
end

N = 2
#ForwardDiff.pickchunksize(length(u0_vec))

adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)

#diffusion pre-exponential factor, diffusion activation energy
guess_params = Float64[log(0.0006951927961775613), 28.718553641885585]
loss(guess_params)

lower_bounds = [log(1e-7), 1.0000]
upper_bounds = [log(1e-1), 100.0000]

optprob = Optimization.OptimizationProblem(optf, guess_params, lb = lower_bounds, ub = upper_bounds)

abstol = 1e-100
#bruh, this whole time I could've just decreased the tolerance to get it to explore more 
#what the fuck, how does that even make any sense

res = Optimization.solve(
    optprob,
    callback = optimization_callback_tracker,
    OptimizationOptimJL.LBFGS(),
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works ok
    f_abstol = abstol,
    g_abstol = abstol,
)

optimized_diffusion_pre_exponential_factor = exp(res.u[1])
optimized_diffusion_activation_energy = res.u[2]

#END OF OPTIMIZATION SOLVING

total_loss, solutions = check_predicted_against_observed([log(optimized_diffusion_pre_exponential_factor), optimized_diffusion_activation_energy])

predicted
observed

plot(trials["T1"].compared_time.mass_fractions, observed["T1"])
plot!(trials["T1"].compared_time.mass_fractions, predicted["T1"])

plot(trials["T2"].compared_time.mass_fractions, observed["T2"])
plot!(trials["T2"].compared_time.mass_fractions, predicted["T2"])

plot(trials["T3"].compared_time.mass_fractions, observed["T3"])
plot!(trials["T3"].compared_time.mass_fractions, predicted["T3"])


using DataFrames

df1 = DataFrame(AA = trials["T1"].compared_time.mass_fractions, AB = observed["T1"], AC = predicted["T1"])
df2 = DataFrame(AA = trials["T2"].compared_time.mass_fractions, AB = observed["T2"], AC = predicted["T2"])
df3 = DataFrame(AA = trials["T3"].compared_time.mass_fractions, AB = observed["T3"], AC = predicted["T3"])

results_output = joinpath(@__DIR__, "dialysis_tubing_parameter_fitting_results.xlsx")

XLSX.writetable(results_output, 
    "T1" => df1, 
    "T2" => df2,
    "T3" => df3
)

#FINAL OPTIMIZED PARAMETERS (DO NOT TOUCH)
    # - A: 0.0006951928076296322
    # - Ea: 28.7185536309222
#

#TODO: check if add_patch! actually creates two connections per facet

loss([log(9.998618464786643e-5), 27.999999950598372])
loss([log(0.00010695283811760264), 24.99999997677348])
loss([log(7.41846435984539e-5), 23.999997648035107])
loss([log(5.129888295561744e-5), 22.9999977703627])
loss([log(3.548090572948086e-5), 21.999997914997068])
loss([log(1.3559637272570382e-5), 19.38750000355438])
loss([log(0.00010707372974861651), 24.97650595937726])
loss([log(0.00040703854247398335), 24.548420545838148])

#START OF INITIAL SEARCHING

total_loss, solutions = check_predicted_against_observed([log(0.00010707372974861651), 24.97650595937726])

predicted
observed

plot(trials["T1"].compared_time.mass_fractions, observed["T1"])
plot!(trials["T1"].compared_time.mass_fractions, predicted["T1"])

plot(trials["T2"].compared_time.mass_fractions, observed["T2"])
plot!(trials["T2"].compared_time.mass_fractions, predicted["T2"])

plot(trials["T3"].compared_time.mass_fractions, observed["T3"])
plot!(trials["T3"].compared_time.mass_fractions, predicted["T3"])


Ea_guess_lb = 4000
Ea_guess_ub = 25000

Ea_range = collect(range(Ea_guess_lb, Ea_guess_ub, 2))

A_lb = 1e-9
A_ub = 5e-6

A_guess = 1e-6

function loss_for_fixed_Ea(ln_A_val, fixed_Ea)
    return loss([ln_A_val[1], fixed_Ea])
end

zoop(x) = loss_for_fixed_Ea(x, 30000)

ForwardDiff.gradient(zoop, [1e-6])
#=
function profile_Ea_scan(Ea_range, tol)
    results = []

    for Ea_val in Ea_range
        prob_1D = Optimization.OptimizationProblem(
            Optimization.OptimizationFunction((x, p) -> loss_for_fixed_Ea(x, Ea_val), 
            Optimization.AutoForwardDiff()),
            [log(A_guess)], # Initial guess for A
            lb = [log(A_lb)], 
            ub = [log(A_ub)],
        )

        sol = Optimization.solve(
            prob_1D, 
            callback = optimization_callback_tracker,
            OptimizationOptimJL.LBFGS(),
            f_abstol = tol,
            g_abstol = tol,
        ) # Brent is optimal for 1D

        push!(results, (Ea=Ea_val, A=exp(sol.u[1]), Loss=sol.objective))
        println("Ea: $Ea_val, Opt_A: $(exp(sol.u[1])), Loss: $(sol.objective)")
    end
    return results
end
=#

tol = 1e-18

results = profile_Ea_scan(Ea_range, tol)

A_range = [results[i].Ea for i in eachindex(results)]
loss_range = [results[i].Loss for i in eachindex(results)]

min_idx = argmin(loss_range)

results[min_idx]
results[min_idx].Ea
results[min_idx].A
results[min_idx].Loss

#using Plots
plot(Ea_range, A_range)
plot(Ea_range, loss_range)
#END OF INITIAL SEARCHING

record_sol = true

sim_file = @__FILE__

u_proto_named = [(; create_views_inline(sol.u[i], system.u_proto_axes)...) for i in eachindex(sol.u)]

#u_named = rebuild_u_named(sol.u, u_proto_named)

cell = grid.cellsets["dialysis_tubing_interior"][1]

mass_fractions_dialysis_tubing_interior = []

for step in eachindex(sol.u)
    push!(mass_fractions_dialysis_tubing_interior, u_proto_named[step].mass_fractions.methylene_blue[cell])
end

mass_fractions_dialysis_tubing_interior

sol.u[1] == sol.u[end]

mass_fractions_beginning = u_proto_named[1].mass_fractions.methylene_blue[cell]

mass_fractions_end = u_proto_named[end].mass_fractions.methylene_blue[cell]

mass_fractions_beginning == mass_fractions_end
mass_fractions_beginning - mass_fractions_end

if record_sol == true
    sol_to_vtk(sol, u_proto_named, grid, sim_file)
end