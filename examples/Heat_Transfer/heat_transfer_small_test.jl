using Revise
using Logging

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

using Unitful

grid_dimensions = (100, 10, 10)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Tetrahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))

n_cells = length(grid.cells)
u_proto = ComponentVector(
    temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto)

n_faces = length(config.geo.cell_neighbor_areas[1])

struct Solid <: AbstractPhysics end

add_setup_syms!(
    config;
    cache_syms_and_units = (
        heat = u"J",
        cache_test = u"s",
        optimial_parameters_test = u"1/1",
    ),
    special_caches = ComponentVector(),
    optimized_parameters = ComponentVector(
        optimial_parameters_test = 1.0
    )
)

function cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    # J/s /= m^3 * kg*m^3 * J/(kg*K)
    # = K/s
    du.temp[cell_id] += du.heat[cell_id] / (vol * u.rho[cell_id] * u.cp[cell_id])
end


add_region!(
    config, "copper";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 270.0u"°C",
    ),
    properties = ComponentVector(
        k = 237.0u"W/(m*K)", 
        rho = 2700.0u"kg/m^3",
        cp = 921.0u"J/(kg*K)",
    ),
    property_update_function = function copper_property_update!(du, u, p, t, cell_id, vol, system)
        
    end,
    region_function =
    function heat_transfer!(du, u, cell_id, vol)
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "steel";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 350.0u"°C",
    ),
    properties = ComponentVector(
        k = 123.0u"W/(m*K)",
        rho = 7800.0u"kg/m^3",
        cp = 450.0u"J/(kg*K)",
    ),
    property_update_function = function steel_property_update!(du, u, p, t, cell_id, vol, system)
        
    end,
    region_function =
    function heat_transfer!(du, u, p, t, cell_id, vol)
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

function heat_diffusion!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / dist

    du.heat[idx_a] -= -k_effective * grad_T * area
end

function solid_solid_flux!(
    du, u, p, t, 
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)

    #hmm, perhaps these physics functions need to be more strictly typed
    #Checking profview, I'm getting some runtime dispatch and GC here, I don't know why 
    heat_diffusion!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function connection_map_function(type_a, type_b)
    typeof(type_a) <: Solid && typeof(type_b) <: Solid && return solid_solid_flux!
end

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

function solve_system!(du, u, p, t, geo, system)
    append_fixed_properties_and_p_to_u!(u, p, system)

    solve_connection_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

test_prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 1000.0), system.p_vec)
sol = solve(test_prob, Tsit5(), tspan = (0.0, 10.0))

t0 = 0.0
tMax = 1000.0
tspan = (t0, tMax)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, system.p_vec, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, system.p_vec)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb, saveat = save_interval)
@time sol = solve(implicit_prob, FBDF(linsolve = KLUFactorization(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF(precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
#728.605 ms (341178 allocations: 1.08 GiB) (non-multithreaded)
#787.708 ms (202070 allocations: 1.07 GiB) (only connections multithreading)
#799.862 ms (219386 allocations: 1.07 GiB) (everything multithreaded)
#864.517 ms (211023 allocations: 1.07 GiB) (only regions multithreading)


VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

record_sol = true

sim_file = @__FILE__

du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, geo, system.p_vec)

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end