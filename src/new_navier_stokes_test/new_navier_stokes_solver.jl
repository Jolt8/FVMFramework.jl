#NOTE:
#This file was written entirely by AI
#the reason I did this is because I wanted to test if my framework was fast enough to actually handle navier stokes and from what I'm seeing, it definitely is!!!
#a problem with 100 x 100 cells (10000 total) only took 114 seconds to solve for a lid-driven cavity simulation which is in the ballpark of OpenFOAM
#thus, while I'm not going to work on navier stokes further personally because I really don't want to sell my soul to the devil of getting navier stokes to work and be physically consistent,
#and also because navier stokes doesn't really seem that useful in chemical engineering for a majority of problems (and is also often too slow even with extremely efficient solvers for parameter optimization), 
#I at least know it's possible if I want to do it in the future.

using Revise
using Logging
using Unitful
using OrdinaryDiffEq
using NonlinearSolve
using Sparspak
using Ferrite
using SparseConnectivityTracer
using ComponentArrays
import ADTypes
using ILUZero
using LinearAlgebra
using StaticArrays
using FVMFramework


grid_dimensions = (10, 10, 10)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

# Grid sets
addcellset!(grid, "fluid", xyz -> true)

n_cells = length(grid.cells)

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3",
    momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_v = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_w = zeros(n_cells)u"kg/(m^2*s)",
    volumetric_energy = zeros(n_cells)u"J/m^3"
)

config = create_fvm_config(grid, u_proto);

add_setup_syms!(
    config;
    cache_syms_and_units = (
        vel_u = u"m/s",
        vel_v = u"m/s",
        vel_w = u"m/s",
        pressure = u"Pa",
        temperature = u"K",
        speed_of_sound = u"m/s",
        #specific_internal_energy = u"J/kg", #not cached right now, just a property
    ),
    special_caches = ComponentVector(),
    second_order_syms = [],
    optimized_parameters = ComponentVector(),
)


struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "flux_functions/navier_stokes_fluid_fluid_flux.jl"))
Revise.includet(joinpath(@__DIR__, "patch_functions/navier_stokes_patch_flux.jl"))

# Mapping connection functions
function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
end

# Fluid parameters (Re = 100)
Re = 100.0
U_lid = 1.0u"m/s"

# Add region and patches
add_region!(
    config, "fluid";
    type = Fluid(),
    initial_conditions = ComponentVector(
        density = 1.0u"kg/m^3",
        momentum_density_u = 0.0u"kg/(m^2*s)",
        momentum_density_v = 0.0u"kg/(m^2*s)",
        momentum_density_w = 0.0u"kg/(m^2*s)",
        volumetric_energy = 0.0u"J/m^3",
    ),
    properties = ComponentVector(
        specific_internal_energy = 0.0u"J/kg",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
    ),
    property_update_function = 
    function update_fluid_properties!(properties, u)
        
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, cell_id, vol)
        du.u_vel[cell_id] += du.u_flow[cell_id] / vol
        du.v_vel[cell_id] += du.v_flow[cell_id] / vol
        du.w_vel[cell_id] += du.w_flow[cell_id] / vol
        du.pressure[cell_id] += du.p_flow[cell_id] / vol
    end,
)

add_patch!(
    config, "left_wall";
    properties = ComponentVector(
        u_bc = 0.0,
        v_bc = 0.0,
        w_bc = 0.0,
    ),
    patch_function = 
    function left_wall_patch_flux!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_face_areas, cell_face_normals, cell_face_distances,
        cell_neighbor_normals, cell_neighbor_distances, 
        cell_volumes
    )
        wall_patch_flux_generic!(
            du, u, p, t,
            idx_a, idx_b, face_idx,
            cell_volumes[idx_a],
            cell_face_areas[idx_a][face_idx], cell_face_normals[idx_a][face_idx], cell_face_distances[idx_a][face_idx],
            cell_neighbor_normals[idx_a], cell_neighbor_distances[idx_a]
        )
    end
)

add_patch!(
    config, "right_wall";
    properties = ComponentVector(
        u_bc = 0.0,
        v_bc = 0.0,
        w_bc = 0.0,
    ),
    patch_function = 
    function right_wall_patch_flux!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_face_areas, cell_face_normals, cell_face_distances,
        cell_neighbor_normals, cell_neighbor_distances, 
        cell_volumes
    )
        wall_patch_flux_generic!(
            du, u, p, t,
            idx_a, idx_b, face_idx,
            cell_volumes[idx_a],
            cell_face_areas[idx_a][face_idx], cell_face_normals[idx_a][face_idx], cell_face_distances[idx_a][face_idx],
            cell_neighbor_normals[idx_a], cell_neighbor_distances[idx_a]
        )
    end
)
add_patch!(
    config, "bottom_wall";
    properties = ComponentVector(
        u_bc = 0.0,
        v_bc = 0.0,
        w_bc = 0.0,
    ),
    patch_function = 
    function bottom_wall_patch_flux!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_face_areas, cell_face_normals, cell_face_distances,
        cell_neighbor_normals, cell_neighbor_distances, 
        cell_volumes
    )
        wall_patch_flux_generic!(
            du, u, p, t,
            idx_a, idx_b, face_idx,
            cell_volumes[idx_a],
            cell_face_areas[idx_a][face_idx], cell_face_normals[idx_a][face_idx], cell_face_distances[idx_a][face_idx],
            cell_neighbor_normals[idx_a], cell_neighbor_distances[idx_a]
        )
    end
)
add_patch!(
    config, "top_lid";
    properties = ComponentVector(
        u_bc = 1.0,
        v_bc = 0.0,
        w_bc = 0.0,
    ),
    patch_function = 
    function top_lid_patch_flux!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_face_areas, cell_face_normals, cell_face_distances,
        cell_neighbor_normals, cell_neighbor_distances, 
        cell_volumes
    )
        wall_patch_flux_generic!(
            du, u, p, t,
            idx_a, idx_b, face_idx,
            cell_volumes[idx_a],
            cell_face_areas[idx_a][face_idx], cell_face_normals[idx_a][face_idx], cell_face_distances[idx_a][face_idx],
            cell_neighbor_normals[idx_a], cell_neighbor_distances[idx_a]
        )
    end
)

# Finish FVM configuration
du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

# System solver function
function solve_system!(du, u, p, t, geo, system)
    solve_connection_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)
p_guess = 0.0
detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

# 1. Direct steady-state solve using NonlinearSolve.jl
f_closure_steady = (du, u, p) -> f_closure_implicit(du, u, p, 0.0)

# Detect Jacobian sparsity for NonlinearProblem
nl_jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_steady(du, u, p_guess), du0_vec, u0_vec, detector
) 
#this scales absolutely abysmally as the number of cells goes up
#for example, a 10x increase in the amount of cells made this take around 100x longer! 

t0 = 0.0
tMax = 100000.0
tspan = (t0, tMax)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(nl_jac_sparsity))
implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

println("Solving the Navier-Stokes ODE system...")
@time sol = solve(
    implicit_prob,
    FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true),
    callback = approximate_time_to_finish_cb,
)
