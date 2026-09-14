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

# Parameters
L = 1.0u"m"
H = 0.1u"m"
Nx = 100
Ny = 1000
Nz = 1

grid_dimensions = (Nx, Ny, Nz)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((ustrip(L), ustrip(L), ustrip(H)))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

# Grid sets
addcellset!(grid, "fluid", xyz -> true)
addfacetset!(grid, "left_wall", xyz -> abs(xyz[1] - 0.0) < 1e-6)
addfacetset!(grid, "right_wall", xyz -> abs(xyz[1] - ustrip(L)) < 1e-6)
addfacetset!(grid, "bottom_wall", xyz -> abs(xyz[2] - 0.0) < 1e-6)
addfacetset!(grid, "top_lid", xyz -> abs(xyz[2] - ustrip(L)) < 1e-6)

n_cells = length(grid.cells)

u_proto = ComponentVector(
    u_vel = zeros(n_cells)u"m/s",
    v_vel = zeros(n_cells)u"m/s",
    w_vel = zeros(n_cells)u"m/s",
    pressure = zeros(n_cells)u"Pa",
)

config = create_fvm_config(grid, u_proto);

add_setup_syms!(
    config;
    cache_syms_and_units = (
        u_flow = u"m^4/s^2",
        v_flow = u"m^4/s^2",
        w_flow = u"m^4/s^2",
        p_flow = u"Pa*m^3/s",
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
        u_vel = 0.0u"m/s",
        v_vel = 0.0u"m/s",
        w_vel = 0.0u"m/s",
        pressure = 0.0u"Pa"
    ),
    properties = ComponentVector(
        rho = 1.0u"kg/m^3",
        nu = (U_lid * L) / Re,
        beta = 10.0u"Pa"
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

#1718 seconds (12.44 M allocations, 55.6 GiB, 1.86 gc time, 0.21 compilation time) for 100 * 1000 cells
#=
@time sol = solve(
    implicit_prob,
    FBDF(linsolve = SparspakFactorization()),
    callback = approximate_time_to_finish_cb,
)
    =#
#for 100000 seconds of sim time:
    #- 10000 cells takes: 114 seconds (500.50 k allocations: 4.071 GiB, 1.55% gc time)
    #- 20000 cells takes: 300.863069 seconds (150.22 M allocations: 15.861 GiB, 1.01% gc time, 81.24% compilation time) (ignore the compilation time and allocations, I don't think those actually had an effect on the solve time)
#algebraicmultigrid sucks ass, it got stuck at a sim time of 5.8 seconds and then it just stalled 
#I also tested out FBDF with the default solver and it just stalls at around 6.8 seconds in, so gradient based optimization with navier stokes would likely be impossible :(


#nl_func = NonlinearFunction(f_closure_steady, jac_prototype = float.(nl_jac_sparsity))
#nl_prob = NonlinearProblem(nl_func, u0_vec, p_guess)

#println("Solving for steady state directly using NonlinearSolve NewtonRaphson...")
#@time nl_sol = solve(nl_prob, NewtonRaphson(concrete_jac = true))

println("ODE solving complete. Rebuilding state for VTK output...")

du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, geo, p_guess, track_progress = true)

u_named[200].u_vel[300]

u_named_test = []

for i in eachindex(sol.u)
    sol_state_named = ComponentVector(sol.u[i], system.state_axes)

    vel_vec = SVector{3, Float64}[
        SVector{3, Float64}(
            #why must we do this cursedness just to get it to look right?
            sol_state_named.w_vel[c],
            -sol_state_named.v_vel[c], 
            sol_state_named.u_vel[c], 
        ) for c in 1:n_cells
    ]

    vel_matrix = reduce(hcat, vel_vec)
    vel_matrix = transpose(vel_matrix)

    if i == 200
        @show vel_vec[300]
        @show vel_matrix[300, :] 
    end
    
    u_named[i] = merge_properties(u_named[i], ComponentVector(velocity = vel_matrix,))

    #=
    step_u = ComponentVector(sol.u[i], system.state_axes)
    vel_vec = [
        SVector{3,Float64}(
            getproperty(step_u, :u_vel)[c],
            getproperty(step_u, :v_vel)[c],
            getproperty(step_u, :w_vel)[c]
        ) for c = 1:n_cells
    ]
    push!(u_named_test, ComponentVector(velocity = vel_vec, pressure = step_u.pressure))
    =#

    #u_named[i] = merge_properties(u_named[i], u_named_test[i])
    
    #merge_properties(u_named[i], (velocity = vel_vec, pressure = step_u.pressure))

    #=
    flow_vec = [SVector{3, Float64}(sol_state_named.u_flow[c], sol_state_named.v_flow[c], sol_state_named.w_flow[c]) for c = 1:n_cells]
    u_named[i] = merge_properties(u_named[i], ComponentVector(flow = flow_vec, ))
    =#
end

u_named_test
u_named[200].velocity

du_named = ComponentVector()
u_named = [ComponentVector(sol.u, system.state_axes) for i in eachindex(sol.u)]

du_named = []
u_named = []

for i in eachindex(sol.u)
    push!(du_named, ComponentVector(sol.u[i], system.state_axes))
    push!(u_named, ComponentVector(sol.u[i], system.state_axes))
end

for i in eachindex(sol.u)
    sol_state_named = ComponentVector(sol.u[i], system.state_axes)
    vel_vec = SVector{3, Float64}[
        SVector{3, Float64}(
            #why must we do this cursedness just to get it to look right?
            sol_state_named.w_vel[c],
            -sol_state_named.v_vel[c], 
            sol_state_named.u_vel[c], 
        ) for c in 1:n_cells
    ]

    vel_matrix = reduce(hcat, vel_vec)
    vel_matrix = transpose(vel_matrix)
    
    u_named[i] = merge_properties(u_named[i], ComponentVector(velocity = vel_matrix,))
end

u_named[1].velocity

root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\julia_cfd_output_files"
println("Saving VTK files to: ", root_dir)
sol_to_vtk(sol, du_named, u_named, grid, geo, @__FILE__, root_dir, track_progress = true)
println("VTK export complete!")
