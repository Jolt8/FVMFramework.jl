#NOTE:
#This file was written entirely by AI
#the reason I did this is because I wanted to test if my framework was fast enough to actually handle navier stokes and from what I'm seeing, it definitely is!!!
#a problem with 100 x 100 cells (10000 total) only took 114 seconds to solve for 100000.0 simulated seconds of a lid-driven cavity simulation which is in the ballpark of OpenFOAM
#thus, while I'm not going to work on navier stokes further personally because I really don't want to sell my soul to the devil of getting navier stokes to work and be physically consistent,
#and also because navier stokes doesn't really seem that useful in chemical engineering for a majority of problems (and is also often too slow even with extremely efficient solvers for parameter optimization), 
#I at least know it's possible if I want to do it in the future.

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
Ny = 10
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


# Fluid parameters (Re = 100)
Re = 100.0
U_lid = 1.0u"m/s"
nu_val = (U_lid * L) / Re
rho_val = 1.0u"kg/m^3"
beta_val = 10.0u"Pa"

dx_val = L / Nx
dy_val = L / Ny
dz_val = H / Nz

fluid_properties = ComponentVector(
    rho = rho_val,
    nu = nu_val,
    beta = beta_val,
    dx = dx_val,
    dy = dy_val,
    dz = dz_val,
)


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

# Connection flux (internal faces)
function fluid_fluid_flux!(
    du,
    u,
    p,
    t,
    idx_a,
    idx_b,
    face_idx,
    cell_neighbor_areas,
    cell_neighbor_normals,
    cell_neighbor_distances,
    cell_volumes
)
    area = cell_neighbor_areas[idx_a][face_idx]
    norm = cell_neighbor_normals[idx_a][face_idx]
    dist = cell_neighbor_distances[idx_a][face_idx]

    # Velocity at cells a and b
    v_a = Ferrite.Vec{3}((u.u_vel[idx_a], u.v_vel[idx_a], u.w_vel[idx_a]))
    v_b = Ferrite.Vec{3}((u.u_vel[idx_b], u.v_vel[idx_b], u.w_vel[idx_b]))

    # Average velocity at face
    v_f = 0.5 * (v_a + v_b)

    rho_eff = 0.5 * (u.rho[idx_a] + u.rho[idx_b])

    # Volumetric flow rate across face from a to b with Rhie-Chow pressure stabilization
    # We subtract a pressure difference term to couple pressure and velocity and prevent checkerboarding
    chi = 0.05 / rho_eff  # stabilization coefficient (can be tuned, typical range is 0.01 - 0.1)
    dp_dn = (u.pressure[idx_b] - u.pressure[idx_a]) / dist
    Q_f = (dot(v_f, norm) - chi * dp_dn) * area

    # Convective flux of velocity using upwind scheme
    if Q_f > 0.0
        F_conv_u = Q_f * u.u_vel[idx_a]
        F_conv_v = Q_f * u.v_vel[idx_a]
        F_conv_w = Q_f * u.w_vel[idx_a]
    else
        F_conv_u = Q_f * u.u_vel[idx_b]
        F_conv_v = Q_f * u.v_vel[idx_b]
        F_conv_w = Q_f * u.w_vel[idx_b]
    end

    # Pressure force on cell a (using physical density)
    p_f = 0.5 * (u.pressure[idx_a] + u.pressure[idx_b])

    F_press_u = (1.0 / rho_eff) * p_f * area * norm[1]
    F_press_v = (1.0 / rho_eff) * p_f * area * norm[2]
    F_press_w = (1.0 / rho_eff) * p_f * area * norm[3]

    # Viscous flux from cell b to cell a
    nu_eff = 0.5 * (u.nu[idx_a] + u.nu[idx_b])

    F_visc_u = nu_eff * area * (u.u_vel[idx_b] - u.u_vel[idx_a]) / dist
    F_visc_v = nu_eff * area * (u.v_vel[idx_b] - u.v_vel[idx_a]) / dist
    F_visc_w = nu_eff * area * (u.w_vel[idx_b] - u.w_vel[idx_a]) / dist

    # Continuity flux (compressibility)
    beta_eff = 0.5 * (u.beta[idx_a] + u.beta[idx_b])
    F_mass = beta_eff * Q_f

    # Accumulate into du for cell a
    du.u_flow[idx_a] -= F_conv_u + F_press_u - F_visc_u
    du.v_flow[idx_a] -= F_conv_v + F_press_v - F_visc_v
    du.w_flow[idx_a] -= F_conv_w + F_press_w - F_visc_w
    du.p_flow[idx_a] -= F_mass
end

# Boundary patch flux
function wall_patch_flux_generic!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals,
    cell_neighbor_distances, cell_volumes,
    u_bc, v_bc, w_bc
)
    area = cell_face_areas[idx_a][face_idx]
    norm = cell_face_normals[idx_a][face_idx]

    # Calculate distance to boundary face based on axis alignment
    dist = 0.0
    if abs(norm[1]) > 0.9
        dist = u.dx[idx_a] / 2.0
    elseif abs(norm[2]) > 0.9
        dist = u.dy[idx_a] / 2.0
    elseif abs(norm[3]) > 0.9
        dist = u.dz[idx_a] / 2.0
    else
        dist = 0.5 * cell_volumes[idx_a]^(1 / 3) # fallback
    end

    # Viscous flux at boundary (drag force)
    F_visc_u = u.nu[idx_a] * area * (u_bc - u.u_vel[idx_a]) / dist
    F_visc_v = u.nu[idx_a] * area * (v_bc - u.v_vel[idx_a]) / dist
    F_visc_w = u.nu[idx_a] * area * (w_bc - u.w_vel[idx_a]) / dist

    # Pressure force from boundary
    F_press_u = (1.0 / u.rho[idx_a]) * u.pressure[idx_a] * area * norm[1]
    F_press_v = (1.0 / u.rho[idx_a]) * u.pressure[idx_a] * area * norm[2]
    F_press_w = (1.0 / u.rho[idx_a]) * u.pressure[idx_a] * area * norm[3]

    # Accumulate
    du.u_flow[idx_a] += F_visc_u - F_press_u
    du.v_flow[idx_a] += F_visc_v - F_press_v
    du.w_flow[idx_a] += F_visc_w - F_press_w
end

# Specialized patch functions (unitless)
stationary_wall_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals,
    cell_neighbor_distances, cell_volumes,
) = wall_patch_flux_generic!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals,
    cell_neighbor_distances, cell_volumes,
    0.0, 0.0, 0.0
)

moving_lid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals,
    cell_neighbor_distances, cell_volumes,
) = wall_patch_flux_generic!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals,
    cell_neighbor_distances, cell_volumes,
    1.0, 0.0, 0.0
)

# Mapping connection functions
function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
end

# Add region and patches
add_region!(
    config,
    "fluid";
    type = Fluid(),
    initial_conditions = ComponentVector(
        u_vel = 0.0u"m/s",
        v_vel = 0.0u"m/s",
        w_vel = 0.0u"m/s",
        pressure = 0.0u"Pa",
    ),
    properties = fluid_properties,
    property_update_function = function update_fluid_properties!(properties, u)
        
    end,
    region_function = function fluid_physics!(du, u, p, t, cell_id, vol)
        du.u_vel[cell_id] += du.u_flow[cell_id] / vol
        du.v_vel[cell_id] += du.v_flow[cell_id] / vol
        du.w_vel[cell_id] += du.w_flow[cell_id] / vol
        du.pressure[cell_id] += du.p_flow[cell_id] / vol
    end,
)

add_patch!(
    config,
    "left_wall";
    properties = ComponentVector(),
    patch_function = stationary_wall_flux!,
)
add_patch!(
    config,
    "right_wall";
    properties = ComponentVector(),
    patch_function = stationary_wall_flux!,
)
add_patch!(
    config,
    "bottom_wall";
    properties = ComponentVector(),
    patch_function = stationary_wall_flux!,
)
add_patch!(
    config,
    "top_lid";
    properties = ComponentVector(),
    patch_function = moving_lid_flux!,
)

# Finish FVM configuration
du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false)

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

@time sol = solve(
    implicit_prob,
    FBDF(linsolve = SparspakFactorization()),
    callback = approximate_time_to_finish_cb,
)
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

u_named = []

du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, geo, p_guess, track_progress = true)

for i in eachindex(u_named)
    step_u = u_named[i]
    vel_vec = [
        [ [step_u.u_vel[c], step_u.v_vel[c], step_u.w_vel[c]] for c = 1:n_cells ]
    ]
    u_named[i] = merge_properties(u_named[i], ComponentVector(velocity = vel_vec, ))

    flow_vec = [
        SVector{3,Float64}(step_u.u_flow[c], step_u.v_flow[c], step_u.w_flow[c]) for
        c = 1:n_cells
    ]
    u_named[i] = merge_properties(u_named[i], ComponentVector(flow = flow_vec, ))
end

root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\julia_cfd_output_files"
println("Saving VTK files to: ", root_dir)
sol_to_vtk(sol, du_named, u_named, grid, geo, @__FILE__, root_dir, track_progress = true)
println("VTK export complete!")
