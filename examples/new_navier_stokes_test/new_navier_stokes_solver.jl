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


grid_x_length = 1.0
grid_y_length = 0.1
grid_z_length = 0.1

grid_dimensions = (100, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((grid_x_length, grid_y_length, grid_z_length))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

n_cells = length(grid.cells)
n_faces = nfacets(grid.cells[1])

# Grid sets
addcellset!(grid, "fluid", xyz -> true)
getcellset(grid, "fluid")

addfacetset!(grid, "supersonic_inlet", xyz -> abs(xyz[1] - 0.0) < 1e-8)
addfacetset!(grid, "supersonic_outlet", xyz -> abs(xyz[1] - grid_x_length) < 1e-8)
addfacetset!(grid, "y_min_wall", xyz -> abs(xyz[2] - 0.0) < 1e-8)
addfacetset!(grid, "y_max_wall", xyz -> abs(xyz[2] - grid_y_length) < 1e-8)
addfacetset!(grid, "z_min_wall", xyz -> abs(xyz[3] - 0.0) < 1e-8)
addfacetset!(grid, "z_max_wall", xyz -> abs(xyz[3] - grid_z_length) < 1e-8)

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

        #flows for later integrating fluxes
        density_flow = u"kg/s",
        momentum_density_u_flow = u"N",
        momentum_density_v_flow = u"N",
        momentum_density_w_flow = u"N",
        volumetric_energy_flow = u"W",
        
        #specific_internal_energy = u"J/kg", #not cached right now, just a property
    ),
    special_caches = ComponentVector(
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(),
)

struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "face_reconstructors/first_order_face_reconstruction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC.jl"))

function HLLC_closure!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals, cell_face_distances,
    cell_neighbor_normals, cell_neighbor_distances, 
    cell_volumes
)
    HLLC!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_face_areas, cell_face_normals, cell_face_distances,
        cell_neighbor_normals, cell_neighbor_distances, 
        cell_volumes,
        first_order_face_reconstruction!
    )
end

function fluid_fluid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals, cell_face_distances,
    cell_neighbor_normals, cell_neighbor_distances, 
    cell_volumes
)
    HLLC_closure!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_face_areas[idx_a][face_idx], cell_face_normals[idx_a][face_idx], cell_face_distances[idx_a][face_idx],
        cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a],
    )
end

# Mapping connection functions
#IMPORTANT: when we switch to evaluating only one flux update per face rather than two, this will determine which idx_a will be and idx_b will be
function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
end

Revise.includet(joinpath(@__DIR__, "navier_stokes_fluid_property_update_functions/fluid_property_update_functions.jl"))
Revise.includet(joinpath(@__DIR__, "sum_and_cap_functions/cap_functions.jl"))

function construct_initial_conditions_from_intuitive_inputs(u)
    momentum_density_u = u.density * u.vel_u
    momentum_density_v = u.density * u.vel_v
    momentum_density_w = u.density * u.vel_w

    internal_energy = u.cv * u.temperature
    kinetic_energy = 0.5 * (u.vel_u^2 + u.vel_v^2 + u.vel_w^2)
    
    volumetric_energy = u.density * (internal_energy + kinetic_energy)

    pressure = (u.R_gas / u.mw) * u.density * u.temperature

    @show speed_of_sound = sqrt((u.cp / u.cv) * pressure / u.density) |> u"m/s"

    return ComponentVector(
        density = u.density,
        momentum_density_u = momentum_density_u,
        momentum_density_v = momentum_density_v,
        momentum_density_w = momentum_density_w,
        volumetric_energy = volumetric_energy
    ), ComponentVector(
        cp = u.cp,
        cv = u.cv,
        R_gas = u.R_gas,
        mw = u.mw
    )
end

fluid_initial_conditions, fluid_properties = construct_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 0.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = 0.0u"m/s",
        density = 1.18u"kg/m^3",
        temperature = 300.0u"K",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol"
    )
)

# Add region and patches
add_region!(
    config, "fluid";
    type = Fluid(),
    initial_conditions = fluid_initial_conditions,
    properties = fluid_properties,
    property_update_function = 
    function update_fluid_properties!(du, u, p, t, system, cell_id, vol, system)
        overall_navier_stokes_property_update!(du, u, p, t, cell_id, vol)
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, cell_id, vol)
        cap_navier_stokes_flow!(du, u, p, t, cell_id, vol)
    end,
)

supersonic_inlet_initial_conditions, supersonic_inlet_properties = construct_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 600.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = 0.0u"m/s",
        density = 1.18u"kg/m^3",
        temperature = 300.0u"K",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol"
    )
)

supersonic_inlet_initial_conditions_stripped = ustrip.(upreferred.(supersonic_inlet_initial_conditions))

add_patch!(
    config, "supersonic_inlet";
    properties = ComponentVector(),
    patch_function = 
    function supersonic_inlet_flux!(
        du, u, p, t, system, 
        idx_a, idx_b, face_idx,
        cell_face_areas, cell_face_normals, cell_face_distances,
        cell_neighbor_normals, cell_neighbor_distances, 
        cell_volumes
    )
        area = cell_face_areas[idx_a][face_idx]
        normal = cell_face_normals[idx_a][face_idx]

        density_bc = supersonic_inlet_initial_conditions_stripped.density
        momentum_density_u_bc = supersonic_inlet_initial_conditions_stripped.momentum_density_u
        momentum_density_v_bc = supersonic_inlet_initial_conditions_stripped.momentum_density_v
        momentum_density_w_bc = supersonic_inlet_initial_conditions_stripped.momentum_density_w
        volumetric_energy_bc = supersonic_inlet_initial_conditions_stripped.volumetric_energy

        gamma_bc = u.cp[idx_a] / u.cv[idx_a]

        _, _, _, pressure_bc, _ =
            primitive_from_conservative(
                density_bc,
                momentum_density_u_bc,
                momentum_density_v_bc,
                momentum_density_w_bc,
                volumetric_energy_bc,
                gamma_bc,
            )

        (
            F_density,
            F_momentum_density_u,
            F_momentum_density_v,
            F_momentum_density_w,
            F_volumetric_energy,
        ) = physical_flux(
            density_bc,
            momentum_density_u_bc,
            momentum_density_v_bc,
            momentum_density_w_bc,
            volumetric_energy_bc,
            pressure_bc,
            normal,
        )

        du.density_flow[idx_a] -= area * F_density

        du.momentum_density_u_flow[idx_a] -= area * F_momentum_density_u
        du.momentum_density_v_flow[idx_a] -= area * F_momentum_density_v
        du.momentum_density_w_flow[idx_a] -= area * F_momentum_density_w

        du.volumetric_energy_flow[idx_a] -= area * F_volumetric_energy
    end
)

add_patch!(
    config, "supersonic_outlet";
    properties = ComponentVector(),
    patch_function = 
    function supersonic_outlet_flux!(
        du, u, p, t, system, 
        idx_a, idx_b, face_idx,
        cell_face_areas, cell_face_normals, cell_face_distances,
        cell_neighbor_normals, cell_neighbor_distances, 
        cell_volumes
    )
        area = cell_face_areas[idx_a][face_idx]
        normal = cell_face_normals[idx_a][face_idx]

        gamma = u.cp[idx_a] / u.cv[idx_a]

        _, _, _, pressure, _ =
            primitive_from_conservative(
                u.density[idx_a],
                u.momentum_density_u[idx_a],
                u.momentum_density_v[idx_a],
                u.momentum_density_w[idx_a],
                u.volumetric_energy[idx_a],
                gamma,
            )

        (
            F_density,
            F_momentum_density_u,
            F_momentum_density_v,
            F_momentum_density_w,
            F_volumetric_energy,
        ) = physical_flux(
            u.density[idx_a],
            u.momentum_density_u[idx_a],
            u.momentum_density_v[idx_a],
            u.momentum_density_w[idx_a],
            u.volumetric_energy[idx_a],
            pressure,
            normal,
        )

        du.density_flow[idx_a] -= area * F_density

        du.momentum_density_u_flow[idx_a] -= area * F_momentum_density_u
        du.momentum_density_v_flow[idx_a] -= area * F_momentum_density_v
        du.momentum_density_w_flow[idx_a] -= area * F_momentum_density_w

        du.volumetric_energy_flow[idx_a] -= area * F_volumetric_energy
    end
)

for name in ["y_min_wall", "y_max_wall", "z_min_wall", "z_max_wall"]
    add_patch!(
        config, name;
        properties = ComponentVector(),
        patch_function = 
        function wall_patch_flux!(
            du, u, p, t, system,
            idx_a, idx_b, face_idx,
            cell_face_areas, cell_face_normals, cell_face_distances,
            cell_neighbor_normals, cell_neighbor_distances, 
            cell_volumes
        )
            gamma = u.cp[idx_a] / u.cv[idx_a]

            _, _, _, pressure, _ = primitive_from_conservative(
                u.density[idx_a],
                u.momentum_density_u[idx_a],
                u.momentum_density_v[idx_a],
                u.momentum_density_w[idx_a],
                u.volumetric_energy[idx_a],
                gamma,
            )

            du.momentum_density_u_flow[idx_a] -= cell_face_areas[idx_a][face_idx] * pressure * cell_face_normals[idx_a][face_idx][1]
            du.momentum_density_v_flow[idx_a] -= cell_face_areas[idx_a][face_idx] * pressure * cell_face_normals[idx_a][face_idx][2]
            du.momentum_density_w_flow[idx_a] -= cell_face_areas[idx_a][face_idx] * pressure * cell_face_normals[idx_a][face_idx][3]
        end
    )
end

#I think we're going to add an additional_data field of the system so that arbitrary data can be passed into any function without allocations
additional_data = (
    #You could put a neural network in here, custom structs, interpolations, just any useful data in general. 
    #This does require that every function have system in its arguments though
)

# Finish FVM configuration
du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, additional_data, check_units = false);

# System solver function
function solve_system!(du, u, p, t, system, geo)
    update_region_groups!(du, u, p, t, system, geo)
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

#=
function f_closure_implicit(du, u, p, t)
    fvm_operator!(du, u, p, t, system, geo, solve_system!)
end
=#

p_guess = 0.0

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

# Detect Jacobian sparsity for NonlinearProblem
jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
) 
#this scales absolutely abysmally as the number of cells goes up
#for example, a 10x increase in the amount of cells made this take around 100x longer! 


function state_is_invalid(u, p, t, system)
    u_named = ComponentVector(u, system.state_axes)

    for i in 1:n_cells
        kinetic_energy_density = 0.5 * (
            u_named.momentum_density_u[i]^2 +
            u_named.momentum_density_v[i]^2 +
            u_named.momentum_density_w[i]^2
        ) / u_named.density[i]
        
        internal_energy_density = u_named.volumetric_energy[i] - kinetic_energy_density

        if !all(
            isfinite, (
                u_named.density[i],
                u_named.momentum_density_u[i],
                u_named.momentum_density_v[i],
                u_named.momentum_density_w[i],
                u_named.volumetric_energy[i],
                internal_energy_density,
            )) || u_named.density[i] <= 0.0 || internal_energy_density <= 0.0
            return true
        end
    end

    return false
end

state_is_invalid_closure = (u, p, t) -> state_is_invalid(u, p, t, system);

#transient
t0 = 0.0
tMax = 1000.0
tspan = (t0, tMax)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

early_dtmax = 100.0
late_dtmax  = 1000.0
switch_time = 400.0

function increase_dtmax!(integrator)
    integrator.opts.dtmax = late_dtmax

    # We changed solver settings, not the state vector.
    u_modified!(integrator, false)

    println("Raised dtmax to $late_dtmax at t = $(integrator.t)")
end

increase_dtmax_cb = DiscreteCallback(
    (u, t, integrator) -> t >= switch_time && integrator.opts.dtmax < late_dtmax,
    increase_dtmax!;
    save_positions = (false, false),
)

callbacks = CallbackSet(
    approximate_time_to_finish_cb,
    increase_dtmax_cb,
)

println("Solving the Navier-Stokes ODE system...")
@time sol = solve(
    implicit_prob,
    #FBDF(linsolve = SparspakFactorization(), nlsolve = NLNewton(relax = 0.5)),
    FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true),
    #FBDF(linsolve = KrylovJL_GMRES(), nlsolve = NLNewton(relax = 0.7), precs = iluzero, concrete_jac = true),
    callback = callbacks,
    isoutofdomain = state_is_invalid_closure,
    #saveat = (tMax / 300)
    #dtmax = 100
    #dtmax = early_dtmax
)

f_closure_steady = (du, u, p) -> f_closure_implicit(du, u, p, 0.0)

nl_func = NonlinearFunction(f_closure_steady, jac_prototype = float.(jac_sparsity))

prob = NonlinearProblem(nl_func, u0_vec, p_guess)

@time sol_steady = solve(prob, NonlinearSolve.NewtonRaphson(concrete_jac = true))

du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, geo, p_guess, track_progress = false);

root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\julia_cfd_output_files"
println("Saving VTK files to: ", root_dir)
sol_to_vtk(sol, du_named, u_named, grid, geo, @__FILE__, root_dir, track_progress = false)
println("VTK export complete!")

