using Revise
using Logging
using Unitful
using OrdinaryDiffEq
using NonlinearSolve
using Sparspak
using Ferrite
using SparseConnectivityTracer
using ForwardDiff
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
    volumetric_energy = zeros(n_cells)u"J/m^3",
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
        k = u"W/(m*K)",
        
        #specific_internal_energy = u"J/kg", #not cached right now, just a property
    ),
    special_caches = ComponentVector(
        grad_density = zeros(n_cells, 3)u"kg/m^4",
        grad_pressure = zeros(n_cells, 3)u"Pa/m",
        grad_vel_u = zeros(n_cells, 3)u"m/(s*m)",
        grad_vel_v = zeros(n_cells, 3)u"m/(s*m)",
        grad_vel_w = zeros(n_cells, 3)u"m/(s*m)",
        grad_temperature = zeros(n_cells, 3)u"K/m",

        
        vel_u_face = zeros(n_cells, n_faces)u"m/s",
        vel_v_face = zeros(n_cells, n_faces)u"m/s",
        vel_w_face = zeros(n_cells, n_faces)u"m/s",
        temperature_face = zeros(n_cells, n_faces)u"J/m^3",
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(),
)

struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "face_reconstructors/first_order_face_reconstruction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC_low_mach_correction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC.jl"))
Revise.includet(joinpath(@__DIR__, "weighted_least_squares/weighted_least_squares.jl"))
Revise.includet(joinpath(@__DIR__, "face_reconstructors/MUSCL_face_reconstruction.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/fluid_viscous_and_diffusive_fluxes.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/wall_viscous_and_diffusive_fluxes.jl"))

const USE_THORNBER = "--thornber" in ARGS
const LOW_MACH_CORRECTION = if USE_THORNBER
    thornber_low_mach_correction
else
    no_low_mach_correction
end
const USE_NO_SLIP_WALLS = "--no-slip-walls" in ARGS

function fluid_fluid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    #=
    populate_weighted_least_squares_face_values!(
        du, u, p, t, system, geo,
        idx_a, face_a,
        idx_b, face_b
    ) #this is to update vel_u_face, vel_v_face, vel_w_face, and temperature_face
    #IMPORTANT: if more face values are needed, this function needs to be changed
    =#
    #NOTE: we only need this if we get tired of calculating grad_u_face values inside different functions
    #right now, all grad_u_face values are only required in fluid_viscous_and_diffusive_flux!

    #I think the difficulties with MUSCL_face_reconsturction! have something to do with the inlet cell, it always seems to be the 
    #one with non-physical values

    HLLC!(
        du, u, p, t, system, geo,
        idx_a, face_a, 
        idx_b, face_b,
        MUSCL_face_reconstruction!,
        #first_order_face_reconstruction!,
        no_low_mach_correction,
    )
    #HLLC!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b, first_order_face_reconstruction!)

    fluid_viscous_and_diffusive_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
end

# Mapping connection functions
#IMPORTANT: when we switch to evaluating only one flux update per face rather than two, this will determine which idx_a will be and idx_b will be
#This makes sense because if we have two fluids that have the same type (flux functions should be the same) they will both be updated randomly
function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
end

Revise.includet(joinpath(@__DIR__, "navier_stokes_fluid_property_update_functions/fluid_property_update_functions.jl"))
Revise.includet(joinpath(@__DIR__, "navier_stokes_fluid_property_update_functions/non_standard_fluid_property_update_functions.jl"))
Revise.includet(joinpath(@__DIR__, "sum_and_cap_functions/cap_functions.jl"))

function construct_initial_conditions_from_intuitive_inputs(u)
    momentum_density_u = u.density * u.vel_u
    momentum_density_v = u.density * u.vel_v
    momentum_density_w = u.density * u.vel_w

    internal_energy = u.cv * u.temperature
    kinetic_energy = 0.5 * (u.vel_u^2 + u.vel_v^2 + u.vel_w^2)
    
    volumetric_energy = u.density * (internal_energy + kinetic_energy)

    R_specific = u.cp - u.cv
    pressure = u.density * R_specific * u.temperature

    #@show speed_of_sound = sqrt((u.cp / u.cv) * pressure / u.density) |> u"m/s"

    return ComponentVector(
        density = u.density,
        momentum_density_u = momentum_density_u,
        momentum_density_v = momentum_density_v,
        momentum_density_w = momentum_density_w,
        volumetric_energy = volumetric_energy,
    ), ComponentVector(
        cp = u.cp,
        cv = u.cv,
        R_gas = u.R_gas,
        mw = u.mw,
        mu = u.mu, 
        prandtl_number = u.prandtl_number,
    )
end

fluid_initial_conditions, fluid_properties = construct_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 1.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = 0.0u"m/s",
        density = 1.18u"kg/m^3",
        temperature = 300.0u"K",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol",
        mu = 1e-5u"Pa*s",
        #k = 0.026u"W/(m*K)",
        prandtl_number = 0.705,
    )
)

# Add region and patches
add_region!(
    config, "fluid";
    type = Fluid(),
    initial_conditions = fluid_initial_conditions,
    properties = fluid_properties,
    property_update_function = 
    function update_fluid_properties!(du, u, p, t, system, geo, cell_id)
        update_k_from_prandtl!(du, u, p, t, system, geo, cell_id)
        overall_navier_stokes_property_update!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, system, geo, cell_id)
        cap_navier_stokes_flow!(du, u, p, t, system, geo, cell_id)
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
        mw = 28.97u"g/mol",
        mu = 1e-5u"Pa*s",
        #k = 0.026u"W/(m*K)",
        prandtl_number = 0.705,
    )
)

supersonic_inlet_initial_conditions_stripped = ustrip.(upreferred.(supersonic_inlet_initial_conditions))

add_patch!(
    config, "supersonic_inlet";
    #type_a = Fluid(),
    #type_b = nothing,
    properties = ComponentVector(),
    patch_function = 
    function supersonic_inlet_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, 
        idx_b, face_b,
    )
        face_area, face_normal, face_distance, vol = boundary_geometry(geo, idx_a, face_a)

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
            face_normal,
        )

        du.density_flow[idx_a] -= face_area * F_density

        du.momentum_density_u_flow[idx_a] -= face_area * F_momentum_density_u
        du.momentum_density_v_flow[idx_a] -= face_area * F_momentum_density_v
        du.momentum_density_w_flow[idx_a] -= face_area * F_momentum_density_w

        du.volumetric_energy_flow[idx_a] -= face_area * F_volumetric_energy
    end
)

add_patch!(
    config, "supersonic_outlet";
    #type_a = Fluid(),
    #type_b = nothing,
    properties = ComponentVector(),
    patch_function = 
    function supersonic_outlet_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, 
        idx_b, face_b,
    )
        face_area, face_normal, face_distance, vol = boundary_geometry(geo, idx_a, face_a)

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
            face_normal
        )

        du.density_flow[idx_a] -= face_area * F_density

        du.momentum_density_u_flow[idx_a] -= face_area * F_momentum_density_u
        du.momentum_density_v_flow[idx_a] -= face_area * F_momentum_density_v
        du.momentum_density_w_flow[idx_a] -= face_area * F_momentum_density_w

        du.volumetric_energy_flow[idx_a] -= face_area * F_volumetric_energy
    end
)

for name in ["y_min_wall", "y_max_wall", "z_min_wall", "z_max_wall"]
    add_patch!(
        config, name;
        #type_a = Fluid(),
        #type_b = nothing,
        properties = ComponentVector(),
        patch_function = 
        function wall_patch_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a, 
            idx_b, face_b,
        )
            face_area, face_normal, face_distance, vol = boundary_geometry(geo, idx_a, face_a)

            gamma = u.cp[idx_a] / u.cv[idx_a]

            _, _, _, pressure, _ = primitive_from_conservative(
                u.density[idx_a],
                u.momentum_density_u[idx_a],
                u.momentum_density_v[idx_a],
                u.momentum_density_w[idx_a],
                u.volumetric_energy[idx_a],
                gamma,
            )

            du.momentum_density_u_flow[idx_a] -= face_area * pressure * face_normal[1]
            du.momentum_density_v_flow[idx_a] -= face_area * pressure * face_normal[2]
            du.momentum_density_w_flow[idx_a] -= face_area * pressure * face_normal[3]

            if USE_NO_SLIP_WALLS
                non_moving_wall_viscous_and_diffusive_flux!(
                    du,
                    u,
                    p,
                    t,
                    system,
                    geo,
                    idx_a,
                    face_a,
                    idx_b,
                    face_b,
                )
            end
        end
    ) 
end

#I think we're going to add an additional_data field of the system so that arbitrary data can be passed into any function without allocations
additional_data = (
    #You could put a neural network in here, custom structs, interpolations, just any useful data in general. 
    #This does require that every function have system in its arguments though
)

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, additional_data, check_units = false);
WLS_STENCIL = build_weighted_least_squares_stencil(geo)

# System solver function
function solve_system!(du, u, p, t, system, geo)
    update_region_groups!(du, u, p, t, system, geo)
    update_weighted_least_squares_gradients!(u, WLS_STENCIL)
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

if "--smoke-test" in ARGS
    smoke_test_derivative = similar(u0_vec)
    f_closure_implicit(smoke_test_derivative, u0_vec, 0.0, 0.0)
    @assert length(u0_vec) == 5 * n_cells
    @assert all(isfinite, u0_vec)
    @assert all(isfinite, smoke_test_derivative)
    println("No-species smoke test passed with $(length(u0_vec)) Euler degrees of freedom.")
    exit(0)
end

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
        density = u_named.density[i]
        if !isfinite(density) || density <= 0.0
            return true
        end
        kinetic_energy_density = 0.5 * (
            u_named.momentum_density_u[i]^2 +
            u_named.momentum_density_v[i]^2 +
            u_named.momentum_density_w[i]^2
        ) / density
        
        internal_energy_density = u_named.volumetric_energy[i] - kinetic_energy_density

        if !all(
            isfinite, (
                density,
                u_named.momentum_density_u[i],
                u_named.momentum_density_v[i],
                u_named.momentum_density_w[i],
                u_named.volumetric_energy[i],
                internal_energy_density,
            )) || internal_energy_density <= 0.0
            return true
        end
    end

    return false
end

state_is_invalid_closure = (u, p, t) -> state_is_invalid(u, p, t, system);

#transient
t0 = 0.0
function command_line_float(prefix, default)
    argument = findfirst(arg -> startswith(arg, prefix), ARGS)
    return isnothing(argument) ? default : parse(Float64, split(ARGS[argument], "=", limit = 2)[2])
end

tMax = command_line_float("--tmax=", 1000.0)
tMax = 3.0
transient_abstol = command_line_float("--abstol=", 1.0e-6)
transient_reltol = command_line_float("--reltol=", 1.0e-4)
tspan = (t0, tMax)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

progress_interval = max((tMax - t0) / 20.0, eps(Float64))
next_progress_time = Ref(t0 + progress_interval)
progress_callback = DiscreteCallback(
    (u, t, integrator) -> t >= next_progress_time[],
    integrator -> begin
        println("sim time: $(integrator.t)s")
        next_progress_time[] += progress_interval
    end;
    save_positions = (false, false),
)

transient_algorithm = Rosenbrock23(
    autodiff = ADTypes.AutoFiniteDiff(),
    linsolve = SparspakFactorization(),
)
save_times = range(t0, tMax; length = 101)
transient_solve_kwargs = (
    #abstol = transient_abstol,
    #reltol = transient_reltol,
    #dt = min(1e-5, tMax - t0),
    isoutofdomain = state_is_invalid_closure,
    #saveat = save_times,
    save_everystep = false,
    maxiters = 1_000_000,
)

@time sol = solve(
    implicit_prob,
    FBDF(linsolve = KrylovJL_GMRES(), autodiff = ADTypes.AutoFiniteDiff()),
    #transient_algorithm;
    #transient_solve_kwargs...,
    callback = approximate_time_to_finish_cb,
)

sol.alg

sol = if "--progress" in ARGS
    @time solve(
        implicit_prob,
        transient_algorithm;
        transient_solve_kwargs...,
        callback = approximate_time_to_finish_cb,
    )
else
    @time solve(implicit_prob, transient_algorithm; transient_solve_kwargs...)
end

println("Transient algorithm: Rosenbrock23(AutoFiniteDiff, Sparspak)")
println("HLLC reconstruction: $(USE_THORNBER ? "Thornber low-Mach" : "ordinary")")
println("Wall model: $(USE_NO_SLIP_WALLS ? "no-slip viscous" : "slip adiabatic")")
println("Transient return code: $(sol.retcode), final time: $(sol.t[end])")
println(
    "Transient steps: accepted=$(sol.destats.naccept), rejected=$(sol.destats.nreject), " *
    "nonlinear failures=$(sol.destats.nnonlinconvfail)",
)
final_transient_state = ComponentVector(sol.u[end], system.state_axes)
final_transient_velocity_u =
    final_transient_state.momentum_density_u ./ final_transient_state.density
println(
    "Final axial velocity range: " *
    "$(extrema(final_transient_velocity_u)) m/s",
)
if !OrdinaryDiffEq.SciMLBase.successful_retcode(sol)
    error("transient solve failed with return code $(sol.retcode)")
end

if !("--steady" in ARGS)
    exit(0)
end

if "--jacobian-smoke-test" in ARGS
    direction = max.(abs.(u0_vec), 1.0)
    jacobian_vector_product = ForwardDiff.jacobian([0.0]) do epsilon
        perturbed_state = u0_vec .+ epsilon[1] .* direction
        derivative = similar(perturbed_state)
        f_closure_implicit(derivative, perturbed_state, 0.0, 0.0)
        return derivative
    end
    @assert all(isfinite, jacobian_vector_product)
    println("No-species ForwardDiff Jacobian-vector smoke test passed.")
    exit(0)
end

f_closure_steady = (du, u, p) -> f_closure_implicit(du, u, p, 0.0)

inlet_density = supersonic_inlet_initial_conditions_stripped.density
inlet_momentum_scale = max(
    abs(supersonic_inlet_initial_conditions_stripped.momentum_density_u),
    inlet_density * sqrt(
        (fluid_properties.cp / fluid_properties.cv) *
        supersonic_inlet_initial_conditions_stripped.volumetric_energy /
        inlet_density,
    ),
)
inlet_energy = supersonic_inlet_initial_conditions_stripped.volumetric_energy
state_scales = similar(u0_vec)
named_state_scales = ComponentVector(state_scales, system.state_axes)
named_state_scales.density .= inlet_density
named_state_scales.momentum_density_u .= inlet_momentum_scale
named_state_scales.momentum_density_v .= inlet_momentum_scale
named_state_scales.momentum_density_w .= inlet_momentum_scale
named_state_scales.volumetric_energy .= inlet_energy

function scaled_steady_residual!(scaled_residual, scaled_state, p)
    physical_state = scaled_state .* state_scales
    if state_is_invalid(physical_state, p, 0.0, system)
        fill!(scaled_residual, 1.0e6)
        return nothing
    end
    physical_residual = similar(physical_state)
    f_closure_steady(physical_residual, physical_state, p)
    scaled_residual .= .-physical_residual ./ state_scales
    return nothing
end

uniform_inlet_state = similar(u0_vec)
named_uniform_inlet_state = ComponentVector(uniform_inlet_state, system.state_axes)
named_uniform_inlet_state.density .=
    supersonic_inlet_initial_conditions_stripped.density
named_uniform_inlet_state.momentum_density_u .=
    supersonic_inlet_initial_conditions_stripped.momentum_density_u
named_uniform_inlet_state.momentum_density_v .=
    supersonic_inlet_initial_conditions_stripped.momentum_density_v
named_uniform_inlet_state.momentum_density_w .=
    supersonic_inlet_initial_conditions_stripped.momentum_density_w
named_uniform_inlet_state.volumetric_energy .=
    supersonic_inlet_initial_conditions_stripped.volumetric_energy

steady_initial_state = if USE_NO_SLIP_WALLS
    sol.u[end]
else
    uniform_inlet_state
end
scaled_steady_initial_state = steady_initial_state ./ state_scales
initial_scaled_residual = similar(scaled_steady_initial_state)
scaled_steady_residual!(
    initial_scaled_residual,
    scaled_steady_initial_state,
    p_guess,
)
println("Initial scaled steady residual (Linf): $(norm(initial_scaled_residual, Inf))")

nl_func = NonlinearFunction(
    scaled_steady_residual!;
    jac_prototype = float.(jac_sparsity),
)
nl_prob = NonlinearLeastSquaresProblem(nl_func, scaled_steady_initial_state, p_guess)
steady_algorithm = LevenbergMarquardt(
    autodiff = ADTypes.AutoFiniteDiff(),
    linsolve = SparspakFactorization(),
)
sol_steady = @time solve(
    nl_prob,
    steady_algorithm;
    abstol = 1e-8,
    reltol = 1e-8,
    maxiters = 500,
)
steady_state = sol_steady.u .* state_scales
final_scaled_residual = similar(sol_steady.u)
scaled_steady_residual!(final_scaled_residual, sol_steady.u, p_guess)
println("Steady return code: $(sol_steady.retcode)")
println("Final scaled steady residual (Linf): $(norm(final_scaled_residual, Inf))")
println("Steady state admissible: $(!state_is_invalid(steady_state, p_guess, 0.0, system))")
if !OrdinaryDiffEq.SciMLBase.successful_retcode(sol_steady)
    error("steady solve failed with return code $(sol_steady.retcode)")
end
#=
du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, geo, p_guess, track_progress = false);

root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\julia_cfd_output_files"

sol_to_vtk(sol, du_named, u_named, grid, geo, @__FILE__, root_dir, track_progress = false)
=#