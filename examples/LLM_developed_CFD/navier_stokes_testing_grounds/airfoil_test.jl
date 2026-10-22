using Revise
using Logging
using Unitful
using OrdinaryDiffEq
using NonlinearSolve
using Sparspak
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using ILUZero
using LinearAlgebra
using StaticArrays
using FVMFramework

mesh_path = joinpath(@__DIR__, "airfoil_mesh/ramjet/ramjet_mesh-Cut.msh")
grid = togrid(mesh_path)

n_cells = length(grid.cells)
n_faces = nfacets(grid.cells[1])

"""

"""

compressible_solver_options = (
    use_sst_turbulence = true,
)

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3",
    momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_v = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_w = zeros(n_cells)u"kg/(m^2*s)",
    volumetric_energy = zeros(n_cells)u"J/m^3",
)

if compressible_solver_options.use_sst_turbulence == true
    u_proto = merge_properties(u_proto, ComponentVector(
        turbulent_kinetic_energy_density = zeros(n_cells)u"J/m^3",
        specific_dissipation_rate_density = zeros(n_cells)u"kg/(m^3*s)",
    ))
end

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
        gamma = u"1",
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

        grad_turbulent_kinetic_energy = zeros(n_cells, 3)u"m/s^2",
        grad_specific_dissipation_rate = zeros(n_cells, 3)u"1/(m*s)",
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(),
)

if compressible_solver_options.use_sst_turbulence == true
    config.cache_syms_and_units = merge(config.cache_syms_and_units, (
        turbulent_kinetic_energy_density = u"J/m^3",
        specific_dissipation_rate_density = u"kg/(m^3*s)",

        turbulent_kinetic_energy = u"m^2/s^2",
        specific_dissipation_rate = u"1/s",
        turbulent_viscosity = u"Pa*s",
        effective_dynamic_viscosity = u"Pa*s",
        effective_thermal_conductivity = u"W/(m*K)",
        wall_distance = u"m",
        strain_rate_squared = u"1/s^2",
        sst_F1 = u"1",
        sst_F2 = u"1",
        turbulent_kinetic_energy_density_flow = u"W",
        specific_dissipation_rate_density_flow = u"kg/s^2",
    ))
    config.special_caches = merge_properties(config.special_caches, ComponentVector(
        grad_turbulent_kinetic_energy = zeros(n_cells, 3)u"m/s^2",
        grad_specific_dissipation_rate = zeros(n_cells, 3)u"1/(m*s)",
    ))
end

struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "face_reconstructors/first_order_face_reconstruction.jl"))
Revise.includet(joinpath(@__DIR__, "turbulence/sst_k_omega.jl"))

Revise.includet(joinpath(@__DIR__, "boundary_flux_helpers/boundary_flux_helpers.jl")) # Boundary flux helpers (`physical_flux` and `primitive_from_conservative`).
Revise.includet(joinpath(@__DIR__, "boundary_conditions/boundary_condition_exports.jl")) # Boundary Conditions

Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC_low_mach_correction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/SLAU2.jl"))

Revise.includet(joinpath(@__DIR__, "weighted_least_squares/weighted_least_squares.jl"))

Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/fluid_viscous_and_diffusive_fluxes.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/wall_viscous_and_diffusive_fluxes.jl"))

Revise.includet(joinpath(@__DIR__, "navier_stokes_fluid_property_update_functions/non_standard_fluid_property_update_functions.jl"))
Revise.includet(joinpath(@__DIR__, "navier_stokes_fluid_property_update_functions/fluid_property_update_functions.jl"))

Revise.includet(joinpath(@__DIR__, "sum_and_cap_functions/cap_functions.jl"))

function fluid_fluid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    SLAU2!(
        du, u, p, t, system, geo,
        idx_a, face_a,
        idx_b, face_b,
        first_order_face_reconstruction!,
    )

    #=
    HLLC!(
        du, u, p, t, system, geo,
        idx_a, face_a,
        idx_b, face_b,
        #MUSCL_face_reconstruction!,
        first_order_face_reconstruction!,
        no_low_mach_correction
        #thornber_low_mach_correction,
    )=#

    fluid_viscous_and_diffusive_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
end

# Mapping connection functions
#IMPORTANT: when we switch to evaluating only one flux update per face rather than two, this will determine which idx_a will be and idx_b will be
#This makes sense because if we have two fluids that have the same type (flux functions should be the same) they will both be updated randomly
function connection_map_function(phys_a, phys_b)
    phys_a isa Fluid && phys_b isa Fluid && return fluid_fluid_flux!
    error("No flux is defined between $(typeof(phys_a)) and $(typeof(phys_b))")
end

function construct_initial_conditions_from_intuitive_inputs(u)
    momentum_density_u = u.density * u.vel_u
    momentum_density_v = u.density * u.vel_v
    momentum_density_w = u.density * u.vel_w

    internal_energy = u.cv * u.temperature
    kinetic_energy = 0.5 * (u.vel_u^2 + u.vel_v^2 + u.vel_w^2)
    
    volumetric_energy = u.density * (internal_energy + kinetic_energy)

    properties_comp_vec = ComponentVector(
        cp = u.cp,
        cv = u.cv,
        R_gas = u.R_gas,
        mw = u.mw,
        mu = u.mu, 
        molecular_viscosity = u.mu,
        prandtl_number = u.prandtl_number,
    )

    if system.compressible_solver_options.use_sst_turbulence == true
        turbulent_kinetic_energy_density = u.density * u.turbulent_kinetic_energy
        specific_dissipation_rate_density = u.density * u.specific_dissipation_rate
        return ComponentVector(
            density = u.density,
            momentum_density_u = momentum_density_u,
            momentum_density_v = momentum_density_v,
            momentum_density_w = momentum_density_w,
            volumetric_energy = volumetric_energy,
            turbulent_kinetic_energy_density = turbulent_kinetic_energy_density,
            specific_dissipation_rate_density = specific_dissipation_rate_density,
        ), properties_comp_vec
    else
        return ComponentVector(
            density = u.density,
            momentum_density_u = momentum_density_u,
            momentum_density_v = momentum_density_v,
            momentum_density_w = momentum_density_w,
            volumetric_energy = volumetric_energy,
        ), properties_comp_vec
    end
end

fluid_initial_conditions, fluid_properties = construct_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 800.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = 0.0u"m/s",
        density = 1.18u"kg/m^3",
        temperature = 300.0u"K",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol",
        mu = 1e-5u"Pa*s",
        turbulent_kinetic_energy = 0.1u"m^2/s^2",
        specific_dissipation_rate = 100.0u"1/s",
        #k = 0.026u"W/(m*K)",
        prandtl_number = 0.705,
    )
)

# Add region and patches
add_region!(
    config, "surrounding_air_volume";
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
    end
)

supersonic_inlet_initial_conditions, supersonic_inlet_properties = construct_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 800.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = 0.0u"m/s",
        density = 1.18u"kg/m^3",
        temperature = 300.0u"K",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol",
        mu = 1e-5u"Pa*s",
        turbulent_kinetic_energy = 54.0u"m^2/s^2",
        specific_dissipation_rate = 1900.0u"1/s",
        #k = 0.026u"W/(m*K)",
        prandtl_number = 0.705,
    )
)

supersonic_inlet_initial_conditions_stripped = ustrip.(upreferred.(supersonic_inlet_initial_conditions))

add_patch!(
    config, "inlet_patch";
    #type_a = Fluid(),
    #type_b = nothing,
    properties = ComponentVector(),
    patch_function = function inlet_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        supersonic_inlet_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a,
            idx_b, face_b,
            inlet_conditions = supersonic_inlet_initial_conditions_stripped,
        )
    end
)

add_patch!(
    config, "outlet_patch";
    #type_a = Fluid(),
    #type_b = nothing,
    properties = ComponentVector(),
    patch_function = supersonic_outlet_flux!
)

for name in ["top_patch", "bottom_patch"]
    add_patch!(
        config, name;
        properties = ComponentVector(),
        patch_function = function top_and_bottom_patch_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a, 
            idx_b, face_b,
        )
            characteristic_farfield_flux!(
                du, u, p, t, system, geo,
                idx_a, face_a,
                idx_b, face_b,
                farfield_conditions = supersonic_inlet_initial_conditions_stripped,
            )
        end
    )
end

add_patch!(
    config, "airfoil_surface_patches";
    properties = ComponentVector(),
    patch_function = function airfoil_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        no_slip_adiabatic_wall_flux!(du, u, p, t, system, geo, idx_a, face_a)
    end,
)

# The mesh is one cell thick. These front/back faces are symmetry planes for
# the nominally two-dimensional airfoil calculation.
add_patch!(
    config, "side_patches";
    properties = ComponentVector(),
    patch_function = symmetry_flux!
)

wls_stencil = build_weighted_least_squares_stencil(config.geo)
sst_wall_distances = create_sst_wall_distances(grid, config.geo, 
    [
        "airfoil_surface_patches",
    ]
)

#I think we're going to add an additional_data field of the system so that arbitrary data can be passed into any function without allocations
additional_data = (
    compressible_solver_options,
    wls_stencil,
    sst_wall_distances
)
#There's an overload that allows you to access the additional data like system.field

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, additional_data, check_units = false);

# System solver function
function solve_system!(du, u, p, t, system, geo)
    update_region_groups!(du, u, p, t, system, geo)

    update_weighted_least_squares_gradients!(u, system.wls_stencil)
    
    if system.compressible_solver_options.use_sst_turbulence == true
        update_sst_closure!(u, system.sst_wall_distances)
    end
    
    #update_MUSCL_gradients!(u, wls_stencil) #REMEMBER TO REACTIVATE THIS!!!
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

if "--smoke-test" in ARGS
    smoke_test_derivative = similar(u0_vec)
    f_closure_implicit(smoke_test_derivative, u0_vec, 0.0, 0.0)
    @assert all(isfinite, u0_vec)
    @assert all(isfinite, smoke_test_derivative)
    println("Airfoil SLAU2 smoke test passed with $(length(u0_vec)) degrees of freedom.")
    exit(0)
end

p_guess = 0.0

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

# Detect Jacobian sparsity for NonlinearProblem
jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
) 
#this scales absolutely abysmally as the number of cells goes up
#for example, a 10x increase in the amount of cells made this take around 100x longer! 

Revise.includet(joinpath(@__DIR__, "suspicious_jacobians.jl"))
#check_suspicious_jacobian(jac_sparsity, du0_vec, u0_vec, p_guess, 0.0, system, geo, f_closure_implicit);
#check_phantom_dependencies(jac_sparsity, u0_vec, system);

#transient
t0 = 0.0
switch_time = 0.002
tMax = 0.019

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
startup_prob = ODEProblem(ode_func, u0_vec, (t0, switch_time), p_guess)

callbacks = CallbackSet(
    approximate_time_to_finish_cb,
    #increase_dtmax_cb,
)

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

@time startup_sol = solve(
    startup_prob,
    Tsit5(),
    #AutoTsit5(FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true)),
    #FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff()),
    #AutoTsit5(FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff())),
    #FBDF(linsolve = KrylovJL_GMRES(), autodiff = AutoFiniteDiff()),
    #AutoTsit5(FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true)),
    #FBDF(linsolve = KrylovJL_GMRES(), nlsolve = NLNewton(relax = 0.7), precs = iluzero, concrete_jac = true),
    callback = callbacks,
    isoutofdomain = state_is_invalid_closure,
)

steady_prob = remake(startup_prob; u0 = startup_sol.u[end], tspan = (switch_time, tMax))

@time steady_sol = solve(
    steady_prob,
    #FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true),
    #FBDF(linsolve = KrylovJL_GMRES()),
    #Tsit5(),
    #AutoTsit5(FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true)),
    FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff()),
    #AutoTsit5(FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff())),
    #FBDF(linsolve = KrylovJL_GMRES(), autodiff = AutoFiniteDiff()),
    #AutoTsit5(FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true)),
    #FBDF(linsolve = KrylovJL_GMRES(), nlsolve = NLNewton(relax = 0.7), precs = iluzero, concrete_jac = true),
    callback = callbacks,
    isoutofdomain = state_is_invalid_closure,
)

sol_combined = (u = startup_sol.u, t = startup_sol.t)
sol_combined = (u = [startup_sol.u..., steady_sol.u...], t = [startup_sol.t..., steady_sol.t...])

#f_closure_steady = (du, u, p) -> f_closure_implicit(du, u, p, 0.0)

#nl_func = NonlinearFunction(f_closure_steady, jac_prototype = float.(jac_sparsity))

#prob = NonlinearProblem(nl_func, u0_vec, p_guess)

#@time sol_steady = solve(prob, NonlinearSolve.NewtonRaphson(concrete_jac = true))

#we do this instead of changing the saveat because saveat is pretty linear and doesn't automatically capture fast and slow transient behaviour well

desired_shown_steps = 250
save_fraction = desired_shown_steps / length(sol_combined.t)

reduced_sol = (
    u = [[sol_combined.u[i] for i in 1:(Int(ceil(1/save_fraction))):length(sol_combined.u)]..., sol_combined.u[end]],
    t = [sol_combined.t[1:(Int(ceil(1/save_fraction))):length(sol_combined.t)]..., sol_combined.t[end]]
)

du_named, u_named = regenerate_fvm_state(reduced_sol, system, solve_system!, geo, p_guess, save_properties_and_caches = true, track_progress = true);

add_xyz_vec_to_u_named!(u_named, :velocity, :vel_u, :vel_v, :vel_w)
#we previously didn't need this when vel used to be a single vector of u, v, and w components, since they're separate named fields, we must now do this explicitly, oh well...

root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\julia_cfd_output_files"

sol_to_vtk(reduced_sol, du_named, u_named, grid, geo, @__FILE__, root_dir, include_zeros_fields = false, track_progress = true)
