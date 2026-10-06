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
    turbulent_kinetic_energy_density = zeros(n_cells)u"J/m^3",
    specific_dissipation_rate_density = zeros(n_cells)u"kg/(m^3*s)",
    species_densities = ComponentVector(
        species_a = zeros(n_cells)u"kg/m^3",
        species_b = zeros(n_cells)u"kg/m^3",
    ),
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
        mass_fractions = (
            species_a = u"kg/kg",
            species_b = u"kg/kg",
        ),
        species_density_flow = (
            species_a = u"kg/s",
            species_b = u"kg/s",
        ),
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
        grad_turbulent_kinetic_energy = zeros(n_cells, 3)u"m/s^2",
        grad_specific_dissipation_rate = zeros(n_cells, 3)u"1/(m*s)",

        
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
Revise.includet(joinpath(@__DIR__, "turbulence/sst_k_omega.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC_low_mach_correction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/SLAU2.jl"))
Revise.includet(joinpath(@__DIR__, "weighted_least_squares/weighted_least_squares.jl"))
Revise.includet(joinpath(@__DIR__, "face_reconstructors/MUSCL_face_reconstruction.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/fluid_viscous_and_diffusive_fluxes.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/wall_viscous_and_diffusive_fluxes.jl"))

const USE_SLAU2 = "--slau2" in ARGS

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

    #if USE_SLAU2
        SLAU2!(
            du, u, p, t, system, geo,
            idx_a, face_a,
            idx_b, face_b,
            first_order_face_reconstruction!,
            #MUSCL_face_reconstruction!,
        )
    #else
    #=
        HLLC!(
            du, u, p, t, system, geo,
            idx_a, face_a,
            idx_b, face_b,
            MUSCL_face_reconstruction!,
            #first_order_face_reconstruction!,
            thornber_low_mach_correction,
        )
            =#
    #end
    #HLLC!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b, first_order_face_reconstruction!)

    fluid_viscous_and_diffusive_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    (
        distance,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)
    
    add_conservative_species_diffusion_flux!(
        du,
        u,
        idx_a,
        idx_b,
        face_area_a,
        distance,
    )
    add_sst_diffusion_flux!(
        du,
        u,
        idx_a,
        idx_b,
        face_area_a,
        face_normal_a,
        distance,
    )
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
    turbulent_kinetic_energy_density =
        u.density * u.turbulent_kinetic_energy
    specific_dissipation_rate_density =
        u.density * u.specific_dissipation_rate
    species_densities = ComponentVector(
        species_a = u.density * u.mass_fractions.species_a,
        species_b = u.density * u.mass_fractions.species_b,
    )

    return ComponentVector(
        density = u.density,
        momentum_density_u = momentum_density_u,
        momentum_density_v = momentum_density_v,
        momentum_density_w = momentum_density_w,
        volumetric_energy = volumetric_energy,
        turbulent_kinetic_energy_density = turbulent_kinetic_energy_density,
        specific_dissipation_rate_density = specific_dissipation_rate_density,
        species_densities = species_densities,
    ), ComponentVector(
        cp = u.cp,
        cv = u.cv,
        R_gas = u.R_gas,
        mw = u.mw,
        mu = u.mu, 
        molecular_viscosity = u.mu,
        prandtl_number = u.prandtl_number,
        diffusion_coefficients = u.diffusion_coefficients,
    )
end

fluid_initial_conditions, fluid_properties = construct_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 1e-4u"m/s",
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
        mass_fractions = ComponentVector(
            species_a = 0.25u"kg/kg",
            species_b = 0.75u"kg/kg",
        ),
        diffusion_coefficients = ComponentVector(
            species_a = 1e-5u"m^2/s",
            species_b = 1e-5u"m^2/s",
        ),
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
        update_sst_primitives!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, system, geo, cell_id)
        cap_navier_stokes_flow!(du, u, p, t, system, geo, cell_id)
        cap_sst_transport!(du, u, geo, cell_id)
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
        turbulent_kinetic_energy = 54.0u"m^2/s^2",
        specific_dissipation_rate = 1900.0u"1/s",
        #k = 0.026u"W/(m*K)",
        prandtl_number = 0.705,
        mass_fractions = ComponentVector(
            species_a = 0.8u"kg/kg",
            species_b = 0.2u"kg/kg",
        ),
        diffusion_coefficients = ComponentVector(
            species_a = 1e-5u"m^2/s",
            species_b = 1e-5u"m^2/s",
        ),
    )
)

supersonic_inlet_initial_conditions_stripped = ustrip.(upreferred.(supersonic_inlet_initial_conditions))
supersonic_inlet_mass_fractions_stripped = ComponentVector(
    species_a = fill(
        supersonic_inlet_initial_conditions_stripped.species_densities.species_a /
            supersonic_inlet_initial_conditions_stripped.density,
        n_cells,
    ),
    species_b = fill(
        supersonic_inlet_initial_conditions_stripped.species_densities.species_b /
            supersonic_inlet_initial_conditions_stripped.density,
        n_cells,
    ),
)

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
        add_prescribed_boundary_species_advection_flux!(
            du,
            idx_a,
            face_area,
            F_density,
            supersonic_inlet_mass_fractions_stripped,
        )
        add_prescribed_boundary_sst_advection_flux!(
            du,
            idx_a,
            face_area,
            F_density,
            supersonic_inlet_initial_conditions_stripped.turbulent_kinetic_energy_density /
                supersonic_inlet_initial_conditions_stripped.density,
            supersonic_inlet_initial_conditions_stripped.specific_dissipation_rate_density /
                supersonic_inlet_initial_conditions_stripped.density,
        )
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
        add_boundary_species_advection_flux!(
            du,
            u,
            idx_a,
            face_area,
            F_density,
        )
        add_boundary_sst_advection_flux!(du, u, idx_a, face_area, F_density)
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

            non_moving_wall_viscous_and_diffusive_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
            add_sst_wall_flux!(du, u, geo, idx_a, face_a)
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
SST_WALL_DISTANCES = [
    min(
        geo.cell_centroids[cell_id][2],
        grid_y_length - geo.cell_centroids[cell_id][2],
        geo.cell_centroids[cell_id][3],
        grid_z_length - geo.cell_centroids[cell_id][3],
    ) for cell_id in eachindex(geo.cell_centroids)
]

# System solver function
function solve_system!(du, u, p, t, system, geo)
    update_region_groups!(du, u, p, t, system, geo)
    update_weighted_least_squares_gradients!(u, WLS_STENCIL)
    update_sst_closure!(u, SST_WALL_DISTANCES)
    #update_MUSCL_gradients!(u, WLS_STENCIL) #REMEMBER TO REACTIVATE THIS!!!
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

if "--setup-only" in ARGS
    setup_derivative = similar(u0_vec)
    f_closure_implicit(setup_derivative, u0_vec, 0.0, 0.0)
    if !all(isfinite, setup_derivative)
        error("SST setup smoke check produced a non-finite residual")
    end
    if USE_SLAU2
        println("Inviscid flux: SLAU2 with conservative species and SST advection")
    else
        println("Inviscid flux: HLLC with conservative species and SST advection")
    end
    println("SST setup smoke check passed with $(length(u0_vec)) conservative degrees of freedom")
    exit()
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
                u_named.turbulent_kinetic_energy_density[i],
                u_named.specific_dissipation_rate_density[i],
            )) || u_named.density[i] <= 0.0 || internal_energy_density <= 0.0
            return true
        end
        if u_named.turbulent_kinetic_energy_density[i] <= 0.0 ||
            u_named.specific_dissipation_rate_density[i] <= 0.0
            return true
        end

        species_density_sum = 0.0
        for species_name in propertynames(u_named.species_densities)
            species_density = getproperty(
                u_named.species_densities,
                species_name,
            )[i]
            if !isfinite(species_density) || species_density < 0.0
                return true
            end
            species_density_sum += species_density
        end
        if !isapprox(species_density_sum, u_named.density[i]; atol = 1e-10, rtol = 1e-10)
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

callbacks = CallbackSet(
    approximate_time_to_finish_cb,
    #increase_dtmax_cb,
)

@time sol = solve(
    implicit_prob,
    #Tsit5(),
    #AutoTsit5(FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true)),
    FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff()),
    #FBDF(linsolve = KrylovJL_GMRES(), autodiff = AutoFiniteDiff()),
    #FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true),
    #FBDF(linsolve = KrylovJL_GMRES(), nlsolve = NLNewton(relax = 0.7), precs = iluzero, concrete_jac = true),
    callback = callbacks,
    #isoutofdomain = state_is_invalid_closure,
)

sol.destats
#=
sol.alg

#f_closure_steady = (du, u, p) -> f_closure_implicit(du, u, p, 0.0)

#nl_func = NonlinearFunction(f_closure_steady, jac_prototype = float.(jac_sparsity))

#prob = NonlinearProblem(nl_func, u0_vec, p_guess)

#@time sol_steady = solve(prob, NonlinearSolve.NewtonRaphson(concrete_jac = true))

#we do this instead of changing the saveat because saveat is pretty linear and doesn't automatically capture fast and slow transient behaviour well
reduced_sol = (
    u = [[sol.u[i] for i in 1:10:length(sol.u)]..., sol.u[end]],
    t = [sol.t[1:10:length(sol.t)]..., sol.t[end]]
)

du_named, u_named = regenerate_fvm_state(reduced_sol, system, solve_system!, geo, p_guess, track_progress = true);

add_xyz_vec_to_u_named!(u_named, :velocity, :vel_u, :vel_v, :vel_w) 
#we previously didn't need this when vel used to be a single vector of u, v, and w components, since they're separate named fields, we must now do this explicitly, oh well...

root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\julia_cfd_output_files"

sol_to_vtk(reduced_sol, du_named, u_named, grid, geo, @__FILE__, root_dir, include_zeros_fields = true, track_progress = true)
