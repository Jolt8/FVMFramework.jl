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
# Boundary flux helpers (`physical_flux` and `primitive_from_conservative`).
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC_low_mach_correction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/SLAU2.jl"))
Revise.includet(joinpath(@__DIR__, "weighted_least_squares/weighted_least_squares.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/fluid_viscous_and_diffusive_fluxes.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/wall_viscous_and_diffusive_fluxes.jl"))

function fluid_fluid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    (
        distance,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

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

    if compressible_solver_options.use_sst_turbulence == true
        add_sst_diffusion_flux!(du, u, idx_a, idx_b, face_area_a, face_normal_a, distance)
    end
end

function airfoil_wall_flux!(du, u, p, t, system, geo, cell_id, face_id)
    face_area, face_normal, _, _ = boundary_geometry(geo, cell_id, face_id)
    gamma = u.cp[cell_id] / u.cv[cell_id]
    _, _, _, pressure, _ = primitive_from_conservative(
        u.density[cell_id],
        u.momentum_density_u[cell_id],
        u.momentum_density_v[cell_id],
        u.momentum_density_w[cell_id],
        u.volumetric_energy[cell_id],
        gamma,
    )

    du.momentum_density_u_flow[cell_id] -= face_area * pressure * face_normal[1]
    du.momentum_density_v_flow[cell_id] -= face_area * pressure * face_normal[2]
    du.momentum_density_w_flow[cell_id] -= face_area * pressure * face_normal[3]
    
    non_moving_wall_viscous_and_diffusive_flux!(du, u, p, t, system, geo, cell_id, face_id, 0, 0)
    
    if compressible_solver_options.use_sst_turbulence
        add_sst_wall_flux!(du, u, geo, cell_id, face_id)
    end
end

# Mapping connection functions
#IMPORTANT: when we switch to evaluating only one flux update per face rather than two, this will determine which idx_a will be and idx_b will be
#This makes sense because if we have two fluids that have the same type (flux functions should be the same) they will both be updated randomly
function connection_map_function(phys_a, phys_b)
    phys_a isa Fluid && phys_b isa Fluid && return fluid_fluid_flux!
    error("No flux is defined between $(typeof(phys_a)) and $(typeof(phys_b))")
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

    properties_comp_vec = ComponentVector(
        cp = u.cp,
        cv = u.cv,
        R_gas = u.R_gas,
        mw = u.mw,
        mu = u.mu, 
        molecular_viscosity = u.mu,
        prandtl_number = u.prandtl_number,
    )

    if compressible_solver_options.use_sst_turbulence == true
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
        if compressible_solver_options.use_sst_turbulence == true
            update_sst_primitives!(du, u, p, t, system, geo, cell_id)
        end
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, system, geo, cell_id)
        cap_navier_stokes_flow!(du, u, p, t, system, geo, cell_id)
        if compressible_solver_options.use_sst_turbulence == true
            cap_sst_transport!(du, u, geo, cell_id)
        end
    end,
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

        if compressible_solver_options.use_sst_turbulence == true
            add_prescribed_boundary_sst_advection_flux!(
                du, idx_a, face_area,
                F_density,
                supersonic_inlet_initial_conditions_stripped.turbulent_kinetic_energy_density /
                    supersonic_inlet_initial_conditions_stripped.density,
                supersonic_inlet_initial_conditions_stripped.specific_dissipation_rate_density /
                    supersonic_inlet_initial_conditions_stripped.density,
            )
        end
    end
)

add_patch!(
    config, "outlet_patch";
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

        if compressible_solver_options.use_sst_turbulence == true
            add_boundary_sst_advection_flux!(du, u, idx_a, face_area, F_density)
        end
    end
)

"""
    characteristic_farfield_flux!(...)

Open far-field boundary for subsonic or supersonic normal flow.  At subsonic
faces the outgoing acoustic invariant is taken from the interior and the
incoming invariant from the prescribed freestream.  Entropy and tangential
velocity are taken from the upwind side.  Supersonic outflow extrapolates the
interior state, while supersonic inflow prescribes the complete freestream
state.
"""
function characteristic_farfield_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b,
)
    face_area, face_normal, _, _ = boundary_geometry(geo, idx_a, face_a)

    gamma = u.cp[idx_a] / u.cv[idx_a]
    density_i = u.density[idx_a]
    momentum_u_i = u.momentum_density_u[idx_a]
    momentum_v_i = u.momentum_density_v[idx_a]
    momentum_w_i = u.momentum_density_w[idx_a]
    energy_i = u.volumetric_energy[idx_a]
    vel_u_i, vel_v_i, vel_w_i, pressure_i, sound_i =
        primitive_from_conservative(
            density_i,
            momentum_u_i,
            momentum_v_i,
            momentum_w_i,
            energy_i,
            gamma,
        )
    normal_vel_i =
        vel_u_i * face_normal[1] +
        vel_v_i * face_normal[2] +
        vel_w_i * face_normal[3]

    density_inf = supersonic_inlet_initial_conditions_stripped.density
    momentum_u_inf = supersonic_inlet_initial_conditions_stripped.momentum_density_u
    momentum_v_inf = supersonic_inlet_initial_conditions_stripped.momentum_density_v
    momentum_w_inf = supersonic_inlet_initial_conditions_stripped.momentum_density_w
    energy_inf = supersonic_inlet_initial_conditions_stripped.volumetric_energy
    vel_u_inf, vel_v_inf, vel_w_inf, pressure_inf, sound_inf =
        primitive_from_conservative(
            density_inf,
            momentum_u_inf,
            momentum_v_inf,
            momentum_w_inf,
            energy_inf,
            gamma,
        )
    normal_vel_inf =
        vel_u_inf * face_normal[1] +
        vel_v_inf * face_normal[2] +
        vel_w_inf * face_normal[3]

    if normal_vel_i >= sound_i
        # Every characteristic leaves the domain.
        density_b = density_i
        momentum_u_b = momentum_u_i
        momentum_v_b = momentum_v_i
        momentum_w_b = momentum_w_i
        energy_b = energy_i
    elseif normal_vel_i <= -sound_i
        # Every characteristic enters the domain.
        density_b = density_inf
        momentum_u_b = momentum_u_inf
        momentum_v_b = momentum_v_inf
        momentum_w_b = momentum_w_inf
        energy_b = energy_inf
    else
        # R+ leaves the domain and R- enters it.
        invariant_plus = normal_vel_i + 2.0 * sound_i / (gamma - 1.0)
        invariant_minus = normal_vel_inf - 2.0 * sound_inf / (gamma - 1.0)
        normal_vel_b = 0.5 * (invariant_plus + invariant_minus)
        sound_b = max(
            0.25 * (gamma - 1.0) * (invariant_plus - invariant_minus),
            1e-8,
        )

        if normal_vel_b >= 0.0
            # Subsonic outflow: entropy and tangential velocity leave the domain.
            entropy_b = pressure_i / density_i^gamma
            tangential_u = vel_u_i - normal_vel_i * face_normal[1]
            tangential_v = vel_v_i - normal_vel_i * face_normal[2]
            tangential_w = vel_w_i - normal_vel_i * face_normal[3]
        else
            # Subsonic inflow: entropy and tangential velocity enter from freestream.
            entropy_b = pressure_inf / density_inf^gamma
            tangential_u = vel_u_inf - normal_vel_inf * face_normal[1]
            tangential_v = vel_v_inf - normal_vel_inf * face_normal[2]
            tangential_w = vel_w_inf - normal_vel_inf * face_normal[3]
        end

        density_b = (sound_b^2 / (gamma * entropy_b))^(1.0 / (gamma - 1.0))
        pressure_b = entropy_b * density_b^gamma
        vel_u_b = tangential_u + normal_vel_b * face_normal[1]
        vel_v_b = tangential_v + normal_vel_b * face_normal[2]
        vel_w_b = tangential_w + normal_vel_b * face_normal[3]
        momentum_u_b = density_b * vel_u_b
        momentum_v_b = density_b * vel_v_b
        momentum_w_b = density_b * vel_w_b
        energy_b = pressure_b / (gamma - 1.0) + 0.5 * density_b * (
            vel_u_b^2 + vel_v_b^2 + vel_w_b^2
        )
    end

    (
        F_density,
        F_momentum_density_u,
        F_momentum_density_v,
        F_momentum_density_w,
        F_volumetric_energy,
    ) = slau2_flux(
        density_i,
        momentum_u_i,
        momentum_v_i,
        momentum_w_i,
        energy_i,
        gamma,
        density_b,
        momentum_u_b,
        momentum_v_b,
        momentum_w_b,
        energy_b,
        gamma,
        face_normal,
    )

    du.density_flow[idx_a] -= face_area * F_density
    du.momentum_density_u_flow[idx_a] -= face_area * F_momentum_density_u
    du.momentum_density_v_flow[idx_a] -= face_area * F_momentum_density_v
    du.momentum_density_w_flow[idx_a] -= face_area * F_momentum_density_w
    du.volumetric_energy_flow[idx_a] -= face_area * F_volumetric_energy

    if compressible_solver_options.use_sst_turbulence == true
        if F_density >= 0.0
            add_boundary_sst_advection_flux!(du, u, idx_a, face_area, F_density)
        else
            add_prescribed_boundary_sst_advection_flux!(
                du, idx_a, face_area, F_density,
                supersonic_inlet_initial_conditions_stripped.turbulent_kinetic_energy_density /
                    supersonic_inlet_initial_conditions_stripped.density,
                supersonic_inlet_initial_conditions_stripped.specific_dissipation_rate_density /
                    supersonic_inlet_initial_conditions_stripped.density,
            )
        end
    end
end

for name in ["top_patch", "bottom_patch"]
    add_patch!(
        config, name;
        properties = ComponentVector(),
        patch_function = characteristic_farfield_flux!,
    )
end

add_patch!(
    config, "airfoil_surface_patches";
    properties = ComponentVector(),
    patch_function = function airfoil_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        airfoil_wall_flux!(du, u, p, t, system, geo, idx_a, face_a)
    end,
)

#=
for name in ["top_patch", "bottom_patch"]
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
            
            if compressible_solver_options.use_sst_turbulence == true
                add_sst_wall_flux!(du, u, geo, idx_a, face_a)
            end
        end
    ) 
end
=#

# The mesh is one cell thick. These front/back faces are symmetry planes for
# the nominally two-dimensional airfoil calculation.
add_patch!(
    config, "side_patches";
    properties = ComponentVector(),
    patch_function = function symmetry_flux!(
        du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b,
    )
        face_area, face_normal, _, _ = boundary_geometry(geo, idx_a, face_a)
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
    end
)

#I think we're going to add an additional_data field of the system so that arbitrary data can be passed into any function without allocations
additional_data = (
    #You could put a neural network in here, custom structs, interpolations, just any useful data in general. 
    #This does require that every function have system in its arguments though
)

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, additional_data, check_units = false);
WLS_STENCIL = build_weighted_least_squares_stencil(geo)

if compressible_solver_options.use_sst_turbulence == true
    airfoil_face_centers = typeof(geo.cell_centroids[1])[]
    for (cell_id, face_id) in getfacetset(grid, "airfoil_surface_patches")
        face_node_ids = Ferrite.faces(grid.cells[cell_id])[face_id]
        face_center = sum(grid.nodes[node_id].x for node_id in face_node_ids) /
            length(face_node_ids)
        push!(airfoil_face_centers, face_center)
    end

    isempty(airfoil_face_centers) && error("No fluid-airfoil interface was found in the mesh")
    SST_WALL_DISTANCES = fill(1.0, n_cells)
    for cell_id in eachindex(geo.cell_centroids)
        SST_WALL_DISTANCES[cell_id] = minimum(
            norm(geo.cell_centroids[cell_id] - face_center)
            for face_center in airfoil_face_centers
        )
    end
end

# System solver function
function solve_system!(du, u, p, t, system, geo)
    update_region_groups!(du, u, p, t, system, geo)

    update_weighted_least_squares_gradients!(u, WLS_STENCIL)
    
    if compressible_solver_options.use_sst_turbulence == true
        update_sst_closure!(u, SST_WALL_DISTANCES)
    end
    
    #update_MUSCL_gradients!(u, WLS_STENCIL) #REMEMBER TO REACTIVATE THIS!!!
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
