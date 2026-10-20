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
using Clapeyron

mesh_path = joinpath(@__DIR__, "meshes/hybrid_rocket_cfd_3d.msh")
grid = togrid(mesh_path)

num_hex = count(c -> isa(c, Hexahedron), grid.cells)
num_wedges = count(c -> isa(c, Wedge), grid.cells)

n_cells = length(grid.cells)
n_faces = nfacets(grid.cells[1])

compressible_solver_options = (
    use_sst_turbulence = true,
)

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3",
    momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_v = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_w = zeros(n_cells)u"kg/(m^2*s)",
    volumetric_energy = zeros(n_cells)u"J/m^3",
    species_densities = ComponentVector(
        nitrous_oxide = zeros(n_cells)u"kg/m^3",
        ethylene = zeros(n_cells)u"kg/m^3",
        nitrogen = zeros(n_cells)u"kg/m^3",
        water = zeros(n_cells)u"kg/m^3",
        carbon_dioxide = zeros(n_cells)u"kg/m^3",
        carbon_monoxide = zeros(n_cells)u"kg/m^3",
        hydrogen = zeros(n_cells)u"kg/m^3",
    ),
    fuel_grain_thickness = zeros(n_cells)u"m"
)

if compressible_solver_options.use_sst_turbulence == true
    u_proto = merge_properties(u_proto, ComponentVector(
        turbulent_kinetic_energy_density = zeros(n_cells)u"J/m^3",
        specific_dissipation_rate_density = zeros(n_cells)u"kg/(m^3*s)",
    ))
end

config = create_fvm_config(grid, u_proto)

Revise.includet(joinpath(@__DIR__, "geometry_editing/geometry_editing.jl"))

config.geo = correct_hybrid_geo(config.geo, grid)

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

        mass_fractions = (
            species_a = u"kg/kg",
            species_b = u"kg/kg",
        ),
        species_density_flow = (
            species_a = u"kg/s",
            species_b = u"kg/s",
        ),
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

path_to_navier_stokes_files = joinpath("..", "LLM_developed_CFD/navier_stokes_testing_grounds")

Revise.includet(joinpath(path_to_navier_stokes_files, "face_reconstructors/first_order_face_reconstruction.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "turbulence/sst_k_omega.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "boundary_flux_helpers/boundary_flux_helpers.jl")) # Boundary flux helpers (`physical_flux` and `primitive_from_conservative`).
Revise.includet(joinpath(path_to_navier_stokes_files, "boundary_conditions/boundary_condition_exports.jl")) # Boundary Conditions

Revise.includet(joinpath(path_to_navier_stokes_files, "riemann_solvers/HLLC_low_mach_correction.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "riemann_solvers/HLLC.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "riemann_solvers/SLAU2.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "weighted_least_squares/weighted_least_squares.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "viscous_and_diffusive_terms/fluid_viscous_and_diffusive_fluxes.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "viscous_and_diffusive_terms/wall_viscous_and_diffusive_fluxes.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "navier_stokes_fluid_property_update_functions/non_standard_fluid_property_update_functions.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "navier_stokes_fluid_property_update_functions/fluid_property_update_functions.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "sum_and_cap_functions/cap_functions.jl"))

"""
Current Regions and Patches:
    - aluminum_housing_to_inner_graphite_nozzle_interface
    - aluminum_housing_to_outer_aluminum_nozzle_casing_interface
    - aluminum_housing_to_phenolic_liner_interface
    - fluid_domain_to_inner_graphite_nozzle_interface
    - fluid_domain_to_outer_aluminum_nozzle_casing_interface
    - fluid_domain_to_phenolic_liner_interface
    - fuel_grain_to_phenolic_liner_interface
    - inner_graphite_nozzle_to_outer_aluminum_nozzle_casing_interface
    - inner_graphite_nozzle_to_phenolic_liner_interface

    - aluminum_housing_wall_patch
    - far_field_patch
    - fluid_domain_wall_patch
    - fuel_wall_patch
    - injector_inlet_patch
    - out_of_plane_patch
    - outer_aluminum_nozzle_casing_wall_patch
    - outlet_patch
    - phenolic_liner_wall_patch
    - symmetry_patch

Current Regions:
    - aluminum_housing_volume
    - fluid_domain
    - fuel_grain_volume
    - inner_graphite_nozzle_volume
    - outer_aluminum_nozzle_casing_volume
    - phenolic_liner_volume
"""

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

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
        #MUSCL_face_reconstruction!,
    )

    fluid_viscous_and_diffusive_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
end

function fluid_solid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
end

function solid_solid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b
)
    heat_diffusion!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Solid && return fluid_fluid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
end

function construct_fluid_initial_conditions_from_intuitive_inputs(u)
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
    
    properties_comp_vec = ComponentVector(
        cp = u.cp,
        cv = u.cv,
        R_gas = u.R_gas,
        mw = u.mw,
        mu = u.mu, 
        molecular_viscosity = u.mu,
        prandtl_number = u.prandtl_number,
        diffusion_coefficients = u.diffusion_coefficients,
    )

    if system.compressible_solver_options.use_sst_turbulence == true
        return ComponentVector(
            density = u.density,
            momentum_density_u = momentum_density_u,
            momentum_density_v = momentum_density_v,
            momentum_density_w = momentum_density_w,
            volumetric_energy = volumetric_energy,
            turbulent_kinetic_energy_density = turbulent_kinetic_energy_density,
            specific_dissipation_rate_density = specific_dissipation_rate_density,
            species_densities = species_densities,
        ), properties_comp_vec
    else
        return ComponentVector(
            density = u.density,
            momentum_density_u = momentum_density_u,
            momentum_density_v = momentum_density_v,
            momentum_density_w = momentum_density_w,
            volumetric_energy = volumetric_energy,
            species_densities = species_densities,
        ), properties_comp_vec
    end
end

function construct_solid_initial_conditions_from_intuitive_inputs(u)
    """
    Requires:
    - density
    - temperature
    - cv
    - cp
    """
    state_comp_vec = ComponentVector(
        density = u.density,
        momentum_density_u = 0.0,
        momentum_density_v = 0.0,
        momentum_density_w = 0.0,
        volumetric_energy = u.density * u.cv * u.temperature,
    )

    properties_comp_vec = ComponentVector(
        cp = u.cp,
        cv = u.cv
    )

    return state_comp_vec, properties_comp_vec
end


