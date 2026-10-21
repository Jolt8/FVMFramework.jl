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

typeof(grid)

n_cells = length(grid.cells)

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
        oxygen = zeros(n_cells)u"kg/m^3",
        carbon_dioxide = zeros(n_cells)u"kg/m^3",
        carbon_monoxide = zeros(n_cells)u"kg/m^3",
        hydrogen = zeros(n_cells)u"kg/m^3",
        water = zeros(n_cells)u"kg/m^3",
    ),
    fuel_grain_thickness = zeros(n_cells)u"m"
)

if compressible_solver_options.use_sst_turbulence == true
    u_proto = merge_properties(u_proto, ComponentVector(
        turbulent_kinetic_energy_density = zeros(n_cells)u"J/m^3",
        specific_dissipation_rate_density = zeros(n_cells)u"kg/(m^3*s)",
    ))
end

config = create_fvm_config(grid, u_proto);

Revise.includet(joinpath(@__DIR__, "geometry_editing/geometry_editing.jl"))

config.geo.cell_neighbor_normals[6, 1]

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
            nitrous_oxide = zeros(n_cells)u"kg/kg",
            ethylene = zeros(n_cells)u"kg/kg",
            nitrogen = zeros(n_cells)u"kg/kg",
            oxygen = zeros(n_cells)u"kg/kg",
            carbon_dioxide = zeros(n_cells)u"kg/kg",
            carbon_monoxide = zeros(n_cells)u"kg/kg",
            hydrogen = zeros(n_cells)u"kg/kg",
            water = zeros(n_cells)u"kg/kg",
        ),
        species_density_flow = (
            nitrous_oxide = zeros(n_cells)u"kg/s",
            ethylene = zeros(n_cells)u"kg/s",
            nitrogen = zeros(n_cells)u"kg/s",
            oxygen = zeros(n_cells)u"kg/s",
            carbon_dioxide = zeros(n_cells)u"kg/s",
            carbon_monoxide = zeros(n_cells)u"kg/s",
            hydrogen = zeros(n_cells)u"kg/s",
            water = zeros(n_cells)u"kg/s",
        ),
    ),
    special_caches = ComponentVector(
        grad_density = zeros(n_cells, 3)u"kg/m^4",
        grad_pressure = zeros(n_cells, 3)u"Pa/m",
        grad_vel_u = zeros(n_cells, 3)u"m/(s*m)",
        grad_vel_v = zeros(n_cells, 3)u"m/(s*m)",
        grad_vel_w = zeros(n_cells, 3)u"m/(s*m)",
        grad_temperature = zeros(n_cells, 3)u"K/m",
        
        vel_u_face = zeros(n_cells, 6)u"m/s",
        vel_v_face = zeros(n_cells, 6)u"m/s",
        vel_w_face = zeros(n_cells, 6)u"m/s",
        temperature_face = zeros(n_cells, 6)u"J/m^3",

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
# Mesh Physical Groups Summary

## Special Volumes
- fluid_domain

## Non-Special Volumes
- aluminum_housing_volume
- fuel_grain_volume
- injector_bulkhead_volume
- injector_disk_volume
- inner_graphite_nozzle_volume
- outer_aluminum_nozzle_casing_volume
- phenolic_liner_volume

## Special Interfaces
- fuel_wall_patch

## Non-Special Interfaces
- aluminum_housing_to_injector_bulkhead_interface
- aluminum_housing_to_inner_graphite_nozzle_interface
- aluminum_housing_to_outer_aluminum_nozzle_casing_interface
- aluminum_housing_to_phenolic_liner_interface
- fluid_domain_to_injector_bulkhead_interface
- fluid_domain_to_injector_disk_interface
- fluid_domain_to_inner_graphite_nozzle_interface
- fluid_domain_to_outer_aluminum_nozzle_casing_interface
- fluid_domain_to_phenolic_liner_interface
- fuel_grain_to_phenolic_liner_interface
- injector_bulkhead_to_injector_disk_interface
- injector_bulkhead_to_phenolic_liner_interface
- inner_graphite_nozzle_to_outer_aluminum_nozzle_casing_interface
- inner_graphite_nozzle_to_phenolic_liner_interface

## Special Wall Patches
- far_field_patch
- out_of_plane_patch
- outlet_patch
- symmetry_patch

## Non-Special Wall Patches
- aluminum_housing_wall_patch
- fluid_domain_wall_patch
- injector_bulkhead_wall_patch
- injector_disk_wall_patch
- outer_aluminum_nozzle_casing_wall_patch
"""

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end
struct FuelGrain <: AbstractPhysics end

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

function heat_diffusion!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b
)
    (
        distance,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / distance

    heat_exchange = -k_effective * grad_T * face_area_a

    du.volumetric_energy[idx_a] -= heat_exchange / volume_a
    du.volumetric_energy[idx_b] += heat_exchange / volume_b
end

function fluid_solid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    heat_diffusion!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) #Don't think this is the most accurate because it neglects convection and radiation
        #TODO: implement convection here
        #Radiation can be saved for later if it's discovered that it's actually needed
    no_slip_adiabatic_wall_flux!(du, u, p, t, system, geo, idx_a, face_a)
end

function no_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    #nothing here since we handle these BCs at the patches
end

function solid_solid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b
)
    heat_diffusion!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Solid && return solid_solid_flux!
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Solid && return fluid_solid_flux!
    typeof(phys_a) <: FuelGrain && typeof(phys_b) <: FuelGrain && return solid_solid_flux!
    typeof(phys_a) <: FuelGrain && typeof(phys_b) <: Solid && return solid_solid_flux!
    typeof(phys_a) <: FuelGrain && typeof(phys_b) <: Fluid && return no_flux! #this one is special
end

Revise.includet(joinpath(@__DIR__, "construct_state_from_intuitive_inputs/construct_state_from_intuitive_inputs.jl"))

fluid_initial_conditions, fluid_properties = construct_initial_conditions_from_intuitive_inputs(
    #TODO: fill out properties of air
    ComponentVector(
        vel_u = 0.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = -1e-8u"m/s",
        density = 1.225u"kg/m^3",
        temperature = 21.13u"°C",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol",
        mu = 1e-5u"Pa*s",
        turbulent_kinetic_energy = 0.1u"m^2/s^2",
        specific_dissipation_rate = 100.0u"1/s",
        #k = 0.026u"W/(m*K)",
        prandtl_number = 0.705,

        mass_fractions = (
            nitrous_oxide = 0.0,
            ethylene = 0.0,
            nitrogen = 0.781,
            oxygen = 0.211,
            carbon_dioxide = 0.0003,
            carbon_monoxide = 0.0,
            hydrogen = 0.0,
            water = 0.0
        )
    )
)

add_region!(
    config, "fluid_domain";
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


fuel_grain_initial_conditions, fuel_grain_properties = construct_solid_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        density = 950.0u"kg/m^3",
        temperature = 21.13u"°C",
        cp = 2000.0u"J/(kg*K)",
        k = 0.45u"W/(m*K)"
    )
)

add_region!(
    config, "fuel_grain_volume";
    type = FuelGrain(), #we're going to have a special flux between this and the fluid defined in an add_patch!() call
    initial_conditions = fuel_grain_initial_conditions,
    properties = fuel_grain_properties,
    property_update_function = 
    function update_fuel_grain_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function fuel_grain_region!(du, u, p, t, system, geo, cell_id)
        fuel_grain_regression!(du, u, p, t, system, geo, cell_id)
    end
)

aluminum_initial_conditions, aluminum_properties = construct_solid_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        density = 2700.0u"kg/m^3",
        temperature = 21.13u"°C",
        cp = 900.0u"J/(kg*K)",
        k = 237.0u"W/(m*K)"
    )
)

aluminum_regions = [
    "aluminum_housing_volume",
    "injector_bulkhead_volume",
    "injector_disk_volume",
    "outer_aluminum_nozzle_casing_volume"
]

for name in aluminum_regions
    add_region!(
        config, name;
        type = Solid(),
        initial_conditions = aluminum_initial_conditions,
        properties = aluminum_properties,
        property_update_function = 
        function update_aluminum_region!(du, u, p, t, system, geo, cell_id)
            update_solid_properties!(du, u, p, t, system, geo, cell_id)
        end,
        region_function = 
        function aluminum_region!(du, u, p, t, system, geo, cell_id)
        end
    )
end

graphite_initial_conditions, graphite_properties = construct_solid_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        density = 1760.0u"kg/m^3",
        temperature = 21.13u"°C",
        cp = 710.0u"J/(kg*K)",
        k = 130.0u"W/(m*K)"
    )
)

#NOTE: I'm making an assumption here
#since we're using extruded graphite, the thermal conductivity in the axial and radial directions are different
#I'm feeling lazy so I'm not going to model it and just go off of the thermal conductivity that the manufactuerers provided directly
add_region!(
    config, "inner_graphite_nozzle_volume";
    type = Solid(), 
    initial_conditions = graphite_initial_conditions,
    properties = graphite_properties,
    property_update_function = 
    function update_graphite_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function graphite_region!(du, u, p, t, system, geo, cell_id)
    end
)

phenolic_liner_initial_conditions, phenolic_liner_properties = construct_solid_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        density = 1250.0u"kg/m^3",
        temperature = 21.13u"°C",
        cp = 1400.0u"J/(kg*K)",
        k = 0.25u"W/(m*K)"
    )
)

add_region!(
    config, "phenolic_liner_volume";
    type = Solid(),
    initial_conditions = phenolic_liner_initial_conditions,
    properties = phenolic_liner_properties,
    property_update_function = 
    function update_phenolic_liner_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function phenolic_liner_region!(du, u, p, t, system, geo, cell_id)
    end
)

#Special Interfaces
add_patch!(
    config, "fuel_wall_patch";
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        heat_diffusion!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) #Don't think this is the most accurate because it neglects convection and radiation
        #TODO: implement convection here
        #Radiation can be saved for later if it's discovered that it's actually needed
        no_slip_adiabatic_wall_flux!(du, u, p, t, system, geo, idx_a, face_a)

        #TODO: put in the regression logic here
    end,
)

#Special Walls
add_patch!(
    config, "far_field_patch";
    patch_function = function farfield_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b, 
    )
        characteristic_farfield_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b, farfield_conditions = fluid_initial_conditions)
        #we reuse the initial conditions for the fluid (air) as our farfield conditions
    end,
)

#TODO: basically reimplement the nitrous oxide injection logic from the earlier design script 
#bro, this is going to be so difficult, oh well...
inlet_conditions = construct_inlet_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 0.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = -50.0u"m/s", #we're going to have this vary with time
        density = 1.225u"kg/m^3", #this is pretty hard to determine
        temperature = 21.13u"°C", #also hard to determine
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol",
        mu = 1e-5u"Pa*s",
        turbulent_kinetic_energy = 0.1u"m^2/s^2",
        specific_dissipation_rate = 100.0u"1/s",
        #k = 0.026u"W/(m*K)",
        prandtl_number = 0.705,

        mass_fractions = (
            nitrous_oxide = 0.0,
            ethylene = 0.0,
            nitrogen = 0.781,
            oxygen = 0.211,
            carbon_dioxide = 0.0003,
            carbon_monoxide = 0.0,
            hydrogen = 0.0,
            water = 0.0
        )
    )
)

add_patch!(
    config, "injector_inlet_patch";
    patch_function = function injector_inlet_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b, 
    )
        subsonic_inlet_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a, 
            idx_b, face_b;
            inlet_conditions = inlet_conditions
        )
        #TODO: we need a new subsonic inlet boundary condition function
    end,
)

add_patch!(
    config, "out_of_plane_patch";
    patch_function = function out_of_plane_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b, 
    )
        #hmm, not sure what the most appropriate BC is to use here
        #TODO: determine what type of bc to use here
    end,
)

add_patch!(
    config, "outlet_patch";
    patch_function = supersonic_outlet_flux! #I think this is appropriate here
)

add_patch!(
    config, "symmetry_patch";
    patch_function = symmetry_flux!
)

function ambient_air_heat_exchange!(
    du, u, p, t, system, geo, 
    idx_a, face_a, 
    idx_b, face_b; 
    ambient_conditions = fluid_initial_conditions
)
    #this is going to be more complex since we should probably use a heat_transfer_coefficient based method so this is somewhat accurate
end

for name in [
    "aluminum_housing_wall_patch",
    "fluid_domain_wall_patch",
    "injector_bulkhead_wall_patch",
    "injector_disk_wall_patch",
    "outer_aluminum_nozzle_casing_wall_patch"
]
    add_patch!(
        config, name;
        patch_function = function ambient_air_heat_exchange_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a, idx_b, face_b
        )
            ambient_air_heat_exchange!(
                du, u, p, t, system, geo,
                idx_a, face_a, idx_b, face_b;
                ambient_conditions = fluid_initial_conditions
            )
        end
    )
end