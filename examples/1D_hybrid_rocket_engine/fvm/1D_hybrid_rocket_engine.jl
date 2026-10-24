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
using PreallocationTools
using DataInterpolations

mesh_path = joinpath(@__DIR__, "meshes/hybrid_rocket_cfd_3d.msh")
grid = togrid(mesh_path)

num_hex = count(c -> isa(c, Hexahedron), grid.cells)
num_wedges = count(c -> isa(c, Wedge), grid.cells)

typeof(grid)

n_cells = length(grid.cells)

compressible_solver_options = (
    use_sst_turbulence = true,
    use_muscl_reconstruction = false,
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


Revise.includet(joinpath(@__DIR__, "..", "design_scripts/overloads/clapeyron_tracer_overloads.jl"))
Revise.includet(joinpath(@__DIR__, "injector_delivery_train_scripts/injector_delivery_train_scripts.jl"))
nitrous_oxide_model = PR(["nitrous oxide"])
injector_u0, injector_properties = get_injector_train_u0_and_properties(nitrous_oxide_model)

u_proto = merge_properties(u_proto, injector_u0)

config = create_fvm_config(grid, u_proto);

Revise.includet(joinpath(@__DIR__, "geometry_editing/geometry_editing.jl"))

correct_hybrid_geo!(config.geo, grid)

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
            nitrous_oxide = u"kg/kg",
            ethylene = u"kg/kg",
            nitrogen = u"kg/kg",
            oxygen = u"kg/kg",
            carbon_dioxide = u"kg/kg",
            carbon_monoxide = u"kg/kg",
            hydrogen = u"kg/kg",
            water = u"kg/kg",
        ),
        species_density_flow = (
            nitrous_oxide = u"kg/s",
            ethylene = u"kg/s",
            nitrogen = u"kg/s",
            oxygen = u"kg/s",
            carbon_dioxide = u"kg/s",
            carbon_monoxide = u"kg/s",
            hydrogen = u"kg/s",
            water = u"kg/s",
        ),
        diffusion_coefficients = (
            nitrous_oxide = u"m^2/s",
            ethylene = u"m^2/s",
            nitrogen = u"m^2/s",
            oxygen = u"m^2/s",
            carbon_dioxide = u"m^2/s",
            carbon_monoxide = u"m^2/s",
            hydrogen = u"m^2/s",
            water = u"m^2/s",
        )
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
        temperature_face = zeros(n_cells, 6)u"K",

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

path_to_navier_stokes_files = joinpath("..", "..", "LLM_developed_CFD/navier_stokes_testing_grounds")

Revise.includet(joinpath(path_to_navier_stokes_files, "face_reconstructors/first_order_face_reconstruction.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "turbulence/sst_k_omega.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "turbulence/sst_k_omega_helper_functions.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "boundary_flux_helpers/boundary_flux_helpers.jl")) # Boundary flux helpers (`physical_flux` and `primitive_from_conservative`).
Revise.includet(joinpath(path_to_navier_stokes_files, "boundary_conditions/boundary_condition_exports.jl")) # Boundary Conditions
Revise.includet(joinpath(@__DIR__, "boundary_condition_functions/subsonic_inlet_flux.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "riemann_solvers/HLLC_low_mach_correction.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "riemann_solvers/HLLC.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "riemann_solvers/SLAU2.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "weighted_least_squares/weighted_least_squares.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "viscous_and_diffusive_terms/fluid_viscous_and_diffusive_fluxes.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "viscous_and_diffusive_terms/wall_viscous_and_diffusive_fluxes.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "navier_stokes_fluid_property_update_functions/non_standard_fluid_property_update_functions.jl"))
Revise.includet(joinpath(path_to_navier_stokes_files, "navier_stokes_fluid_property_update_functions/fluid_property_update_functions.jl"))

Revise.includet(joinpath(path_to_navier_stokes_files, "sum_and_cap_functions/cap_functions.jl"))

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end
struct FuelGrain <: AbstractPhysics end

function fluid_fluid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    (
        dist,
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
        dist,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

    grad_T = (u.temperature[idx_b] - u.temperature[idx_a]) / dist

    heat_exchange = -k_effective * grad_T * face_area_a

    du.volumetric_energy[idx_a] -= heat_exchange / volume_a
    du.volumetric_energy[idx_b] += heat_exchange / volume_b

    #=
    if idx_a == 443
        @show dist
        @show face_area_a
        @show volume_a
        @show volume_b 
        @show u.k[idx_a]
        @show u.k[idx_b]
        @show k_effective
        @show u.temperature[idx_b]
        @show u.temperature[idx_a]
        @show grad_T
        @show heat_exchange
        @show du.volumetric_energy[idx_a]
        @show du.volumetric_energy[idx_b]
        @show u.density[idx_a]
        @show u.cp[idx_a]
    end
    =#
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
    typeof(phys_a) <: Fluid && typeof(phys_b) <: FuelGrain && return no_flux!
    #fluid_solid_flux! #this one is special
end

Revise.includet(joinpath(@__DIR__, "construct_state_from_intuitive_inputs/construct_state_from_intuitive_inputs.jl"))
Revise.includet(joinpath(@__DIR__, "update_solid_properties/update_solid_properties.jl"))

fluid_initial_conditions, fluid_properties = construct_fluid_initial_conditions_from_intuitive_inputs(compressible_solver_options,
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
        ),
        diffusion_coefficients_ambient = (
            nitrous_oxide = 1.4e-5u"m^2/s",
            ethylene = 1.4e-5u"m^2/s",
            nitrogen = 1.8e-5u"m^2/s",
            oxygen = 1.9e-5u"m^2/s",
            carbon_dioxide = 1.6e-5u"m^2/s",
            carbon_monoxide = 2.1e-5u"m^2/s",
            hydrogen = 7.6e-5u"m^2/s",
            water = 2.4e-5u"m^2/s",
        ),
        diffusion_coefficients_ambient_temperature = 21.13u"°C",
        diffusion_coefficients_ambient_pressure = 1.0u"atm"
    )
)

function update_fluid_diffusion_coefficients!(du, u, p, t, system, geo, cell_id)
    for_fields!(u.diffusion_coefficients_ambient, u.diffusion_coefficients) do species, u_diffusion_coefficients_ambient, u_diffusion_coefficients
        u_diffusion_coefficients[species[cell_id]] = u_diffusion_coefficients_ambient[species[cell_id]] * (
            u.temperature[cell_id] / u.diffusion_coefficients_ambient_temperature[cell_id]) ^ 1.8 * 
            (u.diffusion_coefficients_ambient_pressure[cell_id] / u.pressure[cell_id]
        )
    end
end

function add_axisymmetric_inviscid_source!(
    du, u, p, t, system, geo, cell_id
)
    r = geo.cell_centroids[cell_id][1]

    # Radial momentum source for axisymmetric Euler flow without swirl.
    du.momentum_density_u[cell_id] += u.pressure[cell_id] / r
    #Note, this must be (p + ρ u_θ² - τ_θθ) / r in the future

    return nothing
end

add_region!(
    config, "fluid_domain";
    type = Fluid(),
    initial_conditions = fluid_initial_conditions,
    properties = fluid_properties,
    property_update_function = 
    function update_fluid_properties!(du, u, p, t, system, geo, cell_id)
        update_k_from_prandtl!(du, u, p, t, system, geo, cell_id)
        overall_navier_stokes_property_update!(du, u, p, t, system, geo, cell_id)
        update_fluid_diffusion_coefficients!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, system, geo, cell_id)
        cap_navier_stokes_flow!(du, u, p, t, system, geo, cell_id)
        add_axisymmetric_inviscid_source!(du, u, p, t, system, geo, cell_id)
    end
)

fuel_grain_initial_conditions, fuel_grain_properties = construct_solid_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        density = 950.0u"kg/m^3",
        temperature = 21.13u"°C",
        cp = 2000.0u"J/(kg*K)",
        cv = 1500.0u"J/(kg*K)",
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
        #fuel_grain_regression!(du, u, p, t, system, geo, cell_id)
    end
)

aluminum_initial_conditions, aluminum_properties = construct_solid_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        density = 2700.0u"kg/m^3",
        temperature = 21.13u"°C",
        cp = 900.0u"J/(kg*K)",
        cv = 880.0u"J/(kg*K)",
        k = 237.0u"W/(m*K)"
    )
)

aluminum_regions = [
    "aluminum_housing_volume",
    "injector_bulkhead_volume",
    "injector_disk_volume",
    "outer_aluminum_nozzle_casing_volume",
    "flange_fasteners_bottom_volume",
    "flange_fasteners_top_volume",
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
        cv = 700.0u"J/(kg*K)",
        k = 130.0u"W/(m*K)"
    )
)

#NOTE: I'm making an assumption here
#since we're using extruded graphite, the thermal conductivity in the axial and radial directions are different
#I'm feeling lazy so I'm not going to model it and just go off of the thermal conductivity that the manufactuerers provided directly 130.0 [W/(m*K)]
#The radial thermal conductivity is usually around 50-57 [W/(m*K)]
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
        cv = 1350.0u"J/(kg*K)",
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

#TODO: these properties are not particularly accurate, but it doesn't really matter
o_ring_initial_conditions, o_ring_properties = construct_solid_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        density = 1000.0u"kg/m^3",
        temperature = 21.13u"°C",
        cp = 1400.0u"J/(kg*K)",
        cv = 1830.0u"J/(kg*K)",
        k = 0.15u"W/(m*K)"
    )
)

for name in ["housing_o_ring_volume", "injector_disk_o_ring_volume"]
    add_region!(
        config, name;
        type = Solid(),
        initial_conditions = o_ring_initial_conditions,
        properties = o_ring_properties,
        property_update_function = 
        function update_o_ring_region!(du, u, p, t, system, geo, cell_id)
            update_solid_properties!(du, u, p, t, system, geo, cell_id)
        end,
        region_function = 
        function o_ring_region!(du, u, p, t, system, geo, cell_id)
        end
    )
end

Revise.includet(joinpath(@__DIR__, "fuel_grain_regression_physics/get_dominant_fuel_grain_cells_and_subsidiaries.jl"))
Revise.includet(joinpath(@__DIR__, "fuel_grain_regression_physics/geometry_interpolations.jl"))

dominant_fuel_grain_cells,
fuel_grain_subsidiary_cell_ids_to_dominant_cell_id_map,
fuel_grain_dominant_cell_id_to_subsidiaries_map,
fluid_domain_subsidiary_cell_ids_to_dominant_cell_id_map,
fluid_domain_dominant_cell_id_to_subsidiaries_map = 
    get_dominant_fuel_grain_cells_and_subsidiaries(grid, config.geo)

fuel_grain_regression_geometry_interpolations = precompute_fuel_grain_regression_geometry_interpolations(
    grid, config.geo,
    dominant_fuel_grain_cells, 
    fuel_grain_subsidiary_cell_ids_to_dominant_cell_id_map,
    fuel_grain_dominant_cell_id_to_subsidiaries_map,
    fluid_domain_subsidiary_cell_ids_to_dominant_cell_id_map,
    fluid_domain_dominant_cell_id_to_subsidiaries_map
);

function test_udpate_geo_speed(geo, node_coordinates)
    update_fvm_geometry!(geo, node_coordinates)
end

original_node_coordinates = config.geo.node_coordinates

config.geo.node_coordinates

edited_node_coords = [Vector(node_coords) .+= 0.0001 for node_coords in config.geo.node_coordinates]

using BenchmarkTools

@btime update_fvm_geometry!(config.geo, config.geo.node_coordinates);

@time correct_hybrid_geo!(config.geo, grid)

#=
#Special Interfaces
add_patch!(
    config, "fuel_wall_patch";
    owner_region = "fluid_domain",
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

supersonic_inlet_initial_conditions, supersonic_inlet_properties = construct_fluid_initial_conditions_from_intuitive_inputs(compressible_solver_options,
    ComponentVector(
        vel_u = 0.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = -600.0u"m/s",
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
        mass_fractions = (
            nitrous_oxide = 0.0,
            ethylene = 0.0,
            nitrogen = 0.781,
            oxygen = 0.211,
            carbon_dioxide = 0.0003,
            carbon_monoxide = 0.0,
            hydrogen = 0.0,
            water = 0.0
        ),
        diffusion_coefficients_ambient = (
            nitrous_oxide = 1.4e-5u"m^2/s",
            ethylene = 1.4e-5u"m^2/s",
            nitrogen = 1.8e-5u"m^2/s",
            oxygen = 1.9e-5u"m^2/s",
            carbon_dioxide = 1.6e-5u"m^2/s",
            carbon_monoxide = 2.1e-5u"m^2/s",
            hydrogen = 7.6e-5u"m^2/s",
            water = 2.4e-5u"m^2/s",
        ),
        diffusion_coefficients_ambient_temperature = 21.13u"°C",
        diffusion_coefficients_ambient_pressure = 1.0u"atm"
    )
)

supersonic_inlet_initial_conditions_stripped = ustrip.(upreferred.(supersonic_inlet_initial_conditions))

add_patch!(
    config, "injector_inlet_patch";
    patch_function = function injector_inlet_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b, 
    )  
        #=
        (
            injector_mass_flow, 
            post_injector_temperature, 
            post_injector_pressure,
            post_injector_density,
            post_injector_specific_internal_energy,
            injector_velocity,
        )
        =#
        inlet_conditions = update_injector_train!(du, u, p, t, system, geo)
        
        #=
        subsonic_inlet_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a, 
            idx_b, face_b;
            inlet_conditions = inlet_conditions
        )
            =#
        
        #=
        supersonic_inlet_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a,
            idx_b, face_b,
            inlet_conditions = supersonic_inlet_initial_conditions_stripped,
        )
            =#
        
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
    patch_function = function outlet_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b, 
    )
        characteristic_farfield_flux!(
            du, u, p, t, system, geo, 
            idx_a, face_a, idx_b, face_b, 
            farfield_conditions = fluid_initial_conditions
        )
    end
)

add_patch!(
    config, "symmetry_patch";
    patch_function = no_flux!
)

ambient_air_properties = ustrip.(upreferred.(ComponentVector(
    gravity = 9.81u"m/s^2",
    coefficient_of_thermal_expansion  = 3.41e-3u"1/K",
    air_temperature = 21.13u"°C",
    thermal_conductivity = 0.024u"W/(m*K)",
    density = 1.29u"kg/m^3",
    cp = 1006.0u"J/(kg*K)",
    kinematic_viscosity = 1.506e-5u"m^2/s",
    prandtl_number = 0.707
)))

function ambient_air_heat_exchange!(
    du, u, p, t, system, geo, 
    idx_a, face_a, 
    idx_b, face_b; 
    air = nothing
)
    face_area, face_normal, face_distance, vol = boundary_geometry(geo, idx_a, face_a)

    thermal_diffusivity = air.thermal_conductivity / (air.density * air.cp)

    characteristic_length = sqrt(face_area)

    #=
    @show air.gravity
    @show air.coefficient_of_thermal_expansion
    @show u.temperature[idx_a]
    @show air.air_temperature
    @show characteristic_length
    @show air.kinematic_viscosity
    @show thermal_diffusivity
    =#

    rayleigh_number = (
        ((air.gravity * air.coefficient_of_thermal_expansion * (u.temperature[idx_a] - air.air_temperature)) * characteristic_length^3) /
        (air.kinematic_viscosity * thermal_diffusivity)
    )

    if rayleigh_number < 1e-9
        return nothing
    end

    nusselt_number = (
        0.825 + (0.387 * (rayleigh_number)^(1/6)) /
        (1 + (0.492 / air.prandtl_number)^(9/16))^(8/27)
    )^2

    convection_coefficient = nusselt_number * air.thermal_conductivity / characteristic_length

    du.volumetric_energy[idx_a] -= (convection_coefficient * face_area * (u.temperature[idx_a] - air.air_temperature)) / vol

    #NOTE: this does not currently account for Radiation, but I don't think the alunimum would get hot enough for it to matter
end

for name in [
    "aluminum_housing_wall_patch",
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
                air = ambient_air_properties
            )
        end
    )
end

wls_stencil = build_weighted_least_squares_stencil(config.geo)
sst_wall_distances = create_sst_wall_distances(grid, config.geo, 
    [
        "fuel_wall_patch",

        "fluid_domain_to_injector_bulkhead_interface",
        "fluid_domain_to_injector_disk_interface",
        "fluid_domain_to_inner_graphite_nozzle_interface",
        "fluid_domain_to_outer_aluminum_nozzle_casing_interface",
        "fluid_domain_to_phenolic_liner_interface",
        "flange_fasteners_bottom_to_fluid_domain_interface"
    ]
)

nitrous_oxide_model = PR(["nitrous oxide"])

injector_u0, injector_properties = get_injector_train_u0_and_properties(nitrous_oxide_model)

injector_train_properties_vec = Vector(deepcopy(injector_properties))
injector_train_properties_cache = DiffCache(deepcopy(injector_train_properties_vec))
injector_train_properties_axes = getaxes(injector_properties)

valve_opening = 1.0

post_injector_cell_ids = [face.idx[1] for face in getfacetset(grid, "injector_inlet_patch")]

additional_data = (
    #Navier stokes data
    compressible_solver_options = compressible_solver_options,
    wls_stencil = wls_stencil,
    sst_wall_distances = sst_wall_distances,
    
    #Injector train properties
    injector_train_properties_vec = injector_train_properties_vec,
    injector_train_properties_cache = injector_train_properties_cache,
    injector_train_properties_axes = injector_train_properties_axes,
    valve_opening = valve_opening,
    nitrous_oxide_model = nitrous_oxide_model,
    post_injector_cell_ids = post_injector_cell_ids
)

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, additional_data, check_units = false);

u_named = ComponentVector(u0_vec, system.state_axes);

function solve_system!(du, u, p, t, system, geo)
    update_region_groups!(du, u, p, t, system, geo)

    #update_weighted_least_squares_gradients!(u, system.wls_stencil)
    
    if system.compressible_solver_options.use_sst_turbulence == true
        update_sst_closure!(u, system.sst_wall_distances)
    end
    
    #update_MUSCL_gradients!(u, wls_stencil)
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)

    #=
    @show findall(isnan, du.density)
    @show findall(isnan, u.density)
    @show findall(isnan, u.volumetric_energy)
    @show findall(isnan, du.volumetric_energy)
    @show findall(isnan, u.momentum_density_u)
    @show findall(isnan, du.momentum_density_u)
    @show findall(isnan, u.momentum_density_v)
    @show findall(isnan, du.momentum_density_v)
    @show findall(isnan, u.momentum_density_w)
    @show findall(isnan, du.momentum_density_w)
    @show findall(isnan, u.turbulent_kinetic_energy_density)
    @show findall(isnan, du.turbulent_kinetic_energy_density)
    @show findall(isnan, u.specific_dissipation_rate_density)
    @show findall(isnan, du.specific_dissipation_rate_density)
    @show findall(isnan, u.species_densities[1:end])
    @show u.tank
    @show u.mid_section=#
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

p_guess = 0.0

Revise.includet(joinpath(@__DIR__, "state_is_invalid/state_is_invalid.jl"))

state_is_invalid_closure = (u, p, t) -> state_is_invalid(u, p, t, system);

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

t0 = 0.0
switch_time = 0.02
tMax = 10.0

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
startup_prob = ODEProblem(ode_func, u0_vec, (t0, switch_time), p_guess)

callbacks = CallbackSet(
    approximate_time_to_finish_cb,
    #increase_dtmax_cb,
)

Revise.includet(joinpath(@__DIR__, "state_is_invalid/state_is_invalid.jl"))

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
    #FBDF(linsolve = KrylovJL_GMRES(), autodiff = AutoForwardDiff()),
    #AutoTsit5(FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true)),
    #FBDF(linsolve = KrylovJL_GMRES(), nlsolve = NLNewton(relax = 0.7), precs = iluzero, concrete_jac = true),
    callback = callbacks,
    isoutofdomain = state_is_invalid_closure,
)

sol_combined = (u = startup_sol.u, t = startup_sol.t)
#sol_combined = (u = [startup_sol.u..., steady_sol.u...], t = [startup_sol.t..., steady_sol.t...])

desired_shown_steps = 250
save_fraction = desired_shown_steps / length(sol_combined.t)

reduced_sol = (
    u = [[sol_combined.u[i] for i in 1:(Int(ceil(1/save_fraction))):length(sol_combined.u)]..., sol_combined.u[end]],
    t = [sol_combined.t[1:(Int(ceil(1/save_fraction))):length(sol_combined.t)]..., sol_combined.t[end]]
)

du_named, u_named = regenerate_fvm_state(reduced_sol, system, solve_system!, geo, p_guess, save_properties_and_caches = true, track_progress = true);

du_named[end].density[2421]

add_xyz_vec_to_u_named!(u_named, :velocity, :vel_u, :vel_v, :vel_w)

root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\julia_cfd_output_files"

sol_to_vtk(reduced_sol, du_named, u_named, grid, geo, @__FILE__, root_dir, include_zeros_fields = false, track_progress = true)
