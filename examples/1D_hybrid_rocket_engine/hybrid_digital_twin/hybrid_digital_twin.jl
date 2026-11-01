using Revise

using FVMFramework

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
using Clapeyron
using PreallocationTools
using DataInterpolations

object_names = [
    "injector_bulkhead_top",
    "injector_bulkhead_bottom",
    "tank",
    "mid_section",
    "pre_chamber",
    "phenolic_liner_top",
    "aluminum_housing_top",

    #(1-n cells)
    "fluid_domain",
    "fuel_grain",
    "phenolic_liner",
    "aluminum_housing",
    #(1-n cells) end

    #"recirculation_zone", #abstraction; not a physical volume
    #on second thought, I think we're just going to apply some extra physics to the last fluid domain cell

    "nozzle_top",
    "nozzle_bottom",
    "aluminum_housing_near_nozzle_top",
    "aluminum_housing_near_nozzle_bottom",
    "nozzle_retainer"
]

rom_geometry_json_path = joinpath(@__DIR__, "rom_geometry", "hybrid_rom_geometry.json")

rom_data = read_rom_geometry_json(rom_geometry_json_path);

add_virtual_cell!(rom_data;
    name = "nitrous_oxide_tank",
    volume = ustrip(upreferred(1.6u"L")),
    faces = [(
        neighbor_region = "mid_section", 
        neighbor_axial_index = 1,
        area = 0.0, 
        half_cell_distance = 0.0, 
        neighbor_half_cell_distance = 0.0
    )]
)

geo = build_rom_fvm_geometry(rom_data);

n_cells = geo.n_cells

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3",
    temperature = zeros(n_cells)u"K",
    mass_fractions = (
        nitrous_oxide = zeros(n_cells)u"kg/kg",
        ethylene = zeros(n_cells)u"kg/kg",
        nitrogen = zeros(n_cells)u"kg/kg",
        oxygen = zeros(n_cells)u"kg/kg",
        carbon_dioxide = zeros(n_cells)u"kg/kg",
        water = zeros(n_cells)u"kg/kg",
    )
)

config = create_fvm_config_from_geo(geo, u_proto);

nitrous_oxide_model = PR(["nitrous oxide"])
model = ReidIdeal([
    "nitrous oxide",
    "ethylene",
    "nitrogen",
    "oxygen",
    "carbon dioxide",
    "water"
])

#We're going to have to append this to system.additional_data
method = ChemPotDewPressure(
    noncondensables = ["ethylene", "nitrogen", "oxygen", "carbon_dioxide", "water"] #I don't think water will ever condense
)

config.data = merge(config.data, (clapeyron_method = method, nitrous_oxide_model = nitrous_oxide_model, model = model))

injector_u0, injector_properties = get_injector_train_u0_and_properties(nitrous_oxide_model)

u_proto = merge_properties(u_proto, injector_u0)

config = create_fvm_config()

add_setup_syms!(
    config;
    cache_syms_and_units = (
        # Variables that will be capped to a change in a state vector later:
        heat = u"J", 
        nitrous_oxide_vapor_quality = u"kg/kg",
        
        #used for collecting heat fluxes (this is capped to temperature change)

        #updated based on temperature, density, and composition:
        pressure = u"Pa",

        #might be updated based on temp in the future, but for now they're fixed
        #k = u"W/(m*K)",
        #mu = u"Pa*s",
        
        #updated based on temperature, density, and composition:
        cp = u"J/(kg*K)",
        cv = u"J/(kg*K)",
        
        #updated based on composition:
        mw_avg = u"g/mol",

        molar_concentrations = (
            nitrous_oxide = u"mol/m^3",
            ethylene = u"mol/m^3",
            nitrogen = u"mol/m^3",
            oxygen = u"mol/m^3",
            carbon_dioxide = u"mol/m^3",
            water = u"mol/m^3",
        ),
        
        #=
        #These will just stay fixed for now
        diffusion_coefficients = (
            nitrous_oxide = u"m^2/s",
            ethylene = u"m^2/s",
            nitrogen = u"m^2/s",
            oxygen = u"m^2/s",
            carbon_dioxide = u"m^2/s",
            water = u"m^2/s",
        )
        =#

        #Caches specifically for the injector delivery train
        specific_internal_energy = 0.0u"J/kg",
        specific_enthalpy = 0.0u"J/kg",
        injector_velocity = 0.0u"m/s", #this is only used to observe injector velocity over time, it does not affect any physics (yet!)
    ),
    special_caches = (
        
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(
        # Nitrous Oxide Decomposition kinetics:
        nitrous_oxide_decomposition_pre_exponential_factor = 4.2e9u"s^-1",
        nitrous_oxide_decomposition_activation_energy = 222_000u"J/mol",
        nitrous_oxide_decomposition_heat_of_reaction = -82_000u"J/mol",

        # Westbrook kinetics:
        westbrook_dryer_step_1_pre_exponential_factor = 7.5895e7u"m^2.25 * mol^-0.75 * s^-1",
        westbrook_dryer_step_1_activation_energy = 125_520u"J/mol",
        westbrook_dryer_step_1_heat_of_reaction = -757_120u"J/mol",

        westbrook_dryer_step_2_forward_pre_exponential_factor = 1.2589e10u"m^2.25 * mol^-0.75 * s^-1",
        westbrook_dryer_step_2_forward_activation_energy = 167_360u"J/mol",
        westbrook_dryer_step_2_forward_heat_of_reaction = -282_980u"J/mol",

        westbrook_dryer_step_2_reverse_pre_exponential_factor = 5.0000e8u"s^-1",
        westbrook_dryer_step_2_reverse_activation_energy = 167_360u"J/mol",
        westbrook_dryer_step_2_reverse_heat_of_reaction = 282_980u"J/mol",

        # Simple kinetics:
        simple_arrhenius_reaction_model_pre_exponential_factor = 6.3245e7u"m^2.25 * mol^-0.75 * s^-1",
        simple_arrhenius_reaction_model_activation_energy = 125_520u"J/mol",
        simple_arrhenius_reaction_model_heat_of_reaction = -1_323_000u"J/mol",

        # Turbulence constants: (ranges from around 1.0 to 4.0, usually closer from 1.0 to 2.0)
        turbulence_limited_reaction_rate_mix_parameter = 1.5,
        
        # Regression Parameters:
        regression_coefficient_a = 1.1573e-4, #units are m^(2n+1) * kg^(-n) * s^(-1)
        regression_exponent_n = 0.331,
        fuel_heat_of_pyrolysis = 2.3u"MJ/kg",
        min_fuel_regression_temperature = 800.0u"K",

        # Injector System Parameters:
        injector_orifice_area = 4.0u"mm^2",
        injector_discharge_coefficient = 0.7,
        valve_flow_capacity_factor = 10.0u"mm^2",
        valve_opening = 1.0
    ),
)

function cap_heat_flux_to_temp_change!(du, u, p, t, system, geo, cell_id)
    du.temperature[cell_id] += du.heat[cell_id] / (geo.cell_volumes[cell_id] * u.density[cell_id] * u.cp[cell_id])
end

function overall_combustion_model!(du, u, p, t, system, geo, cell_id)
    decompose_nitrous_oxide!(du, u, p, t, system, geo, cell_id)
    combust_ethylene_and_oxygen_westbrook_dryer!(du, u, p, t, system, geo, cell_id) 
    #combust_ethylene_and_oxygen_simple!(du, u, p, t, system, geo, cell_id) #is also an option
end

Revise.includet(joinpath(@__DIR__, "property_update_functions", "system_moles_vec.jl"))
Revise.includet(joinpath(@__DIR__, "property_update_functions", "update_fluid_caches.jl"))


#Don't need this right now since k and mu are fixed
#=
function manual_fluid_cache_updates!(du, u, p, t, system, geo, cell_id)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.mu[cell_id] = properties.mu[cell_id]
end
=#

function update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    #manual_fluid_cache_updates!(du, u, p, t, system, geo, cell_id)
    fluid_property_update_ideal!(du, u, p, t, system, geo, cell_id)
end

function fluid_domain_physics!(du, u, p, t, system, geo, cell_id)
    overall_combustion_model!(du, u, p, t, system, geo, cell_id)
    cap_heat_flux_to_temp_change!(du, u, p, t, system, geo, cell_id)
    cap_fluid_flow!(du, u, p, t, system, geo, cell_id)
end

function normalize!(vec)
    total = sum(vec)
    vec ./= total
end

room_temp = 294.13u"K"

air_initial_conditions = ComponentVector(
    density = 1.1225u"kg/m^3",
    temperature = room_temp,
    mass_fractions = normalize!((
        nitrous_oxide = 0.0,
        ethylene = 0.0,
        nitrogen = 0.781,
        oxygen = 0.211,
        carbon_dioxide = 0.0003,
        carbon_monoxide = 0.0,
        hydrogen = 0.0,
        water = 1e-8
    ))
)

air_properties = ComponentVector(
    #We're assuming these two stay constant which isn't that accurate, but oh well
    k = 0.024u"W/(m*K)", 
    mu = 1.9e-5u"Pa*s", 

    molecular_weights = (
        nitrous_oxide = 44.01u"g/mol", 
        ethylene = 28.05u"g/mol", 
        nitrogen = 28.01u"g/mol", 
        oxygen = 32.00u"g/mol", 
        carbon_dioxide = 44.01u"g/mol", 
        water = 18.02u"g/mol"
    ),
    #these will just stay fixed for now
    diffusion_coefficients = (
        nitrous_oxide = 1.4e-5u"m^2/s",
        ethylene = 1.4e-5u"m^2/s",
        nitrogen = 1.8e-5u"m^2/s",
        oxygen = 1.9e-5u"m^2/s",
        carbon_dioxide = 1.6e-5u"m^2/s",
        carbon_monoxide = 2.1e-5u"m^2/s",
        hydrogen = 7.6e-5u"m^2/s",
        water = 2.4e-5u"m^2/s",
    ),
)

add_region!(
    config, "pre_chamber";
    type = Fluid(),
    initial_conditions = air_initial_conditions,
    properties = air_properties,
    property_update_function = 
    function update_pre_chamber_region!(du, u, p, t, system, geo, cell_id)
        update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function pre_chamber_region!(du, u, p, t, system, geo, cell_id)
        fluid_domain_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "fluid_domain";
    type = Fluid(),
    initial_conditions = air_initial_conditions,
    properties = air_properties,
    property_update_function = 
    function update_fluid_domain_region!(du, u, p, t, system, geo, cell_id)
        update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function fluid_domain_region!(du, u, p, t, system, geo, cell_id)
        fluid_domain_physics!(du, u, p, t, system, geo, cell_id)

        if cell_id == get_cells(geo, "fluid_domain")[end]
            #TODO: do some CFD to create a closure for the recirculation zone here
            #recirculation_zone_closure(du, u, p, t, system, geo, cell_id)
        end
    end
)

add_region!(
    config, "tank",
    type = Fluid(),
    initial_conditions = air_initial_conditions,
    properties = air_properties,
    property_update_function = 
    function update_tank_region!(du, u, p, t, system, geo, cell_id)
        update_tank_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function tank_region!(du, u, p, t, system, geo, cell_id)
        fluid_domain_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "mid_section",
    type = Fluid(),
    initial_conditions = air_initial_conditions,
    properties = air_properties,
    property_update_function = 
    function update_mid_section_region!(du, u, p, t, system, geo, cell_id)
        update_mid_section_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function mid_section_region!(du, u, p, t, system, geo, cell_id)
        fluid_domain_physics!(du, u, p, t, system, geo, cell_id)
    end
)

function update_solid_properties!(du, u, p, t, system, geo, cell_id)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    #u.k[cell_id] = properties.k[cell_id] #not needed right now, but will be in the future
    u.cp[cell_id] = properties.cp[cell_id]

end

function solid_physics!(du, u, p, t, system, geo, cell_id)
    cap_heat_flux_to_temp_change!(du, u, p, t, system, geo, cell_id)
end

aluminum_initial_conditions = ComponentVector(
    density = 2700.0u"kg/m^3",
    temperature = room_temp
)

aluminum_properties = ComponentVector(
    cp = 896.0u"J/(kg*K)",
    k = 167.0u"W/(m*K)"
)

add_region!(
    config, "injector_bulkhead_top";
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = 
    function update_injector_bulkhead_top_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function injector_bulkhead_top_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "injector_bulkhead_bottom";
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = 
    function update_injector_bulkhead_bottom_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function injector_bulkhead_bottom_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "aluminum_housing_top",
)

add_region!(
    config, "aluminum_housing",
)

phenolic_liner_initial_conditions = ComponentVector(
    density = 1250.0u"kg/m^3",
    temperature = room_temp
)

phenolic_liner_properties = ComponentVector(
    cp = 1400.0u"J/(kg*K)",
    k = 0.25u"W/(m*K)"
)

add_region!(
    config, "phenolic_liner_top",
    type = Solid(),
    initial_conditions = phenolic_liner_initial_conditions,
    properties = phenolic_liner_properties,
    property_update_function = 
    function update_phenolic_liner_top_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function phenolic_liner_top_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "phenolic_liner",
    initial_conditions = phenolic_liner_initial_conditions,
    properties = phenolic_liner_properties,
    property_update_function = 
    function update_phenolic_liner_top_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function phenolic_liner_top_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

graphite_initial_conditions = ComponentVector(
    density = 1760.0u"kg/m^3",
    temperature = room_temp
)

graphite_properties = ComponentVector(
    cp = 710.0u"J/(kg*K)",
    k = 130.0u"W/(m*K)"
)

add_region!(
    config, "nozzle_top",
    type = Solid(),
    initial_conditions = graphite_initial_conditions,
    properties = graphite_properties,
    property_update_function = 
    function update_nozzle_top_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function nozzle_top_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "aluminum_housing_near_nozzle_top",
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = 
    function update_aluminum_housing_near_nozzle_top_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function aluminum_housing_near_nozzle_top_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "nozzle_bottom",
    type = Solid(),
    initial_conditions = graphite_initial_conditions,
    properties = graphite_properties,
    property_update_function = 
    function update_nozzle_bottom_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function nozzle_bottom_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "aluminum_housing_near_nozzle_bottom",
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = 
    function update_aluminum_housing_near_nozzle_bottom_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function aluminum_housing_near_nozzle_bottom_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "nozzle_retainer",
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = 
    function update_nozzle_retainer_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function nozzle_retainer_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_patch!(
    owners = "tank",
    connections = "mid_section",
    patch_function = function tank_to_mid_section_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        adjustable_valve_flow!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

#This patch is where the nitrous_oxide injection will happen
add_patch!(
    owners = "mid_section",
    connections = "injector_bulkhead_top",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        #might need a CFD-informed closure due to the convective heat transfer
        #mid_section_to_injector_bulkhead_top_closure_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) 

        #for now, we'll just use regular heat_diffusion
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)

        #this is less pressing than developing a CFD closure for heat transfer from the inner combustion gases 
        #to other components
    end
)

add_patch!(
    owners = "mid_section",
    connections = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        #might need a CFD-informed closure due to the convective heat transfer
        #mid_section_to_injector_bulkhead_top_closure_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) 

        #for now, we'll just use regular heat_diffusion
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)

        #this is less pressing than developing a CFD closure for heat transfer from the inner combustion gases 
        #to other components
    end
)

add_patch!(
    owner = "injector_bulkhead_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "mid_section",
    connection = "pre_chamber",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        injector_valve_flow!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) 
        #this would be the valve equation determining flow between the mid section and first chamber cell
        #purely algebraic
    end
)

#Axial fluxes from top cells to injector_bulkhead_bottom:
#start
add_patch!(
    owner = "pre_chamber",
    connection = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "phenolic_liner_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "aluminum_housing_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)
#end

#Radial fluxes across top cells
#start
add_patch!(
    owner = "pre_chamber",
    connection = "phenolic_liner_top",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        #this one will likely need a CFD closure in the future to account for convection
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "phenolic_liner_top",
    connection = "aluminum_housing_top",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

Revise.includet(joinpath(@__DIR__, "physics", "pressure_driven_mass_flux", "pressure_driven_mass_flux.jl"))

add_patch!(
    owner = "pre_chamber",
    connection = "fluid_domain", #1-n cells,
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        #this is where the most complicated CFD closure will likely need to be to determine the mass flux between these two regions
        pressure_driven_mass_flux!
    end
)

#end

#We could add a flux between the top fluid domain and the fuel grain, but I'm going to ignore it

#Axial fluxes along n cells
#start

#these will be added to solve_system!
function fluid_domain_axial_flux!(du, u, p, t, system, geo)
    fluid_domain_cells = get_cells(geo, "fluid_domain")

    #if we use this, we'll have n-1 axial fuel port id states
    fuel_port_id = 1

    #this requires that fluid_domain cells are ordered so that the downstream cell adjacent to the current cell has an id of upstream_cell_id + 1
    for upstream_cell in fluid_domain_cells[1:end-1]
        downstream_cell = upstream_cell + 1
        #some flux function here

        fuel_port_id = system.fuel_port_ids[upstream_cell, downstream_cell]

        fluid_port_area = pi * u.fuel_port_radius[fuel_port_id]^2
        fuel_grain_area = pi * (geo.max_bounding_box_x[upstream_cell] - geo.centroids[upstream_cell][1])
        interface_area = fuel_grain_area - fluid_port_area

        #use the new geometry to calculate the flux

        pressure_driven_mass_flux!

        fuel_port_id += 1
    end
end

function fuel_grain_axial_flux!(du, u, p, t, system, geo)
    fuel_grain_cells = get_cells(geo, "fuel_grain")

    #this requires that fuel_grain cells are ordered so that the downstream cell adjacent to the current cell has an id of upstream_cell_id + 1
    for upstream_cell in fuel_grain_cells[1:end-1]
        downstream_cell = upstream_cell + 1
        #some flux function here
    end
end

function phenolic_liner_axial_flux!(du, u, p, t, system, geo)
    phenolic_liner_cells = get_cells(geo, "phenolic_liner")

    #this requires that phenolic_liner cells are ordered so that the downstream cell adjacent to the current cell has an id of upstream_cell_id + 1
    for upstream_cell in phenolic_liner_cells[1:end-1]
        downstream_cell = upstream_cell + 1
        #some flux function here
    end
end

function aluminum_housing_axial_flux!(du, u, p, t, system, geo)
    aluminum_housing_cells = get_cells(geo, "aluminum_housing")

    #this requires that aluminum_housing cells are ordered so that the downstream cell adjacent to the current cell has an id of upstream_cell_id + 1
    for upstream_cell in aluminum_housing_cells[1:end-1]
        downstream_cell = upstream_cell + 1
        #some flux function here
    end
end

function axial_flux_functions!(du, u, p, t, system, geo)
    fluid_domain_axial_flux!(du, u, p, t, system, geo)
    fuel_grain_axial_flux!(du, u, p, t, system, geo)
    phenolic_liner_axial_flux!(du, u, p, t, system, geo)
    aluminum_housing_axial_flux!(du, u, p, t, system, geo)
end
#end

#Radial fluxes along cells
#start
add_patch!(
    owner = "fluid_domain",
    connection = "fuel_grain",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        #some flux function 
        #this requires geometry updating
    end
)

add_patch!(
    owner = "fuel_grain",
    connection = "phenolic_liner",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        #some flux function 
        #this requires geometry updating
    end
)

add_patch!(
    owner = "phenolic_liner",
    connection = "aluminum_housing",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        #some flux function 
        #this requires geometry updating
    end
)
#end

add_patch!(
    owner = "fluid_domain",
    connection = "nozzle_top",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        nozzle_outlet_closure!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
        heat_transfer_from_fluid_to_nozzle_top!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

#we could also add a flux between the nozzle top and phenolic liner but I don't think it's necessary

#This is a radial flux between the nozzle graphite and aluminum housing
add_patch!(
    owner = "nozzle_top",
    connection = "aluminum_housing_nozzle",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        nozzle_top_to_aluminum_housing_nozzle_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "nozzle_top",
    connection = "nozzle_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "nozzle_bottom",
    connection = "aluminum_housing_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "nozzle_bottom",
    connection = "nozzle_retainer",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "aluminum_housing_bottom",
    connection = "nozzle_retainer",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

function solve_system!(du, u, p, t, system, geo)
    update_fuel_regression_geometry!(u, system.fuel_regression_model, geo, system.fuel_regression_grid)

    update_region_groups!(du, u, p, t, system, geo)
    
    update_sytem_moles_and_mole_fractions_vec(du, u, p, t, system, geo) #This has to be called after update_region_groups! is called

    axial_flux_functions!(du, u, p, t, system, geo)
    
    #update_MUSCL_gradients!(u, wls_stencil)
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    regress_fuel_layers!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end


system_data = (
    moles_vec = generate_system_moles_vec(config),
    mole_fractions_vec = generate_system_mole_fractions_vec(config)
)


