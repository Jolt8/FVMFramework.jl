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

    "recirculation_zone", #abstraction; not a physical volume

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
    mass = zeros(n_cells)u"kg",
    temperature = zeros(n_cells)u"kg",
    mass_fractions = (
        nitrous_oxide = zeros(n_cells)u"kg/kg",
        ethylene = zeros(n_cells)u"kg/kg",
        nitrogen = zeros(n_cells)u"kg/kg",
        oxygen = zeros(n_cells)u"kg/kg",
        carbon_dioxide = zeros(n_cells)u"kg/kg",
        water = zeros(n_cells)u"kg/kg",
    ),
    nitrous_oxide_vapor_fraction = zeros(n_cells)u"kg",
    nitrous_oxide_liquid_mass = zeros(n_cells)u"kg",
    #not sure if I should keep track of the nitrous oxide vapor fraction or liquid mass
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

room_temp = 294.13u"K"

fluid_domain_volume = geo.cell_volumes[fluid_domain[1]]

air_initial_conditions, air_properties, air_cache_syms_and_units, air_special_caches = get_air_properties(fluid_domain_volume, room_temp)

add_setup_syms!(
    config;
    cache_syms_and_units = (
        heat = u"J", #used for collecting heat fluxes (this is capped to temperature change)

        #updated based on temperature, density, and composition:
        density = u"kg/m^3",
        pressure = u"Pa",
        speed_of_sound = u"m/s",

        #updated based on temp:
        k = u"W/(m*K)",
        mu = u"Pa*s",
        
        #updated based on temperature, density, and composition:
        cp = u"J/(kg*K)",
        cv = u"J/(kg*K)",

        #updated based on cp and cv
        gamma = u"1",
        
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
        vapor_fraction = 0.0u"kg", #we can't use u.nitrous_oxide_vapor_fraction because we can't set a state variable
        injector_velocity = 0.0u"m/s", #this is only used to observe injector velocity over time, it does not affect any physics
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
        min_fuel_regression_temperature = 800.0u"K"

        # Injector System Parameters:
        injector_orifice_area = 4.0u"mm^2",
        injector_discharge_coefficient = 0.7,
        valve_flow_capacity_factor = 10.0u"mm^2",
        valve_opening = 1.0
    ),
)

function overall_combustion_model!(du, u, p, t, system, geo, cell_id)
    decompose_nitrous_oxide!(du, u, p, t, system, geo, cell_id)
    combust_ethylene_and_oxygen_westbrook_dryer!(du, u, p, t, system, geo, cell_id) 
    #combust_ethylene_and_oxygen_simple!(du, u, p, t, system, geo, cell_id) #is also an option
end

function update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    overall_fluid_property_update!(du, u, p, t, system, geo, cell_id)
    update_species_molar_concentrations!(du, u, p, t, system, geo, cell_id)
    update_fluid_diffusion_coefficients!(du, u, p, t, system, geo, cell_id)
end

function fluid_domain_physics!(du, u, p, t, system, geo, cell_id)
    overall_combustion_model!(du, u, p, t, system, geo, cell_id)
    cap_fluid_flow!(du, u, p, t, system, geo, cell_id)
end

add_region!(
    config, "fluid_domain";
    type = Fluid(), #Not sure if we even need a type, my gut feeling says no
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
            recirculation_zone_closure(du, u, p, t, system, geo, cell_id)
        end
    end
)

add_patch!(
    owners = "tank",
    connections = "mid_section",
    patch_function = function tank_to_mid_section_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        mid_section_to_injector_bulkhead_top_closure_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
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
        mid_section_to_injector_bulkhead_top_closure_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owners = "mid_section",
    connections = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        mid_section_to_injector_bulkhead_bottom_closure_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "injector_bulkhead_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(
    owner = "mid_section",
    connection = "pre_chamber",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        mid_section_to_pre_chamber_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) 
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
        #flux function here
    end
)

add_patch!(
    owner = "phenolic_liner_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        #flux function here
    end
)

add_patch!(
    owner = "aluminum_housing_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        #flux function here
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
        #flux function here
    end
)

add_patch!(
    owner = "phenolic_liner_top",
    connection = "aluminum_housing_top",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        #flux function here
    end
)
#end


add_patch!(
    owner = "pre_chamber",
    connection = "fluid_domain", #1-n cells,
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        #flux function here
    end
)

#We could add a flux between the top fluid domain and the fuel grain, but I'm going to ignore it

#Axial fluxes along n cells
#start
#=
add_patch!(
    owner = "fluid_domain",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        axial_fluid_domain_flux!(du, u, p, t, system, interface_area, upstream_cell, downstream_cell)
    end
)
=#

#these will be added to 
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

    axial_flux_functions!(du, u, p, t, system, geo)
    
    #update_MUSCL_gradients!(u, wls_stencil)
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    regress_fuel_layers!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end