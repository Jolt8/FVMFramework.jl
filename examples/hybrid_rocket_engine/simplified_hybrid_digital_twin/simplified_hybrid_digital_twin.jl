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
using JSON

using Roots

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

cad_table_data_json_path = joinpath(@__DIR__, "rom_and_table_inputs", "table_data", "freecad_table_data.json")
Revise.includet(joinpath(@__DIR__, "rom_and_table_inputs", "table_data", "read_freecad_table_data.jl"))
cad_data = ComponentVector(read_freecad_table_data_json(cad_table_data_json_path))

rom_geometry_json_path = joinpath(@__DIR__, "rom_and_table_inputs", "rom_geometry", "hybrid_rom_geometry.json")

rom_data = read_rom_geometry_json(rom_geometry_json_path);

add_virtual_cell!(rom_data;
    name = "tank",
    volume = ustrip(upreferred(1.6u"L")),
    faces = [(
        neighbor_region = "mid_section", 
        neighbor_axial_index = 1,
        area = 0.0, 
        half_cell_distance = 0.0, 
        neighbor_half_cell_distance = 0.0
    )]
)

add_virtual_face!(rom_data;
    region_a = "tank",
    axial_index_a = 1,
    region_b = "pre_chamber",
    axial_index_b = 1,
    area = ustrip(upreferred(cad_data.injector_orifice_area)), 
    half_cell_distance_a = 0.0, 
    half_cell_distance_b = 0.0
)

geo = build_rom_fvm_geometry(rom_data);

n_cells = geo.n_cells

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3", #fluid's still have density
    state_temperature = zeros(n_cells)u"K", #fluid temperatures are now based on NASA CEA
    igniter_mass = [0.0u"kg"]
)

config = create_fvm_config_from_geo(geo, u_proto);

nozzle_inlet_radius = cad_data.nozzle_inlet_diameter / 2
nozzle_throat_radius = cad_data.nozzle_throat_diameter / 2
nozzle_exit_radius = cad_data.nozzle_exit_diameter / 2

nozzle_converging_height = cad_data.nozzle_converging_length
nozzle_diverging_height = cad_data.nozzle_diverging_length

converging_nozzle_surface_area = pi * (nozzle_inlet_radius + nozzle_throat_radius) * sqrt((nozzle_inlet_radius - nozzle_throat_radius)^2 + nozzle_converging_height^2)
diverging_nozzle_surface_area = pi * (nozzle_throat_radius + nozzle_exit_radius) * sqrt((nozzle_exit_radius - nozzle_throat_radius)^2 + nozzle_diverging_height^2)

Revise.includet(joinpath(@__DIR__, "overloads", "clapeyron_tracer_overloads.jl"))

nitrous_oxide_model = PR(["nitrous oxide"])
model = ReidIdeal([
    "nitrous oxide",
    "ethylene",
    "nitrogen",
    "oxygen",
    "carbon dioxide",
    "carbon monoxide",
    "water",
    "hydrogen"
])

#We're going to have to append this to system.additional_data
method = ChemPotDewPressure(
    noncondensables = ["ethylene", "nitrogen", "oxygen", "carbon dioxide", "carbon monoxide", "water", "hydrogen"] #I don't think water will ever condense
)

config.data = merge(config.data, (
    cad_data = ustrip.(upreferred.(cad_data)),
    converging_nozzle_surface_area = ustrip(upreferred(converging_nozzle_surface_area)),
    diverging_nozzle_surface_area = ustrip(upreferred(diverging_nozzle_surface_area)),
    clapeyron_method = method, 
    nitrous_oxide_model = nitrous_oxide_model, 
    model = model,
    tstops = [],
    cold_flow = true,
    hot_fire = true,
    burn_igniter = false
));

Revise.includet(joinpath(@__DIR__, "cea_retrieval", "CEA_lookup_table.jl"));

config.data = merge(config.data, (
    chamber_cp_interp = chamber_cp_interp,
    chamber_mu_interp = chamber_mu_interp,
    chamber_k_interp = chamber_k_interp,
    chamber_prandtl_interp = chamber_prandtl_interp,
    chamber_isp_vac_interp = chamber_isp_vac_interp,
    chamber_cstar_interp = chamber_cstar_interp,
    chamber_T_comb_interp = chamber_T_comb_interp,
    chamber_mw_avg_interp = chamber_mw_avg_interp,
    chamber_gamma_interp = chamber_gamma_interp,

    mass_fractions_interps_vec = mass_fractions_interps_vec
));

add_setup_syms!(
    config;
    cache_syms_and_units = (
        temperature = u"K", #this will now be updated using the state_temperature state variable for solids and NASA CEA tables for fluids
        # Variables that will be capped to a change in a state vector later:
        heat = u"J",
        
        #Fuel regression:
        injector_oxidizer_mass_flow = u"kg/s", #stored at the pre-chamber cell by injector_valve_flow!; drives the regression law
        fuel_regression_rate = u"m/s", #diagnostic, stored at each fuel grain cell
        pyrolyzed_fuel_mass_flow = u"kg/s", #diagnostic, stored at each port gas cell
        nitrous_oxide_vapor_quality = u"kg/kg",
        
        #used for collecting heat fluxes (this is capped to temperature change)

        #updated based on temperature, density, and composition:
        pressure = u"Pa",

        #might be updated based on temp in the future, but for now they're fixed
        #these are now updated with CEA data
        k = u"W/(m*K)",
        mu = u"Pa*s",
        
        #updated based on temperature, density, and composition:
        cp = u"J/(kg*K)",
        cv = u"J/(kg*K)",
        
        #updated based on composition:
        mw_avg = u"g/mol",

        mass_fractions = (
            nitrous_oxide = u"kg/kg",
            ethylene = u"kg/kg",
            nitrogen = u"kg/kg",
            oxygen = u"kg/kg",
            carbon_dioxide = u"kg/kg",
            carbon_monoxide = u"kg/kg",
            water = u"kg/kg",
            hydrogen = u"kg/kg"
        ),

        #Caches specifically for the injector delivery train
        specific_internal_energy = 0.0u"J/kg",
        specific_enthalpy = 0.0u"J/kg",
        outflow_density = u"kg/m^3", #density of what leaves the cell (liquid from the tank bottom); kept separate from the density state
        injector_velocity = 0.0u"m/s", #this is only used to observe injector velocity over time, it does not affect any physics (yet!)

        #TODO: REPLACE THE PRE-CHAMBER TEMPERATURE FLOOR WITH A REAL VAPORIZATION MODEL.
        neglected_injector_cooling_power = u"W", #TODO: latent heat per second the pre-chamber temperature floor is currently not changing. Diagnostic only.

        injector_mass_flow = u"kg/s" #to track injector mass flow over time
    ),
    special_caches = ComponentVector(
        
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(
        # Regression Parameters:
        oxidizer_to_fuel_ratio = 7.0,
        regression_coefficient_a = 1.1573e-4, #units are m^(2n+1) * kg^(-n) * s^(-1)
        regression_exponent_n = 0.331,
        fuel_heat_of_pyrolysis = 2.3u"MJ/kg",
        min_fuel_regression_temperature = 800.0u"K", #gas temperature next to the fuel where regression starts ramping on (full at +100 K)
        fuel_surface_temperature = 750.0u"K", #temperature the pyrolyzed fuel gas leaves the surface at (typical HDPE pyrolysis surface: 700-800 K)

        # Injector System Parameters:
        injector_area = cad_data.injector_orifice_area,
        injector_discharge_coefficient = 0.7,
        valve_flow_capacity_factor = 10.0u"mm^2",
        valve_opening = 1.0,
        minimum_pre_chamber_temperature = 230.0u"K", #temperature threshold where the latent heat of vaporization is ignored in injector_deliver_train_physics
        valve_open_time = 0.1u"s",
        
        # Port Flow Parameters:
        port_darcy_friction_factor = 0.02, #Darcy friction factor of the fuel port wall (smooth tube at Re ~ 1e5); used by pressure_driven_mass_flux

        # Nozzle Parameters:
        nozzle_discharge_coefficient = 0.98,
        ambient_pressure = 1.0u"atm",
        expansion_ratio = 4.0,

        # Igniter Parameters:
        igniter_density = 1.90u"g/cm^3", #based on the density of the sugar rocket recipe
        igniter_heat_of_combustion = 2.30u"kJ/g", #based on the middle value of typical heat of combustion for the sugar rocket recipe
        igniter_burn_time = 0.5u"s",
        igniter_volume = cad_data.igniter_puck_volume |> u"m^3",

        # Other Parameters:
        motor_burn_time = 6.0u"s",
        cea_switch_time = 0.6u"s",
        chamber_pressure = 20.0u"bar",

        # Stop Parameters
        #Largest port radius allowed before the burn stops, leaving the desired web of fuel to protect the liner
        max_fuel_port_radius = (cad_data.phenolic_liner_inner_diameter - 2 * cad_data.fuel_grain_desired_fuel_web_thickness) / 2,
        min_chamber_pressure = 100_000u"Pa", # Solver stops once the pressure drops below this value (also requires 2 seconds to pass before this is allowed to stop the solver)

        # Convection Parameters:
        chamber_convective_heat_transfer_coefficient = 3000.0u"W/(m^2*K)",


        #Very rough heat transfer coefficients right now because we're definitely going to develop a closure for this using CFD
        #We're going slightly high on the converging and diverging numbers because we're neglecting the typical throat heat transfer coefficient

        nozzle_converging_heat_transfer_coefficient = 4000u"W/(m^2*K)", # (ranges between 1000 and 5000)
        nozzle_throat_heat_transfer_coefficient = 12000u"W/(m^2*K)", # (ranges between 10000 and 20000)
        nozzle_diverging_heat_transfer_coefficient = 6000u"W/(m^2*K)", # (ranges between 2000 and 8000)
    ),
)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end
struct FuelGrain <: AbstractPhysics end
struct IgniterPuck <: AbstractPhysics end

push!(config.data.tstops, config.optimized_parameters.valve_open_time)
push!(config.data.tstops, config.optimized_parameters.igniter_burn_time)
push!(config.data.tstops, config.optimized_parameters.cea_switch_time)
push!(config.data.tstops, config.optimized_parameters.motor_burn_time)

Revise.includet(joinpath(@__DIR__, "property_update_functions", "system_moles_vec.jl"))
Revise.includet(joinpath(@__DIR__, "property_update_functions", "pre_cea_update_fluid_caches.jl"))
Revise.includet(joinpath(@__DIR__, "property_update_functions", "update_fluid_caches.jl"))
Revise.includet(joinpath(@__DIR__, "physics", "heat_transfer_physics", "heat_transfer_physics.jl"))
Revise.includet(joinpath(@__DIR__, "physics", "injector_delivery_train_physics", "injector_delivery_train_physics.jl"))

#Don't need this right now since k and mu are fixed
function manual_cache_updates!(du, u, p, t, system, geo, cell_id)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.mu[cell_id] = properties.mu[cell_id]
end

function cap_heat_flux_to_temp_change!(du, u, p, t, system, geo, cell_id)
    du.state_temperature[cell_id] += du.heat[cell_id] / (geo.cell_volumes[cell_id] * u.density[cell_id] * u.cv[cell_id])
end

function update_fluid_properties_non_cea!(du, u, p, t, system, geo, cell_id)
    u.temperature[cell_id] = u.state_temperature[cell_id]

    manual_cache_updates!(du, u, p, t, system, geo, cell_id)

    properties = ComponentVector(system.properties_vec, system.properties_axes)

    for_fields!(u.mass_fractions, properties.mass_fractions) do species, u_mass_fractions, properties_mass_fractions
        u.mass_fractions[species[cell_id]] = properties_mass_fractions[species[cell_id]]
    end

    pre_cea_fluid_property_update_ideal!(du, u, p, t, system, geo, cell_id)
end

function update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    #manual_fluid_cache_updates!(du, u, p, t, system, geo, cell_id)
    if t < p.cea_switch_time || system.hot_fire == false || t > p.motor_burn_time
        update_fluid_properties_non_cea!(du, u, p, t, system, geo, cell_id)
    else
        cea_fluid_property_update!(du, u, p, t, system, geo, cell_id)
        #cea_fluid_property_update_clapeyron!(du, u, p, t, system, geo, cell_id)
    end
    #fluid_property_update_clapeyron!(du, u, p, t, system, geo, cell_id)
end

function fluid_domain_physics!(du, u, p, t, system, geo, cell_id)
    if t < p.cea_switch_time || system.hot_fire == false || t > p.motor_burn_time
        cap_heat_flux_to_temp_change!(du, u, p, t, system, geo, cell_id)
    end
end

function normalize!(vec)
    total = sum(vec)
    vec ./= total
end

#Since we're defining all connections manually, we can just have no_flux! for every region
function no_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    #nothing here since we handle these BCs at the patches
end

function connection_map_function(phys_a, phys_b)
    return no_flux!
end

room_temp = 294.13u"K"

air_initial_conditions = ComponentVector(
    density = 1.225u"kg/m^3",
    state_temperature = room_temp,
)

air_properties = ComponentVector(
    #We're assuming these two stay constant which isn't that accurate, but oh well
    k = 0.024u"W/(m*K)", 
    mu = 1.9e-5u"Pa*s", 

    mass_fractions = normalize!(ComponentVector(
        nitrous_oxide = 0.0,
        ethylene = 0.0,
        nitrogen = 0.755,
        oxygen = 0.232,
        carbon_dioxide = 0.0003,
        carbon_monoxide = 0.0,
        water = 0.0,
        hydrogen = 0.0
    )),

    molecular_weights = (
        nitrous_oxide = 44.01u"g/mol", 
        ethylene = 28.05u"g/mol", 
        nitrogen = 28.01u"g/mol", 
        oxygen = 32.00u"g/mol", 
        carbon_dioxide = 44.01u"g/mol",
        carbon_monoxide = 28.01u"g/mol",
        water = 18.02u"g/mol",
        hydrogen = 2.02u"g/mol"
    ),

    R_gas = 8.314u"J/(mol*K)"
)

#Uses the unitful cad_data: config.data.cad_data has its units stripped, which would turn this mass into a density
initial_igniter_puck_mass = config.optimized_parameters.igniter_density * cad_data.igniter_puck_volume |> u"kg"
config.u_proto.igniter_mass = [initial_igniter_puck_mass]
config.data = merge(config.data, (
    initial_igniter_mass = ustrip(upreferred(initial_igniter_puck_mass)), #stripped like the rest of system data, since the RHS works in plain SI numbers
));

Revise.includet(joinpath(@__DIR__, "physics", "igniter_physics", "igniter_physics.jl"))

function motor_burned_out(du, u, p, t, system, geo)
    return t > p.motor_burn_time
end

add_region!(
    config, "pre_chamber";
    type = Fluid(),
    initial_conditions = air_initial_conditions,
    properties = air_properties,
    property_update_function = 
    function update_pre_chamber_region!(du, u, p, t, system, geo, cell_id)
        update_fluid_properties_non_cea!(du, u, p, t, system, geo, cell_id)
        #update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function pre_chamber_region!(du, u, p, t, system, geo, cell_id)
        #The igniter adds its heat to du.heat, so it has to run before cap_heat_flux_to_temp_change! converts du.heat to temperature
        if system.burn_igniter == true
            burn_igniter_puck!(du, u, p, t, system, geo, cell_id)
        end

        #this isn't very rigorous, we prevent this cell from heating up or cooling down for most of the motor firing
        if t < p.cea_switch_time || system.hot_fire == false || t > p.motor_burn_time 
            cap_heat_flux_to_temp_change!(du, u, p, t, system, geo, cell_id)
        end

        if motor_burned_out(du, u, p, t, system, geo)
            du.density[cell_id] *= 0.0
            du.state_temperature[cell_id] *= 0.0
        end
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

        if motor_burned_out(du, u, p, t, system, geo)
            du.density[cell_id] *= 0.0
            du.state_temperature[cell_id] *= 0.0
        end
    end
)

#=
N2O properties to keep in mind:

vapor_k = 0.017u"W/(m*K)"
vapor_mu = 1.4e-5u"Pa*s"

liquid_k = 0.09u"W/(m*K)"
liquid_mu_lower_bound = 6.0e-5u"Pa*s"
liquid_mu_upper_bound = 3.2e-4u"Pa*s"
=#

#we'll use the vapor properties for now since a large majority of the calculations that involve these values happen in the vapor phase
nitrous_oxide_properties = merge_properties(
    air_properties,
    ComponentVector(
        k = 0.017u"W/(m*K)", 
        mu = 1.4e-5u"Pa*s", 
        mass_fractions = (
            nitrous_oxide = 1.0u"kg/kg",
            ethylene = 0.0u"kg/kg",
            nitrogen = 0.0u"kg/kg",
            oxygen = 0.0u"kg/kg",
            carbon_dioxide = 0.0u"kg/kg",
            carbon_monoxide = 0.0u"kg/kg",
            water = 0.0u"kg/kg",
            hydrogen = 0.0u"kg/kg"
        )
    )
)

#Since we're using a series-valve assumption, the mid_section is no longer needed
#we might use it in the future to allow for heat to transfer to the nitrous oxide from the injector bulkhead
add_region!(
    config, "mid_section";
    type = Fluid(),
    initial_conditions = ComponentVector(
        density = 1.0u"kg/m^3",
        state_temperature = room_temp
    ),
    properties = ComponentVector(),
    property_update_function = 
    function update_mid_section!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function mid_section_region!(du, u, p, t, system, geo, cell_id)
    end
)

add_region!(
    config, "tank";
    type = Fluid(),
    initial_conditions = ComponentVector(
        density = (1.134u"kg" / geo.cell_volumes[get_cellset(geo, "tank")[1]]),
        state_temperature = room_temp,
    ),
    properties = nitrous_oxide_properties,
    property_update_function = 
    function update_tank_region!(du, u, p, t, system, geo, cell_id)
        properties = ComponentVector(system.properties_vec, system.properties_axes)

        for_fields!(u.mass_fractions, properties.mass_fractions) do species, u_mass_fractions, properties_mass_fractions
            u.mass_fractions[species[cell_id]] = properties_mass_fractions[species[cell_id]]
        end

        u.temperature[cell_id] = u.state_temperature[cell_id]

        update_tank_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function tank_region!(du, u, p, t, system, geo, cell_id)
        cap_heat_flux_to_temp_change!(du, u, p, t, system, geo, cell_id)
    end
)

function update_solid_properties!(du, u, p, t, system, geo, cell_id)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    manual_cache_updates!(du, u, p, t, system, geo, cell_id)

    #u.k[cell_id] = properties.k[cell_id] #not needed right now, but will be in the future
    u.cp[cell_id] = properties.cp[cell_id]

    #Solids barely expand, so cp and cv are effectively equal
    u.cv[cell_id] = properties.cp[cell_id]

    u.temperature[cell_id] = u.state_temperature[cell_id]
end

function cap_heat_flux_to_solid_temp_change!(du, u, p, t, system, geo, cell_id)
    du.state_temperature[cell_id] += du.heat[cell_id] / (geo.cell_volumes[cell_id] * u.density[cell_id] * u.cv[cell_id])
end

function solid_physics!(du, u, p, t, system, geo, cell_id)
    cap_heat_flux_to_solid_temp_change!(du, u, p, t, system, geo, cell_id)
end

aluminum_initial_conditions = ComponentVector(
    density = 2700.0u"kg/m^3",
    state_temperature = room_temp
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
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_region!(
    config, "injector_bulkhead_bottom";
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_region!(
    config, "aluminum_housing_top";
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_region!(
    config, "aluminum_housing",
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

phenolic_liner_initial_conditions = ComponentVector(
    density = 1250.0u"kg/m^3",
    state_temperature = room_temp
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
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_region!(
    config, "phenolic_liner",
    type = Solid(),
    initial_conditions = phenolic_liner_initial_conditions,
    properties = phenolic_liner_properties,
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

fuel_grain_initial_conditions = ComponentVector(
    density = 950.0u"kg/m^3",
    state_temperature = room_temp,
)

fuel_grain_properties = ComponentVector(
    cp = 2000.0u"J/(kg*K)",
    k = 0.45u"W/(m*K)"
)

add_region!(
    config, "fuel_grain",
    type = FuelGrain(),
    initial_conditions = fuel_grain_initial_conditions,
    properties = fuel_grain_properties,
    property_update_function = 
    function update_fuel_grain_region!(du, u, p, t, system, geo, cell_id)
        update_solid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function fuel_grain_region!(du, u, p, t, system, geo, cell_id)
        solid_physics!(du, u, p, t, system, geo, cell_id)
    end

)

graphite_initial_conditions = ComponentVector(
    density = 1760.0u"kg/m^3",
    state_temperature = room_temp
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
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_region!(
    config, "aluminum_housing_near_nozzle_top",
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_region!(
    config, "nozzle_bottom",
    type = Solid(),
    initial_conditions = graphite_initial_conditions,
    properties = graphite_properties,
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_region!(
    config, "aluminum_housing_near_nozzle_bottom",
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_region!(
    config, "nozzle_retainer",
    type = Solid(),
    initial_conditions = aluminum_initial_conditions,
    properties = aluminum_properties,
    property_update_function = update_solid_properties!,
    region_function = solid_physics!
)

add_patch!(config;
    owner = "tank",
    connection = "mid_section",
    patch_function = function tank_to_mid_section!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
    end
)

#the mid section is disabled for now 
#in the future, we might want to consider the heat transfered to the nitrous oxide as it travels through the injector bulkhead
add_patch!(config;
    owner = "mid_section",
    connection = "injector_bulkhead_top",
    patch_function = function mid_section_to_injector_bulkhead_top!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
        #might need a CFD-informed closure due to the convective heat transfer
        #mid_section_to_injector_bulkhead_top_closure_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) 
    end
)

add_patch!(config;
    owner = "mid_section",
    connection = "injector_bulkhead_bottom",
    patch_function = function mid_section_to_injector_bulkhead_bottom!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
    end
)

add_patch!(config;
    owner = "injector_bulkhead_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function injector_bulkhead_top_to_injector_bulkhead_bottom!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "tank",
    connection = "pre_chamber",
    patch_function = function mid_section_to_pre_chamber!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        if t > p.valve_open_time && t < p.motor_burn_time
            injector_train_mass_flow!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) 
        end
        #this would be the valve equation determining flow between the mid section and first chamber cell
        #purely algebraic
    end
)

#Axial fluxes from top cells to injector_bulkhead_bottom:
#start
add_patch!(config;
    owner = "pre_chamber",
    connection = "injector_bulkhead_bottom",
    patch_function = function pre_chamber_to_injector_bulkhead_bottom!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        #TODO: WE CANNOT IGNORE CONVECTION (IT WILL LIKELY BE OFF BY A FACTOR OF 1000x)!
        if motor_burned_out(du, u, p, t, system, geo) == false
            heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) #this is pretty negligible
            chamber_convective_heat_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
        end
    end
)

add_patch!(config;
    owner = "phenolic_liner_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function phenolic_liner_top_to_injector_bulkhead_bottom!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "aluminum_housing_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function aluminum_housing_top_to_injector_bulkhead_bottom!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "aluminum_housing_top",
    connection = "aluminum_housing",
    patch_function = function aluminum_housing_top_to_aluminum_housing!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)
#end

#Radial fluxes across top cells

Revise.includet(joinpath(@__DIR__, "physics", "chamber_convective_heat_transfer", "chamber_convective_heat_transfer.jl"))

#start
add_patch!(config;
    owner = "pre_chamber",
    connection = "phenolic_liner_top",
    patch_function = function pre_chamber_to_phenolic_liner_top!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        #this one will likely need a CFD closure in the future to account for convection
        #TODO: WE CANNOT IGNORE CONVECTION (IT WILL LIKELY BE OFF BY A FACTOR OF 1000x)!
        if motor_burned_out(du, u, p, t, system, geo) == false
            heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
            chamber_convective_heat_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
        end
    end
)

add_patch!(config;
    owner = "phenolic_liner_top",
    connection = "aluminum_housing_top",
    patch_function = function phenolic_liner_top_to_aluminum_housing_top!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

Revise.includet(joinpath(@__DIR__, "physics", "pressure_driven_mass_flux", "pressure_driven_mass_flux.jl"))

add_patch!(config;
    owner = "pre_chamber",
    connection = "fluid_domain", #1-n cells,
    patch_function = function pre_chamber_to_fluid_domain!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b
    )
        pressure_driven_mass_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a, idx_b, face_b
        )
    end
)

#end

#We could add a flux between the top fluid domain and the fuel grain, but I'm going to ignore it

#Axial fluxes along n cells
#start

#these are called from solve_system!
#The stations and the faces between them come from build_fuel_port_geometry, ordered by axial position,
#so each axial face connects a lower station to the station directly above it.

function get_axial_flux_locations(geo, axial_region_name)
    axial_flux_locations = NTuple{4, Int}[]

    axial_cells = get_cellset(geo, axial_region_name)

    for cell_id in axial_cells
        for (neighbor_cell_id, local_face, neighbor_local_face) in geo.cell_neighbors[cell_id][2]
            if neighbor_cell_id in axial_cells && cell_id < neighbor_cell_id
                push!(axial_flux_locations, (
                    cell_id, local_face, neighbor_cell_id, neighbor_local_face
                ))
            end
        end
    end

    return axial_flux_locations
end

config.data = merge(
    config.data,
    (
        fluid_domain_axial_flux_locations = get_axial_flux_locations(config.geo, "fluid_domain"),
        fuel_grain_axial_flux_locations = get_axial_flux_locations(config.geo, "fuel_grain"),
        phenolic_liner_axial_flux_locations = get_axial_flux_locations(config.geo, "phenolic_liner"),
        aluminum_housing_axial_flux_locations = get_axial_flux_locations(config.geo, "aluminum_housing")
    )
);

function fluid_domain_axial_flux!(du, u, p, t, system, geo)
    for (idx_a, face_a, idx_b, face_b) in system.fluid_domain_axial_flux_locations
        pressure_driven_mass_flux!(
            du, u, p, t, system, geo,
            idx_a, face_a, idx_b, face_b
        )
    end
end

function fuel_grain_axial_flux!(du, u, p, t, system, geo)
    for (idx_a, face_a, idx_b, face_b) in system.fuel_grain_axial_flux_locations
        heat_diffusion_series_resistance!(
            du, u, p, t, system, geo,
            idx_a, face_a, idx_b, face_b
        )
    end
end

#The liner and housing don't move, so their axial conduction uses the CAD geometry directly.
#Their cell sets are ordered by axial position, so each cell's upper neighbour is the next cell in the set.
function phenolic_liner_axial_flux!(du, u, p, t, system, geo)
    for (idx_a, face_a, idx_b, face_b) in system.phenolic_liner_axial_flux_locations
        heat_diffusion_series_resistance!(
            du, u, p, t, system, geo,
            idx_a, face_a, idx_b, face_b
        )
    end
end

function aluminum_housing_axial_flux!(du, u, p, t, system, geo)
    for (idx_a, face_a, idx_b, face_b) in system.aluminum_housing_axial_flux_locations
        heat_diffusion_series_resistance!(
            du, u, p, t, system, geo,
            idx_a, face_a, idx_b, face_b
        )
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
add_patch!(config;
    owner = "fluid_domain",
    connection = "fuel_grain",
    patch_function = function fluid_domain_to_fuel_grain!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        #we're going to assume that once the gas temperature next to the fuel grain goes above its regression temperature, it stops heating up further
        #since (like water boiling into steam) more heating just results in faster regression and doesn't meaningfully increase the temperature of the fuel grain

        if u.temperature[idx_a] <= p.min_fuel_regression_temperature && motor_burned_out(du, u, p, t, system, geo) == false
            heat_diffusion_series_resistance!(
                du, u, p, t, system, geo,
                idx_a, face_a, idx_b, face_b
            )
            chamber_convective_heat_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
        end
    end
)

add_patch!(config;
    owner = "fuel_grain",
    connection = "phenolic_liner",
    patch_function = function fuel_grain_to_phenolic_liner!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "phenolic_liner",
    connection = "aluminum_housing",
    patch_function = function phenolic_liner_to_aluminum_housing!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        #this one, unlike the others, does not require geometry updating
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)
#end

Revise.includet(joinpath(@__DIR__, "physics", "nozzle_outlet", "nozzle_outlet_closure.jl"))
Revise.includet(joinpath(@__DIR__, "physics", "nozzle_outlet", "nozzle_outlet_simple.jl"))

add_patch!(config, "fluid_domain_to_nozzle_inlet";
    patch_function = function fluid_domain_to_nozzle_inlet!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        #idx_b, and face_b are both 0 here
        nozzle_outlet_simple!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
        #heat_transfer_from_fluid_to_nozzle_top!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
        #geometry information here is not accurate because this is a virtual boundary condition cell defined in CAD, use the CAD geometry data instead

        fluid_cell_id = idx_a
        converging_cell_id = get_cellset(geo, "nozzle_top")[1]
        diverging_cell_id = get_cellset(geo, "nozzle_bottom")[1]

        converging_area = system.converging_nozzle_surface_area
        diverging_area = system.diverging_nozzle_surface_area

        if motor_burned_out(du, u, p, t, system, geo) == false
            du.heat[fluid_cell_id] -= p.nozzle_converging_heat_transfer_coefficient * converging_area * (u.temperature[fluid_cell_id] - u.temperature[converging_cell_id])
            du.heat[converging_cell_id] += p.nozzle_converging_heat_transfer_coefficient * converging_area * (u.temperature[fluid_cell_id] - u.temperature[converging_cell_id])

            du.heat[fluid_cell_id] -= p.nozzle_diverging_heat_transfer_coefficient * diverging_area * (u.temperature[fluid_cell_id] - u.temperature[diverging_cell_id])
            du.heat[diverging_cell_id] += p.nozzle_diverging_heat_transfer_coefficient * diverging_area * (u.temperature[fluid_cell_id] - u.temperature[diverging_cell_id])
        end
    end
)

#we could also add a flux between the nozzle top and phenolic liner but I don't think it's necessary

#This is a radial flux between the nozzle graphite and aluminum housing
add_patch!(config;
    owner = "nozzle_top",
    connection = "aluminum_housing_near_nozzle_top",
    patch_function = function nozzle_top_to_aluminum_housing_near_nozzle_top!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "nozzle_top",
    connection = "nozzle_bottom",
    patch_function = function nozzle_top_to_nozzle_bottom!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "nozzle_bottom",
    connection = "aluminum_housing_near_nozzle_bottom",
    patch_function = function nozzle_bottom_to_aluminum_housing_near_nozzle_bottom!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "nozzle_bottom",
    connection = "nozzle_retainer",
    patch_function = function nozzle_bottom_to_nozzle_retainer!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "aluminum_housing_near_nozzle_bottom",
    connection = "nozzle_retainer",
    patch_function = function aluminum_housing_near_nozzle_bottom_to_nozzle_retainer!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "aluminum_housing_near_nozzle_top",
    connection = "aluminum_housing_near_nozzle_bottom",
    patch_function = function aluminum_housing_near_nozzle_top_to_aluminum_housing_near_nozzle_bottom!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

add_patch!(config;
    owner = "aluminum_housing",
    connection = "aluminum_housing_near_nozzle_top",
    patch_function = function aluminum_housing_to_aluminum_housing_near_nozzle_top!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b
    )
        heat_diffusion_series_resistance!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end
)

n_species = length(propertynames(config.cache_syms_and_units.mass_fractions))

#Same chunk size finish_fvm_config uses for the du/u DiffCaches
jacobian_chunk_size = ForwardDiff.pickchunksize(length(Vector(config.u_proto)))

config.data = merge(config.data, (
    n_species = n_species,
    moles_vec = CellSpeciesVectors(config.geo.n_cells, n_species, jacobian_chunk_size),
    mole_fractions_vec = CellSpeciesVectors(config.geo.n_cells, n_species, jacobian_chunk_size),
));

function solve_system!(du, u, p, t, system, geo)
    #Only does anything during sparsity detection: the tracer get_tmp returns a new #undef array every call, so one
    #buffer is stored here for the whole RHS call. Float64 and Dual calls look up their persistent buffer at each use.
    #Same number type as the per-cell samples (u.density[cell_id] + first(p)) used when reading and writing the buffers.
    select_cell_species_buffer!(system.moles_vec, u.density[1] + first(p))
    select_cell_species_buffer!(system.mole_fractions_vec, u.density[1] + first(p))

    #Moles and mole fractions are already updated per cell inside update_region_groups!
    #(fluid_property_update_ideal! / fluid_property_update_clapeyron!), after mw_avg!
    update_region_groups!(du, u, p, t, system, geo)

    axial_flux_functions!(du, u, p, t, system, geo)

    #update_MUSCL_gradients!(u, wls_stencil)
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)

    solve_region_groups!(du, u, p, t, system, geo)
end

#The leading semicolon keeps the field names; without it, splatting a NamedTuple gives a plain Tuple of values
system_data = (;
    config.data...,
);

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, system_data, check_units = false);

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

p_guess = ustrip.(upreferred.(deepcopy(config.optimized_parameters)))

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

t0 = 0.0
tMax = 10.0

#Dense Jacobian for now. jac_sparsity above comes from LOCAL sparsity detection at t = 0, which only records the
#couplings active at that exact state: the valve is closed (injector flow = 0 * everything) and the nozzle is off
#(chamber at ambient). Once the valve opens, the flowing state has 5525 Jacobian nonzeros and 3406 of them are
#missing from the t = 0 pattern, so the sparse Jacobian was wrong and Newton failed on almost every step
#(~1e-7 s steps). With the dense Jacobian the cold flow runs 0.1 s -> 0.2 s in ~16 s and 168 steps.
#To go back to sparse, the pattern has to cover every regime: e.g. the union of local patterns at t = 0 and at
#states with the valve open, the nozzle flowing, and (for hot fires) the regression and combustion active.
ode_func = ODEFunction(f_closure_implicit)

#ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity)) #old version, sparsity pattern from t = 0 only
startup_prob = ODEProblem(ode_func, u0_vec, (t0, ustrip(upreferred(config.optimized_parameters.motor_burn_time))), p_guess)

Revise.includet(joinpath(@__DIR__, "callbacks", "stop_sol_callbacks.jl"))

callbacks = CallbackSet(
    approximate_time_to_finish_cb,
    min_chamber_pressure_termination_cb
)

Revise.includet(joinpath(@__DIR__, "state_is_invalid", "state_is_invalid.jl"))

state_is_invalid_closure = (u, p, t) -> state_is_invalid(u, p, t, system);

#Implicit from the start. The laminar port flux equalizes pressure between neighbouring port cells on a ~1 ns
#timescale (Jacobian entries ~1e11), so an explicit method like Tsit5 needs nanosecond steps and blows up on its
#first trial step. FBDF steps over that stiffness. The small first dt keeps the initial step from being sized
#by the hundreds of states that barely move.
@time startup_sol = solve(
    startup_prob,
    #FBDF(linsolve = SparspakFactorization(), autodiff = ADTypes.AutoForwardDiff()),
    #dt = 1e-6,
    #Tsit5(), #old version, explicit
    FBDF(autodiff = ADTypes.AutoForwardDiff()),
    dt = 1e-6,
    callback = callbacks,
    isoutofdomain = state_is_invalid_closure,
    tstops = ustrip.(upreferred.(config.data.tstops))
)

function update_fluid_domain_state_temperatures(sol, startup_u_end)
    u0 = deepcopy(startup_u_end)

    du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, config.geo, p_guess, save_properties_and_caches = true, track_progress = false);
    
    fluid_domain_cell_ids = get_cellset(geo, "fluid_domain")

    u0_named = ComponentVector(u0, system.state_axes)

    for cell_id in fluid_domain_cell_ids
        u0_named.state_temperature[cell_id] = u_named[end].temperature[cell_id]
    end

    return Vector(u0_named)
end

post_motor_burnout_prob = remake(
    startup_prob; 
    u0 = update_fluid_domain_state_temperatures(startup_sol, startup_sol.u[end]), 
    tspan = (ustrip(upreferred(config.optimized_parameters.motor_burn_time)), 50.0)
)

@time post_motor_burnout = solve(
    post_motor_burnout_prob,
    FBDF(autodiff = ADTypes.AutoForwardDiff()),
    dt = 1e-6,
    callback = callbacks,
    isoutofdomain = state_is_invalid_closure,
    tstops = ustrip.(upreferred.(config.data.tstops))
)

sol = (u = [startup_sol.u..., post_motor_burnout.u...], t = [startup_sol.t..., post_motor_burnout.t...])

desired_shown_steps = 250
save_fraction = desired_shown_steps / length(sol.t)

reduced_sol = (
    u = [[sol.u[i] for i in 1:(Int(ceil(1/save_fraction))):length(sol.u)]..., sol.u[end]],
    t = [sol.t[1:(Int(ceil(1/save_fraction))):length(sol.t)]..., sol.t[end]]
)

du_named, u_named = regenerate_fvm_state(reduced_sol, system, solve_system!, config.geo, p_guess, save_properties_and_caches = true, track_progress = true);

using Plots

pre_chamber_cell_id = get_cellset(config.geo, "pre_chamber")[1]
plot(reduced_sol.t, [u_named[i].temperature[pre_chamber_cell_id] for i in eachindex(reduced_sol.t)])
#plot(reduced_sol.t .* u"s", [u_named[i].neglected_injector_cooling_power[pre_chamber_cell_id]u"W" |> u"kW" for i in eachindex(reduced_sol.t)])

#plot(reduced_sol.t, [u_named[i].igniter_mass[1] for i in eachindex(reduced_sol.t)])

plot(
    reduced_sol.t, 
    [u_named[i].temperature[get_cellset(config.geo, "injector_bulkhead_top")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], 
    label = "injector_bulkhead_top", 
    title = "Temperaxtures of different regions over time in °C", 
    #legend = :bottomright
    legendfontsize = 7, 
    xticks = (
        [0:5:50; ustrip(upreferred(config.optimized_parameters.motor_burn_time))],
        [string.(collect(0:5:50)); ]
    )
)
plot!(reduced_sol.t, [u_named[i].temperature[get_cellset(config.geo, "injector_bulkhead_bottom")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], label = "injector_bulkhead_bottom")
plot!(reduced_sol.t, [u_named[i].temperature[get_cellset(config.geo, "aluminum_housing_top")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], label = "aluminum_housing_top")
plot!(reduced_sol.t, [u_named[i].temperature[get_cellset(config.geo, "phenolic_liner_top")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], label = "phenolic_liner_top")
plot!(reduced_sol.t, [u_named[i].temperature[get_cellset(config.geo, "aluminum_housing_near_nozzle_top")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], label = "aluminum_housing_converging")
plot!(reduced_sol.t, [u_named[i].temperature[get_cellset(config.geo, "aluminum_housing_near_nozzle_bottom")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], label = "aluminum_housing_diverging")
plot!(reduced_sol.t, [u_named[i].temperature[get_cellset(config.geo, "nozzle_top")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], label = "nozzle_top")
plot!(reduced_sol.t, [u_named[i].temperature[get_cellset(config.geo, "nozzle_bottom")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], label = "nozzle_bottom")
plot!(reduced_sol.t, [u_named[i].temperature[get_cellset(config.geo, "nozzle_retainer")[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)], label = "nozzle_retainer")

function plot_axial_region_temperatures(region_name)
    region_cell_ids = get_cellset(geo, region_name)

    start_idx = 1
    stop_idx = length(reduced_sol.t)

    plt = plot(reduced_sol.t[start_idx:stop_idx], [u_named[i].temperature[region_cell_ids[1]] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)][start_idx:stop_idx])

    for region_cell_id in region_cell_ids[2:end]
        timesteps = reduced_sol.t[start_idx:stop_idx]
        temperature_data = [u_named[i].temperature[region_cell_id] * u"K" |> u"°C" for i in eachindex(reduced_sol.t)][start_idx:stop_idx]

        plt = plot!(timesteps, temperature_data)

        sleep(0.3)

        display(plt)
    end

    display(plt)
end

#pressure, oxidizer/fuel ratio, expansion ratio
chamber_T_comb_interp(ustrip(10u"bar" |> u"Pa"), 7.0, 10.0)
chamber_T_comb_interp(ustrip(20u"bar" |> u"Pa"), 7.0, 10.0)
chamber_T_comb_interp(ustrip(30u"bar" |> u"Pa"), 7.0, 10.0)

chamber_T_comb_interp(ustrip(20u"bar" |> u"Pa"), 1.0, 10.0)
chamber_T_comb_interp(ustrip(20u"bar" |> u"Pa"), 5.0, 10.0)
chamber_T_comb_interp(ustrip(20u"bar" |> u"Pa"), 13.0, 10.0)
chamber_T_comb_interp(ustrip(20u"bar" |> u"Pa"), 15.0, 10.0)

#plot_axial_region_temperatures("fluid_domain")
plot_axial_region_temperatures("fuel_grain")
plot_axial_region_temperatures("phenolic_liner")
plot_axial_region_temperatures("aluminum_housing")

u_named[end].state_temperature[get_cellset(geo, "aluminum_housing")[end]]