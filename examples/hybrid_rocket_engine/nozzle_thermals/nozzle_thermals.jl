using Revise

using FVMFramework

using Logging
using Unitful
#using OrdinaryDiffEq
using NonlinearSolve
#using Sparspak
#using Ferrite
#using FerriteGmsh
#using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
#using ILUZero
using LinearAlgebra
#using StaticArrays
using Clapeyron
using PreallocationTools
using DataInterpolations
using JSON

using OrderedCollections

using CoolProp

using Roots

cad_table_data_json_path = joinpath(@__DIR__, "rom_and_table_inputs", "table_data", "freecad_table_data.json")
Revise.includet(joinpath(@__DIR__, "rom_and_table_inputs", "table_data", "read_freecad_table_data.jl"))
cad_data = ComponentVector(read_freecad_table_data_json(cad_table_data_json_path))

# Variables defined based on cad that are more or less fixed:
cad_data.aluminum_casing_outer_diameter
cad_data.flange_overhang_axial_length

cad_data.nozzle_inlet_diameter
cad_data.nozzle_throat_diameter
cad_data.nozzle_exit_diameter

#=
Variables that are going to be continuously adjusted along the length of the nozzle

    Main design parameters:
        - The wall thickness between cooling channels and the internal flow passage
            - On second thought, we probably want to go as thin as possible that gives some safety factor for both pressure and thermal stresses 
        - Effective diameter of cooling channels (we can increase heat flux by narrowing the channels)
        - Cooling channel width, height and aspect ratio
        - More unique forms of cooling channel geometry (perhaps internal fins if there's an easy way to predict their impact)
        - Pitch of cooling channels / coils of cooling channels per unit length
        - Orifice diameter upstream of the cooling channels (controls pressure experienced in the pipe based on pump volumetric flow)
            - Higher pressure leads to a higher critical heat flux and also a higher saturation temperature
            - 

    Secondary design parameters:
        - The thickness between the cooling channels and the outside air (only really have to size to contain coolant pressure)
        - coolant flow rate
        - coolant pressure
        - coolant inlet temperature
=#

#=
Assumptions and safety factors:
    - I believe that designing the water cooling to deal with twice the expected heat flux into the nozzle would be safe
    - We're going to use a steady state model because the aluminum heat capacity is negligible 
=#

Revise.includet(joinpath(@__DIR__, "nozzle_stations_creator.jl"))

nozzle_stations = create_nozzle_stations(cad_data, n_converging = 100, n_diverging = 100)

Revise.includet(joinpath(@__DIR__, "cea_retrieval", "CEA_lookup_table.jl"));

system = (
    chamber_cp_interp = chamber_cp_interp,
    chamber_mu_interp = chamber_mu_interp,
    chamber_k_interp = chamber_k_interp,
    chamber_prandtl_interp = chamber_prandtl_interp,
    chamber_isp_vac_interp = chamber_isp_vac_interp,
    chamber_cstar_interp = chamber_cstar_interp,
    chamber_T_comb_interp = chamber_T_comb_interp,
    chamber_mw_avg_interp = chamber_mw_avg_interp,
    chamber_gamma_interp = chamber_gamma_interp,

    mass_fractions_interps_vec = mass_fractions_interps_vec,

    cad_data = cad_data,
    nozzle_throat_position = get_nozzle_throat_position(cad_data), #function defined in nozzle_stations_creator.jl

    coolprop_fluid = "water",

    # Model Modes:
    use_frozen_cea_interps = true, #equilibrium overpredicts way too much, we're always going to use frozen cea interps from now on
    forced_convection_mode = :gnielinski_norris, #options are :gnielinski_norris and :gnielinski_rough_pipe
    #thermal_hydraulic_diameter_model = :optimistic, #options are :optimistic or :pessimistic # We're no longer using :pessimistic
    fully_developed_boiling_heat_flux_mode = :jens_lottes, #options are :jens_lottes or :thom
    use_curved_rough_friction_factor = true, 
    use_curved_smooth_friction_factor = true,

    outlet_orifice_area_or_inlet_pressure_adjusted = :outlet_orifice_area, #options are :outlet_orifice_area and :coolant_inlet_pressure
    #:outlet_orifice_area: 
        # outlet_orifice_area is chosen based on the desired inlet pressure
    #:coolant_inlet_pressure:
        # outlet_orifice_area is kept fixed, and the inlet pressure is allowed to vary (uses Roots.jl)
);

#we'll keep application as a parameter only used for clarification where it's needed
#for example, for valve cv, you would say application = :scalar or application = :multiplier
uncertainties = (
    #Since we're uncertain on the rate of recombination in the nozzle, we're going to multiply the bartz correlation results by this value in a parameter sweep
    #Eh, we'll just use frozen cea interps
    #=
    recombination_factor = (
        low = 1.0,
        high = 1.37,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    =#
    #The bartz heat flux correlation is known to overpredict the heat transfer coefficient at the throat
    bart_heat_flux_correlation_uncertainty = (
        low = 0.6,
        high = 1.3,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    radiative_heat_flux_uncertainty = (
        low = 0.5,
        high = 1.5,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    #Represents the uncertainty in the *pressure stress* that the model predicts, not wall thickness
    hot_wall_pressure_stress_uncertainty = (
        low = 0.85,
        high = 1.3,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    #Represents the uncertainty in the *thermal stress* that the model predicts, not wall thickness
    hot_wall_thermal_stress_uncertainty = (
        low = 0.9,
        high = 1.25,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    #When the coolant travels around the nozzle in the helical coolant channels, 
    #I'm pretty sure that the centrifugal force that leads to a lower pressure on the 
    #surface facing into the helix will lead to a lower heat tnrasfer coefficient there 
    #than what would typically be observed.
    #NOTE: there may actually be a way to resolve this by estimating the pressure difference 
    #between the bulk fluid pressure and the closer and further walls and modifying the equations
    #to use that modified pressure for their calculations. For now, we'll just guess.
    coolant_channels_heat_transfer_curvature_modifier = (
        low = 0.6,
        high = 1.0,
        nominal = 0.8,
        type = :model,
        worst_case = :low
    ),
    jens_lottes_fully_developed_boiling_heat_flux_uncertainty = (
        low = 0.7,
        high = 1.3,
        nominal = 1.0,
        type = :model,
        worst_case = :low
    ),
    thom_fully_developed_boiling_heat_flux_uncertainty = (
        low = 0.7,
        high = 1.3,
        nominal = 1.0,
        type = :model,
        worst_case = :low
    ),
    bergles_rohsenow_onset_of_nucleate_boiling_heat_flux_uncertainty = (
        low = 0.6,
        high = 1.4,
        nominal = 1.0,
        type = :model,
        worst_case = :low
    ),

    #Main rocket design parameters:
    #Chamber pressure could definitely go higher or lower than expected
    chamber_pressure = (
        low = 19.0u"bar",
        high = 23.0u"bar",
        nominal = 20.0u"bar",
        type = :parameter,
        worst_case = :high
    ),
    #higher O/F ratios drives chamber temperatures higher, we'll use a higher than normal value for safety
    oxidizer_to_fuel_ratio = (
        low = 6.0,
        high = 10.0,
        nominal = 7.8,
        type = :parameter,
        worst_case = :high
    ),
    #I'm expecting a lower combustion efficiency than most hybrids due to us not having a post-chamber section
    #high is worst case because more complete combustion leads to higher chamber temperatures
    combustion_efficiency = (
        low = 0.8,
        high = 0.9,
        nominal = 0.85,
        type = :parameter,
        worst_case = :high
    ),
    #Not very sure about this one:
    #Bulk 6061 aluminum (non-SLM) has an emissivity of around 0.07-0.09,
    #but SLM aluminum is reported to be somewhere around 0.25-0.45, I'll trust the SLM values for now, but its something to keep in mind
    #Since radiative heat flux is probably going to be small compared to the convective heat flux, this probably won't matter that much
    wall_emissivity = (
        low = 0.25,
        high = 0.45,
        nominal = 0.35,
        type = :parameter,
        worst_case = :high
    ),
    #There are correlations to estimate this one better, but I'm feeling lazy right now so we're just going to sweep it :P
    #typically around 0.05 to 0.20 for hybrid combustion gases
    gas_emissivity = (
        low = 0.05,
        high = 0.20,
        nominal = 0.10,
        type = :parameter,
        worst_case = :high
    ),
    wall_coefficient_of_thermal_expansion = (
        low = 20.0e-6u"1/K",
        high = 24.0e-6u"1/K",
        nominal = 22.0e-6u"1/K",
        type = :parameter,
        worst_case = :high
    ),
    # Parameter Modifiers:
    wall_equivalent_sand_grain_roughness_uncertainty = (
        low = 1.0,
        high = 5.0,
        nominal = 2.0,
        type = :parameter_modifier,
        worst_case = :high
    ),
)

properties = ComponentVector(
    expansion_ratio = cad_data.expansion_ratio,

    #Aluminum's conductivity stays roughly the same throughout the temperatures we're working at
    #This is also taken directly from JLC3DP's specifications for the SLM 6061 aluminum we're using
    wall_thermal_conductivity = 172.0u"W/(m*K)",

    #JLC3DP has a recommended wall thickness of 2.0mm, they say that their minimum wall thickness is 1.5mm, but let's just 
    #use 2.0mm because I'm sure the 0.5mm deviation won't make a difference and likely make it easier to manufacture
    #nevermind, I'm being told that 1.5mm is closer to the absolute maximum of what the wall thickness should be 
    #for regeneratively cooled nozzles, so we'll use that
    min_wall_thickness = 1.5u"mm",

    #JLC3DP cites it as 276 MPa for 6061 SLM aluminum, but at 200degC (about our operating temprerature),
    #it can go as low as 130 MPa, so we'll use that instead
    wall_yield_strength = 130.0u"MPa",

    wall_density = 2700.0u"kg/m^3",
    wall_cp = 896.0u"J/(kg*K)",

    wall_elastic_modulus = 69.5u"GPa",
    wall_poissons_ratio = 0.33,
    wall_roughness = 8.0u"μm",
    wall_equivalent_sand_grain_roughness = 0.0u"m",

    number_of_ribs = 1, 
    #we're not actually going to adjust the number of ribs, but I just put it in some places just 
    #so I know where having more than one rib modifies the geometry functions

    #Outlet orifice properties
    outlet_orifice_discharge_coefficient = 0.61,
    #most sharp edged orifices are within 0.01 of this value, so we'll just use it without making it an uncertain parameter
    #I'm pretty sure that a sharp edged orifice is going to be the most preferrable design in our case because I've heard that it's more 
    #predictable when the fluid flashes after expansion.

    ambient_pressure = 1.0u"atm",

    #CACHED PROPERTIES:
    coolant_inlet_temperature = 21.0u"°C" |> u"K",
    coolant_temperature = 0.0u"K",
    #coolant_inlet_pressure = 140.0u"bar", #this is an optimized parameter 
    coolant_pressure = 0.0u"bar",
    coolant_cp = 0.0u"J/(kg*K)",
    coolant_dynamic_viscosity = 0.0u"Pa*s",
    coolant_thermal_conductivity = 0.0u"W/(m*K)",
    coolant_surface_tension = 0.0u"N/m",
    coolant_saturation_temperature = 0.0u"K",
    coolant_enthalpy = 0.0u"J/kg",
    coolant_density = 0.0u"kg/m^3",
    coolant_prandtl_number = 0.0,
    coolant_saturated_liquid_enthalpy_at_pressure = 0.0u"J/kg",
    latent_heat_of_vaporization_at_pressure = 0.0u"J/kg",

    wall_thickness = 0.0u"mm", #this is updated based on pressure stress, thermal stress, and manufacturing limits

    rough_friction_factor = 0.0,

    parallel_width_between_channel_walls = 0.0u"m", # updated based on parallel_width_between_channel_walls_to_channel_width_ratio
    rib_height = 0.0u"m", # updated based on rib_height_fraction_of_channel_height
    rib_width = 0.0u"m", # updated based on rib_top_width_to_hot_wall_width_ratio
)

Revise.includet(joinpath(@__DIR__, "variable_methods/uncertainties/uncertainty_helper_functions.jl"))

properties = add_uncertainties_to_properties(properties, uncertainties)

Revise.includet(joinpath(@__DIR__, "variable_methods/margins/margin_helper_functions.jl"))

#all margins will have a margin.desired and a margin.found
margins = (
    wall_thermal_stress_margin = (
        desired = 1.5,
        type = :sized_to,
    ),
    wall_pressure_stress_margin = (
        desired = 1.5,
        type = :sized_to,
    ),
    critical_heat_flux_margin = (
        desired = 1.5, #maybe even two, we really don't want to cross this
        type = :rated_against,
    ),
    #slightly different, max_wall_temperature / wall_temperature (so lower == higher margin)
    #We could technically discard of this one and find the yield strength of the aluminum at the operating temperature and compare that to thermal and pressure stresses
    #This is no longer needed because the nozzle_material's yield strength is now evaluated based on temperature!
    #=
    wall_temperature = Margin(
        desired = 0.8,
        type = :rated_against,
    ),
    =#
    #T_sat(station_pressure) - T_coolant
    subcooling_margin = (
        desired = 1.5,
        type = :rated_against
    ),
    #We may also want another margin to check for orifice flashing (IF we're not using an NHNE injector model, although we may just want to avoid flashing entirely)

    #=
    #coolant_heat_flux / nozzle_heat_flux
    #wait, this is not a good metric because we're doing a non-linear solve to find where these two match
    #in that case, I'm not exactly sure how to treat this margin with model uncertainties
    #maybe we just use the wall temperature margin?
    coolant_heat_flux_to_nozzle_heat_flux = Margin(
        desired = 0.0,
        type = :rated_against,
    ),
    =#

)

properties = add_margins_to_properties(properties, margins)

#=
For helix spacing:

In general, at a given radius you can choose any two of 
    - N_s (number of channels)
    - p (i.e., the rib width / pitch) (wall width between channel walls perpendicular to nozzle + channel width)
    - φ (helix angle) 
        - φ is unintuitive, so we're not going to use it
=#

optimized_parameters = (
    channel_width_parallel_to_nozzle = (
        min = 0.5u"mm",
        max = 3.0u"mm",
        guess = 1.0u"mm",
    ),
    channel_height_perpendicular_to_nozzle = (
        min = 2.0u"mm",
        max = 10.0u"mm",
        guess = 5.0u"mm",
    ),
    rib_height_fraction_of_channel_height = (
        min = 0.25,
        max = 0.75,
        guess = 0.5
    ),
    #represents the percentage of the wall along the width of the nozzle the tops of the ribs take up 
    #this determines rib width
    rib_top_width_to_hot_wall_width_ratio = (
        min = 0.4, 
        max = 0.6,
        guess = 0.5,
    ),
    #typically around 0.5 to 1.0 times the channel width (usually chosen as 1.0)
    #holy variable length batman! :)
    #This represents the ratio between the solid wall separating two adjacent channels from other channels
    parallel_width_between_channel_walls_to_channel_width_ratio = (
        min = 0.5,
        max = 1.0,
        guess = 1.0
    ),
    number_of_channels = (
        min = 1,
        max = 30,
        guess = 10
    ),
)

#TRIFECTA (don't know what to call this relationship) between coolant mass flow, outlet_orifice_area, and inlet pressure
coolant_inlet_pressure_optimized_property = (
    min = 50.0u"bar",
    max = 175.0u"bar",
    guess = 117.0u"bar",
)

coolant_mass_flow_optimized_property = (
    min = 0.05u"kg/s",
    max = 0.2u"kg/s",
    guess = 0.1u"kg/s"
)

outlet_orifice_area_optimized_property = (
    min = 1.0u"mm^2",
    max = 10.0u"mm^2",
    guess = 5.0u"mm^2",
)

if system.outlet_orifice_area_or_inlet_pressure_adjusted == :outlet_orifice_area
    #inlet pressure is chosen/fixed and then optimized
    #outlet_orifice_area is continuously adjusted
    #helps choose an outlet_orifice_area if you have a specific pump max pressure you know
    properties = merge_properties(properties, ComponentVector(
        coolant_inlet_pressure = coolant_inlet_pressure_optimized_property.guess,
        coolant_mass_flow = coolant_mass_flow_optimized_property.guess, 
        outlet_orifice_area = 0.0u"mm^2"
    ))
    optimized_parameters = merge(optimized_parameters, (coolant_inlet_pressure = coolant_inlet_pressure_optimized_property, ))
elseif system.outlet_orifice_area_or_inlet_pressure_adjusted == :coolant_inlet_pressure
    #outlet_orifice_area is chosen/fixed and then optimized
    #inlet_pressure is continuously adjusted
    #helps choose an inlet_pressure if you have a fixed outlet_orifice_area
    properties = merge_properties(properties, ComponentVector(
        coolant_inlet_pressure = 0.0u"bar",
        coolant_mass_flow = coolant_mass_flow_optimized_property.guess, 
        outlet_orifice_area = outlet_orifice_area_optimized_property.guess
    ))
    optimized_parameters = merge(optimized_parameters, (outlet_orifice_area = outlet_orifice_area_optimized_property, ))
end

Revise.includet(joinpath(@__DIR__, "variable_methods/optimized_parameters/optimized_parameter_helper_functions.jl"))
properties = add_optimized_parameters_to_properties(properties, optimized_parameters)

function update_properties!(
    p, system, coolant_temperature, coolant_pressure,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature
)
    p.rib_height = p.rib_height_fraction_of_channel_height * p.channel_height_perpendicular_to_nozzle
    p.rib_width = p.rib_top_width_to_hot_wall_width_ratio * p.channel_width_parallel_to_nozzle

    p.parallel_width_between_channel_walls = p.parallel_width_between_channel_walls_to_channel_width_ratio * p.channel_width_parallel_to_nozzle

    p.wall_equivalent_sand_grain_roughness = p.wall_equivalent_sand_grain_roughness_uncertainty * p.wall_roughness

    p.coolant_cp = PropsSI("Cpmass", "T", ustrip(u"K", coolant_temperature), "P", ustrip(u"Pa", coolant_pressure), system.coolprop_fluid) * u"J/(kg*K)"
    p.coolant_dynamic_viscosity = PropsSI("viscosity", "T", ustrip(u"K", coolant_temperature), "P", ustrip(u"Pa", coolant_pressure), system.coolprop_fluid) * u"Pa*s" 
    p.coolant_thermal_conductivity = PropsSI("conductivity", "T", ustrip(u"K", coolant_temperature), "P", ustrip(u"Pa", coolant_pressure), system.coolprop_fluid) * u"W/(m*K)"
    p.coolant_saturation_temperature = PropsSI("T", "P", ustrip(u"Pa", coolant_pressure), "Q", 0, system.coolprop_fluid) * u"K"
    p.coolant_surface_tension = PropsSI("surface_tension", "T", ustrip(u"K", coolant_temperature), "Q", 0, system.coolprop_fluid) * u"N/m"
    p.coolant_density = PropsSI("D", "T", ustrip(u"K", coolant_temperature), "P", ustrip(u"Pa", coolant_pressure), system.coolprop_fluid) * u"kg/m^3"
    p.coolant_enthalpy = PropsSI("H", "T", ustrip(u"K", coolant_temperature), "P", ustrip(u"Pa", coolant_pressure), system.coolprop_fluid) * u"J/kg"
    p.coolant_saturated_liquid_enthalpy_at_pressure = PropsSI("H", "P", ustrip(u"Pa", coolant_pressure), "Q", 0, system.coolprop_fluid) * u"J/kg"
    p.latent_heat_of_vaporization_at_pressure = PropsSI("H", "P", ustrip(u"Pa", coolant_pressure), "Q", 1, system.coolprop_fluid) * u"J/kg" - PropsSI("H", "P", ustrip(u"Pa", coolant_pressure), "Q", 0, system.coolprop_fluid) * u"J/kg"
    p.coolant_prandtl_number = p.coolant_cp * p.coolant_dynamic_viscosity / p.coolant_thermal_conductivity
end

function get_station_area_to_throat_area_based_on_mach_number(
    mach_number,
    p, system,
    nozzle_station
)
    stripped_pressure = ustrip(p.chamber_pressure |> u"Pa")

    cea_inputs = (
        stripped_pressure,
        p.oxidizer_to_fuel_ratio,
        p.expansion_ratio
    )

    cea_kwargs = (frozen = system.use_frozen_cea_interps,)

    gas_gamma = system.chamber_gamma_interp(cea_inputs...; cea_kwargs...)
    
    station_area_to_throat_area = sqrt((1 / mach_number)^2 * (
        (2 / (gas_gamma + 1)) *
        (1 + ((gas_gamma - 1) / 2) * mach_number^2)
    )^((gas_gamma + 1) / (gas_gamma - 1)))

    station_area = station_area_to_throat_area * system.cad_data.nozzle_throat_area

    guessed_station_diameter = 2 * sqrt(station_area / pi)

    return nozzle_station.middle_diameter - guessed_station_diameter
end

function get_bartz_heat_transfer_coefficient(p, system, nozzle_station, gas_side_wall_surface_temperature)
    cad_data = system.cad_data

    stripped_pressure = ustrip(p.chamber_pressure |> u"Pa")

    cea_inputs = (
        stripped_pressure,
        p.oxidizer_to_fuel_ratio,
        p.expansion_ratio,
    )

    cea_kwargs = (frozen = system.use_frozen_cea_interps,)

    if nozzle_station.is_converging
        bounds = (1e-9, 1.0)
    else
        bounds = (1.0, 20.0)
    end

    mach_number = find_zero(mach_number -> get_station_area_to_throat_area_based_on_mach_number(mach_number, p, system, nozzle_station), bounds)

    #standard assumption that the non_ideal_stagnation_temperature is roughly equal to p.combustion_efficiency^2 * ideal_stagnation_temperature 
    gas_stagnation_temperature = p.combustion_efficiency^2 * system.chamber_T_comb_interp(cea_inputs...; cea_kwargs...) * u"K"
    
    gas_prandtl = system.chamber_prandtl_interp(cea_inputs...; cea_kwargs...)
    gas_gamma = system.chamber_gamma_interp(cea_inputs...; cea_kwargs...)

    gas_temperature = gas_stagnation_temperature / (1 + ((gas_gamma - 1) / 2) * mach_number^2)

    #gas_pressure might be needed in the future for nozzle hoop stress
    gas_static_pressure = p.chamber_pressure * (gas_temperature / gas_stagnation_temperature)^(gas_gamma / (gas_gamma - 1))

    recovery_factor = gas_prandtl^(1/3) #this should be accurate for turbulent flow in the nozzle

    adiabatic_gas_at_wall_temperature = (
        gas_stagnation_temperature * (
            (1 + recovery_factor * ((gas_gamma - 1) / 2) * mach_number^2) /
            (1 + ((gas_gamma - 1) / 2) * mach_number^2)
        )   
    )

    boundary_layer_correction_factor = (
        (gas_side_wall_surface_temperature / (2 * gas_stagnation_temperature)) *
        (1 + ((gas_gamma - 1) / 2) * mach_number^2) + 1/2
    )^-0.68 * (
        1 + ((gas_gamma - 1) / 2) * mach_number^2
    )^-0.12

    gas_dynamic_viscosity = system.chamber_mu_interp(cea_inputs...; cea_kwargs...) * u"Pa*s"
    gas_cp = system.chamber_cp_interp(cea_inputs...; cea_kwargs...) * u"J/(kg*K)"
    c_star_effective = p.combustion_efficiency * system.chamber_cstar_interp(cea_inputs...; cea_kwargs...) * u"m/s"

    radius_of_curvature_of_nozzle_wall_at_throat = 1.25 * (cad_data.nozzle_throat_diameter / 2) #this doesn't have a big effect on the predicted heat transfer

    gas_side_convective_heat_transfer_coefficient = (
        (0.026 / cad_data.nozzle_throat_diameter^0.2) *
        ((gas_dynamic_viscosity^0.2 * gas_cp) / gas_prandtl^0.6) *
        ((p.chamber_pressure |> u"Pa") / c_star_effective)^0.8 *
        (cad_data.nozzle_throat_diameter / radius_of_curvature_of_nozzle_wall_at_throat)^0.1 *
        (cad_data.nozzle_throat_area / (pi * (nozzle_station.middle_diameter / 2)^2))^0.9 *
        boundary_layer_correction_factor
    )

    # Model uncertainties:
    gas_side_convective_heat_transfer_coefficient *= p.bart_heat_flux_correlation_uncertainty

    gas_side_convective_heat_transfer_coefficient = ustrip(upreferred(gas_side_convective_heat_transfer_coefficient)) * u"W/(m^2*K)"
    #since we're not using rational exponents, we can convert it to regular units easily, so we'll just strip it

    heat_flux_into_nozzle_walls = gas_side_convective_heat_transfer_coefficient * (adiabatic_gas_at_wall_temperature - gas_side_wall_surface_temperature)

    heat_transfer_wattage_into_nozzle_walls = gas_side_convective_heat_transfer_coefficient * (adiabatic_gas_at_wall_temperature - gas_side_wall_surface_temperature) * nozzle_station.wall_area

    returned_values = ComponentVector(
        wall_area = nozzle_station.wall_area,
        heat_transfer_wattage_into_nozzle_walls = heat_transfer_wattage_into_nozzle_walls,
        heat_flux_into_nozzle_walls = heat_flux_into_nozzle_walls,
        gas_side_convective_heat_transfer_coefficient = gas_side_convective_heat_transfer_coefficient,

        mach_number = mach_number,
        adiabatic_gas_at_wall_temperature = adiabatic_gas_at_wall_temperature,
        gas_temperature = gas_temperature,
        gas_stagnation_temperature = gas_stagnation_temperature, #we'll reuse this in the radiative heat flux function
        boundary_layer_correction_factor = boundary_layer_correction_factor,
        
        gas_static_pressure = gas_static_pressure,
    )

    return returned_values
end

function get_nozzle_radiative_heat_transfer(p, system, nozzle_station, gas_side_wall_surface_temperature, gas_stagnation_temperature)
    stefan_boltzman_constant = 5.670374419e-8 * u"W/(m^2*K^4)"

    radiative_heat_flux = ((p.wall_emissivity + 1) / 2) * p.gas_emissivity * stefan_boltzman_constant * (
        gas_stagnation_temperature^4 - gas_side_wall_surface_temperature^4
    )

    # Model Uncertainties:
    radiative_heat_flux *= p.radiative_heat_flux_uncertainty

    radiative_heat_transfer = radiative_heat_flux * nozzle_station.wall_area

    return radiative_heat_transfer
end

function get_overall_gas_side_heat_transfer(p, system, nozzle_station, gas_side_wall_surface_temperature)
    bartz_correlation_returned_values = get_bartz_heat_transfer_coefficient(p, system, nozzle_station, gas_side_wall_surface_temperature)
    
    radiative_heat_transfer = get_nozzle_radiative_heat_transfer(p, system, nozzle_station, gas_side_wall_surface_temperature, bartz_correlation_returned_values.gas_stagnation_temperature)

    overall_heat_transfer = bartz_correlation_returned_values.heat_transfer_wattage_into_nozzle_walls + radiative_heat_transfer

    return overall_heat_transfer, bartz_correlation_returned_values
end

function get_aluminum_6061_SLM_yield_strength_vs_temp_interp()
    # JLC3DP SLM "Aluminum 6061" — 0.2% yield strength vs. temperature
    # 23 °C: JLC3DP published minimum (≥250 MPa)
    # 100–250 °C: scaled from A6061-RAM2 L-PBF short-time (ASTM E21) retention ratios
    # 300–500 °C: estimated tail, ±30–50 %
    temperatures = [23.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0] .* u"°C"
    yield_strengths = [250.0, 250.0, 230.0, 172.0, 97.0, 56.0, 34.0, 22.0, 14.0, 10.0] .* u"MPa"

    # Anything evaluated at a lower temperature gets the 23°C yield strenght while this function extrapolates linearly above 500°C
    # Nevermind, I would rather have this fail because I'm going to bracket the solver

    return PCHIPInterpolation(
        ustrip.(u"Pa", yield_strengths), 
        ustrip.(u"K", temperatures),
        extrapolation_left = ExtrapolationType.Constant,
        #extrapolation_right = ExtrapolationType.Linear
    )
end

system = merge(system, (wall_yield_strength_at_temp_interp = get_aluminum_6061_SLM_yield_strength_vs_temp_interp(), ));

function get_min_wall_thickness_due_to_hot_wall_pressure_stress(
    p, system, 
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature, 
    coolant_pressure, bartz_correlation_returned_values
)
    gas_static_pressure = bartz_correlation_returned_values.gas_static_pressure

    @show p.hot_wall_pressure_stress_uncertainty
    @show abs(coolant_pressure - gas_static_pressure)
    @show p.channel_width_parallel_to_nozzle^2
    @show p.wall_pressure_stress_margin.desired
    @show p.wall_yield_strength

    temperature_to_evaluate_at = (gas_side_wall_surface_temperature + coolant_side_wall_surface_temperature) / 2
    wall_yield_strength = system.wall_yield_strength_at_temp_interp(ustrip(u"K", temperature_to_evaluate_at)) * u"Pa"

    min_wall_thickness_due_to_pressure_stress = sqrt(
        p.hot_wall_pressure_stress_uncertainty * (abs(coolant_pressure - gas_static_pressure) * p.channel_width_parallel_to_nozzle^2) * p.wall_pressure_stress_margin.desired / 
        (2 * wall_yield_strength)
    )
    
    return min_wall_thickness_due_to_pressure_stress
end

function get_wall_pressure_stress(
    p, system, 
    coolant_pressure, bartz_correlation_returned_values
)
    gas_static_pressure = bartz_correlation_returned_values.gas_static_pressure

    wall_pressure_stress = p.hot_wall_pressure_stress_uncertainty * 
        ((coolant_pressure - gas_static_pressure) * p.channel_width_parallel_to_nozzle^2) * 
        p.wall_pressure_stress_margin.desired / (2 * p.wall_thickness^2)
    
    return wall_pressure_stress
end

function get_max_wall_thickness_due_to_hot_wall_thermal_stress(
    p, system, 
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature, 
    bartz_correlation_returned_values
)
    surface_heat_flux = bartz_correlation_returned_values.heat_flux_into_nozzle_walls

    temperature_to_evaluate_at = gas_side_wall_surface_temperature
    wall_yield_strength = system.wall_yield_strength_at_temp_interp(ustrip(u"K", temperature_to_evaluate_at)) * u"Pa"
    
    # There's a temperature difference version of this as well, but it's not dependent on wall_thickness, so we're not going to use it

    # Thermal stress is self-limiting, so its limit can be 2 * yield_strength instead of just yield_strength.
    max_wall_thickness_due_to_thermal_stress = (
        (2 * (2 * wall_yield_strength / p.wall_thermal_stress_margin.desired) * p.wall_thermal_conductivity * (1 - p.wall_poissons_ratio)) / 
        (p.wall_elastic_modulus * p.wall_coefficient_of_thermal_expansion * surface_heat_flux * p.hot_wall_thermal_stress_uncertainty)
    )

    return max_wall_thickness_due_to_thermal_stress
end

function get_wall_thermal_stress(
    p, system, 
    bartz_correlation_returned_values
)
    surface_heat_flux = bartz_correlation_returned_values.heat_flux_into_nozzle_walls
    
    wall_thermal_stress = (
        (p.wall_elastic_modulus * p.wall_coefficient_of_thermal_expansion * surface_heat_flux * p.wall_thickness * p.hot_wall_thermal_stress_uncertainty) / 
        (2 * p.wall_thermal_conductivity * (1 - p.wall_poissons_ratio))
    )
    wall_thermal_stress /= 2
    
    return wall_thermal_stress
end

function update_min_wall_thickness!(
    p, system,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature, 
    coolant_pressure, bartz_correlation_returned_values
)
    pressure_stress_min_wall_thickness = get_min_wall_thickness_due_to_hot_wall_pressure_stress(p, system, gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature, coolant_pressure, bartz_correlation_returned_values)

    thermal_stress_max_wall_thickness = get_max_wall_thickness_due_to_hot_wall_thermal_stress(p, system, gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature, bartz_correlation_returned_values)

    manufacturing_limits_min_wall_thickness = p.min_wall_thickness

    #min_wall_thickness = max(pressure_stress_min_wall_thickness, manufacturing_limits_min_wall_thickness)
    min_wall_thickness = pressure_stress_min_wall_thickness

    #=
    if min_wall_thickness > thermal_stress_max_wall_thickness
        # I just want to fail loudly right now, when doing actual optimization I will disable this
        error(
            "No wall thickness can satisfy the minimum wall thickness imposeed by pressure stresses and the max wall thickness imposed by thermal stresses.  
            Thermal stress max wall thickness: $(thermal_stress_max_wall_thickness |> u"mm"), 
            pressure stress min wall thickness: $(pressure_stress_min_wall_thickness |> u"mm"), 
            manufacturing limits min wall thickness: $(manufacturing_limits_min_wall_thickness |> u"mm")"
        )
    end
    =#

    p.wall_thickness = min_wall_thickness

    #pressure_stress = get_wall_pressure_stress(p, system, coolant_pressure, bartz_correlation_returned_values)

    #NOTE: remember, to be consistent with the thermal_stress_max_wall_thickness, you have to either: 
        # - multiply the allowable yield strength by 2
        # - divide thermal stress by 2
    #thermal_stress = get_wall_thermal_stress(p, system, bartz_correlation_returned_values)
end

#this has to be called AFTER the wall thickness is determined
function get_wall_heat_transfer(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature
)
    wall_heat_flux = (p.wall_thermal_conductivity / p.wall_thickness) * (gas_side_wall_surface_temperature - coolant_side_wall_surface_temperature)

    wall_heat_transfer = wall_heat_flux * nozzle_station.wall_area

    return wall_heat_transfer
end

#NOTE: I mean, we could probably also consider the top wall that's parallel to the nozzle, but that's very optimistic, so I'm just going to ignore it


#SUMMARY OF HELIX GEOMETRY:
#=
thermal_hydraulic_diameter vs hydraulic_diameter
    - thermal_hydraulic_diameter omits the top wall (channel_width_parallel_to_nozzle)
    - hydraulic_diameter includes it
=#

function get_channel_wetted_perimeter(p, system)
    cold_wall_perimeter = p.channel_width_parallel_to_nozzle #this is the only line that's different from heated perimeter

    prime_hot_wall_perimeter = p.channel_width_parallel_to_nozzle - p.number_of_ribs * p.rib_width

    rib_perimeter = p.number_of_ribs * (2 * p.rib_height + p.rib_width)

    land_perpendicular_wall_perimeter = 2 * p.channel_height_perpendicular_to_nozzle

    return cold_wall_perimeter + prime_hot_wall_perimeter + rib_perimeter + land_perpendicular_wall_perimeter
end

function get_channel_heated_perimeter(p, system)
    prime_hot_wall_perimeter = p.channel_width_parallel_to_nozzle - p.number_of_ribs * p.rib_width

    rib_perimeter = p.number_of_ribs * (2 * p.rib_height + p.rib_width)

    land_perpendicular_wall_perimeter = 2 * p.channel_height_perpendicular_to_nozzle

    return prime_hot_wall_perimeter + rib_perimeter + land_perpendicular_wall_perimeter
end

function get_single_channel_area(p, system)
    channel_area = p.channel_width_parallel_to_nozzle * p.channel_height_perpendicular_to_nozzle - (p.number_of_ribs * p.rib_height * p.rib_width)
    
    return channel_area
end

function get_channel_hydraulic_diameter(p, system)
    single_channel_wetted_perimeter = get_channel_wetted_perimeter(p, system)

    channel_area = get_single_channel_area(p, system)

    return (4 * channel_area) / single_channel_wetted_perimeter
end

function get_channel_thermal_hydrualic_diameter(p, system)
    single_channel_heated_perimeter = get_channel_heated_perimeter(p, system)

    channel_area = get_single_channel_area(p, system)

    return (4 * channel_area) / single_channel_heated_perimeter
end

function get_helical_station(
    p, system, nozzle_station, coolant_temperature, coolant_pressure
)
    pitch = p.channel_width_parallel_to_nozzle + p.parallel_width_between_channel_walls

    channel_centerline_radius = nozzle_station.middle_diameter / 2 + p.wall_thickness + p.channel_height_perpendicular_to_nozzle / 2

    #pitch is measured perpendicular to the channels, so we use asin instead of atan
    if p.number_of_channels <= 2 * pi * channel_centerline_radius / pitch
        channel_helix_angle_relative_to_circumference = asin(
            (p.number_of_channels * pitch) / 
            (2 * pi * channel_centerline_radius)
        )
    else
        error("This channel geometry is not possible, we're probably going to want to add a loss here in the future")
    end

    delta_channel_length_over_delta_slant_height = (2 * pi * channel_centerline_radius) / (p.number_of_channels * pitch)

    channel_length_along_nozzle_station = delta_channel_length_over_delta_slant_height * nozzle_station.slant_length

    helix_curvature_radius = channel_centerline_radius / (cos(channel_helix_angle_relative_to_circumference)^2)

    hydraulic_diameter = get_channel_hydraulic_diameter(p, system)
    thermal_hydraulic_diameter = get_channel_thermal_hydrualic_diameter(p, system)

    single_channel_flow_area = get_single_channel_area(p, system)

    mass_flux = p.coolant_mass_flow / (p.number_of_channels * single_channel_flow_area)

    coolant_velocity = mass_flux / p.coolant_density

    coolant_reynolds_number = (mass_flux * hydraulic_diameter) / p.coolant_dynamic_viscosity

    coolant_dean_number = coolant_reynolds_number * sqrt(
        hydraulic_diameter / (2 * helix_curvature_radius)
    )

    helical_station = ComponentVector(
        pitch = pitch,
        channel_centerline_radius = channel_centerline_radius,
        channel_helix_angle_relative_to_circumference = channel_helix_angle_relative_to_circumference,
        delta_channel_length_over_delta_slant_height = delta_channel_length_over_delta_slant_height,
        channel_length_along_nozzle_station = channel_length_along_nozzle_station,
        helix_curvature_radius = helix_curvature_radius,
        hydraulic_diameter = hydraulic_diameter,
        thermal_hydraulic_diameter = thermal_hydraulic_diameter,
        single_channel_flow_area = single_channel_flow_area,
        mass_flux = mass_flux,
        coolant_velocity = coolant_velocity,
        coolant_reynolds_number = coolant_reynolds_number,
        coolant_dean_number = coolant_dean_number,
    )

    return helical_station
end

function _gnielinski_nusselt_number(
    friction_factor,
    p, system, nozzle_station, helical_station
)
    prandtl_number = p.coolant_dynamic_viscosity * p.coolant_cp / p.coolant_thermal_conductivity #represents prandtl number for the bulk fluid

    nusselt_number = (
        ((friction_factor / 8) * (helical_station.coolant_reynolds_number - 1000) * prandtl_number) /
        (1 + 12.7 * sqrt(friction_factor / 8) * (prandtl_number^(2/3) - 1))
    )

    return nusselt_number
end

function gnielinski_nusselt_number_with_rough_pipe_friction(
    rough_friction_factor, smooth_friction_factor,
    p, system, nozzle_station, helical_station
)

    nusselt_number = _gnielinski_nusselt_number(rough_friction_factor, p, system, nozzle_station, helical_station)

    return nusselt_number
end

function gnielinski_nusselt_number_with_norris_correction(
    rough_friction_factor, smooth_friction_factor,
    p, system, nozzle_station, helical_station
)
    prandtl_number = p.coolant_dynamic_viscosity * p.coolant_cp / p.coolant_thermal_conductivity

    scaling_exponent = 0.68 * prandtl_number^0.215

    rough_to_smooth_ratio = min(rough_friction_factor / smooth_friction_factor, 4.0) #when this ratio is above 4, it's hypothesized that the extra roughness does not increase heat transfer further

    Nu_rough_over_Nu_smooth = (rough_to_smooth_ratio)^scaling_exponent

    nusselt_number = _gnielinski_nusselt_number(
        smooth_friction_factor,
        p, system, nozzle_station, helical_station
    ) * Nu_rough_over_Nu_smooth

    return nusselt_number
end

Revise.includet(joinpath(@__DIR__, "friction_factors.jl"))

function get_forced_convection_coolant_side_heat_flux(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure,
)
    rough_friction_factor = 0.0
    smooth_friction_factor = 0.0

    if system.use_curved_rough_friction_factor
        rough_friction_factor = get_rough_curved_friction_factor(
            p, system, nozzle_station, helical_station
        )
    else
        rough_friction_factor = get_colebrook_rough_straight_friction_factor(
            p, system, nozzle_station, helical_station
        )
    end

    p.rough_friction_factor = rough_friction_factor #cache the rough friction factor for later
    
    if system.use_curved_smooth_friction_factor
        smooth_friction_factor = get_smooth_curved_friction_factor(
            p, system, nozzle_station, helical_station
        )
    else
        smooth_friction_factor = get_smooth_straight_friction_factor(
            p, system, nozzle_station, helical_station
        )
    end

    if system.forced_convection_mode == :gnielinski_norris
        nusselt_number = gnielinski_nusselt_number_with_norris_correction(
            rough_friction_factor, smooth_friction_factor,
            p, system, nozzle_station, helical_station
        )
    elseif system.forced_convection_mode == :gnielinski_rough_pipe
        nusselt_number = gnielinski_nusselt_number_with_rough_pipe_friction(
            rough_friction_factor, smooth_friction_factor,
            p, system, nozzle_station, helical_station
        )
    end

    heat_transfer_coefficient = (nusselt_number * p.coolant_thermal_conductivity) / helical_station.hydraulic_diameter

    if coolant_side_wall_surface_temperature < p.coolant_saturation_temperature - 0.01u"K" #coolprop can be unreliable close to the saturation line
        coolant_side_wall_temperature_cp = PropsSI("Cpmass", "T", ustrip(u"K", coolant_side_wall_surface_temperature), "P", ustrip(u"Pa", coolant_pressure), system.coolprop_fluid)
        coolant_side_wall_temperature_dynamic_viscosity = PropsSI("viscosity", "T", ustrip(u"K", coolant_side_wall_surface_temperature), "P", ustrip(u"Pa", coolant_pressure), system.coolprop_fluid)
        coolant_side_wall_temperature_thermal_conductivity = PropsSI("conductivity", "T", ustrip(u"K", coolant_side_wall_surface_temperature), "P", ustrip(u"Pa", coolant_pressure), system.coolprop_fluid)
    else
        coolant_side_wall_temperature_cp = PropsSI("Cpmass", "P", ustrip(u"Pa", coolant_pressure), "Q", 0, system.coolprop_fluid)
        coolant_side_wall_temperature_dynamic_viscosity = PropsSI("viscosity", "P", ustrip(u"Pa", coolant_pressure), "Q", 0, system.coolprop_fluid)
        coolant_side_wall_temperature_thermal_conductivity = PropsSI("conductivity", "P", ustrip(u"Pa", coolant_pressure), "Q", 0, system.coolprop_fluid)
    end

    coolant_side_wall_temperature_prandtl_number = coolant_side_wall_temperature_dynamic_viscosity * coolant_side_wall_temperature_cp / coolant_side_wall_temperature_thermal_conductivity

    # Gnielinkski's correction for a heated liquid
    if system.forced_convection_mode == :gnielinski_norris || system.forced_convection_mode == :gnielinski_rough_pipe
        heat_transfer_coefficient *= (p.coolant_prandtl_number / coolant_side_wall_temperature_prandtl_number)^0.11
    end

    heat_transfer_coefficient *= p.coolant_channels_heat_transfer_curvature_modifier 
    #I'm not exactly sure if the centrifugal force pushing the fluid against the wall opposite of the hot wall will reduce the heat transfer

    forced_convection_heat_flux = heat_transfer_coefficient * (coolant_side_wall_surface_temperature - coolant_temperature)

    return forced_convection_heat_flux
end

# Onset of nucleate boiling crelations START:
function bergles_rohsenow_onset_of_nucleate_boiling_heat_flux(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    T_sat = PropsSI("T", "P", ustrip(u"Pa", coolant_pressure), "Q", 0, system.coolprop_fluid)

    local_wall_surface_temperature = ustrip(u"K", coolant_side_wall_surface_temperature)

    local_coolant_pressure = ustrip(u"bar", coolant_pressure)

    if (local_wall_surface_temperature - T_sat) >= 0 
        onset_nucleate_boiling_heat_flux = 1082 * local_coolant_pressure^1.156 * (1.8 * (local_wall_surface_temperature - T_sat))^(2.16 / (local_coolant_pressure^0.0234)) * u"W/m^2"
    else
        return 0.0u"W/m^2"
    end

    onset_nucleate_boiling_heat_flux *= p.bergles_rohsenow_onset_of_nucleate_boiling_heat_flux_uncertainty

    return onset_nucleate_boiling_heat_flux
end

function get_onset_boiling_heat_flux(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    
    onset_nucleate_boiling_heat_flux = bergles_rohsenow_onset_of_nucleate_boiling_heat_flux(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
        coolant_temperature, coolant_pressure
    )

    return onset_nucleate_boiling_heat_flux
end
# END

# Fully developed boiling heat transfer correlations START:
function jens_lottes_fully_developed_boiling_heat_flux(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    heat_flux_in_mw_per_m2 = ((ustrip(u"K", coolant_side_wall_surface_temperature) - ustrip(u"K", p.coolant_saturation_temperature)) / 25)^4 * exp(4*ustrip(u"bar", coolant_pressure)/62) * u"MW/m^2"

    fully_developed_boiling_heat_flux = heat_flux_in_mw_per_m2 |> u"W/m^2"

    fully_developed_boiling_heat_flux *= p.jens_lottes_fully_developed_boiling_heat_flux_uncertainty

    return fully_developed_boiling_heat_flux
end

function thom_fully_developed_boiling_heat_flux(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    heat_flux_in_mw_per_m2 = ((ustrip(u"K", coolant_side_wall_surface_temperature) - ustrip(u"K", p.coolant_saturation_temperature)) / 22.65)^2 * exp(2*ustrip(u"bar", coolant_pressure)/87) * u"MW/m^2"

    fully_developed_boiling_heat_flux = heat_flux_in_mw_per_m2 |> u"W/m^2"

    fully_developed_boiling_heat_flux *= p.thom_fully_developed_boiling_heat_flux_uncertainty

    return fully_developed_boiling_heat_flux
end

function get_fully_developed_boiling_heat_flux(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    if system.fully_developed_boiling_heat_flux_mode == :jens_lottes
        fully_developed_boiling_heat_flux = jens_lottes_fully_developed_boiling_heat_flux(
            p, system, nozzle_station, helical_station,
            gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
            coolant_temperature, coolant_pressure
        )
    elseif system.fully_developed_boiling_heat_flux_mode == :thom
        fully_developed_boiling_heat_flux = thom_fully_developed_boiling_heat_flux(
            p, system, nozzle_station, helical_station,
            gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
            coolant_temperature, coolant_pressure
        )
    end

    return fully_developed_boiling_heat_flux
end
#END

# Critical heat flux correlations START:
function get_critical_heat_flux_hall_mudawar_2000(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    #Valid for:
        # only water
        # thermal_hydraulic_diameter: 0.25 - 15mm
        # mass_flux: 300-30_000 kg/(m^2*s)
        # pressure: 1-200 bar
        # local_thermodynamic_equilibrium_quality: -1.0 to -0.05
    mass_flux = p.coolant_mass_flow / (p.number_of_channels * helical_station.single_channel_flow_area)

    liquid_rho = PropsSI("D", "P", ustrip(u"Pa", coolant_pressure), "Q", 0, system.coolprop_fluid) * u"kg/m^3"
    vapor_rho = PropsSI("D", "P", ustrip(u"Pa", coolant_pressure), "Q", 1, system.coolprop_fluid) * u"kg/m^3"

    weber_number = mass_flux^2 * helical_station.thermal_hydraulic_diameter / (liquid_rho * p.coolant_surface_tension)

    local_thermodynamic_equilibrium_quality = (p.coolant_enthalpy - p.coolant_saturated_liquid_enthalpy_at_pressure) / (p.latent_heat_of_vaporization_at_pressure)

    if local_thermodynamic_equilibrium_quality > -0.05 #this correlation is not recommended when this quantity passes this value
        @warn "this correlation is not recommended for this local_thermodynamic_equilibrium_quality"
    end

    C1 = 0.0722; C2 = -0.312; C3 = -0.644; C4 = 0.900; C5 = 0.724;
    
    boiling_number_at_critical_heat_flux = C1 * weber_number^C2 * (liquid_rho / vapor_rho)^C3 * (1 - C4 * (liquid_rho / vapor_rho)^C5 * local_thermodynamic_equilibrium_quality)

    critical_heat_flux = boiling_number_at_critical_heat_flux * mass_flux * p.latent_heat_of_vaporization_at_pressure

    return critical_heat_flux
end

function _nucleate_boiling_onset_wall_temperature(
    onset_wall_temperature,
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    #We need to know the wall temperature where the nucleate boiling heat flux equals the forced convection heat flux
    forced_convection_heat_flux = get_forced_convection_coolant_side_heat_flux(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, onset_wall_temperature,
        coolant_temperature, coolant_pressure
    )

    onset_boiling_heat_flux = get_onset_boiling_heat_flux(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, onset_wall_temperature,
        coolant_temperature, coolant_pressure
    )

    return forced_convection_heat_flux - onset_boiling_heat_flux
end

function find_nucleate_boiling_onset_wall_temperature(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    
    #We need to find the wall temperature where the nucleate boiling heat flux equals the forced convection heat flux

    nucleate_boiling_onset_wall_temperature = find_zero(onset_wall_temperature -> 
        _nucleate_boiling_onset_wall_temperature(
            onset_wall_temperature,
            p, system, nozzle_station, helical_station,
            gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
            coolant_temperature, coolant_pressure
        ), 
        (p.coolant_saturation_temperature, p.coolant_saturation_temperature + 50.0u"K")
    )
    
    return nucleate_boiling_onset_wall_temperature
end

function get_overall_coolant_side_heat_flux(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    forced_convection_heat_flux = get_forced_convection_coolant_side_heat_flux(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
        coolant_temperature, coolant_pressure
    )

    nucleate_boiling_onset_wall_temperature = find_nucleate_boiling_onset_wall_temperature(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
        coolant_temperature, coolant_pressure
    )

    onset_nucleate_boiling_heat_flux = bergles_rohsenow_onset_of_nucleate_boiling_heat_flux(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
        coolant_temperature, coolant_pressure
    )

    fully_developed_boiling_heat_flux = get_fully_developed_boiling_heat_flux(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
        coolant_temperature, coolant_pressure
    )

    fully_developed_boiling_heat_flux_evaluated_at_onset_temperature = get_fully_developed_boiling_heat_flux(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, nucleate_boiling_onset_wall_temperature,
        coolant_temperature, coolant_pressure
    )

    if coolant_side_wall_surface_temperature < nucleate_boiling_onset_wall_temperature
        return forced_convection_heat_flux
    else
        boiling_heat_flux = forced_convection_heat_flux * (
            1 + ((fully_developed_boiling_heat_flux / forced_convection_heat_flux) * 
                (1 - (fully_developed_boiling_heat_flux_evaluated_at_onset_temperature / fully_developed_boiling_heat_flux))
            )^2
        )^0.5

        critical_heat_flux = get_critical_heat_flux_hall_mudawar_2000(
            p, system, nozzle_station, helical_station,
            gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
            coolant_temperature, coolant_pressure
        )

        p.critical_heat_flux_margin.found = critical_heat_flux / boiling_heat_flux #I know that this is overwritten at every station, I'll fix that later

        return boiling_heat_flux
    end
end


function _straight_fin_efficiency(coolant_side_heat_transfer_coefficient, fin_thickness, fin_length, p)
    # Only takes into account rib thickness if fin_length is set to find_length + 1/2 * fin_thickness

    fin_parameter = sqrt((2 * abs(coolant_side_heat_transfer_coefficient)) / (p.wall_thermal_conductivity * fin_thickness)) #units of 1/m

    mL = ustrip(fin_parameter * fin_length)

    if mL == 0
        return 1.0
    end

    fin_efficiency = (tanh(mL)) / (mL)

    return fin_efficiency
end

function get_channel_effective_perimeter(p, system, coolant_side_heat_transfer_coefficient)
    #Notice that we subdivide the problem into 3 parts (prime hot wall, ribs, and land perpendicular walls)
    prime_hot_wall_perimeter = p.channel_width_parallel_to_nozzle - p.number_of_ribs * p.rib_width

    rib_perimeter = p.number_of_ribs * (2 * p.rib_height + p.rib_width)
    rib_efficiency = _straight_fin_efficiency(coolant_side_heat_transfer_coefficient, p.rib_width, p.rib_height + p.rib_width / 2, p)

    #Pessimistically, we could also use:
    #rib_perimeter = p.number_of_ribs * 2 * p.rib_height
    #rib_efficiency = _straight_fin_efficiency(h, p.rib_width, p.rib_height, p)

    land_perpendicular_wall_perimeter = 2 * p.channel_height_perpendicular_to_nozzle
    land_perpendicular_wall_efficiency = _straight_fin_efficiency(coolant_side_heat_transfer_coefficient, p.parallel_width_between_channel_walls, p.channel_height_perpendicular_to_nozzle, p)

    effective_perimeter = prime_hot_wall_perimeter + rib_efficiency * rib_perimeter + land_perpendicular_wall_efficiency * land_perpendicular_wall_perimeter

    return effective_perimeter
end

function get_fin_geometry_modifier_parameters(
    p, system, nozzle_station, helical_station,
    coolant_side_heat_transfer_coefficient
)
    heated_perimeter = get_channel_heated_perimeter(p, system)
    effective_perimeter = get_channel_effective_perimeter(p, system, coolant_side_heat_transfer_coefficient)

    area_ratio = heated_perimeter / helical_station.pitch
    overall_surface_efficiency = effective_perimeter / heated_perimeter

    return area_ratio, overall_surface_efficiency
end

function get_overall_coolant_side_heat_transfer(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure,
)
    unmodified_boiling_heat_flux = get_overall_coolant_side_heat_flux(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
        coolant_temperature, coolant_pressure,
    )

    unmodified_boiling_heat_transfer_coefficient = unmodified_boiling_heat_flux / (coolant_side_wall_surface_temperature - coolant_temperature) 
    #this is probably wrong, no idea how to get the heat transfer coefficient with correlations that return the heat flux directly 

    area_ratio, overall_surface_efficiency = get_fin_geometry_modifier_parameters(
        p, system, nozzle_station, helical_station,
        unmodified_boiling_heat_transfer_coefficient
    )

    boiling_heat_transfer = overall_surface_efficiency * area_ratio * unmodified_boiling_heat_flux * nozzle_station.wall_area

    return boiling_heat_transfer
end

function _converge_wall_surface_temperatures(
    du, u, p, system, nozzle_station,
    coolant_temperature, coolant_pressure
)
    gas_side_wall_surface_temperature = u[1] * u"K"
    coolant_side_wall_surface_temperature = u[2] * u"K"

    update_properties!(
        p, system, coolant_temperature, coolant_pressure,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature
    )

    overall_gas_side_heat_transfer, bartz_correlation_returned_values = get_overall_gas_side_heat_transfer(
        p, system, nozzle_station, 
        gas_side_wall_surface_temperature
    )


    @show gas_side_wall_surface_temperature
    @show coolant_side_wall_surface_temperature
    @show temperature_to_evaluate_at = max(gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature)
    @show p.wall_yield_strength = system.wall_yield_strength_at_temp_interp(ustrip(u"K", temperature_to_evaluate_at)) * u"Pa"

    update_min_wall_thickness!(
        p, system,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
        coolant_pressure, bartz_correlation_returned_values
    )

    helical_station = get_helical_station(
        p, system, nozzle_station, coolant_temperature, coolant_pressure
    )

    wall_heat_transfer = get_wall_heat_transfer(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature
    )

    overall_coolant_side_heat_transfer = get_overall_coolant_side_heat_transfer(
        p, system, nozzle_station, helical_station,
        gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
        coolant_temperature, coolant_pressure
    )

    du[1] = ustrip(u"W", overall_gas_side_heat_transfer - wall_heat_transfer)
    du[2] = ustrip(u"W", wall_heat_transfer - overall_coolant_side_heat_transfer)

    return nothing
end

function march_nozzle_length(
    #we're going to find the coolant_inlet_pressure that causes the orifice upstream of the 
    #cooling channels to discharge to ambient pressure
    #we could make this an ODEProblem, but I just want to get this working
    p, system, nozzle_stations,
)
    p.coolant_temperature = p.coolant_inlet_temperature
    p.coolant_pressure = p.coolant_inlet_pressure

    coolant_temperatures = []
    coolant_pressures = []

    gas_side_wall_surface_temperature = 400.0u"°C" |> u"K"
    coolant_side_wall_surface_temperature = 150.0u"°C" |> u"K"

    for nozzle_station in nozzle_stations
        update_properties!(
            p, system, p.coolant_temperature, p.coolant_pressure,
            gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature
        )

        _converge_wall_surface_temperatures_closure = (du, u, p) -> _converge_wall_surface_temperatures(
            du, u, p, system, nozzle_station,
            p.coolant_temperature, p.coolant_pressure
        )

        lb = ustrip.(uconvert.(u"K", [100.0u"°C", 20.0u"°C"]))
        ub = ustrip.(uconvert.(u"K", [1200.0u"°C", 300.0u"°C"]))
        
        prob = NonlinearProblem(_converge_wall_surface_temperatures_closure, [ustrip(u"K", gas_side_wall_surface_temperature), ustrip(u"K", coolant_side_wall_surface_temperature)], p; lb = lb, ub = ub)

        sol = solve(prob, FastShortcutNonlinearPolyalg(autodiff = AutoFiniteDiff()))

        # Assign the converged wall temperatures to the next guess because they're probably quite close
        gas_side_wall_surface_temperature = sol.u[1] * u"K"
        coolant_side_wall_surface_temperature = sol.u[2] * u"K"

        # We must call the residual to the solution one more time so that we don't get any finite difference values in p
        _converge_wall_surface_temperatures(similar(sol.u), sol.u, p, system, nozzle_station, p.coolant_temperature, p.coolant_pressure)

        overall_gas_side_heat_transfer, bartz_correlation_returned_values = get_overall_gas_side_heat_transfer(
            p, system, nozzle_station, 
            gas_side_wall_surface_temperature
        )

        update_min_wall_thickness!(
            p, system,
            p.coolant_pressure, bartz_correlation_returned_values
        )

        helical_station = get_helical_station(
            p, system, nozzle_station, 
            p.coolant_temperature, p.coolant_pressure
        )

        # Find the change in the next station's enthalpy and pressure
        # Find the coolant side heat transfer once again
        overall_coolant_side_heat_transfer = get_overall_coolant_side_heat_transfer(
            p, system, nozzle_station, helical_station,
            gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
            p.coolant_temperature, p.coolant_pressure
        )

        coolant_pressure_drop = get_helix_frictional_pressure_drop(
            p, system, nozzle_station, helical_station,
            gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
            p.coolant_temperature, p.coolant_pressure
        )

        # update both temperature and pressure
        new_coolant_enthalpy = p.coolant_enthalpy + (overall_coolant_side_heat_transfer / p.coolant_mass_flow)
        p.coolant_temperature = PropsSI("T", "H", ustrip(u"J/kg", new_coolant_enthalpy), "P", ustrip(u"Pa", p.coolant_pressure), system.coolprop_fluid) * u"K"

        p.coolant_pressure -= coolant_pressure_drop

        #update next nozzle station's pressure and temperature
        push!(coolant_temperatures, p.coolant_temperature)
        push!(coolant_pressures, p.coolant_pressure)
    end

    return p, system
end

#OPTION 1 (determines optimial orifice for a given inlet pressure) START
function determine_optimal_orifice!(
    p, system, nozzle_stations,
)
    p, system = march_nozzle_length(p, system, nozzle_stations)

    p.outlet_orifice_area = p.coolant_mass_flow / (p.outlet_orifice_discharge_coefficient * sqrt(2 * p.coolant_density * (p.coolant_pressure - p.ambient_pressure)))

    return nothing
end
#END

#OPTION 2: finds the required pump max outlet pressure for a fixed outlet_orifice_area START
function _outlet_pressure_residual(
    coolant_inlet_pressure, p, system, nozzle_stations
)
    p.coolant_inlet_pressure = coolant_inlet_pressure

    p, system = march_nozzle_length(p, system, nozzle_stations)

    #calculate pressure based on the orifice area
    #NOTE: just for now, we're going to assume incompressible flow and just error if the fluid would flash after expansion
    #in the future, we'll want to use a NHNE injector model so we can be more certain about this
    outlet_pressure = p.coolant_pressure - ((p.coolant_mass_flow / (p.outlet_orifice_discharge_coefficient * p.outlet_orifice_area))^2 / (2 * p.coolant_density))

    return outlet_pressure - p.ambient_pressure
end

function determine_inlet_pressure!(
    p, system, nozzle_stations
)
    p.coolant_inlet_pressure = find_zero(coolant_inlet_pressure -> _outlet_pressure_residual(
        coolant_inlet_pressure, p, system, nozzle_stations), 
        (50.0u"bar", 300.0u"bar")
    )
end
#END

function solve_overall_system!(
    p, system, nozzle_stations
)
    if system.outlet_orifice_area_or_inlet_pressure_adjusted == :outlet_orifice_area
        determine_optimal_orifice!(p, system, nozzle_stations)
    elseif system.outlet_orifice_area_or_inlet_pressure_adjusted == :coolant_inlet_pressure
        determine_inlet_pressure!(p, system, nozzle_stations)
    end
end

nominal_properties = get_nominal_p(properties, uncertainties)

solve_overall_system!(nominal_properties, system, nozzle_stations)

#margin ratio = (capacity − demand) / combined uncertainty of capacity and demand