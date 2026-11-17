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

nozzle_throat_position = cad_data.nozzle_length - cad_data.nozzle_diverging_length

struct NozzlePosition
    top_position::typeof(1.0u"m")
    middle_position::typeof(1.0u"m")
    bottom_position::typeof(1.0u"m")
    top_diameter::typeof(1.0u"m")
    middle_diameter::typeof(1.0u"m")
    bottom_diameter::typeof(1.0u"m")
    is_converging::Bool
    is_diverging::Bool
    half_angle::typeof(1.0u"°")
    wall_area::typeof(1.0u"m^2")
end

function get_diameter_at_nozzle_position(cad_data, nozzle_position, nozzle_throat_position)
    if nozzle_position < nozzle_throat_position
        return cad_data.nozzle_throat_diameter + 2 * (nozzle_throat_position - nozzle_position) * tand(cad_data.converging_half_angle)
    else
        return cad_data.nozzle_throat_diameter + 2 * (nozzle_position - nozzle_throat_position) * tand(cad_data.diverging_half_angle)
    end
end

function frustum_lateral_area(d_top, d_bottom, dz)
    r_top = d_top / 2
    r_bottom = d_bottom / 2
    slant_length = sqrt((r_bottom - r_top)^2 + dz^2)
    return pi * (r_top + r_bottom) * slant_length
end

function create_nozzle_positions(cad_data; n_converging = 100, n_diverging = 100)
    nozzle_throat_position = cad_data.nozzle_length - cad_data.nozzle_diverging_length

    # Separate grids for each cone so a segment edge lands exactly on the throat
    converging_edges = collect(range(0.0u"mm", nozzle_throat_position, n_converging + 1))
    diverging_edges  = collect(range(nozzle_throat_position, cad_data.nozzle_length, n_diverging + 1))
    segment_edges    = vcat(converging_edges, diverging_edges[2:end]) # drop the duplicated throat edge

    nozzle_positions = NozzlePosition[]

    for i in 1:(n_converging + n_diverging)
        top_position = segment_edges[i]
        bottom_position = segment_edges[i+1]
        middle_position = (top_position + bottom_position) / 2

        top_diameter = get_diameter_at_nozzle_position(cad_data, top_position, nozzle_throat_position)
        middle_diameter = get_diameter_at_nozzle_position(cad_data, middle_position, nozzle_throat_position)
        bottom_diameter = get_diameter_at_nozzle_position(cad_data, bottom_position, nozzle_throat_position)

        if middle_position < nozzle_throat_position
            is_converging = true
            is_diverging = false
            half_angle = cad_data.converging_half_angle
        else
            is_converging = false
            is_diverging = false
            half_angle = cad_data.diverging_half_angle
        end

        wall_area = frustum_lateral_area(top_diameter, bottom_diameter, bottom_position - top_position)

        push!(nozzle_positions, NozzlePosition(
            top_position,
            middle_position,
            bottom_position,
            top_diameter,
            middle_diameter,
            bottom_diameter,
            is_converging,
            is_diverging,
            half_angle,
            wall_area
        ))
    end

    return nozzle_positions
end

nozzle_positions = create_nozzle_positions(cad_data)

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
    nozzle_throat_position = cad_data.nozzle_length - cad_data.nozzle_diverging_length

    use_frozen_cea_interps = true, #equilibrium overpredicts way too much, we're always going to use frozen cea interps from now on
);

#we'll keep application as a parameter only used for clarification where it's needed
#for example, for valve cv, you would say application = :scalar or application = :multiplier
uncertainties = (
    #Since we're uncertain on the rate of recombination in the nozzle, we're going to multiply the bartz correlation results by this value in a parameter sweep
    recombination_factor = (
        low = 1.0,
        high = 1.37,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    #The bartz heat flux correlation is known to overpredict the heat transfer coefficient at the throat
    bart_heat_flux_correlation_uncertainty = (
        low = 0.6,
        high = 1.3,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    radiative_heat_flux_factor = (
        low = 0.5,
        high = 1.5,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    hot_wall_pressure_stress_factor = (
        low = 0.85,
        high = 1.3,
        nominal = 1.0,
        type = :model,
        worst_case = :high
    ),
    hot_wall_thermal_stress_factor = (
        low = 0.9,
        high = 1.25,
        nominal = 1.0,
        type = :model,
        worst_case = :high
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
    )
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

    #this is updated based on pressure stress, thermal stress, and manufacturing limits
    wall_thickness = 0.0u"mm",

    #JLC3DP cites it as 276 MPa for 6061 SLM aluminum, but at 200degC (about our operating temprerature),
    #it can go as low as 130 MPa, so we'll use that instead
    wall_yield_strength = 130.0u"MPa",

    wall_density = 2700.0u"kg/m^3",
    wall_cp = 896.0u"J/(kg*K)",

    wall_elastig_modulus = 69.5u"GPa",
    wall_poissons_ratio = 0.33,
    wall_roughness = 8.0u"μm",

    #Other properties that are cached and updated later
    rib_thickness = 0.0u"mm"
)

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
    rib_thickness_fraction_of_channel_height = (
        min = 0.25,
        max = 0.75,
        guess = 0.5,
        type = :parameter
    ),
    #typically around 0.5 to 1.0 times the channel width (usually chosen as 1.0)
    #holy variable length batman! :)
    width_between_channel_walls_perpendicular_to_nozzle_to_channel_width_ratio = (
        min = 0.5,
        max = 1.0,
        guess = 1.0
    ),
    number_of_channels = (
        min = 1,
        max = 30,
        guess = 10
    )
    # Wall thickness is kinda an optimized parameter, but it's determined by the pressure and thermal stresses during the solve 
    #and outputted later
)

function update_properties!(p)
    p.rib_thickness = p.rib_thickness_fraction_of_channel_height * p.channel_height_perpendicular_to_nozzle

    p.width_between_channel_walls_perpendicular_to_nozzle = p.width_between_channel_walls_perpendicular_to_nozzle_to_channel_width_ratio * p.channel_width_parallel_to_nozzle
end

function get_position_to_throat_area_based_on_mach_number(
    mach_number,
    p, system,
    nozzle_position
)
    stripped_pressure = ustrip(p.chamber_pressure |> u"Pa")

    cea_inputs = (
        stripped_pressure,
        p.oxidizer_to_fuel_ratio,
        p.expansion_ratio
    )

    cea_kwargs = (frozen = system.use_frozen_cea_interps,)

    gas_gamma = system.chamber_gamma_interp(cea_inputs...; cea_kwargs...)
    
    position_to_throat_area = sqrt((1 / mach_number)^2 * (
        (2 / (gas_gamma + 1)) *
        (1 + ((gas_gamma - 1) / 2) * mach_number^2)
    )^((gas_gamma + 1) / (gas_gamma - 1)))

    position_area = position_to_throat_area * system.cad_data.nozzle_throat_area

    guessed_position_diameter = 2 * sqrt(position_area / pi)

    return nozzle_position.middle_diameter - guessed_position_diameter
end

function get_bartz_heat_transfer_coefficient(p, system, nozzle_position, gas_side_wall_surface_temperature)
    cad_data = system.cad_data

    stripped_pressure = ustrip(p.chamber_pressure |> u"Pa")

    cea_inputs = (
        stripped_pressure,
        p.oxidizer_to_fuel_ratio,
        p.expansion_ratio,
    )

    cea_kwargs = (frozen = system.use_frozen_cea_interps,)

    if nozzle_position.is_converging
        bounds = (1e-9, 1.0)
    else
        bounds = (1.0, 20.0)
    end

    mach_number = find_zero(mach_number -> get_position_to_throat_area_based_on_mach_number(mach_number, p, system, nozzle_position), bounds)

    #standard assumption that the non_ideal_stagnation_temperature is roughly equal to p.combustion_efficiency^2 * ideal_stagnation_temperature 
    gas_stagnation_temperature = p.combustion_efficiency^2 * system.chamber_T_comb_interp(cea_inputs...; cea_kwargs...) * u"K"
    
    gas_prantdl = system.chamber_prandtl_interp(cea_inputs...; cea_kwargs...)
    gas_gamma = system.chamber_gamma_interp(cea_inputs...; cea_kwargs...)

    gas_temperature = gas_stagnation_temperature / (1 + ((gas_gamma - 1) / 2) * mach_number^2)

    #gas_pressure might be needed in the future for nozzle hoop stress
    gas_pressure = p.chamber_pressure * (gas_temperature / gas_stagnation_temperature)^(gas_gamma / (gas_gamma - 1))

    recovery_factor = gas_prantdl^(1/3) #this should be accurate for turbulent flow in the nozzle

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

    gas_dynamic_viscosity = chamber_mu_interp(cea_inputs...; cea_kwargs...) * u"Pa*s"
    gas_cp = chamber_cp_interp(cea_inputs...; cea_kwargs...) * u"J/(kg*K)"
    c_star_effective = p.combustion_efficiency * chamber_cstar_interp(cea_inputs...; cea_kwargs...) * u"m/s"

    radius_of_curvature_of_nozzle_wall_at_throat = 1.25 * (cad_data.nozzle_throat_diameter / 2) #this doesn't have a big effect on the predicted heat transfer

    gas_side_convective_heat_transfer_coefficient = (
        (0.026 / cad_data.nozzle_throat_diameter^0.2) *
        ((gas_dynamic_viscosity^0.2 * gas_cp) / gas_prantdl^0.6) *
        ((p.chamber_pressure |> u"Pa") / c_star_effective)^0.8 *
        (cad_data.nozzle_throat_diameter / radius_of_curvature_of_nozzle_wall_at_throat)^0.1 *
        (cad_data.nozzle_throat_area / (pi * (nozzle_position.middle_diameter / 2)^2))^0.9 *
        boundary_layer_correction_factor
    )

    # Model uncertainties:
    gas_side_convective_heat_transfer_coefficient *= p.bartz_correlation_safety_factor
    gas_side_convective_heat_transfer_coefficient *= p.bart_heat_flux_correlation_uncertainty

    gas_side_convective_heat_transfer_coefficient = ustrip(upreferred(gas_side_convective_heat_transfer_coefficient)) * u"W/(m^2*K)"
    #since we're not using rational exponents, we can convert it to regular units easily, so we'll just strip it

    heat_flux_into_nozzle_walls = gas_side_convective_heat_transfer_coefficient * (adiabatic_gas_at_wall_temperature - gas_side_wall_surface_temperature)

    heat_transfer_wattage_into_nozzle_walls = gas_side_convective_heat_transfer_coefficient * (adiabatic_gas_at_wall_temperature - gas_side_wall_surface_temperature) * nozzle_position.wall_area

    returned_values = ComponentVector(
        wall_area = nozzle_position.wall_area,
        heat_transfer_wattage_into_nozzle_walls = heat_transfer_wattage_into_nozzle_walls,
        heat_flux_into_nozzle_walls = heat_flux_into_nozzle_walls,
        gas_side_convective_heat_transfer_coefficient = gas_side_convective_heat_transfer_coefficient,

        mach_number = mach_number,
        adiabatic_gas_at_wall_temperature = adiabatic_gas_at_wall_temperature,
        gas_temperature = gas_temperature,
        gas_stagnation_temperature = gas_stagnation_temperature, #we'll reuse this in the radiative heat flux function
        boundary_layer_correction_factor = boundary_layer_correction_factor,
        
        gas_pressure = gas_pressure,
    )

    return returned_values
end

function get_nozzle_radiative_heat_flux(p, system, nozzle_position, gas_side_wall_surface_temperature, gas_stagnation_temperature)
    stefan_boltzman_constant = 5.670374419e-8 * u"W/(m^2*K^4)"

    radiative_heat_flux = ((p.wall_emissivity + 1) / 2) * gas_emissivity * stefan_boltzman_constant * (
        gas_stagnation_temperature^4 - gas_side_wall_surface_temperature^4
    )

    # Model Uncertainties:
    radiative_heat_flux *= p.radiative_heat_flux_factor

    return radiative_heat_flux
end

function overall_gas_side_heat_transfer(p, system, nozzle_position, gas_side_wall_surface_temperature)
    bartz_correlation_returned_values = get_bartz_heat_transfer_coefficient(p, system, nozzle_position, gas_side_wall_surface_temperature)
    
    radiative_heat_transfer = get_nozzle_radiative_heat_flux(p, system, nozzle_position, gas_side_wall_surface_temperature, bartz_correlation_returned_values.gas_stagnation_temperature)

    overall_heat_transfer = bartz_correlation_returned_values.heat_transfer_wattage_into_nozzle_walls + radiative_heat_transfer

    return overall_heat_transfer
end

function get_min_wall_thickness_due_to_hot_wall_pressure_stress(p, system, coolant_pressure, gas_stagnation_pressure)
    #bending_stress = ((coolant_pressure - gas_stagnation_pressure) * p.channel_height_perpendicular_to_nozzle^2) / (2 * )

    min_wall_thickness_due_to_pressure_stress = sqrt(
        ((coolant_pressure - gas_stagnation_pressure) * p.channel_height_perpendicular_to_nozzle^2) / 
        (2 * p.nozzle_material_yield_strength)
    )

    min_wall_thickness_due_to_pressure_stress *= p.pressure_stress_safety_factor 
    #I'm not sure this is accurate now since we're using the inverse of the typical equation
    
    return min_wall_thickness_due_to_pressure_stress
end

function get_min_wall_thickness_due_to_hot_wall_thermal_stress(p, system, bartz_correlation_returned_values)
    wall_thermal_diffusivity = p.wall_thermal_conductivity / (p.wall_density * p.wall_cp)
    surface_heat_flux = bartz_correlation_returned_values.heat_flux_into_nozzle_walls
    
    #There's a temperature difference version of this as well, but it's not dependent on wall_thickness, so we're not going to use it
    min_wall_thickness_due_to_thermal_stress = (
        (2 * p.nozzle_material_yield_strength * p.wall_thermal_conductivity * (1 - p.wall_poissons_ratio)) / 
        (p.wall_elastic_modulus * wall_thermal_diffusivity * surface_heat_flux)
    )

    min_wall_thickness_due_to_thermal_stress *= p.thermal_stress_safety_factor
    #Again, I'm not sure this is accurate now since we're using the inverse of the typical equation

    return min_wall_thickness_due_to_thermal_stress
end

function get_min_wall_thickness(
    p, system,
    coolant_pressure, gas_stagnation_pressure,
    bartz_correlation_returned_values
)
    pressure_stress_min_wall_thickness = get_min_wall_thickness_due_to_hot_wall_pressure_stress(p, system, coolant_pressure, gas_stagnation_pressure)
    thermal_stress_min_wall_thickness = get_min_wall_thickness_due_to_hot_wall_thermal_stress(p, system, bartz_correlation_returned_values)
    manufacturing_limits_min_wall_thickness = p.min_wall_thickness

    min_wall_thickness = max(pressure_stress_min_wall_thickness, thermal_stress_min_wall_thickness, manufacturing_limits_min_wall_thickness)

    p.wall_thickness = min_wall_thickness

    return min_wall_thickness
end

#this has to be called AFTER the wall thickness is determined
function wall_conduction(p, system, nozzle_position, gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature)
    return (p.wall_thermal_conductivity / p.wall_thickness) * (gas_side_wall_surface_temperature - coolant_side_wall_surface_temperature)
end

function get_helical_station_geometry_properties(
    p, system, nozzle_position
)
    channel_centerline_radius = nozzle_radius + p.wall_thickness + p.channel_width_parallel_to_nozzle / 2

    helix_angle
end

function overall_nozzle_position_heat_flux!(
    properties, system, nozzle_position,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,

)

end


using Plots

function plot_important_values(
    properties, system, nozzle_positions, 
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature
)
    returned_values_vec = ComponentVector[]

    for nozzle_position in nozzle_positions
        returned_values = get_bartz_heat_transfer_coefficient(
            properties, system, nozzle_position, upreferred(150.0u"°C")
        )



        push!(returned_values_vec, returned_values)
    end

    for name in propertynames(returned_values_vec[1])
        plt = plot(
            [nozzle_position.middle_position for nozzle_position in nozzle_positions], 
            [getproperty(returned_values, name) for returned_values in returned_values_vec],
            title = String(name)
        )

        display(plt)
    end

    total_heat_flux_into_nozzle = sum(ustrip.(upreferred.([getproperty(returned_values, :heat_transfer_wattage_into_nozzle_walls) for returned_values in returned_values_vec]))) * u"W"

    @show total_heat_flux_into_nozzle

    return nothing
end

properties = add_uncertainties_to_properties(properties, uncertainties)
properties = get_nominal_p(properties, uncertainties)

plot_important_values(
    properties, system, nozzle_positions,
    200.0u"°C", 110.0u"°C"
)

#margin ratio = (capacity − demand) / combined uncertainty of capacity and demand

#NOTE: our main unknowns here that we're going to use a nonlinearsolver for are:
    #- gas_side_wall_temperature
    #- coolant_side_wall_temperature