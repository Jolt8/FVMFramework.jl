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
        - Effective diameter of cooling channels
        - Cooling channel width, height and aspect ratio
        - Pitch of cooling channels / coils of cooling channels per unit length

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

function create_nozzle_positions(cad_data)
    n_nozzle_positions = 1000

    nozzle_positions = NozzlePosition[]

    nozzle_throat_position = cad_data.nozzle_length - cad_data.nozzle_diverging_length

    nozzle_top_and_bottom_positions = collect(range(0.0u"mm", cad_data.nozzle_length, n_nozzle_positions + 1))

    first_middle_position = (nozzle_top_and_bottom_positions[1] + nozzle_top_and_bottom_positions[2]) / 2
    last_middle_position = (nozzle_top_and_bottom_positions[end-1] + nozzle_top_and_bottom_positions[end]) / 2

    nozzle_middle_positions = collect(range(first_middle_position, last_middle_position, n_nozzle_positions))

    for i in 1:n_nozzle_positions
        top_position = nozzle_top_and_bottom_positions[i]
        middle_position = nozzle_middle_positions[i]
        bottom_position = nozzle_top_and_bottom_positions[i+1]

        top_diameter = get_diameter_at_nozzle_position(cad_data, top_position, nozzle_throat_position)
        middle_diameter = get_diameter_at_nozzle_position(cad_data, middle_position, nozzle_throat_position)
        bottom_diameter = get_diameter_at_nozzle_position(cad_data, bottom_position, nozzle_throat_position)

        if middle_position < nozzle_throat_position
            is_converging = true
            half_angle = cad_data.converging_half_angle
        else
            is_converging = false
            half_angle = cad_data.diverging_half_angle
        end

        wall_area = pi * (top_diameter + bottom_diameter) / 2 * abs(bottom_diameter - top_diameter) / (2 * sind(half_angle))

        push!(nozzle_positions, NozzlePosition(
            top_position,
            middle_position,
            bottom_position,
            top_diameter,
            middle_diameter,
            bottom_diameter,
            is_converging,
            half_angle,
            wall_area
        ))
    end

    return nozzle_positions
end

nozzle_positions = create_nozzle_positions(cad_data)

system = (
    use_frozen_cea_interps = true, 
    #equilibrium gives a much higher heat transfer coefficient, so we'll use that to be safe
    #I think if we use equilibrium to get a higher heat transfer coefficient, I would be comfortable dropping the 2x safety factor and swapping it out for like a 1.2x safety factor
)

Revise.includet(joinpath(@__DIR__, "cea_retrieval", "CEA_lookup_table.jl"));

system = merge(system, (
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

    cad_data = cad_data
));

properties = ComponentVector(
    chamber_pressure = 20.0u"bar",
    #higher O/F ratios drives chamber temperatures higher, we'll use a higher than normal value for safety
    oxidizer_to_fuel_ratio = 9.0,
    expansion_ratio = cad_data.expansion_ratio,
    combustion_efficiency = 0.9, #this is pretty generous, but higher combustion efficiency leads to hotter combustion gases
)

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

function get_bartz_heat_transfer_coefficient(
    p, system, nozzle_position, wall_temperature
)
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
    gas_stagnation_temperature = p.combustion_efficiency^2 * chamber_T_comb_interp(cea_inputs...; cea_kwargs...) * u"K"
    
    gas_prantdl = chamber_prandtl_interp(cea_inputs...; cea_kwargs...)
    gas_gamma = chamber_gamma_interp(cea_inputs...; cea_kwargs...)

    gas_temperature = gas_stagnation_temperature / (1 + ((gas_gamma - 1) / 2) * mach_number^2)

    #gas_pressure might be needed in the future for nozzle hoop stress
    gas_pressure = p.chamber_pressure * (gas_temperature / gas_stagnation_temperature)^(gas_gamma / (gas_gamma - 1))

    recovery_factor = gas_prantdl^(1/3) #this should be accurate for turbulent flow in the nozzle

    adiabatic_gas_at_wall_temperature = (
        gas_temperature * (
            (1 + recovery_factor * ((gas_gamma - 1) / 2) * mach_number^2) /
            (1 + ((gas_gamma - 1) / 2) * mach_number^2)
        )   
    )

    boundary_layer_correction_factor = (
        (wall_temperature / (2 * gas_temperature)) *
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

    @show gas_side_convective_heat_transfer_coefficient |> u"W/m^2/K"

    heat_transfer_wattage_into_nozzle_walls = gas_side_convective_heat_transfer_coefficient * (gas_temperature - wall_temperature) * nozzle_position.wall_area

    returned_values = ComponentVector(
        wall_area = nozzle_position.wall_area,
        heat_transfer_wattage_into_nozzle_walls = heat_transfer_wattage_into_nozzle_walls,
        gas_side_convective_heat_transfer_coefficient = gas_side_convective_heat_transfer_coefficient,

        mach_number = mach_number,
        gas_temperature = gas_temperature,
        boundary_layer_correction_factor = boundary_layer_correction_factor,
        
        gas_pressure = gas_pressure,
    )

    return returned_values
end

using Plots

function plot_important_values()
    returned_values_vec = ComponentVector[]

    for nozzle_position in nozzle_positions
        returned_values = get_bartz_heat_transfer_coefficient(
            properties, system, nozzle_position, upreferred(150.0u"°C")
        )

        #@show returned_values

        push!(returned_values_vec, returned_values)
    end

    for name in propertynames(returned_values_vec[1])
        plt = plot(
            [nozzle_position.middle_position for nozzle_position in nozzle_positions], 
            ustrip.(upreferred.([getproperty(returned_values, name) for returned_values in returned_values_vec])),
            title = String(name)
        )

        display(plt)
    end
end

plot_important_values()