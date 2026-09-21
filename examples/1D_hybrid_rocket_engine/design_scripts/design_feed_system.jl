using Revise
using Unitful
using ComponentArrays

Revise.includet(joinpath(@__DIR__, "../common_data/overall_rocket_data.jl"))

overall_rocket_data = get_overall_rocket_data()

tank_pressure = overall_rocket_data.tank_initial_pressure
desired_chamber_pressure = 30u"bar"

overshoot_factor = 1.2

post_adjustable_valve_minimum_pressure = desired_chamber_pressure * overshoot_factor

injector_orifice_area = get_injector_orifice_area_from_valve_cv(tank_pressure, 1.0)


#NOTE: chamber pressure is determined by propellant mass flow, nozzle throat area, and characteristic velocity of the propellants