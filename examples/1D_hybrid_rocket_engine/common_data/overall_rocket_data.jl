function get_overall_rocket_data()
    specific_impulse = 230.0u"s"
    oxidizer_to_fuel_ratio = 7.6

    impulse_per_kg_of_propellants = specific_impulse * 9.81u"m/s^2"

    impulse_per_kg_of_oxidizer = impulse_per_kg_of_propellants * ((1 + oxidizer_to_fuel_ratio) / oxidizer_to_fuel_ratio)
    impulse_per_kg_of_fuel = impulse_per_kg_of_propellants * ((1 + oxidizer_to_fuel_ratio) / 1)

    desired_burn_time = 5.0u"s"
    average_thrust = 1000.0u"N"

    total_impulse = average_thrust * desired_burn_time

    kg_of_oxidizer = total_impulse / impulse_per_kg_of_oxidizer |> u"kg"
    kg_of_fuel = total_impulse / impulse_per_kg_of_fuel |> u"kg"

    oxidizer_mass_flow = kg_of_oxidizer / desired_burn_time

    tank_initial_pressure = 60.0u"bar"
    tank_initial_temperature = 21.0u"°C"
    
    return ComponentVector(
        desired_burn_time = desired_burn_time,
        total_impulse = total_impulse,
        oxidizer_mass_flow = oxidizer_mass_flow,
        required_kg_of_oxidizer = kg_of_oxidizer,
        required_kg_of_fuel = kg_of_fuel,
        tank_initial_pressure = tank_initial_pressure,
        tank_initial_temperature = tank_initial_temperature,
    )
end