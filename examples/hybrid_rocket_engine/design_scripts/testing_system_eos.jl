tank_oxidizer_moles = properties.u0_tank_oxidizer_mass / properties.nitrous_oxide_molecular_weight
tank_oxidizer_internal_energy = Clapeyron.VT0.internal_energy(
    oxidizer_model,
    ustrip(upreferred(properties.tank_volume)),
    ustrip(upreferred(properties.tank_temperature)),
    [ustrip(upreferred(tank_oxidizer_moles))]
)

lower_temperature = 180.0
upper_temperature = 450.0

energy_residual_tank(T) = begin
    result = vt_flash(oxidizer_model, ustrip(upreferred(properties.tank_volume)), T, [ustrip(upreferred(tank_oxidizer_moles))])
    internal_energy(oxidizer_model, result) - tank_oxidizer_internal_energy
end

energy_residual_tank(lower_temperature)
energy_residual_tank(upper_temperature)

#Mid Section
u0_mid_section_density = mass_density(
    oxidizer_model,
    ustrip(upreferred(properties.mid_section_pressure)),
    ustrip(upreferred(properties.mid_section_temperature)),
    [ustrip(upreferred(tank_oxidizer_moles))]
) * u"kg/m^3"
mid_section_mass = u0_mid_section_density * properties.mid_section_volume |> u"kg"

mid_section_oxidizer_moles = mid_section_mass / properties.nitrous_oxide_molecular_weight |> u"mol"
mid_section_internal_energy = Clapeyron.VT0.internal_energy(
    oxidizer_model,
    ustrip(upreferred(properties.mid_section_volume)),
    ustrip(upreferred(properties.mid_section_temperature)),
    [ustrip(upreferred(mid_section_oxidizer_moles))]
)

energy_residual_mid_section(T) = begin
    result = vt_flash(oxidizer_model, ustrip(upreferred(properties.mid_section_volume)), T, [ustrip(upreferred(mid_section_oxidizer_moles))])
    internal_energy(oxidizer_model, result) - mid_section_internal_energy
end

energy_residual_mid_section(lower_temperature)
energy_residual_mid_section(upper_temperature)

#Chamber
u0_chamber_gas_density = mass_density(chamber_model, properties.chamber_pressure, properties.chamber_temperature, [tank_oxidizer_moles, 0.0u"mol"])
chamber_gas_mass = u0_chamber_gas_density * properties.chamber_volume

chamber_gas_oxidizer_moles = (chamber_gas_mass * properties.u0_chamber_nitrous_oxide_mass_fraction) / properties.nitrous_oxide_molecular_weight
chamber_gas_fuel_moles = (chamber_gas_mass * properties.u0_chamber_hdpe_mass_fraction) / properties.ethylene_molecular_weight
chamber_gas_internal_energy = Clapeyron.VT0.internal_energy(
    chamber_model,
    ustrip(upreferred(properties.chamber_volume)),
    ustrip(upreferred(properties.chamber_temperature)),
    [ustrip(upreferred(chamber_gas_oxidizer_moles)), ustrip(upreferred(chamber_gas_fuel_moles))]
)

energy_residual_chamber(T) = begin
    result = vt_flash(chamber_model, ustrip(upreferred(properties.chamber_volume)), T, [ustrip(upreferred(chamber_gas_oxidizer_moles)), ustrip(upreferred(chamber_gas_fuel_moles))])
    internal_energy(chamber_model, result) - chamber_gas_internal_energy
end

energy_residual_chamber(lower_temperature)
energy_residual_chamber(upper_temperature)