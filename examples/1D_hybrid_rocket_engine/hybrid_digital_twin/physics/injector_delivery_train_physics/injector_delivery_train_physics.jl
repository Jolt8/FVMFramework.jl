"""
    pure_component_vt_flash(model, volume, temperature, amount)

Construct the equilibrium `FlashResult` for a subcritical, pure-component
state specified by total volume, temperature, and amount in moles.

This follows the pure-fluid VT flash directly:

- use the EOS pressure outside the saturation dome;
- use saturation pressure and the lever rule inside the dome.

It avoids Clapeyron's generic `vt_flash`/`tx_flash_pure` compilation path,
which is pathologically expensive for this case on Julia 1.12 and Clapeyron
0.6.26. The implementation remains compatible with `ForwardDiff.Dual` and
the tracer overloads used by this project.

This helper is intentionally limited to one subcritical component. A
supercritical or multicomponent feed needs a separate implementation.
"""
function pure_component_vt_flash(model, volume, temperature, amount)
    volume > zero(volume) ||
        throw(DomainError(volume, "flash volume must be positive"))
    amount > zero(amount) ||
        throw(DomainError(amount, "flash amount must be positive"))

    saturation_pressure_value, liquid_molar_volume, vapor_molar_volume =
        saturation_pressure(model, temperature)
    bulk_molar_volume = volume / amount
    pure_composition = SVector(one(amount))

    if bulk_molar_volume < liquid_molar_volume
        state_pressure = pressure(
            model,
            volume,
            temperature,
            SVector(amount),
        )
        return Clapeyron.FlashResult(
            state_pressure,
            temperature,
            [pure_composition],
            [amount],
            [bulk_molar_volume],
            zero(state_pressure);
            sort = false,
            vapour_phase_index = -1,
        )
    elseif bulk_molar_volume > vapor_molar_volume
        state_pressure = pressure(
            model,
            volume,
            temperature,
            SVector(amount),
        )
        return Clapeyron.FlashResult(
            state_pressure,
            temperature,
            [pure_composition],
            [amount],
            [bulk_molar_volume],
            zero(state_pressure);
            sort = false,
            vapour_phase_index = 1,
        )
    else
        vapor_quality =
            (bulk_molar_volume - liquid_molar_volume) /
            (vapor_molar_volume - liquid_molar_volume)
        return Clapeyron.FlashResult(
            saturation_pressure_value,
            temperature,
            [pure_composition, pure_composition],
            [amount * (1 - vapor_quality), amount * vapor_quality],
            [liquid_molar_volume, vapor_molar_volume],
            zero(saturation_pressure_value);
            sort = false,
            vapour_phase_index = 2,
        )
    end
end

function flash_vapor_quality(flash)
    vapor_phase = flash.data.vapour_idx
    total_phase_amount = sum(flash.fractions)

    if vapor_phase == -1
        return zero(total_phase_amount)
    elseif vapor_phase > 0
        return flash.fractions[vapor_phase] / total_phase_amount
    elseif length(flash.fractions) == 1
        # A single supercritical/otherwise unidentified phase has no distinct
        # liquid fraction, so treat it as the homogeneous outlet phase.
        return one(total_phase_amount)
    else
        # Fallback for an unidentified multiphase result. For pure nitrous
        # oxide, the phase with the largest molar volume is the vapor phase.
        vapor_phase = argmax(flash.volumes)
        return flash.fractions[vapor_phase] / total_phase_amount
    end
end

function phase_weighted_mass_property(property_function, model, flash)
    total_phase_amount = sum(flash.fractions)
    weighted_property = zero(property_function(model, flash, firstindex(flash.fractions)))

    # This mole-fraction weighting is also mass-fraction weighting because the
    # feed model contains only one species: nitrous oxide.
    for phase in eachindex(flash.fractions)
        phase_fraction = flash.fractions[phase] / total_phase_amount
        weighted_property += phase_fraction * property_function(model, flash, phase)
    end

    return weighted_property
end

function update_tank_properties!(du, u, p, t, system, geo, cell_id)
    model = system.nitrous_oxide_model

    vol = geo.cell_volumes[cell_id]

    tank_n_moles = u.mass[cell_id] / u.molecular_weights.nitrous_oxide[cell_id]

    tank_flash = pure_component_vt_flash(
        model,
        vol,
        u.temperature[cell_id],
        tank_n_moles,
    )
    u.pressure[cell_id] = pressure(tank_flash)
    u.vapor_quality[cell_id] = flash_vapor_quality(tank_flash)

    # These describe the bulk material stored in the rigid tank. They are
    # intentionally independent of which phase is being withdrawn.
    u.specific_internal_energy[cell_id] = mass_internal_energy(model, tank_flash)
    u.cp[cell_id] = phase_weighted_mass_property(mass_isobaric_heat_capacity, model, tank_flash)
    u.cv[cell_id] = phase_weighted_mass_property(mass_isochoric_heat_capacity, model, tank_flash)

    # These describe the material leaving the bottom of the tank. Withdraw
    # liquid while a distinct liquid inventory remains, then withdraw vapor.
    if u.vapor_quality[cell_id] <= 0.999
        outlet_phase = argmin(tank_flash.volumes)
    else
        vapor_phase = tank_flash.data.vapour_idx
        outlet_phase = vapor_phase > 0 ? vapor_phase : argmax(tank_flash.volumes)
    end
    u.density[cell_id] = mass_density(model, tank_flash, outlet_phase)
    u.specific_enthalpy[cell_id] = mass_enthalpy(model, tank_flash, outlet_phase)
end

function update_mid_section_properties!(du, u, p, t, system, geo, cell_id)
    model = system.nitrous_oxide_model

    vol = geo.cell_volumes[cell_id]

    mid_section_n_moles = u.mass[cell_id] / u.molecular_weights.nitrous_oxide[cell_id]

    mid_section_flash = pure_component_vt_flash(
        model,
        vol,
        u.temperature[cell_id],
        mid_section_n_moles,
    )
    u.pressure[cell_id] = pressure(mid_section_flash)

    u.vapor_quality[cell_id] = flash_vapor_quality(mid_section_flash)
    u.density[cell_id] = mass_density(model, mid_section_flash)
    u.specific_internal_energy[cell_id] = mass_internal_energy(model, mid_section_flash)
    u.specific_enthalpy[cell_id] = mass_enthalpy(model, mid_section_flash)
    u.cp[cell_id] = phase_weighted_mass_property(mass_isobaric_heat_capacity, model, mid_section_flash)
    u.cv[cell_id] = phase_weighted_mass_property(mass_isochoric_heat_capacity, model, mid_section_flash)
end

function adjustable_valve_flow!(
    du, u, p, t, system, geo,
    tank_cell_id, face_a, mid_section_cell_id, face_b,
)
    if u.pressure[tank_cell_id] - u.pressure[mid_section_cell_id] <= 0
        return 0.0
    end
    
    #Notice that valve_opening uses system.valve_opening, this is because we're going to be adjusting this over time
    volumetric_flow = p.valve_opening * p.valve_flow_capacity_factor * sqrt((u.pressure[tank_cell_id] - u.pressure[mid_section_cell_id]) / u.density[tank_cell_id])

    mass_flow = volumetric_flow * u.density[tank_cell_id]

    tank_internal_energy_rate = -mass_flow * u.specific_enthalpy[tank_cell_id]
    mid_section_internal_energy_rate = mass_flow * u.specific_enthalpy[tank_cell_id]

    du.mass[tank_cell_id] += -mass_flow
    du.mass[mid_section_cell_id] += mass_flow

    du.temperature[tank_cell_id] += tank_internal_energy_rate / (u.mass[tank_cell_id] * u.cv[tank_cell_id])
    du.temperature[mid_section_cell_id] += mid_section_internal_energy_rate / (u.mass[mid_section_cell_id] * u.cv[mid_section_cell_id])

    return mass_flow
end

function injector_valve_flow!(
    du, u, p, t, system, geo,
    mid_section_cell_id, face_a, pre_chamber_cell_id, face_b,
)
    post_injector_cell_ids = system.post_injector_cell_ids
    chamber_pressure = sum([u.pressure[cell_id] for cell_id in post_injector_cell_ids]) / length(post_injector_cell_ids)

    if u.pressure[mid_section_cell_id] - chamber_pressure <= 0
        return 0.0
    end

    volumetric_flow = p.injector_discharge_coefficient * p.injector_orifice_area * sqrt((2 * (u.pressure[mid_section_cell_id] - chamber_pressure)) / u.density[mid_section_cell_id])

    oxidizer_mass_flow = volumetric_flow * u.density[mid_section_cell_id]

    mid_section_internal_energy_rate = -oxidizer_mass_flow * u.specific_enthalpy[mid_section_cell_id]
    pre_chamber_internal_energy_rate = oxidizer_mass_flow * u.specific_enthalpy[mid_section_cell_id]

    du.mass[mid_section_cell_id] += -oxidizer_mass_flow
    du.mass[pre_chamber_cell_id] += oxidizer_mass_flow

    du.temperature[mid_section_cell_id] += mid_section_internal_energy_rate / (u.mass[mid_section_cell_id] * u.cv[mid_section_cell_id])
    du.temperature[pre_chamber_cell_id] += pre_chamber_internal_energy_rate / (u.mass[pre_chamber_cell_id] * u.cv[pre_chamber_cell_id])

    u.injector_velocity[pre_chamber_cell_id] = oxidizer_mass_flow / (u.density[mid_section_cell_id] * p.injector_orifice_area)
end