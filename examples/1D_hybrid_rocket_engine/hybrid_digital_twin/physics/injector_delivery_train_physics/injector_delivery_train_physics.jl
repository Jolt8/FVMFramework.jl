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
    volume > zero(volume) || throw(DomainError(volume, "flash volume must be positive"))
    amount > zero(amount) || throw(DomainError(amount, "flash amount must be positive"))

    saturation_pressure_value, liquid_molar_volume, vapor_molar_volume = saturation_pressure(model, temperature)
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
     
"""
True isochoric heat capacity of a rigid cell of pure nitrous oxide, valid on both
sides of the saturation dome:

    cv = (de/dT) at fixed total volume and amount

e: specific internal energy of the cell contents (J/kg); T: temperature (K).

Inside the dome this includes the latent-heat term (e_vapor - e_liquid) * dx/dT,
because warming a saturated mixture at fixed volume also changes its vapor quality x.
The phase-weighted average of single-phase cvs misses that term and is 3.7-4.5x too
low for the 1.134 kg / 1.6 L tank fill between 273 K and 293 K.

Computed with a central difference instead of ForwardDiff: it only uses arithmetic on
the flash results, so it works unchanged when the solver passes Dual numbers or sparsity
tracers into the RHS, without nesting a second ForwardDiff call inside the solver's own.
"""
function two_phase_mass_isochoric_heat_capacity(model, volume, temperature, amount)
    temperature_step = 0.01

    hot_flash = pure_component_vt_flash(model, volume, temperature + temperature_step, amount)
    cold_flash = pure_component_vt_flash(model, volume, temperature - temperature_step, amount)

    hot_specific_internal_energy = mass_internal_energy(model, hot_flash)
    cold_specific_internal_energy = mass_internal_energy(model, cold_flash)

    return (hot_specific_internal_energy - cold_specific_internal_energy) / (2 * temperature_step)
end

function update_tank_properties!(du, u, p, t, system, geo, cell_id)
    model = system.nitrous_oxide_model

    vol = geo.cell_volumes[cell_id]
    
    tank_mass = u.density[cell_id] * vol

    tank_n_moles = tank_mass / u.molecular_weights.nitrous_oxide[cell_id]

    tank_flash = pure_component_vt_flash(
        model,
        vol,
        u.temperature[cell_id],
        tank_n_moles,
    )
    u.pressure[cell_id] = pressure(tank_flash)
    u.nitrous_oxide_vapor_quality[cell_id] = flash_vapor_quality(tank_flash)

    # These describe the bulk material stored in the rigid tank. They are
    # intentionally independent of which phase is being withdrawn.
    u.specific_internal_energy[cell_id] = mass_internal_energy(model, tank_flash)

    # cv includes the latent-heat term inside the dome. cp has no finite two-phase value
    # (at constant pressure, heat only boils liquid), so the phase-weighted cp is kept only as
    # a placeholder; the temperature equations for this cell use cv.
    u.cp[cell_id] = phase_weighted_mass_property(mass_isobaric_heat_capacity, model, tank_flash)
    u.cv[cell_id] = two_phase_mass_isochoric_heat_capacity(model, vol, u.temperature[cell_id], tank_n_moles)

    #u.cv[cell_id] = phase_weighted_mass_property(mass_isochoric_heat_capacity, model, tank_flash) #old version, misses the latent-heat term

    # These describe the material leaving the bottom of the tank. Withdraw
    # liquid while a distinct liquid inventory remains, then withdraw vapor.
    if u.nitrous_oxide_vapor_quality[cell_id] <= 0.999
        outlet_phase = argmin(tank_flash.volumes)
    else
        vapor_phase = tank_flash.data.vapour_idx
        outlet_phase = vapor_phase > 0 ? vapor_phase : argmax(tank_flash.volumes)
    end

    # The outlet-phase density goes in its own cache. Writing it to u.density would overwrite
    # the density state (total tank inventory) behind the solver's back.
    u.outflow_density[cell_id] = mass_density(model, tank_flash, outlet_phase)
    u.specific_enthalpy[cell_id] = mass_enthalpy(model, tank_flash, outlet_phase)
end

function update_mid_section_properties!(du, u, p, t, system, geo, cell_id)
    model = system.nitrous_oxide_model

    vol = geo.cell_volumes[cell_id]

    mid_section_mass = u.density[cell_id] * vol

    mid_section_n_moles = mid_section_mass / u.molecular_weights.nitrous_oxide[cell_id]

    mid_section_flash = pure_component_vt_flash(
        model,
        vol,
        u.temperature[cell_id],
        mid_section_n_moles,
    )
    u.pressure[cell_id] = pressure(mid_section_flash)

    u.nitrous_oxide_vapor_quality[cell_id] = flash_vapor_quality(mid_section_flash)

    # The mid section is small and treated as well mixed, so what leaves is the bulk mixture
    u.outflow_density[cell_id] = mass_density(model, mid_section_flash)
    
    u.specific_internal_energy[cell_id] = mass_internal_energy(model, mid_section_flash)
    u.specific_enthalpy[cell_id] = mass_enthalpy(model, mid_section_flash)

    # Same treatment as the tank: cv includes the latent-heat term, cp is a placeholder inside the dome
    u.cp[cell_id] = phase_weighted_mass_property(mass_isobaric_heat_capacity, model, mid_section_flash)
    u.cv[cell_id] = two_phase_mass_isochoric_heat_capacity(model, vol, u.temperature[cell_id], mid_section_n_moles)

    #u.cv[cell_id] = phase_weighted_mass_property(mass_isochoric_heat_capacity, model, mid_section_flash) #old version, misses the latent-heat term
end

function adjustable_valve_flow!(
    du, u, p, t, system, geo,
    tank_cell_id, face_a, mid_section_cell_id, face_b,
)
    if u.pressure[tank_cell_id] - u.pressure[mid_section_cell_id] <= 0
        return 0.0
    end
    
    #Notice that valve_opening uses system.valve_opening, this is because we're going to be adjusting this over time
    #The flow is set by the phase actually leaving the tank (liquid while any remains), not the tank's bulk density

    volumetric_flow = p.valve_opening * p.valve_flow_capacity_factor * sqrt((u.pressure[tank_cell_id] - u.pressure[mid_section_cell_id]) / u.outflow_density[tank_cell_id])

    mass_flow = volumetric_flow * u.outflow_density[tank_cell_id]

    tank_volume = geo.cell_volumes[tank_cell_id]
    mid_section_volume = geo.cell_volumes[mid_section_cell_id]

    tank_mass = u.density[tank_cell_id] * tank_volume
    mid_section_mass = u.density[mid_section_cell_id] * mid_section_volume

    # Mass conservation, written for the density states of two fixed-volume cells
    du.density[tank_cell_id] += -mass_flow / tank_volume
    du.density[mid_section_cell_id] += mass_flow / mid_section_volume

    # Energy. Both cells use Clapeyron's absolute energies from the same nitrous model, so the reference state cancels.
    # Donor (tank):     d(M e)/dt = -ṁ h_out, dM/dt = -ṁ   ->   M de/dt = -ṁ (h_out - e_tank)
    # Receiver (mid):   d(M e)/dt = +ṁ h_out, dM/dt = +ṁ   ->   M de/dt = +ṁ (h_out - e_mid)
    # with de = cv dT, where cv is the two-phase cv from the property update

    tank_internal_energy_rate = -mass_flow * (u.specific_enthalpy[tank_cell_id] - u.specific_internal_energy[tank_cell_id])
    mid_section_internal_energy_rate = mass_flow * (u.specific_enthalpy[tank_cell_id] - u.specific_internal_energy[mid_section_cell_id])

    du.temperature[tank_cell_id] += tank_internal_energy_rate / (tank_mass * u.cv[tank_cell_id])
    du.temperature[mid_section_cell_id] += mid_section_internal_energy_rate / (mid_section_mass * u.cv[mid_section_cell_id])

    return mass_flow
end

function open_injector_valve!(du, u, p, t, system, geo)
    if t > p.valve_open_time
        p.valve_opening = 1.0
    else
        p.valve_opening = 0.0
    end
end

function injector_valve_flow!(
    du, u, p, t, system, geo,
    mid_section_cell_id, face_a, pre_chamber_cell_id, face_b,
)
    open_injector_valve!(du, u, p, t, system, geo)
    
    # Reset the diagnostic first, so it reads zero (not a stale value) whenever there is no injection
    u.neglected_injector_cooling_power[pre_chamber_cell_id] = 0.0

    post_injector_cell_ids = get_cellset(geo, "pre_chamber")[1]
    chamber_pressure = sum([u.pressure[cell_id] for cell_id in post_injector_cell_ids]) / length(post_injector_cell_ids)

    if u.pressure[mid_section_cell_id] - chamber_pressure <= 0
        return 0.0
    end

    mid_section_outflow_density = u.outflow_density[mid_section_cell_id]

    injector_orifice_area = system.cad_data.injector_orifice_area

    volumetric_flow = p.injector_discharge_coefficient * injector_orifice_area * sqrt((2 * (u.pressure[mid_section_cell_id] - chamber_pressure)) / mid_section_outflow_density)

    oxidizer_mass_flow = volumetric_flow * mid_section_outflow_density

    # The fuel regression law needs the injected oxidizer flow (G_ox = ṁ_ox / port area), so it's stored here.
    # It stays zero, as the cache is reset every call, whenever the early return above stops the injection.
    u.injector_oxidizer_mass_flow[pre_chamber_cell_id] = oxidizer_mass_flow

    mid_section_volume = geo.cell_volumes[mid_section_cell_id]
    pre_chamber_volume = geo.cell_volumes[pre_chamber_cell_id]

    mid_section_mass = u.density[mid_section_cell_id] * mid_section_volume
    pre_chamber_mass = u.density[pre_chamber_cell_id] * pre_chamber_volume

    # Mass conservation, written for the density states of two fixed-volume cells
    du.density[mid_section_cell_id] += -oxidizer_mass_flow / mid_section_volume
    du.density[pre_chamber_cell_id] += oxidizer_mass_flow / pre_chamber_volume


    # Enthalpy of the injected nitrous on the chamber's sensible basis (h = cp T for gas).
    # The pre-chamber works on the ideal-gas sensible basis used by the chamber fluxes (h = cp T, e = cv T),
    # while the mid section uses Clapeyron's absolute values. The two references don't match, so the
    # injected enthalpy is moved onto the chamber basis by splitting it into an ideal-gas part and a residual:
    #     h_injected = cp_N2O,ideal(T_mid) * T_mid + (h_mid - h_N2O,ideal(T_mid))
    # The residual (h_mid - h_ideal) doesn't depend on the reference state. For liquid nitrous it is strongly
    # negative (about -236 kJ/kg for saturated liquid at 293 K with PR): it is the energy the injected
    # liquid has to absorb to become gas.
    mid_section_specific_enthalpy = u.specific_enthalpy[mid_section_cell_id]

    mid_section_temperature = u.temperature[mid_section_cell_id]
    mid_section_pressure = u.pressure[mid_section_cell_id]

    nitrous_oxide_molar_mass = u.molecular_weights.nitrous_oxide[mid_section_cell_id]

    nitrous_oxide_ideal_gas_model = Clapeyron.idealmodel(system.nitrous_oxide_model)

    clapeyron_ideal_gas_specific_enthalpy = enthalpy(nitrous_oxide_ideal_gas_model, mid_section_pressure, mid_section_temperature) / nitrous_oxide_molar_mass

    injected_residual_specific_enthalpy = mid_section_specific_enthalpy - clapeyron_ideal_gas_specific_enthalpy

    # Pure nitrous oxide is the first species in system.model (ReidIdeal); 1e5 Pa is a placeholder since ideal-gas cp doesn't depend on pressure.
    # One entry per species in system.model (7 since carbon monoxide was added); a short vector reads past its end inside Clapeyron.
    pure_nitrous_oxide_composition = SVector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    nitrous_oxide_ideal_gas_cp = mass_isobaric_heat_capacity(system.model, 1e5, mid_section_temperature, pure_nitrous_oxide_composition)

    injected_sensible_specific_enthalpy = nitrous_oxide_ideal_gas_cp * mid_section_temperature + injected_residual_specific_enthalpy


    # Species. All injected nitrous oxide joins the pre-chamber gas immediately (instant vaporization):
    #     dY_i/dt = (Y_i,injected - Y_i,pre_chamber) * ṁ / M_pre_chamber
    # Y_injected is 1 for nitrous oxide and 0 for everything else. This is applied as -Y_i * ṁ / M
    # for every species, plus ṁ / M for nitrous oxide.
    for_fields!(du.mass_fractions, u.mass_fractions) do species, du_mass_fractions, u_mass_fractions
        du_mass_fractions[species[pre_chamber_cell_id]] -= u_mass_fractions[species[pre_chamber_cell_id]] * oxidizer_mass_flow / pre_chamber_mass
    end

    du.mass_fractions.nitrous_oxide[pre_chamber_cell_id] += oxidizer_mass_flow / pre_chamber_mass


    # Energy of the donor (mid section). Same derivation as the valve:
    #     M de/dt = -ṁ (h_mid - e_mid)
    mid_section_specific_internal_energy = u.specific_internal_energy[mid_section_cell_id]

    mid_section_internal_energy_rate = -oxidizer_mass_flow * (mid_section_specific_enthalpy - mid_section_specific_internal_energy)

    du.temperature[mid_section_cell_id] += mid_section_internal_energy_rate / (mid_section_mass * u.cv[mid_section_cell_id])


    # Energy of the receiver (pre-chamber):
    #     M de/dt = +ṁ (h_injected - e_pre_chamber)
    # The liquid is vaporized instantly, so its whole latent heat comes out of the pre-chamber gas.
    pre_chamber_temperature = u.temperature[pre_chamber_cell_id]

    pre_chamber_sensible_specific_internal_energy = u.cv[pre_chamber_cell_id] * pre_chamber_temperature

    pre_chamber_internal_energy_rate = oxidizer_mass_flow * (injected_sensible_specific_enthalpy - pre_chamber_sensible_specific_internal_energy)

    # Temperature floor. The pre-chamber is a single ideal-gas cell with no fuel, so instant vaporization alone
    # would drag it far below anything physical (a real head end holds a cold spray near T_sat, roughly
    # 230-260 K at chamber pressure). Below 240 K the cooling is ramped off, reaching zero at 230 K.
    #
    # This is a deliberate energy shortcut: the latent heat that is cut off is simply not charged. Once the
    # floor is active the model is treating some liquid as if it vaporized for free, which biases the chamber
    # slightly hot (the full latent heat is about 4% of the burned propellant's energy). Heating is never limited.
    minimum_pre_chamber_temperature = 230.0
    floor_ramp_width = 10.0

    if pre_chamber_internal_energy_rate < 0
        ramp_position = clamp((pre_chamber_temperature - minimum_pre_chamber_temperature) / floor_ramp_width, 0.0, 1.0)

        # Smoothstep: 0 at the floor, 1 at floor + ramp width, with zero slope at both ends so the solver sees no kink
        cooling_fraction = ramp_position^2 * (3.0 - 2.0 * ramp_position)

        # TODO: this is energy the model is not accounting for. Track it so it stays visible until the floor is replaced.
        # The power is what the full cooling would have been minus what was actually applied (W); 
        neglected_cooling_power = -(1.0 - cooling_fraction) * pre_chamber_internal_energy_rate

        u.neglected_injector_cooling_power[pre_chamber_cell_id] = neglected_cooling_power

        pre_chamber_internal_energy_rate = cooling_fraction * pre_chamber_internal_energy_rate
    end

    du.temperature[pre_chamber_cell_id] += pre_chamber_internal_energy_rate / (pre_chamber_mass * u.cv[pre_chamber_cell_id])

    injector_orifice_area = system.cad_data.injector_orifice_area

    u.injector_velocity[pre_chamber_cell_id] = oxidizer_mass_flow / (mid_section_outflow_density * injector_orifice_area)
end