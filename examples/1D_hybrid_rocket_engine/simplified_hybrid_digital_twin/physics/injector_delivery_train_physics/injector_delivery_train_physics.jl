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

"""
Mark a nitrous oxide cell's properties as invalid (NaN) instead of flashing a non-physical state.

The solver evaluates trial states it may later reject. If a trial step drains a cell to zero or negative mass,
pure_component_vt_flash would throw and end the whole solve. Writing NaN instead makes every flux that uses this
cell NaN, so the step's error estimate is NaN and the solver rejects it and retries with a smaller step.
Nothing is clamped, so an invalid state can never be carried forward silently.
"""
function set_nitrous_oxide_properties_invalid!(u, cell_id)
    not_a_number = oftype(u.density[cell_id], NaN)

    u.pressure[cell_id] = not_a_number
    u.nitrous_oxide_vapor_quality[cell_id] = not_a_number
    u.outflow_density[cell_id] = not_a_number
    u.specific_internal_energy[cell_id] = not_a_number
    u.specific_enthalpy[cell_id] = not_a_number
    u.cp[cell_id] = not_a_number
    u.cv[cell_id] = not_a_number
end

function update_tank_properties!(du, u, p, t, system, geo, cell_id)
    model = system.nitrous_oxide_model

    vol = geo.cell_volumes[cell_id]

    tank_mass = u.density[cell_id] * vol

    tank_n_moles = tank_mass / u.molecular_weights.nitrous_oxide[cell_id]

    # A rejected trial step can drain the cell; flag it instead of letting the flash throw.
    # Written as !(x > 0) so NaN is caught too (NaN <= 0 is false).
    if !(tank_n_moles > 0) || !(u.temperature[cell_id] > 0)
        set_nitrous_oxide_properties_invalid!(u, cell_id)
        return nothing
    end

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

    # The mid section holds under a milligram of vapor before the valve opens, so a trial step that's too long can
    # empty it. Flag that state instead of letting the flash throw, so the solver rejects the step and shortens it.
    # Written as !(x > 0) so NaN is caught too (NaN <= 0 is false).
    if !(mid_section_n_moles > 0) || !(u.temperature[cell_id] > 0)
        set_nitrous_oxide_properties_invalid!(u, cell_id)
        return nothing
    end

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

"""
Square root of an orifice pressure drop, with a finite slope at zero:

    sqrt_regularized(ΔP) = ΔP / (ΔP^2 + ΔP_reg^2)^(1/4)

ΔP: pressure drop across the orifice (Pa); ΔP_reg: regularization pressure (Pa).

Orifice laws give ṁ ∝ sqrt(ΔP), whose slope is infinite at ΔP = 0. That is exactly where the feed train sits before
the valve opens and wherever two cells equalize, and an infinite slope in the Jacobian stalls the implicit solver's
Newton iterations. This form equals sqrt(ΔP) to within 0.5% once ΔP > 7 ΔP_reg, and is linear below ΔP_reg.
With ΔP_reg = 1 kPa it is exact for every real operating point (pressure drops of bar to tens of bar).
"""
function regularized_square_root_of_pressure_drop(pressure_drop)
    regularization_pressure = 1000.0

    return pressure_drop / (pressure_drop^2 + regularization_pressure^2)^0.25
end

"""
Fraction the tank valve is open at time t: closed until valve_open_time, then valve_opening.

Computed from t on every call rather than stored in p. Writing it into p from inside the RHS (the old
open_injector_valve!) overwrote the valve_opening parameter, and because the valve patch runs before the injector
patch, the valve always saw the value set during the previous RHS call. The jump at valve_open_time is in tstops.
"""
function tank_valve_opening_fraction(p, t)
    if t > p.valve_open_time
        return p.valve_opening
    else
        return zero(p.valve_opening)
    end
end

function adjustable_valve_flow!(
    du, u, p, t, system, geo,
    tank_cell_id, face_a, mid_section_cell_id, face_b,
)
    valve_pressure_drop = u.pressure[tank_cell_id] - u.pressure[mid_section_cell_id]

    # No backflow into the tank: the feed train only ever carries nitrous from the tank toward the chamber
    if valve_pressure_drop <= 0
        return 0.0
    end

    #The flow is set by the phase actually leaving the tank (liquid while any remains), not the tank's bulk density
    #    ṁ = opening * K_v * sqrt(ρ_out * ΔP)       (from Q = opening * K_v * sqrt(ΔP / ρ_out) and ṁ = ρ_out Q)
    #ṁ: mass flow (kg/s); K_v: valve_flow_capacity_factor (m^2); ρ_out: density of the phase leaving the tank (kg/m^3)
    tank_outflow_density = u.outflow_density[tank_cell_id]

    valve_opening = tank_valve_opening_fraction(p, t)

    mass_flow = valve_opening * p.valve_flow_capacity_factor * sqrt(tank_outflow_density) * regularized_square_root_of_pressure_drop(valve_pressure_drop)

    #Old version (exact sqrt has an infinite slope at ΔP = 0, and the opening came from a value written into p):
    #volumetric_flow = p.valve_opening * p.valve_flow_capacity_factor * sqrt((u.pressure[tank_cell_id] - u.pressure[mid_section_cell_id]) / u.outflow_density[tank_cell_id])
    #mass_flow = volumetric_flow * u.outflow_density[tank_cell_id]

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

"""
Nitrous flow from the tank, through the valve and the injector in series, into the pre-chamber.

The mid section between the valve and the injector holds only 0.47 cm^3. Modelled as its own cell, it filled in
about half a millisecond after the valve opened and then sat at the edge of the vapor dome, forcing 1e-8 to 1e-7 s
steps with repeated Newton failures. Its dynamics last microseconds, the burn lasts seconds, and none of it shows
up in thrust, chamber pressure or O/F. So the two orifices are treated as quasi-steady and in series instead:
the same mass flow passes both, and the mid-section pressure is whatever makes that true.

    valve:      ΔP_valve    = ṁ^2 / (K_v^2 ρ)               (from ṁ = K_v sqrt(ρ ΔP_valve))
    injector:   ΔP_injector = ṁ^2 / (2 ρ (C_d A)^2)         (from ṁ = C_d A sqrt(2 ρ ΔP_injector))
    together:   ΔP_total = ΔP_valve + ΔP_injector
         →      ṁ = K_v * sqrt(2) * C_d A * sqrt(ρ ΔP_total) / sqrt(2 (C_d A)^2 + K_v^2)

ṁ: oxidizer mass flow (kg/s); K_v: valve_flow_capacity_factor times the opening fraction (m^2), zero while the valve
is closed, which gives ṁ = 0 without dividing by zero; C_d A: injector discharge coefficient times orifice area (m^2);
ρ: density of the phase leaving the tank (kg/m^3), assumed not to change through the line (incompressible orifices);
ΔP_total: tank pressure minus pre-chamber pressure (Pa).

The flow is isenthalpic (throttling), so the stream arrives in the pre-chamber with the tank outflow's enthalpy.
The mid-section cell stays in the geometry with frozen states; its pressure cache is overwritten with the series
pressure P_mid = P_pre-chamber + ΔP_injector so it can still be plotted and compared with a pressure transducer.
"""
function injector_valve_flow!(
    du, u, p, t, system, geo,
    mid_section_cell_id, face_a, pre_chamber_cell_id, face_b,
)

    # Reset the diagnostic first, so it reads zero (not a stale value) whenever there is no injection
    u.neglected_injector_cooling_power[pre_chamber_cell_id] = 0.0

    tank_cell_id = get_cellset(geo, "tank")[1]

    chamber_pressure = u.pressure[pre_chamber_cell_id]

    feed_pressure_drop = u.pressure[tank_cell_id] - chamber_pressure

    # The mid section sits at chamber pressure while nothing flows (until the valve flow below overwrites it)
    u.pressure[mid_section_cell_id] = chamber_pressure

    # No backflow: chamber gas never flows back into the pure-nitrous feed train
    if feed_pressure_drop <= 0
        return 0.0
    end

    tank_outflow_density = u.outflow_density[tank_cell_id]

    injector_orifice_area = system.cad_data.injector_orifice_area

    injector_effective_area = p.injector_discharge_coefficient * injector_orifice_area

    valve_effective_area = tank_valve_opening_fraction(p, t) * p.valve_flow_capacity_factor

    oxidizer_mass_flow = valve_effective_area * sqrt(2.0) * injector_effective_area * sqrt(tank_outflow_density) * regularized_square_root_of_pressure_drop(feed_pressure_drop) / sqrt(2 * injector_effective_area^2 + valve_effective_area^2)

    # Mid-section pressure diagnostic: chamber pressure plus the injector's share of the pressure drop
    injector_pressure_drop = oxidizer_mass_flow^2 / (2 * tank_outflow_density * injector_effective_area^2)

    u.pressure[mid_section_cell_id] = chamber_pressure + injector_pressure_drop

    #Old version (the mid section was its own dynamic cell, fed by adjustable_valve_flow! and drained here):
    #injector_pressure_drop = u.pressure[mid_section_cell_id] - chamber_pressure
    #oxidizer_mass_flow = p.injector_discharge_coefficient * injector_orifice_area * sqrt(2 * mid_section_outflow_density) * regularized_square_root_of_pressure_drop(injector_pressure_drop)

    # The fuel regression law needs the injected oxidizer flow (G_ox = ṁ_ox / port area), so it's stored here.
    # It stays zero, as the cache is reset every call, whenever the early return above stops the injection.
    u.injector_oxidizer_mass_flow[pre_chamber_cell_id] = oxidizer_mass_flow

    tank_volume = geo.cell_volumes[tank_cell_id]
    pre_chamber_volume = geo.cell_volumes[pre_chamber_cell_id]

    tank_mass = u.density[tank_cell_id] * tank_volume
    pre_chamber_mass = u.density[pre_chamber_cell_id] * pre_chamber_volume

    # Mass conservation, written for the density states of two fixed-volume cells
    du.density[tank_cell_id] += -oxidizer_mass_flow / tank_volume
    du.density[pre_chamber_cell_id] += oxidizer_mass_flow / pre_chamber_volume


    # Enthalpy of the injected nitrous on the chamber's sensible basis (h = cp T for gas).
    # The pre-chamber works on the ideal-gas sensible basis used by the chamber fluxes (h = cp T, e = cv T),
    # while the mid section uses Clapeyron's absolute values. The two references don't match, so the
    # injected enthalpy is moved onto the chamber basis by splitting it into an ideal-gas part and a residual:
    #     h_injected = cp_N2O,ideal(T_tank) * T_tank + (h_tank,out - h_N2O,ideal(T_tank))
    # The residual (h - h_ideal) doesn't depend on the reference state. For liquid nitrous it is strongly
    # negative (about -236 kJ/kg for saturated liquid at 293 K with PR): it is the energy the injected
    # liquid has to absorb to become gas. Throttling through the valve and injector keeps the enthalpy
    # constant, so the tank outflow's enthalpy is what arrives.
    tank_outflow_specific_enthalpy = u.specific_enthalpy[tank_cell_id]

    tank_temperature = u.temperature[tank_cell_id]
    tank_pressure = u.pressure[tank_cell_id]

    nitrous_oxide_molar_mass = u.molecular_weights.nitrous_oxide[tank_cell_id]

    nitrous_oxide_ideal_gas_model = Clapeyron.idealmodel(system.nitrous_oxide_model)

    clapeyron_ideal_gas_specific_enthalpy = enthalpy(nitrous_oxide_ideal_gas_model, tank_pressure, tank_temperature) / nitrous_oxide_molar_mass

    injected_residual_specific_enthalpy = tank_outflow_specific_enthalpy - clapeyron_ideal_gas_specific_enthalpy

    # Pure nitrous oxide is the first species in system.model (ReidIdeal); 1e5 Pa is a placeholder since ideal-gas cp doesn't depend on pressure.
    # One entry per species in system.model (7 since carbon monoxide was added); a short vector reads past its end inside Clapeyron.
    pure_nitrous_oxide_composition = SVector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    nitrous_oxide_ideal_gas_cp = mass_isobaric_heat_capacity(system.model, 1e5, tank_temperature, pure_nitrous_oxide_composition)

    injected_sensible_specific_enthalpy = nitrous_oxide_ideal_gas_cp * tank_temperature + injected_residual_specific_enthalpy


    # Species. All injected nitrous oxide joins the pre-chamber gas immediately (instant vaporization):
    #     dY_i/dt = (Y_i,injected - Y_i,pre_chamber) * ṁ / M_pre_chamber
    # Y_injected is 1 for nitrous oxide and 0 for everything else. This is applied as -Y_i * ṁ / M
    # for every species, plus ṁ / M for nitrous oxide.
    for_fields!(du.mass_fractions, u.mass_fractions) do species, du_mass_fractions, u_mass_fractions
        du_mass_fractions[species[pre_chamber_cell_id]] -= u_mass_fractions[species[pre_chamber_cell_id]] * oxidizer_mass_flow / pre_chamber_mass
    end

    du.mass_fractions.nitrous_oxide[pre_chamber_cell_id] += oxidizer_mass_flow / pre_chamber_mass


    # Energy of the donor (tank):
    #     M de/dt = -ṁ (h_out - e_tank)
    # with de = cv dT, where cv is the two-phase cv from the tank property update
    tank_specific_internal_energy = u.specific_internal_energy[tank_cell_id]

    tank_internal_energy_rate = -oxidizer_mass_flow * (tank_outflow_specific_enthalpy - tank_specific_internal_energy)

    du.temperature[tank_cell_id] += tank_internal_energy_rate / (tank_mass * u.cv[tank_cell_id])


    # Energy of the receiver (pre-chamber):
    #     M de/dt = +ṁ (h_injected - e_pre_chamber)
    # The liquid is vaporized instantly, so its whole latent heat comes out of the pre-chamber gas.
    pre_chamber_temperature = u.temperature[pre_chamber_cell_id]

    pre_chamber_sensible_specific_internal_energy = u.cv[pre_chamber_cell_id] * pre_chamber_temperature

    # The injected enthalpy is split into the stream arriving as vapor at the tank temperature and the residual
    # (latent heat plus the real-gas term) that it takes to get there:
    #     ṁ (h_injected - e_pre) = ṁ (cp_N2O T_tank - e_pre)   +   ṁ h_residual
    #                              vapor inflow (always applied)   latent cooling (floored below)
    vapor_inflow_power = oxidizer_mass_flow * (nitrous_oxide_ideal_gas_cp * tank_temperature - pre_chamber_sensible_specific_internal_energy)

    latent_cooling_power = oxidizer_mass_flow * injected_residual_specific_enthalpy

    # Temperature floor. The pre-chamber is a single ideal-gas cell with no fuel, so instant vaporization alone
    # would drag it far below anything physical (a real head end holds a cold spray near T_sat, roughly
    # 230-260 K at chamber pressure). Below 240 K the latent cooling is ramped off, reaching zero at 230 K.
    #
    # Only the latent part is floored. The earlier version cut the whole injector energy term, which also removed
    # the heating from nitrous being pushed into the cell; the cell then kept venting into the port with nothing
    # balancing the flow-work cooling and fell to ~47 K, and the solver crawled. With the vapor inflow kept, steady
    # injection with no latent cooling would hold the cell near the tank temperature, so the floor can hold.
    #
    # This is a deliberate energy shortcut: the latent heat that is cut off is simply not charged. Once the
    # floor is active the model is treating some liquid as if it vaporized for free, which biases the chamber
    # slightly hot (the full latent heat is about 4% of the burned propellant's energy). Heating is never limited.
    minimum_pre_chamber_temperature = 230.0
    floor_ramp_width = 10.0

    if latent_cooling_power < 0
        ramp_position = clamp((pre_chamber_temperature - minimum_pre_chamber_temperature) / floor_ramp_width, 0.0, 1.0)

        # Smoothstep: 0 at the floor, 1 at floor + ramp width, with zero slope at both ends so the solver sees no kink
        cooling_fraction = ramp_position^2 * (3.0 - 2.0 * ramp_position)

        # TODO: this is energy the model is not accounting for. Track it so it stays visible until the floor is replaced.
        # The power is what the full latent cooling would have been minus what was actually applied (W)
        neglected_cooling_power = -(1.0 - cooling_fraction) * latent_cooling_power

        u.neglected_injector_cooling_power[pre_chamber_cell_id] = neglected_cooling_power

        latent_cooling_power = cooling_fraction * latent_cooling_power
    end

    pre_chamber_internal_energy_rate = vapor_inflow_power + latent_cooling_power

    #Old version (floored the whole injector energy term, including the vapor inflow's heating):
    #pre_chamber_internal_energy_rate = oxidizer_mass_flow * (injected_sensible_specific_enthalpy - pre_chamber_sensible_specific_internal_energy)
    #then pre_chamber_internal_energy_rate = cooling_fraction * pre_chamber_internal_energy_rate whenever it was negative

    du.temperature[pre_chamber_cell_id] += pre_chamber_internal_energy_rate / (pre_chamber_mass * u.cv[pre_chamber_cell_id])

    u.injector_velocity[pre_chamber_cell_id] = oxidizer_mass_flow / (tank_outflow_density * injector_orifice_area)
end