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

    if isnan(saturation_pressure_value) || bulk_molar_volume < liquid_molar_volume || bulk_molar_volume > vapor_molar_volume
        state_pressure = pressure(model, volume, temperature, SVector(amount))
        is_liquid = !isnan(saturation_pressure_value) && bulk_molar_volume < liquid_molar_volume
        return Clapeyron.FlashResult(
            state_pressure, temperature,
            [pure_composition], [amount], [bulk_molar_volume], zero(state_pressure);
            sort = false,
            vapour_phase_index = is_liquid ? -1 : 1,
        )
    else
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
end

"""
    pure_component_ps_flash(model, pressure_value, molar_entropy)

Construct the equilibrium `FlashResult` for a subcritical, pure-component
state specified by pressure and molar entropy (J/(mol*K), one mole basis).

Same reasoning as `pure_component_vt_flash`: Clapeyron's generic `ps_flash`
has to compile its whole flash machinery for every new number type (Float64,
ForwardDiff.Dual, sparsity tracers), which takes practically forever. For a
pure fluid the PS flash reduces to:

- saturation temperature at the given pressure;
- the lever rule on entropy inside the dome;
- a 1D Newton solve on T (ds/dT = cp/T at fixed p) outside the dome.
"""
function pure_component_ps_flash(model, pressure_value, molar_entropy)
    saturation_temperature_value, liquid_molar_volume, vapor_molar_volume = saturation_temperature(model, pressure_value)
    pure_composition = SVector(one(saturation_temperature_value))

    liquid_molar_entropy = Clapeyron.VT_entropy(model, liquid_molar_volume, saturation_temperature_value)
    vapor_molar_entropy = Clapeyron.VT_entropy(model, vapor_molar_volume, saturation_temperature_value)

    if liquid_molar_entropy <= molar_entropy <= vapor_molar_entropy
        vapor_quality = (molar_entropy - liquid_molar_entropy) / (vapor_molar_entropy - liquid_molar_entropy)
        amount = one(vapor_quality)
        return Clapeyron.FlashResult(
            pressure_value,
            saturation_temperature_value,
            [pure_composition, pure_composition],
            [amount * (1 - vapor_quality), amount * vapor_quality],
            [liquid_molar_volume, vapor_molar_volume],
            zero(pressure_value);
            sort = false,
            vapour_phase_index = 2,
        )
    end

    phase, phase_entropy, vapour_phase_index = molar_entropy < liquid_molar_entropy ?
        (:liquid, liquid_molar_entropy, -1) :
        (:vapour, vapor_molar_entropy, 1)

    #Start from the saturation boundary, assuming constant cp: s - s_sat = cp * log(T / T_sat)
    saturation_cp = isobaric_heat_capacity(model, pressure_value, saturation_temperature_value; phase = phase)
    temperature = saturation_temperature_value * exp((molar_entropy - phase_entropy) / saturation_cp)

    for _ in 1:20
        entropy_error = entropy(model, pressure_value, temperature; phase = phase) - molar_entropy
        temperature -= entropy_error * temperature / isobaric_heat_capacity(model, pressure_value, temperature; phase = phase)
        abs(entropy_error) < 1e-9 && break
    end

    molar_volume = volume(model, pressure_value, temperature; phase = phase)
    return Clapeyron.FlashResult(
        pressure_value,
        temperature,
        [pure_composition],
        [one(molar_volume)],
        [molar_volume],
        zero(pressure_value);
        sort = false,
        vapour_phase_index = vapour_phase_index,
    )
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

    u.specific_internal_energy[cell_id] = mass_internal_energy(model, tank_flash)

    u.cp[cell_id] = phase_weighted_mass_property(mass_isobaric_heat_capacity, model, tank_flash)
    u.cv[cell_id] = two_phase_mass_isochoric_heat_capacity(model, vol, u.temperature[cell_id], tank_n_moles)
    
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

function tank_valve_opening_fraction(p, t)
    if t > p.valve_open_time
        return p.valve_opening
    else
        return zero(p.valve_opening)
    end
end

function adjustable_valve_flow(
    mid_section_pressure,
    du, u, p, t, system, geo,
    tank_cell_id, pre_chamber_cell_id
)
    valve_pressure_drop = u.pressure[tank_cell_id] - mid_section_pressure

    if valve_pressure_drop <= 0
        return 0.0
    end

    volumetric_flow = tank_valve_opening_fraction(p, t) * p.valve_flow_capacity_factor * sqrt((valve_pressure_drop) / u.outflow_density[tank_cell_id])

    mass_flow = volumetric_flow * u.outflow_density[tank_cell_id]

    return mass_flow
end

function get_spi_mass_flow(
    mid_section_pressure,
    du, u, p, t, system, geo,
    tank_cell_id, pre_chamber_cell_id
)
    #nitrous_oxide_liquid_density = mass_density(system.nitrous_oxide_model, mid_section_pressure, u.temperature[tank_cell_id])

    spi_mass_flow = p.injector_discharge_coefficient * p.injector_area * sqrt(2 * u.outflow_density[tank_cell_id] * (mid_section_pressure - u.pressure[pre_chamber_cell_id]))

    return spi_mass_flow
end

function get_hem_mass_flow(
    mid_section_pressure,
    du, u, p, t, system, geo,
    tank_cell_id, pre_chamber_cell_id
)
    pre_injector_enthalpy = u.specific_enthalpy[tank_cell_id]
    injector_inlet_state = ph_flash(system.nitrous_oxide_model, mid_section_pressure, pre_injector_enthalpy * u.molecular_weights.nitrous_oxide[tank_cell_id])
    starting_molar_entropy = entropy(system.nitrous_oxide_model, injector_inlet_state)

    #The isentropic expansion goes from the injector inlet down towards the chamber pressure; the mass flux peaks where it chokes
    n_pressure_steps = 30
    largest_mass_flux = zero(pre_injector_enthalpy)

    for i in 1:n_pressure_steps
        chamber_pressure_test = mid_section_pressure + (u.pressure[pre_chamber_cell_id] - mid_section_pressure) * (i - 1) / (n_pressure_steps - 1)

        post_injector_result = pure_component_ps_flash(system.nitrous_oxide_model, chamber_pressure_test, starting_molar_entropy)
        post_injector_density = mass_density(system.nitrous_oxide_model, post_injector_result)
        post_injector_enthalpy = mass_enthalpy(system.nitrous_oxide_model, post_injector_result)

        #max guards the first step, where the enthalpy drop is zero up to rounding
        mass_flux = post_injector_density * sqrt(max(2 * (pre_injector_enthalpy - post_injector_enthalpy), zero(pre_injector_enthalpy)))

        if mass_flux < largest_mass_flux
            break
        end

        largest_mass_flux = mass_flux
    end

    return p.injector_discharge_coefficient * p.injector_area * largest_mass_flux
end

function get_nhne_injector_flow(
    mid_section_pressure,
    du, u, p, t, system, geo,
    tank_cell_id, pre_chamber_cell_id
)   
    spi_mass_flow = get_spi_mass_flow(
        mid_section_pressure,
        du, u, p, t, system, geo,
        tank_cell_id, pre_chamber_cell_id
    )

    hem_mass_flow = get_hem_mass_flow(
        mid_section_pressure,
        du, u, p, t, system, geo,
        tank_cell_id, pre_chamber_cell_id
    )

    non_eq_parameter = sqrt((mid_section_pressure - u.pressure[pre_chamber_cell_id]) / 
        (saturation_pressure(system.nitrous_oxide_model, u.temperature[tank_cell_id])[1] - u.pressure[pre_chamber_cell_id]))

    nhne_mass_flow = (non_eq_parameter / (1 + non_eq_parameter)) * spi_mass_flow + (1 / (1 + non_eq_parameter)) * hem_mass_flow

    return nhne_mass_flow
end


function series_valve_mass_flow_resid(
    mid_section_pressure,
    du, u, p, t, system, geo,
    tank_cell_id, pre_chamber_cell_id
)
    valve_mass_flow = adjustable_valve_flow(
        mid_section_pressure,
        du, u, p, t, system, geo,
        tank_cell_id, pre_chamber_cell_id
    )

    nhne_injector_mass_flow = get_nhne_injector_flow(
        mid_section_pressure,
        du, u, p, t, system, geo,
        tank_cell_id, pre_chamber_cell_id
    )

    return valve_mass_flow - nhne_injector_mass_flow
end

function injector_train_mass_flow!(
    du, u, p, t, system, geo,
    tank_cell_id, face_a, pre_chamber_cell_id, face_b
)
    tank_cell_id = get_cellset(geo, "tank")[1]
    pre_chamber_cell_id = get_cellset(geo, "pre_chamber")[1]

    @show u.temperature[pre_chamber_cell_id]
    @show u.pressure[pre_chamber_cell_id]
    @show u.density[pre_chamber_cell_id]
    @show u.pressure[tank_cell_id]

    true_mid_section_pressure = find_zero(
        mid_section_pressure -> series_valve_mass_flow_resid(
            mid_section_pressure,
            du, u, p, t, system, geo,
            tank_cell_id, pre_chamber_cell_id
        ),
        (u.pressure[pre_chamber_cell_id], u.pressure[tank_cell_id])
    )

    #We can now use the derived pressure in the cheapest valve mass flow equation
    true_mass_flow = adjustable_valve_flow(
        true_mid_section_pressure,
        du, u, p, t, system, geo,
        tank_cell_id, pre_chamber_cell_id
    )

    du.density[tank_cell_id] += -true_mass_flow / geo.cell_volumes[tank_cell_id]
    du.density[pre_chamber_cell_id] += true_mass_flow / geo.cell_volumes[pre_chamber_cell_id]

    tank_internal_energy_rate = -true_mass_flow * (u.specific_enthalpy[tank_cell_id] - u.specific_internal_energy[tank_cell_id])

    tank_mass = u.density[tank_cell_id] * geo.cell_volumes[tank_cell_id]
    pre_chamber_mass = u.density[pre_chamber_cell_id] * geo.cell_volumes[pre_chamber_cell_id]

    du.state_temperature[tank_cell_id] += tank_internal_energy_rate / (tank_mass * u.cv[tank_cell_id])
    du.state_temperature[pre_chamber_cell_id] -= tank_internal_energy_rate / (pre_chamber_mass * u.cv[pre_chamber_cell_id])

    #I really don't want to separate latent heat cooling to decrease the pre_chamber temperature properly, so I'm just going to neglect it :P
end







