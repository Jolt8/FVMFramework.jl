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
        vapor_fraction =
            (bulk_molar_volume - liquid_molar_volume) /
            (vapor_molar_volume - liquid_molar_volume)
        return Clapeyron.FlashResult(
            saturation_pressure_value,
            temperature,
            [pure_composition, pure_composition],
            [amount * (1 - vapor_fraction), amount * vapor_fraction],
            [liquid_molar_volume, vapor_molar_volume],
            zero(saturation_pressure_value);
            sort = false,
            vapour_phase_index = 2,
        )
    end
end

function flash_vapor_fraction(flash)
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
        weighted_property +=
            phase_fraction * property_function(model, flash, phase)
    end

    return weighted_property
end

function get_injector_train_u0_and_properties(nitrous_oxide_model)
    properties = ComponentVector(
        tank = (
            volume = 1.6u"L",
            density = 0.0u"kg/m^3",
            pressure = 0.0u"bar",
            cp = 0.0u"J/(kg*K)",
            cv = 0.0u"J/(kg*K)",
            specific_internal_energy = 0.0u"J/kg",
            specific_enthalpy = 0.0u"J/kg",
            vapor_fraction = 0.0u"kg/kg",
        ),
        mid_section = (
            volume = 474.05u"mm^3",
            density = 0.0u"kg/m^3",
            pressure = 0.0u"bar",
            cp = 0.0u"J/(kg*K)",
            cv = 0.0u"J/(kg*K)",
            specific_internal_energy = 0.0u"J/kg",
            specific_enthalpy = 0.0u"J/kg",
            vapor_fraction = 0.0u"kg/kg",
        ),

        nitrous_oxide_molecular_weight = 44.013u"g/mol",
        injector_orifice_area = 4.00u"mm^2",
        injector_discharge_coefficient = 0.7,
        valve_flow_capacity_factor = 10.0u"mm^2",
    )

    u0_mid_section_mass = mass_density(nitrous_oxide_model, 1.0u"atm", 298.0u"K", [1.0]) * properties.mid_section.volume |> u"kg"

    injector_u0 = ComponentVector(
        tank = (
            mass = 1.134u"kg",
            temperature = 298.0u"K",
        ),
        mid_section = (
            mass = u0_mid_section_mass,
            temperature = 298.0u"K"
        )
    )

    injector_u0 = ustrip.(upreferred.(injector_u0))
    properties = ustrip.(upreferred.(properties))

    return injector_u0, properties
end

function update_injector_train_properties!(du, u, p, t, system, geo, props)
    # props = system.injector_train_properties
    model = system.nitrous_oxide_model

    tank_n_moles = u.tank.mass / props.nitrous_oxide_molecular_weight

    tank_flash = pure_component_vt_flash(
        model,
        props.tank.volume,
        u.tank.temperature,
        tank_n_moles,
    )
    props.tank.pressure = pressure(tank_flash)
    props.tank.vapor_fraction = flash_vapor_fraction(tank_flash)

    # These describe the bulk material stored in the rigid tank. They are
    # intentionally independent of which phase is being withdrawn.
    props.tank.specific_internal_energy = mass_internal_energy(model, tank_flash)
    props.tank.cp = phase_weighted_mass_property(mass_isobaric_heat_capacity, model, tank_flash)
    props.tank.cv = phase_weighted_mass_property(mass_isochoric_heat_capacity, model, tank_flash)

    # These describe the material leaving the bottom of the tank. Withdraw
    # liquid while a distinct liquid inventory remains, then withdraw vapor.
    if props.tank.vapor_fraction <= 0.999
        outlet_phase = argmin(tank_flash.volumes)
    else
        vapor_phase = tank_flash.data.vapour_idx
        outlet_phase = vapor_phase > 0 ? vapor_phase : argmax(tank_flash.volumes)
    end
    props.tank.density = mass_density(model, tank_flash, outlet_phase)
    props.tank.specific_enthalpy = mass_enthalpy(model, tank_flash, outlet_phase)

    mid_section_n_moles = u.mid_section.mass / props.nitrous_oxide_molecular_weight

    mid_section_flash = pure_component_vt_flash(
        model,
        props.mid_section.volume,
        u.mid_section.temperature,
        mid_section_n_moles,
    )
    props.mid_section.pressure = pressure(mid_section_flash)

    props.mid_section.vapor_fraction = flash_vapor_fraction(mid_section_flash)
    props.mid_section.density = mass_density(model, mid_section_flash)
    props.mid_section.specific_internal_energy = mass_internal_energy(model, mid_section_flash)
    props.mid_section.specific_enthalpy = mass_enthalpy(model, mid_section_flash)
    props.mid_section.cp = phase_weighted_mass_property(mass_isobaric_heat_capacity, model, mid_section_flash)
    props.mid_section.cv = phase_weighted_mass_property(mass_isochoric_heat_capacity, model, mid_section_flash)

    return nothing
end


function adjustable_valve_flow(du, u, p, t, system, geo, props)

    if props.tank.pressure - props.mid_section.pressure <= 0
        return 0.0
    end
    
    #Notice that valve_opening uses system.valve_opening, this is because we're going to be adjusting this over time
    volumetric_flow = system.valve_opening * props.valve_flow_capacity_factor * sqrt((props.tank.pressure - props.mid_section.pressure) / props.tank.density)

    mass_flow = volumetric_flow * props.tank.density

    return mass_flow
end

function injector_valve_flow(du, u, p, t, system, geo, props)

    post_injector_cell_ids = system.post_injector_cell_ids
    chamber_pressure = sum([u.pressure[cell_id] for cell_id in post_injector_cell_ids]) / length(post_injector_cell_ids)

    if props.mid_section.pressure - chamber_pressure <= 0
        return 0.0
    end

    volumetric_flow = props.injector_discharge_coefficient * props.injector_orifice_area * sqrt((2 * (props.mid_section.pressure - chamber_pressure)) / props.mid_section.density)

    oxidizer_mass_flow = volumetric_flow * props.mid_section.density

    return oxidizer_mass_flow
end

function unpack_props_cache(u, system)
    props_vec = get_tmp(system.injector_train_properties_cache, u.tank.mass)
    props_vec .= 0.0

    props_values = system.injector_train_properties_vec

    props_vec .+= props_values

    props = ComponentVector(props_vec, system.injector_train_properties_axes)

    return props
end

function update_injector_train!(du, u, p, t, system, geo)
    props = unpack_props_cache(u, system)

    u.mid_section.temperature = 300.0

    #pressure(oxidizer_model, 100000, 300, [1.0])
    #vt_flash(oxidizer_model, 1.0, 300, 1.0)
    #vt_flash(oxidizer_model, props.tank.volume, u.tank.temperature, u.tank.mass / props.nitrous_oxide_molecular_weight)

    update_injector_train_properties!(du, u, p, t, system, geo, props)

    tank_to_mid_section_mass_flow = adjustable_valve_flow(du, u, p, t, system, geo, props)
    injector_mass_flow = injector_valve_flow(du, u, p, t, system, geo, props)

    #Tank
    tank_mass_rate = -tank_to_mid_section_mass_flow
    tank_internal_energy_rate =
        -tank_to_mid_section_mass_flow * props.tank.specific_enthalpy
    du.tank.mass = tank_mass_rate
    du.tank.temperature =
        (tank_internal_energy_rate -
         props.tank.specific_internal_energy * tank_mass_rate) /
        (u.tank.mass * props.tank.cv)

    #Mid Section
    mid_section_mass_rate =
        tank_to_mid_section_mass_flow - injector_mass_flow
    mid_section_internal_energy_rate =
        tank_to_mid_section_mass_flow * props.tank.specific_enthalpy -
        injector_mass_flow * props.mid_section.specific_enthalpy
    du.mid_section.mass = mid_section_mass_rate
    du.mid_section.temperature =
        (mid_section_internal_energy_rate -
         props.mid_section.specific_internal_energy * mid_section_mass_rate) /
        (u.mid_section.mass * props.mid_section.cv)

    injector_velocity = injector_mass_flow / (props.mid_section.density * props.injector_orifice_area)

    return ComponentVector(
        mass_flow = injector_mass_flow,
        mid_section_temperature = u.mid_section.temperature,
        mid_section_pressure = props.mid_section.pressure,
        mid_section_density = props.mid_section.density,
        mid_section_specific_internal_energy = props.mid_section.specific_internal_energy,
        injector_velocity = injector_velocity,
        injector_area = props.injector_orifice_area,
    )
end

    

    



    
