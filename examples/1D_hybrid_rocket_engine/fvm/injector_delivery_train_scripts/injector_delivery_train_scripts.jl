function get_injector_train_u_proto()
    injector_u_proto = ComponentVector(
        tank = (
            mass = 0.0u"kg",
            temperature = 0.0u"K",
        ),
        mid_section = (
            mass = 0.0u"kg",
            temperature = 0.0u"K"
        )
    )

    return injector_u_proto
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

function update_injector_train_properties!(du, u, p, t, system, geo)
    props = system.injector_train_properties
    model = system.nitrous_oxide_model

    tank_n_moles = u.tank.mass / props.nitrous_oxide_molecular_weight

    props.tank.pressure = pressure(model, props.tank.volume, u.tank.temperature, [tank_n_moles])

    tank_flash = tp_flash(model, props.tank.pressure, u.tank.temperature, [tank_n_moles])

    props.tank.specific_internal_energy = mass_internal_energy(model, tank_flash)
    props.tank.enthalpy = mass_enthalpy(model, tank_flash)

    p.tank_vapor_fraction = tank_flash.fractions[vapor_phase] / sum(tank_flash.fractions)

    if p.tank_vapor_fraction <= 0.999
        props.tank.density = mass_density(oxidizer_model, tank_flash, 1) #get liquid density because we drawing from the bottom of the tank
        props.tank.enthalpy = mass_enthalpy(oxidizer_model, tank_flash, 1)
        props.tank.cp = mass_isobaric_heat_capacity(model, tank_flash, 1)
        props.tank.cv = mass_isochoric_heat_capacity(model, tank_flash, 1)
    else
        props.tank.density = mass_density(oxidizer_model, tank_flash) #otherwise, we will be drawing from the remaining vapor in the tank
        props.tank.enthalpy = mass_enthalpy(oxidizer_model, tank_flash)
        props.tank.cp = mass_isobaric_heat_capacity(model, tank_flash, 2)
        props.tank.cv = mass_isochoric_heat_capacity(model, tank_flash, 2)
    end

    mid_section_n_moles = u.mid_section.mass / props.nitrous_oxide_molecular_weight

    props.mid_section.pressure = pressure(model, props.mid_section.volume, u.mid_section.temperature, [mid_section_n_moles])

    mid_section_flash = tp_flash(model, props.mid_section.pressure, u.mid_section.temperature, [mid_section_n_moles])

    props.mid_section.vapor_fraction = mid_section_flash.fractions[vapor_phase] / sum(mid_section_flash.fractions)

    props.mid_section.density = mass_density(oxidizer_model, mid_section_flash) #otherwise, we will be drawing from the remaining vapor in the tank
    props.mid_section.enthalpy = mass_enthalpy(oxidizer_model, mid_section_flash)
    props.mid_section.cp = mass_isobaric_heat_capacity(model, mid_section_flash, 2)
    props.mid_section.cv = mass_isochoric_heat_capacity(model, mid_section_flash, 2)
end


function adjustable_valve_flow(du, u, p, t, system, geo)
    props = system.injector_train_properties

    if props.tank.pressure - props.mid_section.pressure <= 0
        return 0.0
    end
    
    #Notice that valve_opening uses system.valve_opening, this is because we're going to be adjusting this over time
    volumetric_flow = system.valve_opening * props.valve_flow_capacity_factor * sqrt((props.tank.pressure - props.mid_section.pressure) / props.tank.density)

    mass_flow = volumetric_flow * props.tank.density

    return mass_flow
end

function injector_valve_flow(du, u, p, t, system, geo)
    props = system.injector_train_properties

    post_injector_cell_ids = system.post_injector_cell_ids
    chamber_pressure = sum([u.pressure[cell_id] for cell_id in post_injector_cell_ids]) / length(post_injector_cell_ids)

    if props.mid_section.pressure - chamber_pressure <= 0
        return 0.0
    end

    volumetric_flow = p.injector_discharge_coefficient * p.injector_orifice_area * sqrt((2 * (props.mid_section.pressure - chamber_pressure)) / props.mid_section.density)

    oxidizer_mass_flow = volumetric_flow * props.mid_section.density

    return oxidizer_mass_flow
end

function update_injector_train!(du, u, p, t, system, geo)
    props = system.injector_train_properties
    model = system.nitrous_oxide_model

    tank_to_mid_section_mass_flow = adjustable_valve_flow(du, u, p, t, system, geo)
    injector_mass_flow = injector_valve_flow(du, u, p, t, system, geo)

    #Tank
    du.tank.mass = -tank_to_mid_section_mass_flow

    #du_tank_internal_energy = tank_to_mid_section_mass_flow * props.tank.specific_enthalpy
    #du.tank.temperature = (du_tank_internal_energy - props.tank.specific_internal_energy * tank_to_mid_section_mass_flow) / (u.tank.mass * props.tank.cv)

    #Mid Section
    du.mid_section.mass = tank_to_mid_section_mass_flow - injector_mass_flow

    #du_mid_section_internal_energy = tank_to_mid_section_mass_flow * props.tank.specific_enthalpy - injector_mass_flow * props.mid_section.specific_enthalpy
    #du.mid_section.temperature = (du_mid_section_internal_energy - props.mid_section.specific_internal_energy * (tank_to_mid_section_mass_flow - injector_mass_flow)) / (u.mid_section.mass * props.mid_section.cv)

    injector_velocity = injector_mass_flow / (props.mid_section.density * props.injector_orifice_area)

    return ComponentVector(
        mass_flow = tank_to_mid_section_mass_flow,
        mid_section_temperature = u.mid_section.temperature,
        mid_section_pressure = props.mid_section.pressure,
        mid_section_density = props.mid_section.density,
        mid_section_specific_internal_energy = props.mid_section.specific_internal_energy,
        injector_velocity = injector_velocity,
    )
end

    

    



    
