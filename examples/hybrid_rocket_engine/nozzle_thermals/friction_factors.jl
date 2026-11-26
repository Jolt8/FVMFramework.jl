#These were starting to take up too much space, so I decided to move them

function _colebrook_friction_factor(
    friction_factor,
    p, system, nozzle_position, helical_station_properties
)
    rhs_friction_factor = (
        -2.0 * log10(
            ((p.wall_equivalent_sand_grain_roughness / helical_station_properties.hydraulic_diameter) / 3.7) + 
            (2.51 / (helical_station_properties.coolant_reynolds_number * sqrt(friction_factor)))
        )
    )

    lhs_friction_factor = (1 / sqrt(friction_factor))

    return lhs_friction_factor - rhs_friction_factor
end

function get_colebrook_friction_factor(
    p, system, nozzle_position, helical_station_properties
)
    return find_zero(friction_factor -> _colebrook_friction_factor(
        friction_factor, p, system, nozzle_position, helical_station_properties
    ),
        (0.001, 100_000)
    )
end

function get_smooth_friction_factor(
    p, system, nozzle_position, helical_station_properties
)
    smooth_friction_factor = (0.790 * log(helical_station_properties.coolant_reynolds_number) - 1.64)^(-2)

    return smooth_friction_factor
end

function get_smooth_curved_friction_factor(
    p, system, nozzle_position, helical_station_properties
)
    curvature_ratio = (helical_station_properties.hydraulic_diameter / 2) / helical_station_properties.helix_curvature_radius

    smooth_curved_friction_factor = 0.304 * helical_station_properties.coolant_reynolds_number^(-0.25) * 0.029 * curvature_ratio^(0.5)

    return smooth_curved_friction_factor
end

function curved_rough_friction_factor(
    p, system, nozzle_position, helical_station_properties
)
    rough_friction_factor = get_colebrook_friction_factor(
        p, system, nozzle_position, helical_station_properties
    )

    smooth_friction_factor = get_smooth_friction_factor(
        p, system, nozzle_position, helical_station_properties
    )

    smooth_curved_friction_factor = get_smooth_curved_friction_factor(
        p, system, nozzle_position, helical_station_properties
    )

    return rough_friction_factor * (smooth_friction_factor / smooth_curved_friction_factor)
end

# Pressure drop equation:
function get_helix_frictional_pressure_drop(
    p, system, nozzle_station, helical_station,
    gas_side_wall_surface_temperature, coolant_side_wall_surface_temperature,
    coolant_temperature, coolant_pressure
)
    mass_flux = p.coolant_mass_flow / (p.number_of_channels * single_channel_flow_area)

    coolant_velocity = mass_flux / p.coolant_density

    coolant_pressure_drop = p.rough_friction_factor * (helical_station.channel_length_along_nozzle_station / helical_station.hydraulic_diameter) * ((p.coolant_density * coolant_velocity^2) / 2)

    return coolant_pressure_drop
end