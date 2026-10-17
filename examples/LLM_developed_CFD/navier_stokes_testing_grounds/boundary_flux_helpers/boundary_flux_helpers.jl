function physical_flux(
    density,
    momentum_density_u,
    momentum_density_v,
    momentum_density_w,
    volumetric_energy,
    pressure,
    normal
)
    vel_u = momentum_density_u / density
    vel_v = momentum_density_v / density
    vel_w = momentum_density_w / density

    normal_velocity =
        vel_u * normal[1] +
        vel_v * normal[2] +
        vel_w * normal[3]

    return (
        density * normal_velocity,

        momentum_density_u * normal_velocity +
            pressure * normal[1],

        momentum_density_v * normal_velocity +
            pressure * normal[2],

        momentum_density_w * normal_velocity +
            pressure * normal[3],

        (volumetric_energy + pressure) * normal_velocity
    )
end

function primitive_from_conservative(
    density,
    momentum_density_u,
    momentum_density_v,
    momentum_density_w,
    volumetric_energy,
    specific_heat_ratio
)
    density = max(density, 1e-10)

    vel_u = momentum_density_u / density
    vel_v = momentum_density_v / density
    vel_w = momentum_density_w / density

    kinetic_energy_density = 0.5 * (
        momentum_density_u^2 +
        momentum_density_v^2 +
        momentum_density_w^2
    ) / density

    internal_energy_density = max((volumetric_energy - kinetic_energy_density), 1e-10) 

    pressure = (specific_heat_ratio - 1) * internal_energy_density

    speed_of_sound = sqrt(specific_heat_ratio * pressure / density)

    return vel_u, vel_v, vel_w, pressure, speed_of_sound
end