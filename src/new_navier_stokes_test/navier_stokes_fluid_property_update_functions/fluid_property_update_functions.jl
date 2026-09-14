function update_velocities!(du, u, p, t, cell_id, vol)
    u.vel_u[cell_id] = u.momentum_density_u[cell_id] / u.density[cell_id]
    u.vel_v[cell_id] = u.momentum_density_v[cell_id] / u.density[cell_id]
    u.vel_w[cell_id] = u.momentum_density_w[cell_id] / u.density[cell_id]
end

function update_specific_energy!(du, u, p, t, cell_id, vol)
    velocity_magnitude_squared = u.vel_u[cell_id]^2 + u.vel_v[cell_id]^2 + u.vel_w[cell_id]^2

    u.volumetric_energy[cell_id] = u.specific_internal_energy[cell_id] + 0.5 * velocity_magnitude_squared
end

function update_pressure_ideal!(du, u, p, t, cell_id, vol)
    specific_heat_ratio = u.cp[cell_id] / u.cv[cell_id]

    momentum_magnitude_squared = u.momentum_density_u[cell_id]^2 + u.momentum_density_v[cell_id]^2 + u.momentum_density_w[cell_id]^2

    u.pressure[cell_id] = (specific_heat_ratio - 1) * (
        u.volumetric_energy[cell_id] - momentum_magnitude_squared / (2 * u.density[cell_id])
    )
end

function update_speed_of_sound_ideal(du, u, p, t, cell_id, vol)
    specific_heat_ratio = u.cp[cell_id] / u.cv[cell_id]

    u.speed_of_sound[cell_id] = sqrt(specific_heat_ratio * u.pressure[cell_id] / u.density[cell_id])
end

function update_temperature_ideal(du, u, p, t, cell_id, vol)
    u.temperature[cell_id] = u.volumetric_energy[cell_id] / u.cv[cell_id]
end