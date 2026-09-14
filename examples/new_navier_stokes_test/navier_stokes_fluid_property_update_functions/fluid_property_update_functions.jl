#velocity
function get_velocity(u, cell_id)
    return (
        u.momentum_density_u[cell_id] / u.density[cell_id],
        u.momentum_density_v[cell_id] / u.density[cell_id],
        u.momentum_density_w[cell_id] / u.density[cell_id]
    )
end

function update_velocities!(du, u, p, t, cell_id, vol)
    u.vel_u[cell_id], u.vel_v[cell_id], u.vel_w[cell_id] = get_velocity(u, cell_id)
end

#specific energy
#=
function get_specific_energy(u, cell_id)
    return 0.5 * (u.vel_u[cell_id]^2 + u.vel_v[cell_id]^2 + u.vel_w[cell_id]^2)
end

function update_specific_energy!(du, u, p, t, cell_id, vol)
    u.specific_energy[cell_id] = get_specific_energy(u, cell_id)
end
=#

#temperature
function get_temperature_ideal(u, cell_id)
    KE_per_vol = 0.5 * (
        u.momentum_density_u[cell_id]^2 +
        u.momentum_density_v[cell_id]^2 +
        u.momentum_density_w[cell_id]^2
    ) / u.density[cell_id]

    internal_energy_density = u.volumetric_energy[cell_id] - KE_per_vol

    return internal_energy_density / (u.density[cell_id] * u.cv[cell_id])
end

function update_temperature_ideal!(du, u, p, t, cell_id, vol)
    u.temperature[cell_id] = get_temperature_ideal(u, cell_id)
end

#pressure
function get_pressure_ideal(u, cell_id)
    return u.density[cell_id] * (u.R_gas[cell_id] / u.mw[cell_id]) * u.temperature[cell_id]
end

function update_pressure_ideal!(du, u, p, t, cell_id, vol)
    u.pressure[cell_id] = get_pressure_ideal(u, cell_id)
end


#speed of sound
function get_speed_of_sound_ideal(u, cell_id)
    specific_heat_ratio = u.cp[cell_id] / u.cv[cell_id]

    return sqrt(specific_heat_ratio * u.pressure[cell_id] / u.density[cell_id])
end

function update_speed_of_sound_ideal!(du, u, p, t, cell_id, vol)
    u.speed_of_sound[cell_id] = get_speed_of_sound_ideal(u, cell_id)
end


function overall_navier_stokes_property_update!(du, u, p, t, cell_id, vol)
    update_velocities!(du, u, p, t, cell_id, vol)
    #update_specific_energy!(du, u, p, t, cell_id, vol)
    update_temperature_ideal!(du, u, p, t, cell_id, vol)
    update_pressure_ideal!(du, u, p, t, cell_id, vol)
    update_speed_of_sound_ideal!(du, u, p, t, cell_id, vol)
end