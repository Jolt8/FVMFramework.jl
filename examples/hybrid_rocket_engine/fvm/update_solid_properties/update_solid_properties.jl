#temperature
function get_solid_temperature_ideal(u, cell_id)
    density = max(u.density[cell_id], 1e-10)

    internal_energy_density = max(u.volumetric_energy[cell_id], 1e-10)

    return internal_energy_density / (density * u.cp[cell_id])
end

function update_solid_temperature_ideal!(du, u, p, t, system, geo, cell_id)
    u.temperature[cell_id] = get_solid_temperature_ideal(u, cell_id)
end

function update_solid_properties!(du, u, p, t, system, geo, cell_id)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    update_solid_temperature_ideal!(du, u, p, t, system, geo, cell_id)
end
