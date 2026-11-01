function generate_system_moles_vec(config)
    return [zeros(length(config.u_proto.mass_fractions)) for _ in 1:config.geo.n_cells]
end

function update_cell_moles_vec!(du, u, p, t, system, geo, cell_id)
    species_idx = 1

    cell_mass = u.density[cell_id] * geo.cell_volumes[cell_id]
    
    for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
        system.moles_vec[cell_id][species_idx] = (u_mass_fractions[species[cell_id]] * cell_mass) / u_molecular_weights[species[cell_id]]

        species_idx += 1
    end
end

function generate_system_mole_fractions_vec(config)
    return [zeros(length(config.u_proto.mass_fractions)) for _ in 1:config.geo.n_cells]
end

function update_cell_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    species_idx = 1
    
    for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
        system.mole_fractions_vec[cell_id][species_idx] = u_mass_fractions[species[cell_id]] * (u.mw_avg[cell_id] / u_molecular_weights[species[cell_id]])

        species_idx += 1
    end
end

function update_cell_moles_and_mole_fractions_vec(du, u, p, t, system, geo, cell_id)
    update_cell_moles_vec!(du, u, p, t, system, geo, cell_id)
    update_cell_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
end

