function populate_system_mole_fractions_vec!(du, u, p, t, system, geo)
    for cell_id in 1:n_cells
        species_idx = 1
        
        for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
            system.mole_fractions_vec[cell_id][species_idx] = u_mass_fractions[species[cell_id]] * (u.mw_avg[cell_id] / u_molecular_weights[species[cell_id]])

            species_idx += 1
        end
    end
end