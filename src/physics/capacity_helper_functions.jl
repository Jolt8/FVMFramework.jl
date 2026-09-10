
"""
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)

Applies a change in a cells heat to a change in the cell's temperature.

Parameters:
    - du: State derivatives vector
    - u: State vector
    - cell_id: ID of the cell we are calculating the capacity for
    - vol: Volume of the cell
"""
function cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    # J/s /= m^3 * kg*m^3 * J/(kg*K)
    # = K/s
    du.temp[cell_id] += du.heat[cell_id] / (vol * u.rho[cell_id] * u.cp[cell_id])
end

"""
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)

Applies a change in a cells mass to a change in the cell's pressure.

Parameters:
    - du: State derivatives vector
    - u: State vector
    - cell_id: ID of the cell we are calculating the capacity for
    - vol: Volume of the cell
"""
function cap_mass_flux_to_pressure_change_ideal!(du, u, cell_id, vol)
    # kg/s /= (m^3 / (J/(mol*K) * K))
    #remember: J = Pa*m^3
    # = Pa/s
    du_moles = du.mass[cell_id] / u.mw_avg[cell_id]
    du.pressure[cell_id] += (du_moles * u.R_gas[cell_id] * u.temp[cell_id]) / vol
end

"""
    cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)

Applies a change in the mass of each species entering a

Parameters:
    - du: State derivatives vector
    - u: State vector
    - cell_id: ID of the cell we are calculating the capacity for
    - vol: Volume of the cell
"""
function cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    total_mass = vol * u.rho[cell_id]

    for_fields!(du.mass_fractions, u.mass_fractions, du.species_masses) do species, du_mass_fractions, u_mass_fractions, species_masses
        du_mass_fractions[species[cell_id]] += (species_masses[species[cell_id]] - u_mass_fractions[species[cell_id]] * du.mass[cell_id]) / total_mass
    end
end

