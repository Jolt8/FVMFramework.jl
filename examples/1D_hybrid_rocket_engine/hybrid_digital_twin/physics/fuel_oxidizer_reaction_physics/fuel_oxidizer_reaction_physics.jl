function decompose_nitrous_oxide!(du, u, p, t, system, geo, cell_id)
    if u.mass_fractions.nitrous_oxide[cell_id] <= 1e-10
        return
    end

    rate_constant = p.nitrous_oxide_decomposition_pre_exponential_factor * exp(-p.nitrous_oxide_decomposition_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))
    
    rate = rate_constant * u.molar_concentrations.nitrous_oxide[cell_id] #[mol/(m^3*s)]

    du.mass_fractions.nitrous_oxide[cell_id] -= 2 * rate * u.molecular_weights.nitrous_oxide[cell_id] / u.density[cell_id]
    du.mass_fractions.oxygen[cell_id] += rate * u.molecular_weights.oxygen[cell_id] / u.density[cell_id]
    du.mass_fractions.nitrogen[cell_id] += 2 * rate * u.molecular_weights.nitrogen[cell_id] / u.density[cell_id]

    mol_per_s = rate * geo.cell_volumes[cell_id]

    du.heat[cell_id] -= mol_per_s * p.nitrous_oxide_decomposition_heat_of_reaction
end

"""
Calculate the turbulence-limited consumption rate of a reactant-equivalent mass.

`limiting_reactant_mass_density` is expressed as the mass density of the reference
reactant for the reaction. The returned rate therefore has units of kg/(m^3*s);
each reaction converts it to mol/(m^3*s) with that reference reactant's molar mass.
"""
function turbulence_limited_reactant_mass_rate_density(u, p, cell_id, limiting_reactant_mass_density)
    turbulence_mixing_timescale = u.density[cell_id] / u.specific_dissipation_rate_density[cell_id]

    turbulence_limited_mass_rate_density = p.turbulence_limited_reaction_rate_mix_parameter * limiting_reactant_mass_density / turbulence_mixing_timescale

    return turbulence_limited_mass_rate_density
end
