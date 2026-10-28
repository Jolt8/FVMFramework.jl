#TODO: In the future we MUST model incomplete combustion 

function update_species_molar_concentrations!(du, u, p, t, system, geo, cell_id)
    for_fields!(u.molar_concentrations, u.species_densities, u.molar_masses) do species, u_molar_concentrations, u_species_densities, update_species_mole_fractions
        u_molar_concentrations[species[cell_id]] = u_species_densities[species[cell_id]] / u.molar_masses[species[cell_id]]
    end
end

function decompose_nitrous_oxide!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.nitrous_oxide_decomposition_pre_exponential_factor * exp(-p.nitrous_oxide_decomposition_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    nitrous_oxide_molar_concentration = u.species_densities.nitrous_oxide[cell_id] / u.molar_masses.nitrous_oxide[cell_id]

    rate = rate_constant * nitrous_oxide_molar_concentration #[mol/(m^3*s)]

    du.species_densities.nitrous_oxide[cell_id] += -2 * rate * u.molar_masses.nitrous_oxide[cell_id]
    du.species_densities.oxygen[cell_id] += rate * u.molar_masses.oxygen[cell_id]
    du.species_densities.nitrogen[cell_id] += 2 * rate * u.molar_masses.nitrogen[cell_id]
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

Revise.includet("combustion_reaction_model/simple_arrhenius_reaction_model.jl")
export combust_ethylene_and_oxygen_simple!
Revise.includet("combustion_reaction_model/westbrook_dryer_reaction_mechanism.jl")
export combust_ethylene_and_oxygen_westbrook_dryer!
