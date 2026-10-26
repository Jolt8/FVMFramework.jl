#TODO: In the future we MUST model incomplete combustion 

function update_species_molar_concentrations!(du, u, p, t, system, geo, cell_id)
    for_fields!(u.molar_concentrations, u.species_densities, u.molar_masses) do species, u_molar_concentrations, u_species_densities, update_species_mole_fractions
        u_molar_concentrations[species[cell_id]] = u_species_densities[species[cell_id]] / u.molar_masses[species[cell_id]]
    end
end

function decompose_nitrous_oxide!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.nitrous_oxide_decomposition_pre_exponential_factor * exp(-p.nitrous_oxide_decomposition_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    nitrous_oxide_mole_fraction = u.species_densities.nitrous_oxide[cell_id] / u.molar_masses.nitrous_oxide[cell_id]

    reaction_rate_mol_per_s = rate_constant * (nitrous_oxide_mole_fraction)^2

    change_in_nitrous_oxide_mass = -2 * reaction_rate_mol_per_s * u.molar_masses.nitrous_oxide[cell_id]
    change_in_oxygen_mass = reaction_rate_mol_per_s * u.molar_masses.oxygen[cell_id]
    change_in_nitrogen_mass = 2 * reaction_rate_mol_per_s * u.molar_masses.nitrogen[cell_id]

    du.species_densities.nitrous_oxide[cell_id] += change_in_nitrous_oxide_mass / geo.cell_volumes[cell_id]
    du.species_densities.oxygen[cell_id] += change_in_oxygen_mass / geo.cell_volumes[cell_id]
    du.species_densities.nitrogen[cell_id] += change_in_nitrogen_mass / geo.cell_volumes[cell_id]
end

function turbulence_limited_reaction_rate(du, u, p, t, system, geo, cell_id)
    turbulence_mixing_timescale = u.density[cell_id] / u.specific_dissipation_rate_density[cell_id]

    limiting_species_term = min(u.species_densities.ethylene[cell_id], u.species_densities.oxygen[cell_id] / 3.42) #3.42 is the stoichiometric ratio of oxygen to ethylene required by mass for complete combustion

    turbulence_preducted_reaction_rate = p.turbulence_limited_reaction_rate_mix_parameter * (u.density[cell_id] / turbulence_mixing_timescale) * limiting_species_term

    return turbulence_preducted_reaction_rate
end

include("combustion_reaction_model/simple_arrhenius_reaction_model.jl")
export combust_ethylene_and_oxygen_simple!
include("combustion_reaction_model/westbrook_dryer_reaction_mechanism.jl")
export combust_ethylene_and_oxygen_westbrook_dryer!