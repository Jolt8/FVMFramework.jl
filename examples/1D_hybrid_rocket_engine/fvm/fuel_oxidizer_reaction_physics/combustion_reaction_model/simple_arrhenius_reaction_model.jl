

"""
C2H4 + 3O2 --> 2CO2 + 2H2O
simple_arrhenius_reaction_model_pre_exponential_factor = 6.3245e7 [m^3/(mol*s)]
simple_arrhenius_reaction_model_activation_energy = 125_520 [J/mol]
simple_arrhenius_reaction_model_heat_of_reaction = -1.323e6 [J/mol]
"""
function combust_ethylene_and_oxygen_simple!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.simple_arrhenius_reaction_model_pre_exponential_factor * exp(-p.simple_arrhenius_reaction_model_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    arrhenius_predicted_reaction_rate = rate_constant * u.molar_concentrations.ethylene[cell_id]^0.1 * u.molar_concentrations.oxygen[cell_id]^1.65

    oxygen_mass_required_per_ethylene_mass = 3 * u.molar_masses.oxygen[cell_id] / u.molar_masses.ethylene[cell_id]
    oxygen_as_ethylene_equivalent_mass_density = u.species_densities.oxygen[cell_id] / oxygen_mass_required_per_ethylene_mass
    limiting_ethylene_equivalent_mass_density = min(u.species_densities.ethylene[cell_id], oxygen_as_ethylene_equivalent_mass_density)

    turbulence_limited_ethylene_mass_rate_density = turbulence_limited_reactant_mass_rate_density(u, p, cell_id, limiting_ethylene_equivalent_mass_density)
    turbulence_limited_reaction_rate = turbulence_limited_ethylene_mass_rate_density / u.molar_masses.ethylene[cell_id]
    
    reaction_rate = min(arrhenius_predicted_reaction_rate, turbulence_limited_reaction_rate) #[mol/(m^3*s)]
    
    du.species_densities.ethylene[cell_id] += -1 * reaction_rate * u.molar_masses.ethylene[cell_id]
    du.species_densities.oxygen[cell_id] += -3 * reaction_rate * u.molar_masses.oxygen[cell_id]
    du.species_densities.carbon_dioxide[cell_id] += 2 * reaction_rate * u.molar_masses.carbon_dioxide[cell_id]
    du.species_densities.water[cell_id] += 2 * reaction_rate * u.molar_masses.water[cell_id]

    du.volumetric_energy[cell_id] += reaction_rate * -p.simple_arrhenius_reaction_model_heat_of_reaction
end
