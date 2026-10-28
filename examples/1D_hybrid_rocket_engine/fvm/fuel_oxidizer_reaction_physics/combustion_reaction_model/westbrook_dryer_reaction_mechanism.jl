"""
C2H4 + 2O2 --> 2CO + 2H2O
westbrook_dryer_step_1_pre_exponential_factor = 7.5895e7 [m^2.25 * mol^-0.75 * s^-1]
westbrook_dryer_step_1_activation_energy = 125520 [J/mol]
"""
function westbrook_dryer_step_1!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.westbrook_dryer_step_1_pre_exponential_factor * exp(-p.westbrook_dryer_step_1_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.ethylene[cell_id]^0.1 * u.molar_concentrations.oxygen[cell_id]^1.65

    oxygen_mass_required_per_ethylene_mass = 2 * u.molar_masses.oxygen[cell_id] / u.molar_masses.ethylene[cell_id]
    oxygen_as_ethylene_equivalent_mass_density = u.species_densities.oxygen[cell_id] / oxygen_mass_required_per_ethylene_mass
    limiting_ethylene_equivalent_mass_density = min(u.species_densities.ethylene[cell_id], oxygen_as_ethylene_equivalent_mass_density)

    turbulence_limited_ethylene_mass_rate_density = turbulence_limited_reactant_mass_rate_density(u, p, cell_id, limiting_ethylene_equivalent_mass_density)
    turbulence_limited_rate = turbulence_limited_ethylene_mass_rate_density / u.molar_masses.ethylene[cell_id]

    rate = min(chemical_rate, turbulence_limited_rate) #[mol/(m^3*s)]

    du.species_densities.ethylene[cell_id] -= rate * u.molar_masses.ethylene[cell_id]
    du.species_densities.oxygen[cell_id] -= 2 * rate * u.molar_masses.oxygen[cell_id]
    du.species_densities.carbon_monoxide[cell_id] += 2 * rate * u.molar_masses.carbon_monoxide[cell_id]
    du.species_densities.water[cell_id] += 2 * rate * u.molar_masses.water[cell_id]

    mol_per_s = rate * geo.cell_volumes[cell_id]
    
    du.volumetric_energy[cell_id] -= (mol_per_s * p.westbrook_dryer_step_1_heat_of_reaction) / geo.cell_volumes[cell_id]
end

"""
CO + 0.5O2 --> CO2
westbrook_dryer_step_2_forward_pre_exponential_factor = 1.2589e10 [m^2.25 * mol^-0.75 * s^-1]
westbrook_dryer_step_2_forward_activation_energy = 167360 [J/mol]
"""
function westbrook_dryer_step_2_forward!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.westbrook_dryer_step_2_forward_pre_exponential_factor * exp(-p.westbrook_dryer_step_2_forward_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.carbon_monoxide[cell_id] * u.molar_concentrations.water[cell_id]^0.5 * u.molar_concentrations.oxygen[cell_id]^0.25

    oxygen_mass_required_per_carbon_monoxide_mass = 0.5 * u.molar_masses.oxygen[cell_id] / u.molar_masses.carbon_monoxide[cell_id]
    oxygen_as_carbon_monoxide_equivalent_mass_density = u.species_densities.oxygen[cell_id] / oxygen_mass_required_per_carbon_monoxide_mass
    limiting_carbon_monoxide_equivalent_mass_density = min(u.species_densities.carbon_monoxide[cell_id], oxygen_as_carbon_monoxide_equivalent_mass_density)

    turbulence_limited_carbon_monoxide_mass_rate_density = turbulence_limited_reactant_mass_rate_density(u, p, cell_id, limiting_carbon_monoxide_equivalent_mass_density)
    turbulence_limited_rate = turbulence_limited_carbon_monoxide_mass_rate_density / u.molar_masses.carbon_monoxide[cell_id]

    rate = min(chemical_rate, turbulence_limited_rate) #[mol/(m^3*s)]

    du.species_densities.carbon_monoxide[cell_id] -= rate * u.molar_masses.carbon_monoxide[cell_id]
    du.species_densities.oxygen[cell_id] -= 0.5 * rate * u.molar_masses.oxygen[cell_id]
    du.species_densities.carbon_dioxide[cell_id] += rate * u.molar_masses.carbon_dioxide[cell_id]

    mol_per_s = rate * geo.cell_volumes[cell_id]
    
    du.volumetric_energy[cell_id] -= (mol_per_s * p.westbrook_dryer_step_2_forward_heat_of_reaction) / geo.cell_volumes[cell_id]
end

"""
CO2 --> CO + 0.5O2
westbrook_dryer_step_2_reverse_pre_exponential_factor = 5.0000e8 [s^-1]
westbrook_dryer_step_2_reverse_activation_energy = 167360 [J/mol]
"""
function westbrook_dryer_step_2_reverse!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.westbrook_dryer_step_2_reverse_pre_exponential_factor * exp(-p.westbrook_dryer_step_2_reverse_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.carbon_dioxide[cell_id]

    # CO2 dissociation is a homogeneous unimolecular step. It does not require
    # two reactants to mix, so the chemical rate is not turbulence-limited.
    rate = chemical_rate #[mol/(m^3*s)]

    du.species_densities.carbon_dioxide[cell_id] -= rate * u.molar_masses.carbon_dioxide[cell_id]
    du.species_densities.carbon_monoxide[cell_id] += rate * u.molar_masses.carbon_monoxide[cell_id]
    du.species_densities.oxygen[cell_id] += 0.5 * rate * u.molar_masses.oxygen[cell_id]

    mol_per_s = rate * geo.cell_volumes[cell_id]
    
    du.volumetric_energy[cell_id] -= (mol_per_s * p.westbrook_dryer_step_2_reverse_heat_of_reaction) / geo.cell_volumes[cell_id]
end

function combust_ethylene_and_oxygen_westbrook_dryer!(du, u, p, t, system, geo, cell_id)
    westbrook_dryer_step_1!(du, u, p, t, system, geo, cell_id)
    westbrook_dryer_step_2_forward!(du, u, p, t, system, geo, cell_id)
    westbrook_dryer_step_2_reverse!(du, u, p, t, system, geo, cell_id)
end

