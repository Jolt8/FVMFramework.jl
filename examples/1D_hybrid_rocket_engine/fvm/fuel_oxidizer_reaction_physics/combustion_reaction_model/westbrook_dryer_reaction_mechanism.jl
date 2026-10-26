"""
C2H4 + 2O2 --> 2CO + 2H2O
westbrook_dryer_step_1_pre_exponential_factor = 7.5895e7 [m^2.25 * mol^-0.75 * s^-1]
westbrook_dryer_step_1_activation_energy = 125520 [J/mol]
"""
function westbrook_dryer_step_1!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.westbrook_dryer_step_1_pre_exponential_factor * exp(-p.westbrook_dryer_step_1_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.ethylene[cell_id]^0.1 * u.molar_concentrations.oxygen[cell_id]^1.65

    turbulence_limited_rate = turbulence_limited_reaction_rate(du, u, p, t, system, geo, cell_id)

    rate = min(chemical_rate, turbulence_limited_rate) #[mol/(m^3*s)]

    du.species_densities.ethylene[cell_id] -= rate * u.species_molar_masses.ethylene[cell_id]
    du.species_densities.oxygen[cell_id] -= 2 * rate * u.species_molar_masses.oxygen[cell_id]
    du.species_densities.carbon_monoxide[cell_id] += 2 * rate * u.species_molar_masses.carbon_monoxide[cell_id]
    du.species_densities.water[cell_id] += 2 * rate * u.species_molar_masses.water[cell_id]
end

"""
CO + 0.5O2 --> CO2
westbrook_dryer_step_2_forward_pre_exponential_factor = 1.2589e10 [m^2.25 * mol^-0.75 * s^-1]
westbrook_dryer_step_2_forward_activation_energy = 167360 [J/mol]
"""
function westbrook_dryer_step_2_forward!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.westbrook_dryer_step_2_forward_pre_exponential_factor * exp(-p.westbrook_dryer_step_2_forward_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.carbon_monoxide[cell_id] * u.molar_concentrations.water[cell_id]^0.5 * u.molar_concentrations.oxygen[cell_id]^0.25

    turbulence_limited_rate = turbulence_limited_reaction_rate(du, u, p, t, system, geo, cell_id)

    rate = min(chemical_rate, turbulence_limited_rate)

    du.species_densities.carbon_monoxide[cell_id] -= rate * u.species_molar_masses.carbon_monoxide[cell_id]
    du.species_densities.oxygen[cell_id] -= 0.5 * rate * u.species_molar_masses.oxygen[cell_id]
    du.species_densities.carbon_dioxide[cell_id] += rate * u.species_molar_masses.carbon_dioxide[cell_id]
end

"""
CO2 --> CO + 0.5O2
westbrook_dryer_step_2_reverse_pre_exponential_factor = 5.0000e8 [s^-1]
westbrook_dryer_step_2_reverse_activation_energy = 167360 [J/mol]
"""
function westbrook_dryer_step_2_reverse!(du, u, p, t, system, geo, cell_id)
    rate_constant = p.westbrook_dryer_step_2_reverse_pre_exponential_factor * exp(-p.westbrook_dryer_step_2_reverse_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.carbon_dioxide[cell_id]

    turbulence_limited_rate = turbulence_limited_reaction_rate(du, u, p, t, system, geo, cell_id)

    rate = min(chemical_rate, turbulence_limited_rate)

    du.species_densities.carbon_dioxide[cell_id] -= rate * u.species_molar_masses.carbon_dioxide[cell_id]
    du.species_densities.carbon_monoxide[cell_id] += rate * u.species_molar_masses.carbon_monoxide[cell_id]
    du.species_densities.oxygen[cell_id] += 0.5 * rate * u.species_molar_masses.oxygen[cell_id]
end

