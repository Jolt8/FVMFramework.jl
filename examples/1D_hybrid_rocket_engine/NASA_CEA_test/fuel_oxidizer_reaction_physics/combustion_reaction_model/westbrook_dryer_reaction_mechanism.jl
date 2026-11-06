"""
C2H4 + 2O2 --> 2CO + 2H2O
westbrook_dryer_step_1_pre_exponential_factor = 7.5895e7 [m^2.25 * mol^-0.75 * s^-1]
westbrook_dryer_step_1_activation_energy = 125520 [J/mol]
"""
function westbrook_dryer_step_1!(du, u, p, t, system, geo, cell_id)
    if system.use_smooth_reaction_ramp
        if u.mass_fractions.ethylene[cell_id] <= 0.0 || u.mass_fractions.oxygen[cell_id] <= 0.0
            return
        end

        ethylene_ramp = smooth_ramp(u.mass_fractions.ethylene[cell_id], p.species_min_mass_fraction_ramp_activation, p.species_max_mass_fraction_ramp_max)
        oxygen_ramp = smooth_ramp(u.mass_fractions.oxygen[cell_id], p.species_min_mass_fraction_ramp_activation, p.species_max_mass_fraction_ramp_max)

        min_ramp = min(ethylene_ramp, oxygen_ramp)
    else
        if u.mass_fractions.ethylene[cell_id] <= 1e-10 || u.mass_fractions.oxygen[cell_id] <= 1e-10
            return
        end
    end

    rate_constant = p.westbrook_dryer_step_1_pre_exponential_factor * exp(-p.westbrook_dryer_step_1_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.ethylene[cell_id]^0.1 * u.molar_concentrations.oxygen[cell_id]^1.65

    #=
    oxygen_mass_required_per_ethylene_mass = 2 * u.molecular_weights.oxygen[cell_id] / u.molecular_weights.ethylene[cell_id]
    oxygen_as_ethylene_equivalent_mass_density = u.species_densities.oxygen[cell_id] / oxygen_mass_required_per_ethylene_mass
    limiting_ethylene_equivalent_mass_density = min(u.species_densities.ethylene[cell_id], oxygen_as_ethylene_equivalent_mass_density)

    turbulence_limited_ethylene_mass_rate_density = turbulence_limited_reactant_mass_rate_density(u, p, cell_id, limiting_ethylene_equivalent_mass_density)
    turbulence_limited_rate = turbulence_limited_ethylene_mass_rate_density / u.molecular_weights.ethylene[cell_id]

    rate = min(chemical_rate, turbulence_limited_rate) #[mol/(m^3*s)]
    =#

    rate = chemical_rate

    if system.use_smooth_reaction_ramp
        rate *= 3*min_ramp^2 - 2*min_ramp^3
    end

    du.mass_fractions.ethylene[cell_id] -= rate * u.molecular_weights.ethylene[cell_id] / u.density[cell_id]
    du.mass_fractions.oxygen[cell_id] -= 2 * rate * u.molecular_weights.oxygen[cell_id] / u.density[cell_id]
    du.mass_fractions.carbon_monoxide[cell_id] += 2 * rate * u.molecular_weights.carbon_monoxide[cell_id] / u.density[cell_id]
    du.mass_fractions.water[cell_id] += 2 * rate * u.molecular_weights.water[cell_id] / u.density[cell_id]

    mol_per_s = rate * geo.cell_volumes[cell_id]

    du.internal_energy_density[cell_id] -= (mol_per_s * p.westbrook_dryer_step_1_heat_of_reaction) / geo.cell_volumes[cell_id]
end

"""
CO + 0.5O2 --> CO2
westbrook_dryer_step_2_forward_pre_exponential_factor = 1.2589e10 [m^2.25 * mol^-0.75 * s^-1]
westbrook_dryer_step_2_forward_activation_energy = 167360 [J/mol]
"""
function westbrook_dryer_step_2_forward!(du, u, p, t, system, geo, cell_id)
    if system.use_smooth_reaction_ramp
        if u.mass_fractions.carbon_monoxide[cell_id] <= 0.0 || u.mass_fractions.oxygen[cell_id] <= 0.0 || u.mass_fractions.water[cell_id] <= 0.0
            return
        end
        
        carbon_monoxide_ramp = smooth_ramp(u.mass_fractions.carbon_monoxide[cell_id], p.species_min_mass_fraction_ramp_activation, p.species_max_mass_fraction_ramp_max)
        oxygen_ramp = smooth_ramp(u.mass_fractions.oxygen[cell_id], p.species_min_mass_fraction_ramp_activation, p.species_max_mass_fraction_ramp_max)
        water_ramp = smooth_ramp(u.mass_fractions.water[cell_id], p.species_min_mass_fraction_ramp_activation, p.species_max_mass_fraction_ramp_max)

        min_ramp = min(carbon_monoxide_ramp, oxygen_ramp, water_ramp)
    else
        if u.mass_fractions.carbon_monoxide[cell_id] <= 1e-10 || u.mass_fractions.oxygen[cell_id] <= 1e-10 || u.mass_fractions.water[cell_id] <= 1e-10
            return
        end
    end

    rate_constant = p.westbrook_dryer_step_2_forward_pre_exponential_factor * exp(-p.westbrook_dryer_step_2_forward_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.carbon_monoxide[cell_id] * u.molar_concentrations.water[cell_id]^0.5 * u.molar_concentrations.oxygen[cell_id]^0.25

    #=
    oxygen_mass_required_per_carbon_monoxide_mass = 0.5 * u.molecular_weights.oxygen[cell_id] / u.molecular_weights.carbon_monoxide[cell_id]
    oxygen_as_carbon_monoxide_equivalent_mass_density = u.species_densities.oxygen[cell_id] / oxygen_mass_required_per_carbon_monoxide_mass
    limiting_carbon_monoxide_equivalent_mass_density = min(u.species_densities.carbon_monoxide[cell_id], oxygen_as_carbon_monoxide_equivalent_mass_density)

    turbulence_limited_carbon_monoxide_mass_rate_density = turbulence_limited_reactant_mass_rate_density(u, p, cell_id, limiting_carbon_monoxide_equivalent_mass_density)
    turbulence_limited_rate = turbulence_limited_carbon_monoxide_mass_rate_density / u.molecular_weights.carbon_monoxide[cell_id]

    rate = min(chemical_rate, turbulence_limited_rate) #[mol/(m^3*s)]
    =#

    rate = chemical_rate

    if system.use_smooth_reaction_ramp
        rate *= 3*min_ramp^2 - 2*min_ramp^3
    end

    du.mass_fractions.carbon_monoxide[cell_id] -= rate * u.molecular_weights.carbon_monoxide[cell_id] / u.density[cell_id]
    du.mass_fractions.oxygen[cell_id] -= 0.5 * rate * u.molecular_weights.oxygen[cell_id] / u.density[cell_id]
    du.mass_fractions.carbon_dioxide[cell_id] += rate * u.molecular_weights.carbon_dioxide[cell_id] / u.density[cell_id]

    mol_per_s = rate * geo.cell_volumes[cell_id]
    
    du.internal_energy_density[cell_id] -= (mol_per_s * p.westbrook_dryer_step_2_forward_heat_of_reaction) / geo.cell_volumes[cell_id]
end

"""
CO2 --> CO + 0.5O2
westbrook_dryer_step_2_reverse_pre_exponential_factor = 5.0000e8 [s^-1]
westbrook_dryer_step_2_reverse_activation_energy = 167360 [J/mol]
"""
function westbrook_dryer_step_2_reverse!(du, u, p, t, system, geo, cell_id)
    if system.use_smooth_reaction_ramp
        if u.mass_fractions.carbon_dioxide[cell_id] <= 0.0
            return
        end

        carbon_dioxide_ramp = smooth_ramp(u.mass_fractions.carbon_dioxide[cell_id], p.species_min_mass_fraction_ramp_activation, p.species_max_mass_fraction_ramp_max)

        min_ramp = min(carbon_dioxide_ramp)
    else
        if u.mass_fractions.carbon_dioxide[cell_id] <= 1e-10
            return
        end
    end

    rate_constant = p.westbrook_dryer_step_2_reverse_pre_exponential_factor * exp(-p.westbrook_dryer_step_2_reverse_activation_energy / (u.R_gas[cell_id] * u.temperature[cell_id]))

    chemical_rate = rate_constant * u.molar_concentrations.carbon_dioxide[cell_id]

    # CO2 dissociation is a homogeneous unimolecular step. It does not require two reactants to mix, so the chemical rate is not turbulence-limited.
    rate = chemical_rate #[mol/(m^3*s)]

    if system.use_smooth_reaction_ramp
        rate *= 3*min_ramp^2 - 2*min_ramp^3
    end

    du.mass_fractions.carbon_dioxide[cell_id] -= rate * u.molecular_weights.carbon_dioxide[cell_id] / u.density[cell_id]
    du.mass_fractions.carbon_monoxide[cell_id] += rate * u.molecular_weights.carbon_monoxide[cell_id] / u.density[cell_id]
    du.mass_fractions.oxygen[cell_id] += 0.5 * rate * u.molecular_weights.oxygen[cell_id] / u.density[cell_id]

    mol_per_s = rate * geo.cell_volumes[cell_id]
    
    du.internal_energy_density[cell_id] -= (mol_per_s * p.westbrook_dryer_step_2_reverse_heat_of_reaction) / geo.cell_volumes[cell_id]
end

function combust_ethylene_and_oxygen_westbrook_dryer!(du, u, p, t, system, geo, cell_id)
    westbrook_dryer_step_1!(du, u, p, t, system, geo, cell_id)
    westbrook_dryer_step_2_forward!(du, u, p, t, system, geo, cell_id)
    westbrook_dryer_step_2_reverse!(du, u, p, t, system, geo, cell_id)
end
