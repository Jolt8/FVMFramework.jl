#we should probably move to using cell instead of cell_id inside these functions and cell_a cell_b 
#however in the main loop inside the op function, we should still use cell because cell can sometimes be used
function net_reaction_rate!(du, u, cell, reaction, reaction_name, kr_A, kr_Ea) #reaction is u.reactions[reaction_id]
    kf = (reaction.kf_A[1] * exp(-reaction.kf_Ea[1] / (R_gas * u.temp[cell])))
    kr = (kr_A * exp(-kr_Ea / (R_gas * u.temp[cell])))

    forward_term = 1.0
    map(keys(reaction.reactant_stoich_coeffs)) do species_name
        concentration = u.molar_concentrations[species_name][cell]
        stoich_coeff = reaction.reactant_stoich_coeffs[species_name][1]
        forward_term *= concentration^stoich_coeff
    end

    reverse_term = 1.0
    map(keys(reaction.product_stoich_coeffs)) do species_name
        concentration = u.molar_concentrations[species_name][cell]
        stoich_coeff = reaction.product_stoich_coeffs[species_name][1]
        reverse_term *= concentration^stoich_coeff
    end

    u.net_rates[reaction_name][1] = ((kf * forward_term) - (kr * reverse_term))
end

#I think I know what is causing the slowdown here, the reaciton dynamic dispatch seems to be causing some issues
#whenever a reactin property is grabbed using reaciton.something, it causes GC


#after rewriting this, I think it would be better to have u.species = ComponentVector(molecular_weight = 10.0, mass_fraction = 0.4, molar_concentration = 0.3, diffusion_coefficient = 1e-5, etc.)
#also, I'm absolutely loving writing physics functions while only referencing du and u
#it's absolutely amazing to write out the equations almost in their base form and then only have to keep track where each variable is stored 
#however in the future, I'm sure this will become more and more concrete and it will only become more and more intuitive
function power_law_react_cell!(du, u, cell_id, reaction, reaction_name, vol)
    #find reverse pre exponential_factor
    K_ref = K_gibbs_free(u, cell_id, reaction) 

    kr_A = (reaction.kf_A[1] / K_ref) * exp(-reaction.heat_of_reaction[1] / (R_gas * u.temp[cell_id])) 

    #find reverse Ea
    kr_Ea = reaction.kf_Ea[1] - reaction.heat_of_reaction[1] 

    #mutates reaction.net_rate
    net_reaction_rate!(du, u, cell_id, reaction, reaction_name, kr_A, kr_Ea) * u.reactions_kg_cat[reaction_name][cell_id] 
    #rate returned by net_reactions_rate is in mol / (kg_cat * s), so we times by the cell's kg_cat / m3
    #rate is now in mol / (m3 * s)

    map(keys(du.molar_concentrations)) do species_name
        du.molar_concentrations[species_name][cell_id] += u.net_rates[reaction_name][1] * reaction.stoich_coeffs[species_name][1]
        # rate (mol/(m3*s)) * MW (g/mol) / rho (g/m3) = unitless/s
        du.mass_fractions[species_name][cell_id] += (u.net_rates[reaction_name][1] * reaction.stoich_coeffs[species_name][1] * u.molecular_weights[species_name][1]) / u.rho[cell_id]
    end

    du.heat[cell_id] += u.net_rates[reaction_name][1] * (-reaction.heat_of_reaction[1]) * vol 
    # rate (mol/(m3*s)) * H (J/mol) * vol (m3) = J/s = Watts
end