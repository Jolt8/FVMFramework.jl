function unit_independent_max(val, min_val)
    return max(ustrip(val), min_val) * unit(val)
end

function van_t_hoff(A, dH, T, R)
    #K = A * exp(dH/RT)
    return A * exp(-dH / (R * T))
end

#this is for Peppley-Amphlett methanol steam reforming (MSR) kinetics

"""
    PAM_reforming_react_cell!(du, u, cell_id, vol)

    Calculates the net reaction rates using the Peppley-Amphlett methanol steam reforming kinetics.
    Note that this is likely not accurate beacuse I think some unit discrepancies between this package and the original source causes it to be innacurate.
    Thus, this likely needs more testing and more validation, but I'm not working on methanol reforming right now, so I'll just leave it for now.
    This does serve as a good example of the kind of complicated physics that are possible in this package.

    
    Parameters:
    - du - The residual array to be updated
    - u - The state array
    - cell_id - The ID of the cell to calculate the reaction rates for
    - vol - The volume of the cell
"""
function PAM_reforming_react_cell!(du, u, cell_id, vol)
    mw_avg!(u, cell_id)
    rho_ideal!(u, cell_id)
    molar_concentrations!(u, cell_id)
    
    #we divide by 1e-5 because PAM parameters are in bar
    #we add a tiny epsilon (1e-12) to prevent division by zero in the denominators/roots
    #during the start-up phase when concentrations are close to 0.0

    conversion_factor = ((u.R_gas[cell_id] * u.temp[cell_id]) * 1e-5)
    P_CH3OH = unit_independent_max(u.molar_concentrations.methanol[cell_id] * conversion_factor, 1e-20)
    P_H2O = unit_independent_max(u.molar_concentrations.water[cell_id] * conversion_factor, 1e-20)
    P_CO = unit_independent_max(u.molar_concentrations.carbon_monoxide[cell_id] * conversion_factor, 1e-20)
    P_H2 = unit_independent_max(u.molar_concentrations.hydrogen[cell_id] * conversion_factor, 1e-20)
    P_CO2 = unit_independent_max(u.molar_concentrations.carbon_dioxide[cell_id] * conversion_factor, 1e-20)

    #P_CH3OH = u.molar_concentrations.methanol[cell_id] * conversion_factor
    #P_H2O = u.molar_concentrations.water[cell_id] * conversion_factor
    #P_CO = u.molar_concentrations.carbon_monoxide[cell_id] * conversion_factor
    #P_H2 = u.molar_concentrations.hydrogen[cell_id] * conversion_factor
    #P_CO2 = u.molar_concentrations.carbon_dioxide[cell_id] * conversion_factor

    #this is weird, I get different results depending on whether or not I use max or not
    #if I don't use it, for some reason the hydrogen conversion goes down a lot and the carbon monoxide concentration 
    #decreases along the length of the reactor, and so does hydrogen concentration
    #hmm, I think it's from the WGS reaction going in the reverse direction, but for some reason carbon dioxide concentration does not seem to increase
    #I think we should really look into returning the extent of each reaction per cell to see what's going on

    #if ustrip(P_H2) > 1e-8 && ustrip(P_CH3OH) > 1e-8 #ok this was causing a massive slowdown
        #=
        println("P_CH3OH: $(P_CH3OH)")
        println("P_H2O: $(P_H2O)")
        println("P_CO: $(P_CO)")
        println("P_H2: $(P_H2)")
        println("P_CO2: $(P_CO2)")
        println("")
        =#

        K_CH3O = van_t_hoff(
            u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_A.CH3O[cell_id], 
            u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_dH.CH3O[cell_id], 
            u.temp[cell_id],
            u.R_gas[cell_id]
        )
        K_HCOO = van_t_hoff(
            u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_A.HCOO[cell_id], 
            u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_dH.HCOO[cell_id], 
            u.temp[cell_id],
            u.R_gas[cell_id]
        )
        K_OH = van_t_hoff(
            u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_A.OH[cell_id], 
            u.reactions.reforming_reactions.MSR_rxn.van_t_hoff_dH.OH[cell_id], 
            u.temp[cell_id],
            u.R_gas[cell_id]
        )

        #=
        println("K_CH3O: $(K_CH3O)")
        println("K_HCOO: $(K_HCOO)")
        println("K_OH: $(K_OH)")
        println("")
        =#

        term_CH3O = K_CH3O * (P_CH3OH / sqrt(P_H2))
        term_HCOO = K_HCOO * (P_CO2 * sqrt(P_H2))
        term_OH = K_OH * (P_H2O / sqrt(P_H2))

        #=
        println("term_CH3O: $(term_CH3O)")
        println("term_HCOO: $(term_HCOO)")
        println("term_OH: $(term_OH)")
        println("")
        =#

        DEN = 1.0 + term_CH3O + term_HCOO + term_OH

        #println("DEN: $(DEN)")
        #println("")

        for_fields!(u.net_rates.reforming_reactions) do reaction_name, net_rates
            net_rates[reaction_name[cell_id]] = 0.0
        end

        #MSR
        k_MSR = u.reactions.reforming_reactions.MSR_rxn.kf_A[cell_id] * exp(-u.reactions.reforming_reactions.MSR_rxn.kf_Ea[cell_id] / (u.R_gas[cell_id] * u.temp[cell_id]))
        K_eq_MSR = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.MSR_rxn)

        driving_force = 1.0 - ((P_CO2 * P_H2^3) / (K_eq_MSR * P_CH3OH * P_H2O))

        u.net_rates.reforming_reactions.MSR_rxn[cell_id] = (k_MSR * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2)

        #=
        println("temp: $(u.temp[cell_id])")
        println("kf_Ea: $(u.reactions.reforming_reactions.MSR_rxn.kf_Ea[cell_id])")
        println("exp(-kf_Ea / (R_gas * temp)): $(exp(-u.reactions.reforming_reactions.MSR_rxn.kf_Ea[cell_id] / (u.R_gas[cell_id] * u.temp[cell_id])))")
        println("kf_A: $(u.reactions.reforming_reactions.MSR_rxn.kf_A[cell_id])")
        println("k_MSR: $(k_MSR)")    
        
        println("K_eq_MSR: $(K_eq_MSR)")
        println("driving_force: $(driving_force)")
        println("u.net_rates.reforming_reactions.MSR_rxn[1]: $(u.net_rates.reforming_reactions.MSR_rxn[1])")
        println("")
        =#

        #MD
        k_MD = u.reactions.reforming_reactions.MD_rxn.kf_A[cell_id] * exp(-u.reactions.reforming_reactions.MD_rxn.kf_Ea[cell_id] / (u.R_gas[cell_id] * u.temp[cell_id]))
        K_eq_MD = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.MD_rxn)
        driving_force = 1.0 - ((P_CO * P_H2^2) / (K_eq_MD * P_CH3OH))

        u.net_rates.reforming_reactions.MD_rxn[cell_id] = (k_MD * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2)

        #=
        println("k_MD: $(k_MD)")
        println("K_eq_MD: $(K_eq_MD)")
        println("driving_force: $(driving_force)")
        println("u.net_rates.reforming_reactions.MD_rxn[1]: $(u.net_rates.reforming_reactions.MD_rxn[1])")
        println("")
        =#
        

        #WGS
        k_WGS = u.reactions.reforming_reactions.WGS_rxn.kf_A[cell_id] * exp(-u.reactions.reforming_reactions.WGS_rxn.kf_Ea[cell_id] / (u.R_gas[cell_id] * u.temp[cell_id]))
        K_eq_WGS = K_gibbs_free(u, cell_id, u.reactions.reforming_reactions.WGS_rxn)
        driving_force = 1.0 - ((P_CO2 * P_H2) / (K_eq_WGS * P_CO * P_H2O))

        u.net_rates.reforming_reactions.WGS_rxn[cell_id] = (k_WGS * K_OH * ((P_CO * P_H2O) / sqrt(P_H2)) * driving_force) / (DEN^2)

        #=
        println("k_WGS: $(k_WGS)")
        println("K_eq_WGS: $(K_eq_WGS)")
        println("driving_force: $(driving_force)")
        println("u.net_rates.reforming_reactions.WGS_rxn[1]: $(u.net_rates.reforming_reactions.WGS_rxn[1])")
        println("")
        =#

        for_fields!(du.molar_concentrations) do species, du_molar_concentrations
            du_molar_concentrations[species[cell_id]] = 0.0
        end

        #MSR_rxn
        for_fields!(u.reactions.reforming_reactions.MSR_rxn.stoich_coeffs, du.molar_concentrations) do species_name, stoich_coeffs, du_molar_concentrations
            du_molar_concentrations[species_name[cell_id]] += u.net_rates.reforming_reactions.MSR_rxn[cell_id] * stoich_coeffs[species_name[cell_id]] * u.reactions_kg_cat.reforming_reactions.MSR_rxn[cell_id]
        end

        #MD_rxn
        for_fields!(u.reactions.reforming_reactions.MD_rxn.stoich_coeffs, du.molar_concentrations) do species_name, stoich_coeffs, du_molar_concentrations
            du_molar_concentrations[species_name[cell_id]] += u.net_rates.reforming_reactions.MD_rxn[cell_id] * stoich_coeffs[species_name[cell_id]] * u.reactions_kg_cat.reforming_reactions.MD_rxn[cell_id]
        end

        #WGS_rxn
        for_fields!(u.reactions.reforming_reactions.WGS_rxn.stoich_coeffs, du.molar_concentrations) do species_name, stoich_coeffs, du_molar_concentrations
            du_molar_concentrations[species_name[cell_id]] += u.net_rates.reforming_reactions.WGS_rxn[cell_id] * stoich_coeffs[species_name[cell_id]] * u.reactions_kg_cat.reforming_reactions.WGS_rxn[cell_id]
        end

        for_fields!(du.mass_fractions, du.molar_concentrations, u.molecular_weights) do species_name, du_mass_fractions, du_molar_concentrations, u_molecular_weights
            du_mass_fractions[species_name[cell_id]] += (du_molar_concentrations[species_name[cell_id]] * u_molecular_weights[species_name[cell_id]] / u.rho[cell_id])
            # rate (mol/(m3*s)) * MW (g/mol) / rho (g/m3) = unitless/s
        end

        du.heat[cell_id] += u.net_rates.reforming_reactions.MSR_rxn[cell_id] * (-u.reactions.reforming_reactions.MSR_rxn.heat_of_reaction[cell_id]) * u.reactions_kg_cat.reforming_reactions.MSR_rxn[cell_id] * vol
        du.heat[cell_id] += u.net_rates.reforming_reactions.MD_rxn[cell_id] * (-u.reactions.reforming_reactions.MD_rxn.heat_of_reaction[cell_id]) * u.reactions_kg_cat.reforming_reactions.MD_rxn[cell_id] * vol
        du.heat[cell_id] += u.net_rates.reforming_reactions.WGS_rxn[cell_id] * (-u.reactions.reforming_reactions.WGS_rxn.heat_of_reaction[cell_id]) * u.reactions_kg_cat.reforming_reactions.WGS_rxn[cell_id] * vol

        #for_fields!(u.net_rates.reforming_reactions, u.reactions.reforming_reactions.heat_of_reaction) do reforming_reaction, net_rates, heat_of_reactions
        #    du.heat[cell_id] += net_rates[reforming_reaction[cell_id]] * (-heat_of_reactions[reforming_reaction[cell_id]]) * vol 
        #    # rate (mol/(m3*s)) * H (J/mol) * vol (m3) = J/s = Watts
        #end
    #end
end