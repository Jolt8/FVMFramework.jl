function mw_avg!(du, u, p, t, system, geo, cell_id)
    # A Ref, like species_idx in the moles vector: a plain variable reassigned inside the do-block closure gets
    # wrapped in an untyped Core.Box, which allocates and dispatches dynamically on every species.
    # zero(u.mw_avg[cell_id]) starts it at the right number type (Float64, Dual, or tracer).
    inverse_mw_sum = Ref(zero(u.mw_avg[cell_id]))

    for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
        inverse_mw_sum[] += u_mass_fractions[species[cell_id]] / u_molecular_weights[species[cell_id]]
    end

    u.mw_avg[cell_id] = 1.0 / inverse_mw_sum[]
end

function cea_ideal_eos!(du, u, p, t, system, geo, cell_id)
    universal_gas_constant = 8.314462618

    u.pressure[cell_id] = u.density[cell_id] * universal_gas_constant * u.temperature[cell_id] / u.mw_avg[cell_id]
end

function cea_clapeyron_eos!(du, u, p, t, system, geo, cell_id)
    result = vt_flash(system.model, geo.cell_volumes[cell_id], u.temperature[cell_id], system.moles_vec[cell_id])

    u.pressure[cell_id] = pressure(system.model, result)
end

function converge_cea_pressure(cea_pressure, u, p, system, geo, cell_id)
    cea_parameters = (cea_pressure, p.oxidizer_to_fuel_ratio, p.expansion_ratio)

    interp_idx = Ref(1)

    for_fields!(u.mass_fractions) do species, u_mass_fractions
        u_mass_fractions[species[cell_id]] = system.mass_fractions_interps_vec[interp_idx](cea_parameters)
        
        interp_idx[] += 1
    end

    update_cell_moles_and_mole_fractions_vec!(0, u, p, 0, system, geo, cell_id)
    
    result = vt_flash(system.model, geo.cell_volumes[cell_id], u.temperature[cell_id], system.moles_vec[cell_id])

    calculated_pressure = pressure(system.model, result)

    return cea_pressure - calculated_pressure
end

function cea_property_update!(du, u, p, t, system, geo, cell_id)
    converged_pressure = find_zero(cea_pressure -> converge_cea_pressure(cea_pressure, u, p, system, geo, cell_id), (100_000, 5_000_000))

    u.pressure[cell_id] = converged_pressure

    cea_parameters = (converged_pressure, p.oxidizer_to_fuel_ratio, p.expansion_ratio)
    
    u.cp[cell_id] = system.chamber_cp_interp(cea_parameters)
    gamma = system.chamber_gamma_interp(cea_parameters)
    u.cv[cell_id] = u.cp[cell_id] / gamma

    u.mu[cell_id] = system.chamber_mu_interp(cea_parameters)
    u.k[cell_id] = system.chamber_k_interp(cea_parameters)
    u.temperature[cell_id] = system.chamber_T_comb_interp(cea_parameters)

    interp_idx = Ref(1)

    for_fields!(u.mass_fractions) do species, u_mass_fractions
        u_mass_fractions[species[cell_id]] = mass_fractions_interps_vec[interp_idx](cea_parameters)
        
        interp_idx[] += 1
    end
end

function cea_fluid_property_update_ideal!(du, u, p, t, system, geo, cell_id)
    cea_property_update!(du, u, p, t, system, geo, cell_id)
    mw_avg!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_and_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    #cea_ideal_eos!(du, u, p, t, system, geo, cell_id)
end

function cea_fluid_property_update_clapeyron!(du, u, p, t, system, geo, cell_id)
    cea_property_update!(du, u, p, t, system, geo, cell_id)
    mw_avg!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_and_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    #cea_clapeyron_eos!(du, u, p, t, system, geo, cell_id)
end