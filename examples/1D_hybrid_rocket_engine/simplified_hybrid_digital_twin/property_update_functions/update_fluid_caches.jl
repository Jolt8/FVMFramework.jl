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

function ideal_eos!(du, u, p, t, system, geo, cell_id)
    universal_gas_constant = 8.314462618

    u.pressure[cell_id] = u.density[cell_id] * universal_gas_constant * u.temperature[cell_id] / u.mw_avg[cell_id]
end

function clapeyron_eos!(du, u, p, t, system, geo, cell_id)
    result = vt_flash(system.model, geo.cell_volumes[cell_id], u.temperature[cell_id], system.moles_vec[cell_id], system.clapeyron_method)

    u.pressure[cell_id] = pressure(system.model, result)
end

function cea_property_update!(du, u, p, t, system, geo, cell_id)
    #=
    NOTE:
        - Since we require species mass fractions to get pressure using a Peng-Robinson EOS and we require pressure to get species mass fractions through the closure
        - I'm going to just assume that an ideal gas pressure is right
        - This isn't very rigorous and will probably lead to innaccuracies, but it's good enough for now without having to introduce a NonlinearSolver
    =#

    universal_gas_constant = 8.314462618

    ideal_gas_pressure = u.density[cell_id] * universal_gas_constant * u.temperature[cell_id] / u.mw_avg[cell_id]
    
    cea_parameters = (ideal_gas_pressure, p.oxidizer_to_fuel_ratio, p.expansion_ratio)

    u.cp[cell_id] = chamber_cp_interp(cea_parameters)
    gamma = chamber_gamma_interp(cea_parameters)
    u.cv[cell_id] = u.cp[cell_id] / gamma

    u.mu[cell_id] = chamber_mu_interp(cea_parameters)
    u.k[cell_id] = chamber_k_interp(cea_parameters)
    u.temperature[cell_id] = chamber_T_comb_interp(cea_parameters)

    u.peng_robinson_pressure_to_ideal_gas_pressure_ratio[cell_id] = u.pressure[cell_id] / ideal_gas_pressure
end

function fluid_property_update_ideal!(du, u, p, t, system, geo, cell_id)
    mw_avg!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_and_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    cea_property_update!(du, u, p, t, system, geo, cell_id)
    ideal_eos!(du, u, p, t, system, geo, cell_id)
end

function fluid_property_update_clapeyron!(du, u, p, t, system, geo, cell_id)
    mw_avg!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_and_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    cea_property_update!(du, u, p, t, system, geo, cell_id)
    clapeyron_eos!(du, u, p, t, system, geo, cell_id)
end