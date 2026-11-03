function mw_avg!(du, u, p, t, system, geo, cell_id)
    inv_mw = 0.0

    for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
        inv_mw += u_mass_fractions[species[cell_id]] / u_molecular_weights[species[cell_id]]
    end

    u.mw_avg[cell_id] = 1.0 / inv_mw
end


function molar_concentrations!(du, u, p, t, system, geo, cell_id)
    for_fields!(u.mass_fractions, u.molecular_weights, u.molar_concentrations) do species, u_mass_fractions, u_molecular_weights, u_molar_concentrations
        u_molar_concentrations[species[cell_id]] = (u.density[cell_id] * u_mass_fractions[species[cell_id]]) / u_molecular_weights[species[cell_id]]
    end
end

function ideal_eos!(du, u, p, t, system, geo, cell_id)
    universal_gas_constant = 8.314462618

    u.pressure[cell_id] = u.density[cell_id] * universal_gas_constant * u.temperature[cell_id] / u.mw_avg[cell_id]
end

function ideal_gas_heat_capacities!(du, u, p, t, system, geo, cell_id)
    # Mixture cp and cv from the ReidIdeal polynomials in system.model, using the cell's mole fractions.
    # The Reid polynomials are only fitted up to roughly 1000-1500 K. Above that they diverge
    # (N2 cp is ~0 at 2000 K and negative at 3000 K), so the property temperature is clamped.
    # At 3000 K this underestimates cp by roughly 10% for N2, 14% for CO2 and 23% for H2O.
    # Replace with NASA 7-coefficient polynomials before trusting chamber temperatures.
    maximum_reid_polynomial_temperature = 1000.0

    property_temperature = min(u.temperature[cell_id], maximum_reid_polynomial_temperature)

    cell_mole_fractions = system.mole_fractions_vec[cell_id]

    u.cp[cell_id] = mass_isobaric_heat_capacity(system.model, u.pressure[cell_id], property_temperature, cell_mole_fractions)
    u.cv[cell_id] = mass_isochoric_heat_capacity(system.model, u.pressure[cell_id], property_temperature, cell_mole_fractions)
end

function clapeyron_eos!(du, u, p, t, system, geo, cell_id)
    result = vt_flash(system.model, geo.cell_volumes[cell_id], u.temperature[cell_id], system.moles_vec[cell_id], system.clapeyron_method)

    #Clapeyron's property functions take (model, result). The mass_ versions return J/(kg*K) to match the cp and cv caches.
    #The aggregate cp/cv calls error when more than one phase is present, so the phase-weighted helper is used.
    #This is only correct for single-phase results (always the case for ReidIdeal). A two-phase cell needs the
    #latent-heat term; see two_phase_mass_isochoric_heat_capacity in the injector delivery train physics.
    u.pressure[cell_id] = pressure(system.model, result)
    u.cp[cell_id] = phase_weighted_mass_property(mass_isobaric_heat_capacity, system.model, result)
    u.cv[cell_id] = phase_weighted_mass_property(mass_isochoric_heat_capacity, system.model, result)

    u.nitrous_oxide_vapor_quality[cell_id] = flash_vapor_quality(result)
end


function fluid_property_update_ideal!(du, u, p, t, system, geo, cell_id)
    mw_avg!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_and_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    molar_concentrations!(du, u, p, t, system, geo, cell_id)
    ideal_eos!(du, u, p, t, system, geo, cell_id)
    ideal_gas_heat_capacities!(du, u, p, t, system, geo, cell_id) #cp and cv were never set for ideal-gas cells before
end

function fluid_property_update_clapeyron!(du, u, p, t, system, geo, cell_id)
    mw_avg!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_and_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    molar_concentrations!(du, u, p, t, system, geo, cell_id)
    clapeyron_eos!(du, u, p, t, system, geo, cell_id)
end