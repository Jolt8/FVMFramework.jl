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
    cell_mass = u.density[cell_id] * geo.cell_volumes[cell_id]

    u.pressure[cell_id] = ((cell_mass * u.mw_avg[cell_id]) * u.R_gas[cell_id] * u.temperature[cell_id]) / geo.cell_volumes[cell_id]
end

function clapeyron_eos!(du, u, p, t, system, geo, cell_id)
    result = vt_flash(mode, geo.cell_volumes[cell_id], u.temperature[cell_id], view(system.moles_vec, cell_id))

    u.pressure[cell_id] = pressure(result, system.model)
    u.cp[cell_id] = isobaric_heat_capacity(result, system.model)
    u.cv[cell_id] = isochoric_heat_capacity(result, system.model)
    
    u.vapor_quality[cell_id] = result.fractions[2]
end


function fluid_property_update_ideal!(du, u, p, t, system, geo, cell_id)
    mw_avg!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_and_mole_fractions_vec(du, u, p, t, system, geo, cell_id)
    molar_concentrations!(du, u, p, t, system, geo, cell_id)
    ideal_eos!(du, u, p, t, system, geo, cell_id)
end

function fluid_property_update_clapeyron!(du, u, p, t, system, geo, cell_id)
    mw_avg!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_and_mole_fractions_vec(du, u, p, t, system, geo, cell_id)
    molar_concentrations!(du, u, p, t, system, geo, cell_id)
    clapeyron_eos!(du, u, p, t, system, geo, cell_id)
end