function igniter_mass_burn_rate_at_t(du, u, p, t, system, geo)
    initial_burnable_mass = system.data.initial_igniter_mass
    burn_time = p.igniter_burn_time

    if initial_burnable_mass <= zero(initial_burnable_mass) || burn_time <= zero(burn_time) ||
       t < zero(t) || t >= burn_time ||
       u.igniter_mass[1] <= zero(u.igniter_mass[1])
        return 0.0
    end

    nominal_mass_burn_rate = initial_burnable_mass / burn_time

    return nominal_mass_burn_rate
end

function update_igniter_combustion!(du, u, p, t, system, geo)
    du.igniter_mass[1] = -igniter_mass_burn_rate_at_t(du, u, p, t, system, geo)
end

function update_igniter_puck_combustion_source!(du, u, p, t, system, geo, pre_chamber_cell_id)
    if igniter_burnt(du, u, p, t, system, geo)
        return nothing
    end

    vol = geo.cell_volumes[pre_chamber_cell_id]    

    du.density[pre_chamber_cell_id] += du.igniter_mass[1] / vol
    du.heat[pre_chamber_cell_id] += du.igniter_mass[1] * p.igniter_heat_of_combustion

    #the fractions are based on the typical solid and gaseous product yield of combusting a typical sugar rocket mixture
    #the combustion products are typically made up of ~54% solid products and ~46% gaseous products by mass 
    #NOTE: the mass of solid products is ignored and all the reaction heat is added to the gas products
    #this means that this might produce more gas than theoretically possible, but it conserves the igniter's mass
    pre_chamber_mass = u.density[pre_chamber_cell_id] * vol

    du.mass_fractions.water[pre_chamber_cell_id] += (0.17 - u.mass_fractions.water[pre_chamber_cell_id]) * (du.igniter_mass[1] / pre_chamber_mass)
    du.mass_fractions.carbon_dioxide[pre_chamber_cell_id] += (0.10 - u.mass_fractions.carbon_dioxide[pre_chamber_cell_id]) * (du.igniter_mass[1] / pre_chamber_mass)
    du.mass_fractions.carbon_monoxide[pre_chamber_cell_id] += (0.18 - u.mass_fractions.carbon_monoxide[pre_chamber_cell_id]) * (du.igniter_mass[1] / pre_chamber_mass)
    du.mass_fractions.nitrogen[pre_chamber_cell_id] += (0.11 - u.mass_fractions.nitrogen[pre_chamber_cell_id]) * (du.igniter_mass[1] / pre_chamber_mass)
    #TODO: if we add hydrogen as a product later, we must include it here (but it's contribution is tiny ~ 0.7% (0.007))

    return nothing
end

function igniter_burnt(du, u, p, t, system, geo)
    burnout_mass_tolerance = 1e-9 * system.igniter.initial_burnable_mass
    return t >= system.igniter.burn_time || u.igniter_mass[1] <= burnout_mass_tolerance
end
