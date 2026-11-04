# Igniter puck burning in the pre-chamber.
#
# The puck burns at a constant rate over igniter_burn_time, and only in a hot fire. Its products are a sugar-rocket
# mixture: about 54% solids by mass and 46% gas. The solids take up no volume and add no pressure, so only the gas
# products enter the pre-chamber gas. All of the heat of combustion still goes into the gas, because the hot solid
# particles give up their heat to it almost immediately (their own heat capacity is neglected, which makes the gas
# slightly hotter than reality).


"""
Mass burn rate of the igniter puck (kg/s, positive while burning):

    ṁ_igniter = m_igniter,0 / t_burn        while t < t_burn and the puck has mass left, otherwise 0

m_igniter,0: initial puck mass (kg), system.initial_igniter_mass; t_burn: igniter_burn_time (s).
The switch-off at t_burn is a jump, so igniter_burn_time is in tstops to land a step exactly on it.
"""
function igniter_mass_burn_rate(u, p, t, system)
    if !system.hot_fire
        return zero(u.igniter_mass[1])
    end

    if t >= p.igniter_burn_time || u.igniter_mass[1] <= 0
        return zero(u.igniter_mass[1])
    end

    return system.initial_igniter_mass / p.igniter_burn_time + zero(u.igniter_mass[1])
end


"""
Burn the igniter puck into the pre-chamber gas.

Gas products, as mass yields per kg of igniter burned (they sum to the gas yield, about 0.56; the rest is solids):
    water 0.17, carbon dioxide 0.10, carbon monoxide 0.18, nitrogen 0.11
    (hydrogen, ~0.007, isn't a species yet; add it here if it's added to the model)

    ṁ_gas = ṁ_igniter * gas_yield
    mass:     dρ/dt   = ṁ_gas / V
    species:  dY_i/dt = (y_i - Y_i) ṁ_gas / M        y_i = yield_i / gas_yield, the gas product composition
    energy:   M de/dt = ṁ_igniter ΔH_c + ṁ_gas (h_products - e)

ṁ_igniter: igniter burn rate (kg/s); ṁ_gas: gas product mass flow (kg/s); V, M: pre-chamber volume (m^3) and gas
mass (kg); ΔH_c: igniter_heat_of_combustion (J/kg), the heat released per kg of igniter burned; h_products: the
products' sensible enthalpy at the 298 K reference ΔH_c is quoted at, cp * 298.15 K (J/kg); e = cv T: the cell's
specific internal energy (J/kg). The last term is the same mixing form as every other inflow into a cell: the new
gas arrives at the reference temperature and has to be brought up to the cell's.
"""
function burn_igniter_puck!(du, u, p, t, system, geo, pre_chamber_cell_id)
    igniter_burn_rate = igniter_mass_burn_rate(u, p, t, system)

    du.igniter_mass[1] -= igniter_burn_rate


    # Gas product yields per kg of igniter burned
    water_yield = 0.17
    carbon_dioxide_yield = 0.10
    carbon_monoxide_yield = 0.18
    nitrogen_yield = 0.11

    gas_yield = water_yield + carbon_dioxide_yield + carbon_monoxide_yield + nitrogen_yield

    gas_product_mass_flow = igniter_burn_rate * gas_yield


    # Mass: only the gas products enter the gas
    pre_chamber_volume = u.cell_volume[pre_chamber_cell_id]

    pre_chamber_gas_mass = u.density[pre_chamber_cell_id] * pre_chamber_volume

    du.density[pre_chamber_cell_id] += gas_product_mass_flow / pre_chamber_volume


    # Species: dilute every species by the incoming gas, then add each product's share of it
    for_fields!(du.mass_fractions, u.mass_fractions) do species, du_mass_fractions, u_mass_fractions
        du_mass_fractions[species[pre_chamber_cell_id]] -= u_mass_fractions[species[pre_chamber_cell_id]] * gas_product_mass_flow / pre_chamber_gas_mass
    end

    du.mass_fractions.water[pre_chamber_cell_id] += (water_yield / gas_yield) * gas_product_mass_flow / pre_chamber_gas_mass
    du.mass_fractions.carbon_dioxide[pre_chamber_cell_id] += (carbon_dioxide_yield / gas_yield) * gas_product_mass_flow / pre_chamber_gas_mass
    du.mass_fractions.carbon_monoxide[pre_chamber_cell_id] += (carbon_monoxide_yield / gas_yield) * gas_product_mass_flow / pre_chamber_gas_mass
    du.mass_fractions.nitrogen[pre_chamber_cell_id] += (nitrogen_yield / gas_yield) * gas_product_mass_flow / pre_chamber_gas_mass


    # Energy, added to du.heat (W) and converted to temperature by cap_heat_flux_to_temp_change!
    heat_of_combustion_power = igniter_burn_rate * p.igniter_heat_of_combustion

    reference_temperature = 298.15

    gas_product_specific_enthalpy = u.cp[pre_chamber_cell_id] * reference_temperature

    pre_chamber_specific_internal_energy = u.cv[pre_chamber_cell_id] * u.temperature[pre_chamber_cell_id]

    gas_product_inflow_power = gas_product_mass_flow * (gas_product_specific_enthalpy - pre_chamber_specific_internal_energy)

    du.heat[pre_chamber_cell_id] += heat_of_combustion_power + gas_product_inflow_power
end
