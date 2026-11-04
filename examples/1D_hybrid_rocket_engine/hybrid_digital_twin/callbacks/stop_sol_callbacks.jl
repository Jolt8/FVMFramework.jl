# Stop the burn when the fuel port reaches the outside diameter of the grain.
# Only the positive crossing is active because regression increases port_diameter.
function port_diameter_limit_expanded(u, t, integrator, system)
    u_named = ComponentVector(u, system.state_axes)
    p_named = ComponentVector(integrator.p, system.p_axes)

    # Stop the solver once the largest port radius at any station reaches the radius that still leaves the desired
    # web of fuel against the liner. Only the fuel grain cells hold a port radius; every other cell's entry is zero,
    # so the stations are looped over instead of the whole vector.
    largest_fuel_port_radius = 0.0

    for station in system.fuel_port_geometry.stations
        largest_fuel_port_radius = max(largest_fuel_port_radius, u_named.fuel_port_radius[station.fuel_grain_cell_id])
    end

    # Negative while there's fuel to spare, crossing zero upward when the web limit is reached
    # system.p_axes stores every parameter as a 1-element vector, hence the [1]
    return largest_fuel_port_radius - p_named.max_fuel_port_radius[1]
end

port_diameter_limit = (u, t, integrator) -> port_diameter_limit_expanded(u, t, integrator, system)

port_diameter_termination_cb = ContinuousCallback(
    port_diameter_limit,
    terminate!,
    nothing;
    save_positions = (true, false),
)

function min_chamber_pressure_limit_expanded(u, t, integrator, system, geo)
    u_named = ComponentVector(u, system.state_axes)
    p_named = ComponentVector(integrator.p, system.p_axes)

    pre_chamber_cell_id = get_cellset(geo, "pre_chamber")[1]

    # The callback only sees the states, not the caches, so the pre-chamber pressure is recomputed here with the
    # same ideal gas law as ideal_eos!:
    #     P = ρ R_u T / MW_mix,    1 / MW_mix = Σ Y_i / MW_i
    # ρ: pre-chamber gas density (kg/m^3); R_u: universal gas constant (J/(mol*K)); T: temperature (K);
    # Y_i: species mass fractions (kg/kg); MW_i: species molecular weights (kg/mol), read from the region properties.
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    inverse_mixture_molecular_weight = Ref(0.0)

    for_fields!(u_named.mass_fractions, properties.molecular_weights) do species, u_mass_fractions, molecular_weights
        inverse_mixture_molecular_weight[] += u_mass_fractions[species[pre_chamber_cell_id]] / molecular_weights[species[pre_chamber_cell_id]]
    end

    universal_gas_constant = 8.314462618

    pre_chamber_pressure = u_named.density[pre_chamber_cell_id] * universal_gas_constant * u_named.temperature[pre_chamber_cell_id] * inverse_mixture_molecular_weight[]

    if t < 2.0 #we make sure that all rockets can burn for at least 2 seconds before applying this
        return 1.0
    else
        # system.p_axes stores every parameter as a 1-element vector, hence the [1]
        return p_named.min_chamber_pressure[1] - pre_chamber_pressure
    end
end

min_chamber_pressure_limit = (u, t, integrator) -> min_chamber_pressure_limit_expanded(u, t, integrator, system, geo)

min_chamber_pressure_termination_cb = ContinuousCallback(
    min_chamber_pressure_limit,
    terminate!,
    nothing;
    save_positions = (true, false),
)