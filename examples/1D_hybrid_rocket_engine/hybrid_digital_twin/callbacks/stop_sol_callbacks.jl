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
    return largest_fuel_port_radius - p_named.max_fuel_port_radius
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

    result = vt_flash(system.model, geo.cell_volumes[pre_chamber_cell_id], u_named.temperature[pre_chamber_cell_id], system.moles_vec[pre_chamber_cell_id], system.clapeyron_method)

    pre_chamber_pressure = pressure(system.model, result)

    if t < 2.0 #we make sure that all rockets can burn for at least 2 seconds before applying this
        return 1.0
    else
        return p_named.min_chamber_pressure - pre_chamber_pressure
    end
end

min_chamber_pressure_limit = (u, t, integrator) -> min_chamber_pressure_limit_expanded(u, t, integrator, system, geo)

min_chamber_pressure_termination_cb = ContinuousCallback(
    min_chamber_pressure_limit,
    terminate!,
    nothing;
    save_positions = (true, false),
)