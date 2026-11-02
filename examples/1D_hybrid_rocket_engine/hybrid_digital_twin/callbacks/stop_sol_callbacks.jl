# Stop the burn when the fuel port reaches the outside diameter of the grain.
# Only the positive crossing is active because regression increases port_diameter.
function port_diameter_limit_expanded(u, t, integrator, system)
    u_named = ComponentVector(u, system.state_axes)
    p_named = ComponentVector(integrator.p, system.p_axes)

    # Stop the solver if any fuel port radius is below the desired residual thickness to protect the motor
    min_fuel_port_radius = 0.0

    for i in eachindex(u_named.fuel_port_radius)
        min_fuel_port_radius = min(min_fuel_port_radius, u_named.fuel_port_radius[i])
    end

    return min_fuel_port_radius - p_named.min_fuel_port_radius
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

    pre_chamber_cell_id = get_cells(geo, "pre_chamber")[1]

    result = vt_flash(system.model, geo.cell_volumes[pre_chamber_cell_id], u_named.temperature[pre_chamber_cell_id], system.moles_vec[pre_chamber_cell_id], system.clapeyron_method)

    pre_chamber_pressure = pressure(result, system.model)

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