# Stop the burn when the fuel port reaches the outside diameter of the grain.
# Only the positive crossing is active because regression increases port_diameter.
function port_diameter_limit_expanded(u, t, integrator, u_axes, p_axes, local_update_state!, local_oxidizer_model, local_chamber_model)
    u_named = ComponentVector(u, u_axes)
    p_named = ComponentVector(integrator.p, p_axes)

    local_update_state!([0.0], u_named, p_named, t, local_oxidizer_model, local_chamber_model)

    #@show u_named.port_diameter - p_named.final_fuel_grain_void_diameter

    #As soon as this hits zero, the simulation is stopped
    u_named.port_diameter - p_named.final_fuel_grain_void_diameter
end

port_diameter_limit = let
    u_axes_local = u_axes
    p_axes_local = p_axes
    update_state_local! = update_state!
    oxidizer_model_local = oxidizer_model
    chamber_model_local = chamber_model

    (u, t, integrator) -> port_diameter_limit_expanded(u, t, integrator, u_axes_local, p_axes_local, update_state_local!, oxidizer_model_local, chamber_model_local)
end

port_diameter_termination_cb = ContinuousCallback(
    port_diameter_limit,
    terminate!,
    nothing;
    save_positions = (true, false),
)

#Chamber Pressure Callback
function chamber_pressure_limit_expanded(u, t, integrator, u_axes, p_axes, local_update_state!, local_oxidizer_model, local_chamber_model)
    u_named = ComponentVector(u, u_axes)
    p_named = ComponentVector(integrator.p, p_axes)

    local_update_state!([0.0], u_named, p_named, t, local_oxidizer_model, local_chamber_model)

    #@show p_named.chamber_pressure - 500_000.0
    #@show 1e6 - p_named.chamber_pressure

    #As soon as this hits zero, the simulation is stopped
    #We stop simulating after the chamber pressure reaches this value
    100_000 - p_named.chamber_pressure #Pa
end

chamber_pressure_limit = let
    u_axes_local = u_axes
    p_axes_local = p_axes
    update_state_local! = update_state!
    oxidizer_model_local = oxidizer_model
    chamber_model_local = chamber_model

    (u, t, integrator) -> chamber_pressure_limit_expanded(u, t, integrator, u_axes_local, p_axes_local, update_state_local!, oxidizer_model_local, chamber_model_local)
end

chamber_pressure_termination_cb = ContinuousCallback(
    chamber_pressure_limit,
    terminate!,
    nothing;
    save_positions = (true, false),
)

fuel_burnout_condition = let
    u_axes_local = u_axes
    p_axes_local = p_axes

    function (u, t, integrator)
        u_named = ComponentVector(u, u_axes_local)
        p_named = ComponentVector(integrator.p, p_axes_local)

        final_diameter = p_named.u0_fuel_grain_void_diameter + p_named.additional_fuel_grain_void_diameter

        # Disable this event after it has fired.
        if p_named.burned_out == 1.0
            return 1.0
        end

        return u_named.port_diameter - final_diameter
    end
end

fuel_burnout_affect! = let
    u_axes_local = u_axes
    p_axes_local = p_axes

    function (integrator)
        u_named = ComponentVector(integrator.u, u_axes_local)
        p_named = ComponentVector(integrator.p, p_axes_local)

        p_named.burned_out = 1.0
    end
end

fuel_burnout_cb = ContinuousCallback(
    fuel_burnout_condition,
    fuel_burnout_affect!,
    nothing;
    save_positions = (true, true),
)