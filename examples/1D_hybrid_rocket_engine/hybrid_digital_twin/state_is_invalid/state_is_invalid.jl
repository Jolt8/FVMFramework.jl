# Passed to the solver as isoutofdomain: a step whose end state fails this check is rejected and retried with a
# smaller step, instead of being accepted and carried forward.
#
# The checks are physical: every cell must hold positive mass at a positive temperature, and nothing may be NaN.
# NaN appears when a property update met a non-physical trial state (see set_nitrous_oxide_properties_invalid!).
function state_is_invalid(u, p, t, system)
    u_named = ComponentVector(u, system.state_axes)

    for cell_id in eachindex(u_named.density)
        cell_density = u_named.density[cell_id]
        cell_temperature = u_named.temperature[cell_id]

        if !isfinite(cell_density) || !isfinite(cell_temperature)
            return true
        end

        if cell_density <= 0 || cell_temperature <= 0
            return true
        end
    end

    # Any other non-finite state (mass fractions, port radius, igniter mass) also means the step went wrong
    if !all(isfinite, u)
        return true
    end

    return false
end
