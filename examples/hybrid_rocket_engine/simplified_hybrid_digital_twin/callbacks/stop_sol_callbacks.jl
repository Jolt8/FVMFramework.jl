function min_chamber_pressure_limit_expanded(u, t, integrator, system, geo)
    u_named = ComponentVector(u, system.state_axes)
    p_named = ComponentVector(integrator.p, system.p_axes)

    #using the pre_chamber seems like the right choice because we want to stop when all the internals of the motor 
    #go below ambient pressure
    cell_id = get_cellset(geo, "pre_chamber")[1]
    #cell_id = get_cellset(geo, "fluid_domain")[1]

    #NOTE: if this callback is used when the density that's measured is already below the threshold, the solver will not stop

    if t < 2.0 #we make sure that all rockets can burn for at least 2 seconds before applying this
        return 1.0
    else
        # system.p_axes stores every parameter as a 1-element vector, hence the [1]
        return 2.0 - u_named.density[cell_id] #we'll just terminate when the density goes below a certain value :)
    end
end

min_chamber_pressure_limit = (u, t, integrator) -> min_chamber_pressure_limit_expanded(u, t, integrator, system, geo)

min_chamber_pressure_termination_cb = ContinuousCallback(
    min_chamber_pressure_limit,
    terminate!,
    nothing;
    save_positions = (true, false),
)