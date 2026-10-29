if do_plots
    using Plots

    plot(reduced_sol.t, [u_timestep.tank.mass for u_timestep in u_named], label = "Tank itrous Oxide Mass")
    plot!(reduced_sol.t, [u_timestep.tank.temperature for u_timestep in u_named], label = "Tank Temperature")

    plot(reduced_sol.t, [u_timestep.igniter_mass[1] for u_timestep in u_named], label = "Igniter Mass")

    #=
    n_regression_points_observed = 10

    for i in 1:(length(u_named[1].fuel_grain_thickness))
    =#
end