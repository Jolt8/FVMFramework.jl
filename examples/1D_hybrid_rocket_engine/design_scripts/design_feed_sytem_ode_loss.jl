
function trainsient_system_design_loss(theta, u0, p, theta_axes, u_axes, p_axes, oxidizer_model, chamber_model, theta_to_u_map, theta_to_p_map, p_to_u_map, append_optimized_parameters!, update_u0!, isoutofdomain_feedsystem, optimized_cb_set, problem_template)
    theta = ComponentVector(theta, theta_axes)
    
    # Promote p to the type of theta (which will be Dual during ForwardDiff) so it can accept Duals
    u0 = eltype(theta).(u0)
    p = ComponentVector(eltype(theta).(p), p_axes)

    append_optimized_parameters!(theta, u0, p, theta_to_u_map, theta_to_p_map, p_to_u_map)
    update_u0!(u0, p, 0.0, oxidizer_model, chamber_model)

    #=
    prob = get!(task_local_storage(), :implicit_prob) do    
        # Build a new closure bound to the thread-isolated copies
        f_closure = (du, u, p, t) -> system_ode!(du, u, p, t, oxidizer_model, chamber_model, p_axes)

        ode_func = ODEFunction(f_closure)#, jac_prototype = float.(jac_sparsity))

        tspan = (0.0, p.simulation_time)

        prob = ODEProblem(ode_func, u0, tspan, p)
    end
    =#
    prob = remake(problem_template; u0 = Vector(u0), p = Vector(p), tspan = (0.0, p.simulation_time))

    sol = 0
    
    try
    sol = solve(prob,
        saveat = (p.simulation_time / 200),
        callback = optimized_cb_set,
        isoutofdomain = isoutofdomain_feedsystem,
        maxiters = 1000,
        #dtmin = 1e-5
    )
    catch e
        println(e)
        return theta, p, 1e10
    end

    did_not_finish_loss = 0.0

    if !(sol.retcode == SciMLBase.ReturnCode.Success || sol.retcode == SciMLBase.ReturnCode.Terminated)
        #@show p.final_fuel_grain_void_diameter
        u_named = ComponentVector(sol.u[end], u_axes)
        did_not_finish_loss += 1e10 + 1e8 * p.final_fuel_grain_void_diameter - u_named.port_diameter
    end

    #Losses updated every iteration
    injector_velocity_loss = 0.0
    pressure_drop_ratio_loss = 0.0
    oxidizer_to_fuel_ratio_loss = 0.0
    thrust_loss = 0.0
    
    #Loss updated once
    unburned_fuel_loss = 0.0
    unutilized_oxidizer_loss = 0.0

    #Performance Metrics
    cummulative_impulse = 0.0
    depletion_time = 0.0

    show_loss_function_compositions = false

    du_temporary = ComponentVector(similar(sol.u[1]), u_axes)

    found_depletion_time = false

    last_fuel_mass = p.fuel_mass

    thrust_over_time = []
    pressures_over_time = []
    tank_pressure_over_time = []
    tank_oxidizer_mass_over_time = []
    mid_section_pressure_over_time = []
    isp_over_time = []
    injector_velocity_over_time = []
    pressure_drop_ratio_over_time = []
    port_diameter_over_time = []

    for i in eachindex(sol.u)
        curr_t = sol.t[i]
        u_named = ComponentVector(sol.u[i], u_axes)

        update_state!(du_temporary, u_named, p, curr_t, oxidizer_model, chamber_model)

        push!(pressures_over_time, p.chamber_pressure)
        push!(tank_pressure_over_time, p.tank_pressure)
        push!(mid_section_pressure_over_time, p.mid_section_pressure)
        push!(port_diameter_over_time, u_named.port_diameter)
        push!(tank_oxidizer_mass_over_time, u_named.tank_oxidizer_mass)

        if valve_is_closed(u_named, p, curr_t) == true && found_depletion_time == false
            depletion_time = curr_t
            found_depletion_time = true
            #unutilized_oxidizer_loss = 0.001 * abs2(u_named.tank_oxidizer_mass)
            #unburned_fuel_loss = 0.001 * abs2(p.fuel_mass)
            #break #stop evaluating loss if the fuel has burned out
        end

        adjustable_valve_pressure_drop = adjustable_valve_flow!(du_temporary, u_named, p, curr_t)

        injector_valve_pressure_drop, oxidizer_mass_flow = injector_valve_flow!(du_temporary, u_named, p, curr_t)

        #@show p.target_injector_velocity
        #@show oxidizer_mass_flow / (p.mid_section_density * p.injector_orifice_area)
        if show_loss_function_compositions
            @show "injector_velocity_loss"
            @show p.target_injector_velocity
            @show oxidizer_mass_flow / (p.mid_section_density * p.injector_orifice_area)
        end
        injector_velocity = oxidizer_mass_flow / (p.mid_section_density * p.injector_orifice_area)
        push!(injector_velocity_over_time, injector_velocity)
        injector_velocity_loss += 0.00001 * (1 / length(sol.u)) * abs2(p.target_injector_velocity - injector_velocity)
        
        if show_loss_function_compositions
            @show "pressure_drop_ratio_loss"
            @show p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio
            @show injector_valve_pressure_drop / adjustable_valve_pressure_drop
        end

        if adjustable_valve_pressure_drop > 1e-9 #sometimes the very last step of the solve can result in the adjustable valve pressure drop being zero
            pressure_drop_ratio = injector_valve_pressure_drop / adjustable_valve_pressure_drop
            push!(pressure_drop_ratio_over_time, pressure_drop_ratio)
            pressure_drop_ratio_loss += 0.00001 * (1 / length(sol.u)) * abs2(p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio - pressure_drop_ratio)
        else
            push!(pressure_drop_ratio_over_time, 0)
        end

        #OBSERVATION: it seems like we're going to have to weigh the importance of injector velocity against adjustable valve authority
        fuel_mass_flow = 0.0
        if found_depletion_time == false
            fuel_mass_flow = regression_rate!(du_temporary, u_named, p, curr_t, oxidizer_mass_flow)
        end

        #=
        update_mass_fractions!(du_temporary, u_named, p, curr_t, oxidizer_mass_flow, fuel_mass_flow)

        combustion_zone_energy_conservation!(du_temporary, u_named, p, curr_t, chamber_model, oxidizer_mass_flow, fuel_mass_flow, chamber_gas_mass_flow_out)
        =#

        #if we want, we can get the derivative of the cummulative_impulse over time to plot the thrust profile of the engine!
        if i > 1
            u_named_prev = ComponentVector(sol.u[i-1], u_axes)
            dt = sol.t[i] - sol.t[i-1]

            #Cummulative impulse loss calcs
            oxidizer_used = u_named_prev.tank_oxidizer_mass - u_named.tank_oxidizer_mass
            fuel_used = p.fuel_density * (π / 4) * (u_named.port_diameter^2 - u_named_prev.port_diameter^2) * p.fuel_grain_length
            
            oxidizer_mass_flow = oxidizer_used / dt
            fuel_mass_flow = fuel_used / dt

            propellant_used = oxidizer_used + fuel_used
            propellant_mass_flow = oxidizer_mass_flow + fuel_mass_flow

            oxidizer_to_fuel_ratio = oxidizer_mass_flow / max(fuel_mass_flow, 1e-9)

            p.propellant_isp = isp_interpolator_Pa(p.chamber_pressure, oxidizer_to_fuel_ratio, p.expansion_ratio)

            push!(isp_over_time, p.propellant_isp)

            p.propellant_characteristic_velocity = cstar_interpolator_Pa(p.chamber_pressure, oxidizer_to_fuel_ratio, p.expansion_ratio)

            chamber_gas_mass_flow_out = p.nozzle_discharge_coefficient * ((p.chamber_pressure * p.nozzle_throat_area) / p.propellant_characteristic_velocity)

            thrust_produced = p.propellant_isp * p.gravity * chamber_gas_mass_flow_out

            push!(thrust_over_time, thrust_produced)

            cummulative_impulse += thrust_produced * dt

            # O/F ratio is now handled by integrating total mass used at the end of the simulation
            
            p.average_thrust += (1 / (length(sol.t) - 1)) * thrust_produced
            thrust_loss += (1 / (length(sol.t) - 1)) * 0.001 * abs2(p.desired_average_thrust - thrust_produced)
            if show_loss_function_compositions
                @show "thrust_loss"
                @show p.desired_average_thrust
                @show thrust_produced
            end
        end

        if i == length(sol.u)
            if depletion_time == 0.0 && found_depletion_time == false
                @warn "the simulation did no run long enough to completely burn out the fuel grain"
                @show remaining_fuel_grain = (u_named.port_diameter - p.final_fuel_grain_void_diameter)

                @show u_named.port_diameter
                @show p.final_fuel_grain_void_diameter
            end
        end
    end

    if true == true
        injector_velocity_plt = plot(sol.t, injector_velocity_over_time, label = "Injector Velocity", xlabel = "Time [s]", ylabel = "Injector Velocity [m/s]")
        display(injector_velocity_plt)

        pressure_drop_ratio_plt = plot(sol.t, pressure_drop_ratio_over_time, label = "Injector to Adjustable Valve Pressure Drop Ratio", xlabel = "Time [s]", ylabel = "Pressure Drop Ratio")
        display(pressure_drop_ratio_plt)

        isp_plt = plot(sol.t, [0.0, isp_over_time...], label = "ISP", xlabel = "Time [s]", ylabel = "ISP [s]")
        display(isp_plt)

        @show maximum(isp_over_time)

        tank_mass_plt = plot(sol.t, tank_oxidizer_mass_over_time, label = "Tank Oxidizer Mass", xlabel = "Time [s]", ylabel = "Oxidizer Mass [kg]")
        display(tank_mass_plt)

        port_diameter_plot = plot(sol.t, port_diameter_over_time, label = "Port Diameter", xlabel = "Time [s]", ylabel = "Port Diameter [m]")
        display(port_diameter_plot)

        pressure_plot = plot(sol.t, tank_pressure_over_time, label = "Tank Pressure", xlabel = "Time [s]", ylabel = "Pressure [Pa]")
        plot!(pressure_plot, sol.t, mid_section_pressure_over_time, label = "Mid Section Pressure")
        plot!(pressure_plot, sol.t, pressures_over_time, label = "Chamber Pressure")
        display(pressure_plot)

        thrust_plt = plot(sol.t, [0.0, thrust_over_time...], label = "Thrust", xlabel = "Time [s]", ylabel = "Thrust [N]")
        display(thrust_plt)

        @show ComponentVector(sol.u[end], u_axes).tank_oxidizer_mass
    end

    p.cummulative_impulse = cummulative_impulse
    impulse_loss = 0.0000001 * abs2(p.desired_impulse - cummulative_impulse)
    if show_loss_function_compositions
        @show "impulse_loss"
        @show p.desired_impulse
        @show cummulative_impulse
    end

    p.burn_time = depletion_time
    burn_time_loss = 0.01 * abs2(p.desired_burn_time - depletion_time)
    if show_loss_function_compositions
        @show "burn_time_loss"
        @show p.desired_burn_time
        @show depletion_time
    end

    u_start = ComponentVector(sol.u[1], u_axes)
    u_end = ComponentVector(sol.u[end], u_axes)
    
    total_fuel_burned = p.fuel_density * (pi * (u_end.port_diameter / 2)^2 - pi * (u_start.port_diameter / 2)^2) * p.fuel_grain_length
    #@show total_fuel_burned

    total_oxidizer_used = u_start.tank_oxidizer_mass - u_end.tank_oxidizer_mass
    #@show total_oxidizer_used
    
    p.oxidizer_to_fuel_ratio = total_oxidizer_used / max(total_fuel_burned, 1e-9)
    oxidizer_to_fuel_ratio_loss = 0.1 * abs2(p.desired_oxidizer_to_fuel_ratio - p.oxidizer_to_fuel_ratio)

    #above_max_fuel_grain_diameter_loss = 0.01 * abs2(p.u0_fuel_grain_void_diameter + p.additional_fuel_grain_void_diameter - p.fuel_grain_max_diameter)
    #we'll just enforce a bound on final_fuel_grain_void_diameter

    all_losses = ComponentVector(
        #Misc Losses
        did_not_finish_loss = did_not_finish_loss,

        #Updated every solver iteration
        injector_velocity_loss = injector_velocity_loss,
        pressure_drop_ratio_loss = pressure_drop_ratio_loss,
        oxidizer_to_fuel_ratio_loss = oxidizer_to_fuel_ratio_loss,
        #thrust_loss = thrust_loss, #this seems like a bad metric, it allows short burn times with better average thrusts

        #Updated once at the end of the simulation
        impulse_loss = impulse_loss,
        burn_time_loss = burn_time_loss,
        unburned_fuel_loss = unburned_fuel_loss,
        unutilized_oxidizer_loss = unutilized_oxidizer_loss
    )

    #@show sum(all_losses)

    return theta, p, all_losses
end
