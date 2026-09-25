
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

    @show "before ode solve"
    #try
    sol = solve(prob,
        saveat = (p.simulation_time / 200),
        callback = optimized_cb_set,
        isoutofdomain = isoutofdomain_feedsystem,
        maxiters = 1000,
        dtmin = 1e-5
    )
    #catch e
        #println(e)
        #return theta, p, 1e10
    #end
    @show "after ode solve"

    @show sol.retcode

    if !(sol.retcode == SciMLBase.ReturnCode.Success || sol.retcode == SciMLBase.ReturnCode.Terminated)
        #@show p.final_fuel_grain_void_diameter
        u_named = ComponentVector(sol.u[end], u_axes)
        @show u_named.port_diameter
        return theta, p, (1e10 + p.final_fuel_grain_void_diameter - u_named.port_diameter) * 1e8 #we want the solver to prioritize runs that got closer to expending all the fuel even if they failed
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

    for i in eachindex(sol.u)
        curr_t = sol.t[i]
        u_named = ComponentVector(sol.u[i], u_axes)

        update_state!(du_temporary, u_named, p, curr_t, oxidizer_model, chamber_model)

        if (u_named.port_diameter >= p.final_fuel_grain_void_diameter) && found_depletion_time == false
            depletion_time = curr_t
            found_depletion_time = true
            unutilized_oxidizer_loss = 0.001 * abs2(u_named.tank_oxidizer_mass)
            break #stop evaluating loss if the fuel has burned out
        end

        if (u_named.tank_oxidizer_mass <= 1e-6) && found_depletion_time == false
            depletion_time = curr_t
            found_depletion_time = true
            unburned_fuel_loss = 0.001 * abs2(p.fuel_mass)
            break #stop evaluating loss if the oxidizer has burned out
        end

        adjustable_valve_pressure_drop = adjustable_valve_flow!(du_temporary, u_named, p, curr_t)

        injector_valve_pressure_drop, oxidizer_mass_flow = injector_valve_flow!(du_temporary, u_named, p, curr_t)

        #@show p.target_injector_velocity
        #@show oxidizer_mass_flow / (p.mid_section_density * p.injector_orifice_area)
        if show_loss_function_compositions
            @info "injector_velocity_loss"
            @info p.target_injector_velocity
            @info oxidizer_mass_flow / (p.mid_section_density * p.injector_orifice_area)
        end
        injector_velocity_loss += 0.00001 * (1 / length(sol.u)) * abs2(p.target_injector_velocity - oxidizer_mass_flow / (p.mid_section_density * p.injector_orifice_area))
        
        if show_loss_function_compositions
            @info "pressure_drop_ratio_loss"
            @info p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio
            @info injector_valve_pressure_drop / adjustable_valve_pressure_drop
        end
        pressure_drop_ratio_loss += 0.00001 * (1 / length(sol.u)) * abs2(p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio - (injector_valve_pressure_drop / adjustable_valve_pressure_drop))

        #OBSERVATION: it seems like we're going to have to weigh the importance of injector velocity against adjustable valve authority

        fuel_mass_flow = regression_rate!(du_temporary, u_named, p, curr_t, oxidizer_mass_flow)

        #=
        chamber_gas_mass_flow_out = nozzle_outlet!(du_temporary, u_named, p, curr_t)

        update_mass_fractions!(du_temporary, u_named, p, curr_t, oxidizer_mass_flow, fuel_mass_flow)

        combustion_zone_energy_conservation!(du_temporary, u_named, p, curr_t, chamber_model, oxidizer_mass_flow, fuel_mass_flow, chamber_gas_mass_flow_out)
        =#

        #if we want, we can get the derivative of the cummulative_impulse over time to plot the thrust profile of the engine!
        if i > 1
            u_named_prev = ComponentVector(sol.u[i-1], u_axes)
            dt = sol.t[i] - sol.t[i-1]

            #Cummulative impulse loss calcs
            oxidizer_used = -(u_named.tank_oxidizer_mass - u_named_prev.tank_oxidizer_mass)
            fuel_used = -(p.fuel_mass - last_fuel_mass)
            last_fuel_mass = p.fuel_mass

            propellant_used = oxidizer_used + fuel_used

            cummulative_impulse += p.propellant_isp * p.gravity * propellant_used

            #Oxidizer to fuel ratio loss
            p.oxidizer_to_fuel_ratio += (1 / length(sol.u)) * (oxidizer_used / fuel_used)
            oxidizer_to_fuel_ratio_loss += (1 / length(sol.u)) * 0.001 * abs2(p.desired_oxidizer_to_fuel_ratio - (oxidizer_used / fuel_used))

            #Thrust loss calcs
            oxidizer_mass_flow = oxidizer_used / dt
            fuel_mass_flow = fuel_used / dt

            propallant_mass_flow = oxidizer_mass_flow + fuel_mass_flow

            thrust_produced = p.propellant_isp * p.gravity * propallant_mass_flow
            
            p.average_thrust += (1 / (length(sol.t) - 1)) * thrust_produced
            thrust_loss += (1 / (length(sol.t) - 1)) * 0.001 * abs2(p.desired_average_thrust - thrust_produced)
            if show_loss_function_compositions
                @info "thrust_loss"
                @info p.desired_average_thrust
                @info thrust_produced
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

    p.cummulative_impulse = cummulative_impulse
    impulse_loss = 0.0001 * abs2(p.desired_impulse - cummulative_impulse)
    if show_loss_function_compositions
        @info "impulse_loss"
        @info p.desired_impulse
        @info cummulative_impulse
    end

    p.burn_time = depletion_time
    burn_time_loss = 0.01 * abs2(p.desired_burn_time - depletion_time)
    if show_loss_function_compositions
        @info "burn_time_loss"
        @info p.desired_burn_time
        @info depletion_time
    end

    #above_max_fuel_grain_diameter_loss = 0.01 * abs2(p.u0_fuel_grain_void_diameter + p.additional_fuel_grain_void_diameter - p.fuel_grain_max_diameter)
    #we'll just enforce a bound on final_fuel_grain_void_diameter

    all_losses = ComponentVector(
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

    @show sum(all_losses)

    return theta, p, all_losses
end