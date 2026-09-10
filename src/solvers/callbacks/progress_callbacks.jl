function every_step(u, t, integrator) # Event when condition(u,t,integrator) <= 0
    return true
end

function affect_simple!(integrator)
    println("sim time: $(integrator.t)s")
end

show_t_progress = DiscreteCallback(every_step, affect_simple!, save_positions=(false, false))

global last_print_time_FVMFramework_approx_time_to_finish_cb = time() #this is bad, we shouldn't do this

#we have to add an average of the previous n real times left 
function approx_time_left_affect!(integrator)
    curr_time = time()

    dt_sim_time = integrator.t - integrator.tprev

    sim_tspan = first(integrator.opts.tstops) #tstops is a bitset

    approximate_steps_left = (sim_tspan - integrator.t) / dt_sim_time

    dt_real_time = curr_time - last_print_time_FVMFramework_approx_time_to_finish_cb

    approximate_real_time_left = approximate_steps_left * dt_real_time

    println("sim time: $(integrator.t)s")
    println("approximate time until finished: $(approximate_real_time_left)s")
    #we could also convert seconds to minutes, hours, or days depending on how large the time span is

    #furthermore, I don't know how this would be achieved, but it would be interesting to see an actively updated time since last step value
    #sometimes the simulation hangs on one step and after I check back, I don't know if that's expected
    global last_print_time_FVMFramework_approx_time_to_finish_cb = time() #really long variable name so this never gets interacted with 
end

#I wonder if Extrapolations.jl would be cheap enough to run here

#for some reason this affect causes the simulation to end early

approx_time_left_affect_closure! = (integrator) -> approx_time_left_affect!(integrator) #, last_print_time)
#we could maybe use a closure to avoid the global

#I mean, we could make save_positions = (true, false) to log an additional capture of sol.u for almost no additional cost
#even though we could only call the callback at saveat t values, having the callback trigger at very small timespans in the beginning is useful


"""
    approximate_time_to_finish_cb(x)

This is a discrete callback that can be used to approximate the time unti a simulation is finished

Note that this uses a global variable (last_print_time_FVMFramework_approx_time_to_finish_cb)
"""
function approximate_time_to_finish_cb(x) end
approximate_time_to_finish_cb = DiscreteCallback(every_step, approx_time_left_affect_closure!, save_positions=(false, false))

#Update: this is so much more useful than I ever would've guessed, getting an approximate time until finished is amazing


#another thing that would be nice to implement is a method to send my simulation to my computer at home and then run it
#another thing that would be nice is a method to save the state of my simulation and then resume it later. 
#   - I guess we could log the current state of all u values and then use that. 