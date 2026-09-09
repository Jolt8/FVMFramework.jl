function regenerate_fvm_state(sol, system, solve_system!, geo, p_guess; u_additional_information = ComponentVector(), track_progress = false)
    #we don't need an du_additional_information because we're not updating any new fields each time and because we can just put derivatives in u_additional_information
    du_list = ComponentVector[]
    u_list = ComponentVector[]

    u_named = [ComponentVector(sol.u[i], system.state_axes) for i in eachindex(sol.u)]
    du_named = [ComponentVector(deepcopy(sol.u[i]), system.state_axes) for i in eachindex(sol.u)]

    temporary_cache = ComponentVector(system.cache_vec, system.cache_axes)
    temporary_cache .= 0.0

    t_last = sol.t[1]

    for (i, t) in enumerate(sol.t)
        u = merge_properties(u_named[i], deepcopy(ComponentVector(system.properties_vec, system.properties_axes)))
        u = merge_properties(u, deepcopy(ComponentVector(system.cache_vec, system.cache_axes)))
        u = merge_properties(u, u_additional_information)

        du = merge_properties(du_named[i], deepcopy(ComponentVector(system.cache_vec, system.cache_axes)))
        du .= 0.0

        solve_system!(du, u, p_guess, t, geo, system) #solve_sytem! is extremely cheap to run once we've already solved it

        #however, since it updates u_named, we need to set it back to the original
        #if this breaks in the future, we can just fallback to only running update_fluid_properties! each iteration instead of solve_system!
        #if we don't run the entire solve_system! stuff like u.mass_face or u.heat will not get updated
        for propertyname in propertynames(u_named[i])
            getproperty(u, propertyname)[1:end] = getproperty(u_named[i], propertyname)[1:end] #we have to use getproperty() because just doing u[propertyname] doesn't update it
        end

        #TODO: we may also want to log the experimental transducer data into the vtk file so that we can observe it
        #this would be very useful for stuff like cell_velocity_uncertainties
        push!(du_list, du)
        push!(u_list, u)

        #TODO: get the caches finite differences to work with nested fields
        if i == 1
            #for the first time step we can't use finite differences so we will just set it to 0.0
            for propertyname in propertynames(temporary_cache) #we use finite differences for apprroximating the time derivative of the cached variables who derivatives are not cached
                if iszero(du_list[i][propertyname]) #we don't want to overwrite cached variables that already have a value
                    getproperty(du_list[i], propertyname)[1:end] .= 0.0
                end
            end
        else
            for propertyname in propertynames(temporary_cache) #we use finite differences for apprroximating the time derivative of the cached variables who derivatives are not cached
                if iszero(du_list[i][propertyname])  #we don't want to overwrite cached variables that already have a value
                    for cell_id in eachindex(temporary_cache[propertyname])
                        getproperty(du_list[i], propertyname)[cell_id] = (u_list[i][propertyname][cell_id] - u_list[i-1][propertyname][cell_id]) / (t - t_last)
                    end
                end
            end
        end

        t_last = t

        if track_progress == true
            println("regenerating timestep $t of $(sol.t[end]), remaining steps = $(length(sol.t)-i)")
        end
    end

    return du_list, u_list
end