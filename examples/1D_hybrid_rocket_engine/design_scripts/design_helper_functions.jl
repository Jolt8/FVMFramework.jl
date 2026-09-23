function append_optimized_parameters!(theta, u, p, theta_to_u_map, theta_to_p_map, p_to_u_map)
    for (u_idx, theta_idx) in theta_to_u_map
        u[u_idx] = theta[theta_idx]
    end

    for (p_idx, theta_idx) in theta_to_p_map
        p[p_idx] = theta[theta_idx]
    end

    for (u_idx, p_idx) in p_to_u_map
        u[u_idx] = p[p_idx]
    end
end

function create_theta_guess(optimized_properties, u, properties)
    theta_guess_dict = OrderedDict{Symbol, Number}()
    theta_lb_dict = OrderedDict{Symbol, Number}()
    theta_ub_dict = OrderedDict{Symbol, Number}()

    theta_to_u_map = Tuple{Int, Int}[]
    theta_to_p_map = Tuple{Int, Int}[]

    for optimized_property in optimized_properties
        if optimized_property.name in propertynames(u)
            theta_guess_dict[optimized_property.name] = ustrip(upreferred(getproperty(u, optimized_property.name)))
            theta_lb_dict[optimized_property.name] = ustrip(upreferred(optimized_property.lb))
            theta_ub_dict[optimized_property.name] = ustrip(upreferred(optimized_property.ub))
        end

        if optimized_property.name in propertynames(properties)
            theta_guess_dict[optimized_property.name] = ustrip(upreferred(getproperty(properties, optimized_property.name)))
            theta_lb_dict[optimized_property.name] = ustrip(upreferred(optimized_property.lb))
            theta_ub_dict[optimized_property.name] = ustrip(upreferred(optimized_property.ub))
        end
    end

    theta_guess = ComponentVector(theta_guess_dict)
    theta_lb = ComponentVector(theta_lb_dict)
    theta_ub = ComponentVector(theta_ub_dict)

    for optimized_property in optimized_properties
        if optimized_property.name in propertynames(u)
            theta_guess_dict[optimized_property.name] = ustrip(upreferred(getproperty(u, optimized_property.name)))
            theta_lb_dict[optimized_property.name] = ustrip(upreferred(optimized_property.lb))
            theta_ub_dict[optimized_property.name] = ustrip(upreferred(optimized_property.ub))

            #Push these into maps that will append the indexes of the guessed parameters to u
            push!(theta_to_u_map, (label2index(u, string(optimized_property.name))[1], label2index(theta_guess, string(optimized_property.name))[1]))
            #remember, this logs (updated_idx, pulled_from_idx)
        elseif optimized_property.name in propertynames(properties)
            theta_guess_dict[optimized_property.name] = ustrip(upreferred(getproperty(properties, optimized_property.name)))
            theta_lb_dict[optimized_property.name] = ustrip(upreferred(optimized_property.lb))
            theta_ub_dict[optimized_property.name] = ustrip(upreferred(optimized_property.ub))

            #Push these into maps that will append the indexes of the guessed parameters to p
            push!(theta_to_p_map, (label2index(properties, string(optimized_property.name))[1], label2index(theta_guess, string(optimized_property.name))[1]))   
            #remember, this logs (updated_idx, pulled_from_idx)
        else
            @warn "Property $optimized_property not found in properties or u"
        end
    end

    p_to_u_map = Tuple{Int, Int}[]

    for name in propertynames(u)
        new_name = Symbol("u0_" * string(name))
        if new_name in propertynames(properties)
            push!(p_to_u_map, (label2index(u, string(name))[1], label2index(properties, string(new_name))[1]))
            #remember, this logs (updated_idx, pulled_from_idx)
        end
    end

    theta_axes = getaxes(theta_guess)
    u_axes = getaxes(u)
    p_axes = getaxes(properties)

    return (
        theta_guess, theta_lb, theta_ub, 
        theta_axes, u_axes, p_axes,
        theta_to_u_map, theta_to_p_map, p_to_u_map
    )
end

function build_implicit_prob(f_closure, du0_vec, u0_vec, thermocouple_data, p_guess)
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
    )

    ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

    t0 = 0.0
    tMax = ustrip(upreferred(thermocouple_data.timestamps[end]))
    tspan = (t0, tMax)

    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

    return implicit_prob
end

function uv_flash_via_vt(
    model,
    target_internal_energy,
    volume,
    moles;
    temperature_bounds = (180.0, 450.0),
)
    energy_residual(T) = begin
        result = vt_flash(model, volume, T, moles)
        internal_energy(model, result) - target_internal_energy
    end

    lower_temperature, upper_temperature = temperature_bounds

    lower_residual = energy_residual(lower_temperature)
    upper_residual = energy_residual(upper_temperature)

    lower_residual * upper_residual <= 0 ||
        error("Temperature bounds do not bracket the requested internal energy")

    temperature = find_zero(
        energy_residual,
        temperature_bounds,
        Roots.Brent(),
    )

    #=
    energy_residual = (temperature, _) -> Clapeyron.VT0.internal_energy(model, volume, temperature, moles) - target_internal_energy

    temperature_problem = NonlinearProblem(energy_residual, temperature_bounds[1])
    temperature_solution = solve(temperature_problem, NewtonRaphson(); abstol = 1e-8, reltol = 1e-8)
    =#

    return vt_flash(model, volume, temperature, moles)
end