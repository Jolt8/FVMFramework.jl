function _riemann_pressure_function(pressure, density, reference_pressure, sound_speed, gamma)
    if pressure > reference_pressure
        coefficient_a = 2.0 / ((gamma + 1.0) * density)
        coefficient_b = (gamma - 1.0) * reference_pressure / (gamma + 1.0)
        square_root = sqrt(coefficient_a / (pressure + coefficient_b))
        value = (pressure - reference_pressure) * square_root
        derivative = square_root * (
            1.0 - 0.5 * (pressure - reference_pressure) / (pressure + coefficient_b)
        )
    else
        pressure_ratio = pressure / reference_pressure
        exponent = (gamma - 1.0) / (2.0 * gamma)
        value = 2.0 * sound_speed / (gamma - 1.0) * (pressure_ratio^exponent - 1.0)
        derivative = pressure_ratio^(-(gamma + 1.0) / (2.0 * gamma)) / (
            density * sound_speed
        )
    end
    return value, derivative
end

function _sod_star_state(gamma)
    density_left = 1.0
    velocity_left = 0.0
    pressure_left = 1.0
    density_right = 0.125
    velocity_right = 0.0
    pressure_right = 0.1
    sound_speed_left = sqrt(gamma * pressure_left / density_left)
    sound_speed_right = sqrt(gamma * pressure_right / density_right)

    pressure = max(
        1e-12,
        0.5 * (pressure_left + pressure_right) -
        0.125 * (velocity_right - velocity_left) *
        (density_left + density_right) * (sound_speed_left + sound_speed_right),
    )
    for iteration in 1:50
        function_left, derivative_left = _riemann_pressure_function(
            pressure,
            density_left,
            pressure_left,
            sound_speed_left,
            gamma,
        )
        function_right, derivative_right = _riemann_pressure_function(
            pressure,
            density_right,
            pressure_right,
            sound_speed_right,
            gamma,
        )
        updated_pressure = pressure - (
            function_left + function_right + velocity_right - velocity_left
        ) / (derivative_left + derivative_right)
        updated_pressure = max(updated_pressure, 1e-12)
        if abs(updated_pressure - pressure) <= 1e-13 * max(updated_pressure, pressure, 1.0)
            pressure = updated_pressure
            break
        end
        pressure = updated_pressure
    end

    function_left, derivative_left = _riemann_pressure_function(
        pressure,
        density_left,
        pressure_left,
        sound_speed_left,
        gamma,
    )
    function_right, derivative_right = _riemann_pressure_function(
        pressure,
        density_right,
        pressure_right,
        sound_speed_right,
        gamma,
    )
    velocity = 0.5 * (
        velocity_left + velocity_right + function_right - function_left
    )
    return pressure, velocity
end

function sod_exact_primitive(x, time; gamma = 1.4, discontinuity = 0.5)
    if time <= 0.0
        if x < discontinuity
            return (density = 1.0, velocity = 0.0, pressure = 1.0)
        end
        return (density = 0.125, velocity = 0.0, pressure = 0.1)
    end

    density_left = 1.0
    velocity_left = 0.0
    pressure_left = 1.0
    density_right = 0.125
    velocity_right = 0.0
    pressure_right = 0.1
    sound_speed_left = sqrt(gamma * pressure_left / density_left)
    sound_speed_right = sqrt(gamma * pressure_right / density_right)
    pressure_star, velocity_star = _sod_star_state(gamma)
    similarity_coordinate = (x - discontinuity) / time

    if similarity_coordinate <= velocity_star
        if pressure_star > pressure_left
            shock_speed = velocity_left - sound_speed_left * sqrt(
                (gamma + 1.0) * pressure_star / (2.0 * gamma * pressure_left) +
                (gamma - 1.0) / (2.0 * gamma)
            )
            if similarity_coordinate <= shock_speed
                return (density = density_left, velocity = velocity_left, pressure = pressure_left)
            end
            pressure_ratio = pressure_star / pressure_left
            density_star = density_left * (
                pressure_ratio + (gamma - 1.0) / (gamma + 1.0)
            ) / (
                (gamma - 1.0) * pressure_ratio / (gamma + 1.0) + 1.0
            )
            return (density = density_star, velocity = velocity_star, pressure = pressure_star)
        end

        sound_speed_star = sound_speed_left * (
            pressure_star / pressure_left
        )^((gamma - 1.0) / (2.0 * gamma))
        head_speed = velocity_left - sound_speed_left
        tail_speed = velocity_star - sound_speed_star
        if similarity_coordinate <= head_speed
            return (density = density_left, velocity = velocity_left, pressure = pressure_left)
        elseif similarity_coordinate >= tail_speed
            density_star = density_left * (pressure_star / pressure_left)^(1.0 / gamma)
            return (density = density_star, velocity = velocity_star, pressure = pressure_star)
        end

        velocity = 2.0 / (gamma + 1.0) * (
            sound_speed_left + 0.5 * (gamma - 1.0) * velocity_left + similarity_coordinate
        )
        sound_speed = 2.0 / (gamma + 1.0) * (
            sound_speed_left + 0.5 * (gamma - 1.0) * (velocity_left - similarity_coordinate)
        )
        density = density_left * (sound_speed / sound_speed_left)^(2.0 / (gamma - 1.0))
        pressure = pressure_left * (sound_speed / sound_speed_left)^(2.0 * gamma / (gamma - 1.0))
        return (density = density, velocity = velocity, pressure = pressure)
    end

    if pressure_star > pressure_right
        shock_speed = velocity_right + sound_speed_right * sqrt(
            (gamma + 1.0) * pressure_star / (2.0 * gamma * pressure_right) +
            (gamma - 1.0) / (2.0 * gamma)
        )
        if similarity_coordinate >= shock_speed
            return (density = density_right, velocity = velocity_right, pressure = pressure_right)
        end
        pressure_ratio = pressure_star / pressure_right
        density_star = density_right * (
            pressure_ratio + (gamma - 1.0) / (gamma + 1.0)
        ) / (
            (gamma - 1.0) * pressure_ratio / (gamma + 1.0) + 1.0
        )
        return (density = density_star, velocity = velocity_star, pressure = pressure_star)
    end

    sound_speed_star = sound_speed_right * (
        pressure_star / pressure_right
    )^((gamma - 1.0) / (2.0 * gamma))
    head_speed = velocity_right + sound_speed_right
    tail_speed = velocity_star + sound_speed_star
    if similarity_coordinate >= head_speed
        return (density = density_right, velocity = velocity_right, pressure = pressure_right)
    elseif similarity_coordinate <= tail_speed
        density_star = density_right * (pressure_star / pressure_right)^(1.0 / gamma)
        return (density = density_star, velocity = velocity_star, pressure = pressure_star)
    end

    velocity = 2.0 / (gamma + 1.0) * (
        -sound_speed_right + 0.5 * (gamma - 1.0) * velocity_right + similarity_coordinate
    )
    sound_speed = 2.0 / (gamma + 1.0) * (
        sound_speed_right - 0.5 * (gamma - 1.0) * (velocity_right - similarity_coordinate)
    )
    density = density_right * (sound_speed / sound_speed_right)^(2.0 / (gamma - 1.0))
    pressure = pressure_right * (sound_speed / sound_speed_right)^(2.0 * gamma / (gamma - 1.0))
    return (density = density, velocity = velocity, pressure = pressure)
end

function check_sod_benchmark(case, solution, final_time)
    final_state = ComponentVector(solution.u[end], case.system.state_axes)
    density_errors = Float64[]
    pressure_errors = Float64[]
    exact_density = Float64[]
    exact_pressure = Float64[]

    for cell_id in eachindex(case.geo.cell_centroids)
        x = case.geo.cell_centroids[cell_id][1]
        exact = sod_exact_primitive(x, final_time; gamma = case.gamma)
        density = final_state.density[cell_id]
        kinetic_energy_density = 0.5 * (
            final_state.momentum_density_u[cell_id]^2 +
            final_state.momentum_density_v[cell_id]^2 +
            final_state.momentum_density_w[cell_id]^2
        ) / density
        pressure = (case.gamma - 1.0) * (
            final_state.volumetric_energy[cell_id] - kinetic_energy_density
        )
        push!(density_errors, abs(density - exact.density))
        push!(pressure_errors, abs(pressure - exact.pressure))
        push!(exact_density, exact.density)
        push!(exact_pressure, exact.pressure)
    end

    density_l1_error = sum(case.geo.cell_volumes .* density_errors) / sum(case.geo.cell_volumes)
    pressure_l1_error = sum(case.geo.cell_volumes .* pressure_errors) / sum(case.geo.cell_volumes)
    # A first-order 32-cell scheme necessarily smears the shock and contact.
    # These bounds are more than twice the measured coarse-grid error while
    # still rejecting a stationary or seriously misplaced discontinuity.
    density_tolerance = 0.08
    pressure_tolerance = 0.08
    passed = density_l1_error <= density_tolerance && pressure_l1_error <= pressure_tolerance
    return (
        passed = passed,
        summary = "Sod shock tube agrees with the exact Euler Riemann solution",
        metrics = Dict{String, Any}(
            "density_L1_error" => density_l1_error,
            "pressure_L1_error" => pressure_l1_error,
            "cells" => length(case.geo.cell_volumes),
            "final_time" => final_time,
        ),
        expected = Dict{String, Any}(
            "density_L1_maximum" => density_tolerance,
            "pressure_L1_maximum" => pressure_tolerance,
            "oracle" => "independent exact ideal-gas Riemann solution sampled at cell centers",
        ),
        diagnostics = Dict{String, Any}(
            "numerical_density" => collect(final_state.density),
            "exact_density" => exact_density,
            "exact_pressure" => exact_pressure,
        ),
    )
end
