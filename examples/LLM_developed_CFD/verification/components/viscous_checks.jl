function check_weighted_least_squares_linear_gradient(n_cells)
    case = build_compressible_case(n_cells)
    intercept = 1.7
    slope = -0.45
    values = [intercept + slope * center[1] for center in case.geo.cell_centroids]
    gradients = zeros(n_cells, 3)
    stencil = build_weighted_least_squares_stencil(case.geo)
    populate_weighted_least_squares_gradient!(gradients, values, stencil)

    expected = repeat(reshape([slope, 0.0, 0.0], 1, 3), n_cells, 1)
    maximum_error = maximum(abs.(gradients .- expected))
    tolerance = 200.0 * eps(Float64)
    return (
        passed = maximum_error <= tolerance,
        summary = "weighted least squares recovers a resolved linear gradient",
        metrics = Dict{String, Any}(
            "maximum_absolute_error" => maximum_error,
            "cells_checked" => n_cells,
        ),
        expected = Dict{String, Any}(
            "gradient" => [slope, 0.0, 0.0],
            "absolute_tolerance" => tolerance,
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_corrected_face_gradient_projection()
    inverse_root_three = inv(sqrt(3.0))
    normal = [inverse_root_three, inverse_root_three, inverse_root_three]
    distance = 0.37
    value_a = -0.8
    prescribed_normal_derivative = 1.25
    value_b = value_a + distance * prescribed_normal_derivative
    gradient_a = [0.4, -1.1, 0.7]
    gradient_b = [-0.2, 0.3, 1.5]
    average_gradient = 0.5 .* (gradient_a .+ gradient_b)
    expected_gradient = average_gradient .+ (
        prescribed_normal_derivative - dot(average_gradient, normal)
    ) .* normal

    calculated_gradient = collect(corrected_face_gradient(
        value_a,
        value_b,
        gradient_a[1], gradient_a[2], gradient_a[3],
        gradient_b[1], gradient_b[2], gradient_b[3],
        normal,
        distance,
    ))
    maximum_error = maximum(abs.(calculated_gradient .- expected_gradient))
    normal_error = abs(dot(calculated_gradient, normal) - prescribed_normal_derivative)
    tangential_error = norm(
        (calculated_gradient .- dot(calculated_gradient, normal) .* normal) .-
        (average_gradient .- dot(average_gradient, normal) .* normal),
    )
    tolerance = 100.0 * eps(Float64)
    return (
        passed = maximum_error <= tolerance && normal_error <= tolerance && tangential_error <= tolerance,
        summary = "face correction replaces only the full normal-gradient component",
        metrics = Dict{String, Any}(
            "maximum_component_error" => maximum_error,
            "normal_derivative_error" => normal_error,
            "tangential_gradient_error" => tangential_error,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "normal_derivative" => prescribed_normal_derivative,
            "tangential_component" => "unchanged from the interpolated cell gradient",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_all_transport_gradients_are_corrected()
    normal = normalize([0.6, -0.3, 0.7])
    distance = 0.52
    field_data = (
        (name = "u", value_a = 0.1, derivative = 0.8, grad_a = [0.4, 0.2, -0.5], grad_b = [-0.1, 0.7, 0.3]),
        (name = "v", value_a = -0.4, derivative = -1.2, grad_a = [0.9, -0.8, 0.2], grad_b = [0.5, 0.6, -0.4]),
        (name = "w", value_a = 1.1, derivative = 0.35, grad_a = [-0.7, 0.1, 0.8], grad_b = [0.2, -0.9, 0.4]),
        (name = "temperature", value_a = 290.0, derivative = 17.0, grad_a = [8.0, -2.0, 3.0], grad_b = [-4.0, 5.0, 9.0]),
    )

    errors = Dict{String, Float64}()
    for field in field_data
        value_b = field.value_a + distance * field.derivative
        corrected = collect(corrected_face_gradient(
            field.value_a,
            value_b,
            field.grad_a[1], field.grad_a[2], field.grad_a[3],
            field.grad_b[1], field.grad_b[2], field.grad_b[3],
            normal,
            distance,
        ))
        errors[field.name] = abs(dot(corrected, normal) - field.derivative)
    end

    maximum_error = maximum(values(errors))
    tolerance = 1e-12
    return (
        passed = maximum_error <= tolerance,
        summary = "u, v, w, and temperature use the same corrected normal gradient",
        metrics = Dict{String, Any}(
            "normal_derivative_errors" => errors,
            "maximum_error" => maximum_error,
        ),
        expected = Dict{String, Any}("absolute_tolerance" => tolerance),
        diagnostics = Dict{String, Any}(),
    )
end

function check_sst_admissibility_guards()
    state = ComponentVector(
        density = [1.0, 1.0],
        momentum_density_u = zeros(2),
        momentum_density_v = zeros(2),
        momentum_density_w = zeros(2),
        volumetric_energy = fill(2.5, 2),
        turbulent_kinetic_energy_density = [0.1, -0.01],
        specific_dissipation_rate_density = [2.0, 0.0],
    )
    violations = state_violations(
        collect(state),
        getaxes(state);
        gamma = 1.4,
        cv = 717.5,
    )
    variables = Set(violation["variable"] for violation in violations)
    expected_variables = Set([
        "turbulent_kinetic_energy_density",
        "specific_dissipation_rate_density",
    ])
    return (
        passed = variables == expected_variables,
        summary = "SST k and omega admissibility guards reject non-positive states",
        metrics = Dict{String, Any}(
            "reported_variables" => sort!(collect(variables)),
            "violation_count" => length(violations),
        ),
        expected = Dict{String, Any}(
            "reported_variables" => sort!(collect(expected_variables)),
            "constraint" => "k > 0 and omega > 0",
        ),
        diagnostics = Dict{String, Any}("violations" => violations),
    )
end

function _first_internal_face(case)
    for connection_group in case.system.connection_groups
        for (idx_a, neighbor_list) in connection_group.cell_neighbors
            for (idx_b, face_a, face_b) in neighbor_list
                return idx_a, face_a, idx_b, face_b
            end
        end
    end
    error("verification case has no internal face")
end

function _independent_face_gradient(value_a, value_b, gradient_a, gradient_b, normal, distance)
    average_gradient = 0.5 .* (gradient_a .+ gradient_b)
    normal_derivative = (value_b - value_a) / distance
    return average_gradient .+ (
        normal_derivative - dot(average_gradient, normal)
    ) .* normal
end

function check_complete_viscous_face_flux()
    case = build_compressible_case(2)
    idx_a, face_a, idx_b, face_b = _first_internal_face(case)
    distance, area_a, normal_a, _, _, area_b, _, _, _ = interface_geometry(
        case.geo,
        idx_a,
        face_a,
        idx_b,
        face_b,
    )
    state = ComponentVector(
        vel_u = [0.2, 0.5],
        vel_v = [-0.3, 0.4],
        vel_w = [0.8, -0.2],
        temperature = [300.0, 340.0],
        mu = [0.01, 0.02],
        k = [0.2, 0.4],
        grad_vel_u = [0.4 -0.2 0.7; -0.1 0.5 0.3],
        grad_vel_v = [0.8 0.1 -0.4; -0.3 0.6 0.2],
        grad_vel_w = [-0.5 0.9 0.1; 0.4 -0.2 0.8],
        grad_temperature = [12.0 -4.0 7.0; -5.0 9.0 3.0],
    )
    derivative = ComponentVector(
        momentum_density_u_flow = zeros(2),
        momentum_density_v_flow = zeros(2),
        momentum_density_w_flow = zeros(2),
        volumetric_energy_flow = zeros(2),
    )
    fluid_viscous_and_diffusive_flux!(
        derivative,
        state,
        nothing,
        0.0,
        nothing,
        case.geo,
        idx_a,
        face_a,
        idx_b,
        face_b,
    )

    gradient_u = _independent_face_gradient(
        state.vel_u[idx_a], state.vel_u[idx_b],
        collect(state.grad_vel_u[idx_a, :]), collect(state.grad_vel_u[idx_b, :]),
        normal_a, distance,
    )
    gradient_v = _independent_face_gradient(
        state.vel_v[idx_a], state.vel_v[idx_b],
        collect(state.grad_vel_v[idx_a, :]), collect(state.grad_vel_v[idx_b, :]),
        normal_a, distance,
    )
    gradient_w = _independent_face_gradient(
        state.vel_w[idx_a], state.vel_w[idx_b],
        collect(state.grad_vel_w[idx_a, :]), collect(state.grad_vel_w[idx_b, :]),
        normal_a, distance,
    )
    gradient_temperature = _independent_face_gradient(
        state.temperature[idx_a], state.temperature[idx_b],
        collect(state.grad_temperature[idx_a, :]), collect(state.grad_temperature[idx_b, :]),
        normal_a, distance,
    )
    dynamic_viscosity = 2.0 * state.mu[idx_a] * state.mu[idx_b] /
        (state.mu[idx_a] + state.mu[idx_b])
    conductivity = 2.0 * state.k[idx_a] * state.k[idx_b] /
        (state.k[idx_a] + state.k[idx_b])
    velocity_divergence = gradient_u[1] + gradient_v[2] + gradient_w[3]
    stress = [
        2.0 * dynamic_viscosity * gradient_u[1] - (2.0 / 3.0) * dynamic_viscosity * velocity_divergence  dynamic_viscosity * (gradient_u[2] + gradient_v[1])  dynamic_viscosity * (gradient_u[3] + gradient_w[1]);
        dynamic_viscosity * (gradient_u[2] + gradient_v[1])  2.0 * dynamic_viscosity * gradient_v[2] - (2.0 / 3.0) * dynamic_viscosity * velocity_divergence  dynamic_viscosity * (gradient_v[3] + gradient_w[2]);
        dynamic_viscosity * (gradient_u[3] + gradient_w[1])  dynamic_viscosity * (gradient_v[3] + gradient_w[2])  2.0 * dynamic_viscosity * gradient_w[3] - (2.0 / 3.0) * dynamic_viscosity * velocity_divergence
    ]
    traction = stress * collect(normal_a)
    face_velocity = 0.5 .* [
        state.vel_u[idx_a] + state.vel_u[idx_b],
        state.vel_v[idx_a] + state.vel_v[idx_b],
        state.vel_w[idx_a] + state.vel_w[idx_b],
    ]
    energy_flux = dot(face_velocity, traction) +
        conductivity * dot(gradient_temperature, normal_a)
    expected = [
        area_a .* traction;
        area_a * energy_flux;
        -area_b .* traction;
        -area_b * energy_flux;
    ]
    observed = [
        derivative.momentum_density_u_flow[idx_a],
        derivative.momentum_density_v_flow[idx_a],
        derivative.momentum_density_w_flow[idx_a],
        derivative.volumetric_energy_flow[idx_a],
        derivative.momentum_density_u_flow[idx_b],
        derivative.momentum_density_v_flow[idx_b],
        derivative.momentum_density_w_flow[idx_b],
        derivative.volumetric_energy_flow[idx_b],
    ]
    maximum_error = maximum(abs.(observed .- expected))
    conservation_error = maximum(abs.(observed[1:4] .+ observed[5:8]))
    tolerance = 2e-13
    return (
        passed = maximum_error <= tolerance && conservation_error <= tolerance,
        summary = "the complete Newtonian stress and Fourier face flux match an independent oracle",
        metrics = Dict{String, Any}(
            "maximum_flux_error" => maximum_error,
            "face_conservation_error" => conservation_error,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "oracle" => "independent tensor stress, traction, viscous work, and heat flux",
        ),
        diagnostics = Dict{String, Any}(),
    )
end
