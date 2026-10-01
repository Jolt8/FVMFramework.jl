function _one_dimensional_transport_residual(
    field,
    diffusivity,
    source,
    left_value,
    right_value,
)
    n_cells = length(field)
    cell_width = 1.0 / n_cells
    if source isa Number
        residual_values = fill(source, n_cells)
    else
        if length(source) != n_cells
            throw(DimensionMismatch("the source vector must have one value per cell"))
        end
        residual_values = collect(source)
    end
    normal = [0.0, 1.0, 0.0]
    zero_gradient = (0.0, 0.0, 0.0)

    left_gradient = (field[1] - left_value) / (0.5 * cell_width)
    residual_values[1] -= diffusivity * left_gradient / cell_width

    for left_cell in 1:(n_cells - 1)
        right_cell = left_cell + 1
        face_gradient = corrected_face_gradient(
            field[left_cell],
            field[right_cell],
            zero_gradient...,
            zero_gradient...,
            normal,
            cell_width,
        )[2]
        face_flux = diffusivity * face_gradient
        residual_values[left_cell] += face_flux / cell_width
        residual_values[right_cell] -= face_flux / cell_width
    end

    right_gradient = (right_value - field[end]) / (0.5 * cell_width)
    residual_values[end] += diffusivity * right_gradient / cell_width
    return residual_values
end

function _steady_transport_solution(
    n_cells;
    diffusivity,
    source,
    left_value,
    right_value,
)
    zero_field = zeros(n_cells)
    affine_residual = _one_dimensional_transport_residual(
        zero_field,
        diffusivity,
        source,
        left_value,
        right_value,
    )
    operator_matrix = zeros(n_cells, n_cells)
    for column in 1:n_cells
        basis_field = zeros(n_cells)
        basis_field[column] = 1.0
        operator_matrix[:, column] .= _one_dimensional_transport_residual(
            basis_field,
            diffusivity,
            source,
            left_value,
            right_value,
        ) .- affine_residual
    end
    field = -(operator_matrix \ affine_residual)
    cell_centers = [(cell_id - 0.5) / n_cells for cell_id in 1:n_cells]
    residual_norm = norm(_one_dimensional_transport_residual(
        field,
        diffusivity,
        source,
        left_value,
        right_value,
    ), Inf)
    return field, cell_centers, residual_norm
end

function _analytical_transport_errors(calculated, cell_centers, analytical)
    errors = calculated .- analytical.(cell_centers)
    cell_width = 1.0 / length(calculated)
    return (
        L1 = sum(abs.(errors)) * cell_width,
        L2 = sqrt(sum(abs2, errors) * cell_width),
        Linf = maximum(abs.(errors)),
    )
end

function check_couette_flow(n_cells)
    kinematic_viscosity = 0.08
    lower_wall_velocity = 0.0
    upper_wall_velocity = 1.5
    velocity, centers, residual_norm = _steady_transport_solution(
        n_cells;
        diffusivity = kinematic_viscosity,
        source = 0.0,
        left_value = lower_wall_velocity,
        right_value = upper_wall_velocity,
    )
    analytical(y) = lower_wall_velocity +
        (upper_wall_velocity - lower_wall_velocity) * y
    errors = _analytical_transport_errors(velocity, centers, analytical)
    expected_shear_rate = upper_wall_velocity - lower_wall_velocity
    calculated_shear_rate = (velocity[2] - velocity[1]) * n_cells
    shear_error = abs(calculated_shear_rate - expected_shear_rate)
    tolerance = 2e-12
    return (
        passed = errors.Linf <= tolerance && shear_error <= tolerance && residual_norm <= tolerance,
        summary = "plane Couette flow reproduces the linear velocity and constant shear",
        metrics = Dict{String, Any}(
            "L1_error" => errors.L1,
            "L2_error" => errors.L2,
            "Linf_error" => errors.Linf,
            "shear_rate_error" => shear_error,
            "steady_residual_Linf" => residual_norm,
            "cells" => n_cells,
        ),
        expected = Dict{String, Any}(
            "Linf_tolerance" => tolerance,
            "analytical_profile" => "u(y) = U*y/H",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_poiseuille_flow(grid_sizes)
    kinematic_viscosity = 0.08
    acceleration = 0.6
    studies = NamedTuple[]
    for n_cells in grid_sizes
        velocity, centers, residual_norm = _steady_transport_solution(
            n_cells;
            diffusivity = kinematic_viscosity,
            source = acceleration,
            left_value = 0.0,
            right_value = 0.0,
        )
        analytical(y) = acceleration * y * (1.0 - y) / (2.0 * kinematic_viscosity)
        errors = _analytical_transport_errors(velocity, centers, analytical)
        push!(studies, (
            cells = n_cells,
            L2 = errors.L2,
            Linf = errors.Linf,
            residual = residual_norm,
        ))
    end
    observed_orders = [
        log(studies[index - 1].L2 / studies[index].L2) /
        log(studies[index].cells / studies[index - 1].cells)
        for index in 2:length(studies)
    ]
    asymptotic_order = minimum(last(observed_orders, min(2, length(observed_orders))))
    minimum_order = 1.8
    finest_tolerance = 1e-3
    passed =
        asymptotic_order >= minimum_order &&
        studies[end].Linf <= finest_tolerance &&
        studies[end].residual <= 1e-10
    return (
        passed = passed,
        summary = "plane Poiseuille flow converges to the pressure-driven parabola",
        metrics = Dict{String, Any}(
            "grid_results" => [Dict(
                "cells" => item.cells,
                "L2_error" => item.L2,
                "Linf_error" => item.Linf,
                "steady_residual_Linf" => item.residual,
            ) for item in studies],
            "observed_L2_orders" => observed_orders,
            "asymptotic_minimum_order" => asymptotic_order,
        ),
        expected = Dict{String, Any}(
            "minimum_asymptotic_order" => minimum_order,
            "finest_Linf_tolerance" => finest_tolerance,
            "analytical_profile" => "u(y) = g*y*(H-y)/(2*nu)",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_one_dimensional_heat_conduction(n_cells)
    thermal_diffusivity = 0.12
    cold_temperature = 280.0
    hot_temperature = 430.0
    temperature, centers, residual_norm = _steady_transport_solution(
        n_cells;
        diffusivity = thermal_diffusivity,
        source = 0.0,
        left_value = cold_temperature,
        right_value = hot_temperature,
    )
    analytical(y) = cold_temperature + (hot_temperature - cold_temperature) * y
    errors = _analytical_transport_errors(temperature, centers, analytical)
    expected_gradient = hot_temperature - cold_temperature
    calculated_gradient = (temperature[2] - temperature[1]) * n_cells
    gradient_error = abs(calculated_gradient - expected_gradient)
    tolerance = 2e-10
    residual_scale = thermal_diffusivity * abs(expected_gradient)
    scaled_residual = residual_norm / residual_scale
    scaled_residual_tolerance = 1e-10
    return (
        passed =
            errors.Linf <= tolerance &&
            gradient_error <= tolerance &&
            scaled_residual <= scaled_residual_tolerance,
        summary = "one-dimensional conduction reproduces a linear temperature field and heat flux",
        metrics = Dict{String, Any}(
            "L1_error_K" => errors.L1,
            "L2_error_K" => errors.L2,
            "Linf_error_K" => errors.Linf,
            "temperature_gradient_error_K_per_m" => gradient_error,
            "steady_residual_Linf" => residual_norm,
            "scaled_steady_residual_Linf" => scaled_residual,
            "cells" => n_cells,
        ),
        expected = Dict{String, Any}(
            "Linf_tolerance_K" => tolerance,
            "scaled_residual_tolerance" => scaled_residual_tolerance,
            "analytical_profile" => "T(y) = T_cold + (T_hot-T_cold)*y/H",
        ),
        diagnostics = Dict{String, Any}(),
    )
end
