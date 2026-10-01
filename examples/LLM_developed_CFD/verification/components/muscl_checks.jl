function _muscl_scalar_reconstructions(case, values)
    stencil = build_weighted_least_squares_stencil(case.geo)
    gradients = zeros(length(values), 3)
    populate_weighted_least_squares_gradient!(gradients, values, stencil)
    reconstructions = NamedTuple[]

    for connection_group in case.system.connection_groups
        for (idx_a, neighbor_list) in connection_group.cell_neighbors
            for (idx_b, face_a, face_b) in neighbor_list
                (
                    distance,
                    face_area_a, face_normal_a, face_distance_a, volume_a,
                    face_area_b, face_normal_b, face_distance_b, volume_b,
                ) = interface_geometry(case.geo, idx_a, face_a, idx_b, face_b)
                value_a = _MUSCL_limited_face_value(
                    values,
                    gradients,
                    idx_a,
                    face_normal_a,
                    face_distance_a,
                    case.geo,
                )
                value_b = _MUSCL_limited_face_value(
                    values,
                    gradients,
                    idx_b,
                    face_normal_b,
                    face_distance_b,
                    case.geo,
                )
                x_face_a = case.geo.cell_centroids[idx_a][1] + face_distance_a * face_normal_a[1]
                x_face_b = case.geo.cell_centroids[idx_b][1] + face_distance_b * face_normal_b[1]
                push!(reconstructions, (
                    cell = idx_a,
                    neighbor = idx_b,
                    value = value_a,
                    x_face = x_face_a,
                ))
                push!(reconstructions, (
                    cell = idx_b,
                    neighbor = idx_a,
                    value = value_b,
                    x_face = x_face_b,
                ))
            end
        end
    end
    return reconstructions
end

function check_muscl_constant_reconstruction(n_cells)
    case = build_compressible_case(n_cells)
    constant_value = 3.25
    reconstructed = _muscl_scalar_reconstructions(case, fill(constant_value, n_cells))
    maximum_error = maximum(abs(item.value - constant_value) for item in reconstructed)
    tolerance = 20.0 * eps(constant_value)
    return (
        passed = maximum_error <= tolerance,
        summary = "MUSCL exactly preserves a constant scalar field",
        metrics = Dict{String, Any}(
            "maximum_absolute_error" => maximum_error,
            "faces_checked" => length(reconstructed),
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "oracle" => "constant-field identity",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_muscl_linear_reconstruction(n_cells)
    case = build_compressible_case(n_cells)
    intercept = 2.0
    slope = 0.3
    values = [intercept + slope * centroid[1] for centroid in case.geo.cell_centroids]
    reconstructed = _muscl_scalar_reconstructions(case, values)
    errors = [abs(item.value - (intercept + slope * item.x_face)) for item in reconstructed]
    maximum_error = maximum(errors)

    # Venkatakrishnan regularization approaches the unlimited exact-linear
    # reconstruction with refinement, so finite-grid equality is not exact.
    tolerance = 2e-4 / n_cells
    return (
        passed = maximum_error <= tolerance,
        summary = "MUSCL reconstructs a linear field to its regularized finite-grid accuracy",
        metrics = Dict{String, Any}(
            "maximum_absolute_error" => maximum_error,
            "grid_cells" => n_cells,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "oracle" => "analytical linear value at each face",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function _muscl_smooth_error(n_cells, norm_name)
    case = build_compressible_case(n_cells)
    smooth_field(x) = 1.0 + 0.2 * sinpi(2.0 * x)
    values = [smooth_field(centroid[1]) for centroid in case.geo.cell_centroids]
    reconstructed = _muscl_scalar_reconstructions(case, values)
    interior = filter(
        item -> 2 < item.cell < n_cells - 1,
        reconstructed,
    )
    errors = [item.value - smooth_field(item.x_face) for item in interior]
    weights = ones(length(errors))
    return error_norm(errors, weights, norm_name)
end

function check_muscl_smooth_convergence(grid_sizes)
    study = run_convergence_study(
        grid_sizes,
        _muscl_smooth_error;
        norm_name = :L2,
    )
    measured_orders = filter(!isnothing, study.observed_orders)
    asymptotic_order = minimum(last(measured_orders, min(2, length(measured_orders))))
    minimum_order = 1.8
    return (
        passed = asymptotic_order >= minimum_order,
        summary = "smooth MUSCL face reconstruction approaches second order",
        metrics = Dict{String, Any}(
            "asymptotic_minimum_order" => asymptotic_order,
            "convergence_table" => convergence_table(study),
        ),
        expected = Dict{String, Any}(
            "minimum_asymptotic_order" => minimum_order,
            "norm" => "L2",
            "oracle" => "analytical periodic sine field at internal faces",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_muscl_discontinuity(n_cells)
    case = build_compressible_case(n_cells)
    lower_value = 0.5
    upper_value = 1.5
    values = zeros(n_cells)
    for cell_id in eachindex(values)
        if case.geo.cell_centroids[cell_id][1] < 0.5
            values[cell_id] = lower_value
        else
            values[cell_id] = upper_value
        end
    end
    reconstructed = _muscl_scalar_reconstructions(case, values)
    reconstructed_values = [item.value for item in reconstructed]
    all_finite = all(isfinite, reconstructed_values)
    minimum_value = minimum(reconstructed_values)
    maximum_value = maximum(reconstructed_values)
    allowed_overshoot = 0.05 * (upper_value - lower_value)
    within_robust_bounds =
        minimum_value >= lower_value - allowed_overshoot &&
        maximum_value <= upper_value + allowed_overshoot
    positive = minimum_value > 0.0

    return (
        passed = all_finite && within_robust_bounds && positive,
        summary = "limited MUSCL remains finite, positive, and bounded near a jump",
        metrics = Dict{String, Any}(
            "minimum_reconstructed_value" => minimum_value,
            "maximum_reconstructed_value" => maximum_value,
            "all_finite" => all_finite,
        ),
        expected = Dict{String, Any}(
            "minimum_allowed" => lower_value - allowed_overshoot,
            "maximum_allowed" => upper_value + allowed_overshoot,
            "strict_positivity" => true,
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function _maximum_dependency_radius(pattern, metadata)
    maximum_radius = 0
    for input_index in axes(pattern, 2)
        input_cell = metadata[input_index].cell
        for output_index in axes(pattern, 1)
            if pattern[output_index, input_index]
                radius = abs(metadata[output_index].cell - input_cell)
                maximum_radius = max(maximum_radius, radius)
            end
        end
    end
    return maximum_radius
end

function check_muscl_stencil_expansion(n_cells)
    first_order_case = build_compressible_case(
        n_cells;
        profile = :jacobian,
        spatial_method = :first_order,
    )
    muscl_case = build_compressible_case(
        n_cells;
        profile = :jacobian,
        spatial_method = :muscl,
    )
    first_pattern, first_thresholds = empirical_jacobian_sparsity(first_order_case)
    muscl_pattern, muscl_thresholds = empirical_jacobian_sparsity(muscl_case)
    first_metadata = state_index_metadata(first_order_case)
    muscl_metadata = state_index_metadata(muscl_case)
    first_radius = _maximum_dependency_radius(first_pattern, first_metadata)
    muscl_radius = _maximum_dependency_radius(muscl_pattern, muscl_metadata)
    muscl_declared = declared_jacobian_sparsity(muscl_case)
    missing_declared_dependencies = count(
        muscl_pattern[output_index, input_index] &&
        !muscl_declared[output_index, input_index]
        for output_index in axes(muscl_pattern, 1), input_index in axes(muscl_pattern, 2)
    )

    passed =
        first_radius <= 1 &&
        muscl_radius >= 2 &&
        nnz(muscl_pattern) > nnz(first_pattern) &&
        missing_declared_dependencies == 0
    return (
        passed = passed,
        summary = "MUSCL expands the empirical stencil and tracer sparsity contains it",
        metrics = Dict{String, Any}(
            "first_order_observed_nonzeros" => nnz(first_pattern),
            "muscl_observed_nonzeros" => nnz(muscl_pattern),
            "first_order_cell_radius" => first_radius,
            "muscl_cell_radius" => muscl_radius,
            "muscl_declared_nonzeros" => nnz(muscl_declared),
            "missing_declared_dependencies" => missing_declared_dependencies,
        ),
        expected = Dict{String, Any}(
            "first_order_maximum_radius" => 1,
            "muscl_minimum_radius" => 2,
            "missing_declared_dependencies" => 0,
        ),
        diagnostics = Dict{String, Any}(),
    )
end
