function declared_jacobian_sparsity(case, state = case.u0)
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
    derivative = zeros(length(state))
    sparsity = ADTypes.jacobian_sparsity(
        (output, input) -> case.rhs!(output, input, case.p, 0.0),
        derivative,
        state,
        detector,
    )
    return sparse(Bool.(sparsity))
end

function _state_scales(case, state)
    named_state = ComponentVector(state, case.system.state_axes)
    scales = similar(state)
    named_scales = ComponentVector(scales, case.system.state_axes)
    for variable in propertynames(named_state)
        values = getproperty(named_state, variable)
        variable_scale = max(norm(values, Inf), 1.0)
        getproperty(named_scales, variable) .= variable_scale
    end
    return scales
end

function empirical_jacobian_sparsity(case, state = case.u0)
    baseline = residual(case, state)
    scales = _state_scales(case, state)
    n_variables = length(state)
    observed = falses(n_variables, n_variables)
    column_thresholds = zeros(n_variables)

    for input_index in 1:n_variables
        step = sqrt(eps(Float64)) * max(abs(state[input_index]), scales[input_index])
        positive_state = copy(state)
        negative_state = copy(state)
        positive_state[input_index] += step
        negative_state[input_index] -= step

        positive_residual = residual(case, positive_state)
        negative_residual = residual(case, negative_state)
        derivative_column = (positive_residual .- negative_residual) ./ (2.0 * step)
        column_scale = max(norm(derivative_column, Inf), norm(baseline, Inf), 1.0)
        threshold = 2e-7 * column_scale
        column_thresholds[input_index] = threshold
        observed[:, input_index] .= abs.(derivative_column) .> threshold
    end
    return sparse(observed), column_thresholds
end

function _dependency_cells(pattern, column, metadata)
    return sort!(unique(
        metadata[row].cell for row in axes(pattern, 1) if pattern[row, column]
    ))
end

function check_empirical_sparsity(case)
    state = copy(case.u0)
    declared = declared_jacobian_sparsity(case, state)
    observed, thresholds = empirical_jacobian_sparsity(case, state)
    metadata = state_index_metadata(case)

    unexpected = Dict{String, Any}[]
    inactive_declared_count = 0
    for input_index in axes(observed, 2)
        for output_index in axes(observed, 1)
            if observed[output_index, input_index] && !declared[output_index, input_index]
                if length(unexpected) < 25
                    push!(unexpected, Dict{String, Any}(
                        "perturbed_index" => input_index,
                        "perturbed_variable" => string(metadata[input_index].variable),
                        "perturbed_cell" => metadata[input_index].cell,
                        "residual_index" => output_index,
                        "residual_variable" => string(metadata[output_index].variable),
                        "residual_cell" => metadata[output_index].cell,
                        "observed_dependency_cells" => _dependency_cells(observed, input_index, metadata),
                        "declared_dependency_cells" => _dependency_cells(declared, input_index, metadata),
                        "finite_difference_threshold" => thresholds[input_index],
                    ))
                end
            elseif declared[output_index, input_index] && !observed[output_index, input_index]
                inactive_declared_count += 1
            end
        end
    end

    unexpected_count = count(
        observed[output_index, input_index] && !declared[output_index, input_index]
        for output_index in axes(observed, 1), input_index in axes(observed, 2)
    )
    return (
        passed = unexpected_count == 0,
        summary = "finite-difference dependencies are contained in tracer-declared sparsity",
        metrics = Dict{String, Any}(
            "degrees_of_freedom" => length(state),
            "observed_nonzeros" => nnz(observed),
            "declared_nonzeros" => nnz(declared),
            "unexpected_dependency_count" => unexpected_count,
            "state_inactive_declared_count" => inactive_declared_count,
        ),
        expected = Dict{String, Any}(
            "unexpected_dependency_count" => 0,
            "comparison_policy" => "declared supersets are allowed because one state cannot activate every branch",
        ),
        diagnostics = Dict{String, Any}(
            "unexpected_dependencies" => unexpected,
        ),
    )
end

function check_ad_jvp(case, random_seed)
    state = copy(case.u0)
    scales = _state_scales(case, state)
    jacobian_ad = ForwardDiff.jacobian(input -> residual(case, input), state)
    random_number_generator = MersenneTwister(random_seed)
    direction_reports = Dict{String, Any}[]
    maximum_relative_error = 0.0
    classifications = String[]
    tolerance = 2e-5
    base_step = cbrt(eps(Float64))
    step_multipliers = (4.0, 1.0, 0.25)

    for direction_id in 1:3
        dimensionless_direction = randn(random_number_generator, length(state))
        dimensionless_direction ./= norm(dimensionless_direction)
        direction = scales .* dimensionless_direction
        ad_product = jacobian_ad * direction
        step_sweep = Dict{String, Any}[]
        best_relative_error = Inf
        best_absolute_error = Inf
        best_maximum_component_error = Inf
        best_maximum_error_index = 0
        best_step = base_step

        for step_multiplier in step_multipliers
            step = base_step * step_multiplier
            positive_residual = residual(case, state .+ step .* direction)
            negative_residual = residual(case, state .- step .* direction)
            finite_difference_product = (positive_residual .- negative_residual) ./ (2.0 * step)
            difference = ad_product .- finite_difference_product
            absolute_error = norm(difference)
            relative_error = absolute_error / max(
                norm(ad_product),
                norm(finite_difference_product),
                eps(Float64),
            )
            maximum_component_error, maximum_error_index = findmax(abs.(difference))
            push!(step_sweep, Dict{String, Any}(
                "step" => step,
                "absolute_error_norm" => absolute_error,
                "relative_error" => relative_error,
                "maximum_component_error" => maximum_component_error,
                "maximum_error_index" => maximum_error_index,
            ))
            if relative_error < best_relative_error
                best_relative_error = relative_error
                best_absolute_error = absolute_error
                best_maximum_component_error = maximum_component_error
                best_maximum_error_index = maximum_error_index
                best_step = step
            end
        end

        if best_relative_error <= tolerance
            classification = "good agreement"
        elseif best_relative_error <= 10.0 * tolerance
            classification = "probable finite-difference truncation/roundoff sensitivity"
        else
            classification = "definite derivative mismatch"
        end
        push!(classifications, classification)
        maximum_relative_error = max(maximum_relative_error, best_relative_error)
        push!(direction_reports, Dict{String, Any}(
            "direction" => direction_id,
            "absolute_error_norm" => best_absolute_error,
            "relative_error" => best_relative_error,
            "maximum_component_error" => best_maximum_component_error,
            "maximum_error_index" => best_maximum_error_index,
            "selected_step" => best_step,
            "step_sweep" => step_sweep,
            "classification" => classification,
        ))
    end

    return (
        passed = maximum_relative_error <= tolerance,
        summary = "ForwardDiff Jacobian-vector products agree with central finite differences",
        metrics = Dict{String, Any}(
            "maximum_relative_error" => maximum_relative_error,
            "directions" => direction_reports,
            "base_finite_difference_step" => base_step,
        ),
        expected = Dict{String, Any}(
            "relative_tolerance" => tolerance,
            "random_seed" => random_seed,
            "finite_difference_scheme" => "centered",
        ),
        diagnostics = Dict{String, Any}(
            "classifications" => classifications,
        ),
    )
end
