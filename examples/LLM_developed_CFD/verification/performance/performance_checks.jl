const PERFORMANCE_BASELINE_PATH = joinpath(@__DIR__, "performance_baselines.toml")

function _performance_baseline_name(case, name_suffix, time_method)
    spatial_name = isempty(name_suffix) ? "first_order" : name_suffix
    n_cells = length(ComponentVector(case.u0, case.system.state_axes).density)
    return "$(n_cells)_cells_$(spatial_name)_$(time_method)"
end

function _performance_limit(reference_value, relative_allowance, absolute_allowance)
    return ceil(Int, reference_value * (1.0 + relative_allowance) + absolute_allowance)
end

function check_performance_baseline(baseline_name, statistics)
    baseline_document = TOML.parsefile(PERFORMANCE_BASELINE_PATH)
    baselines = baseline_document["baselines"]
    if !haskey(baselines, baseline_name)
        return _verification_result(
            "performance_$(baseline_name)",
            "performance",
            false,
            "no persistent solver-performance baseline exists for this case";
            diagnostics = Dict{String, Any}(
                "baseline_path" => PERFORMANCE_BASELINE_PATH,
                "missing_baseline" => baseline_name,
            ),
        )
    end

    baseline = baselines[baseline_name]
    count_metrics = (
        "accepted_steps",
        "rejected_steps",
        "rhs_evaluations",
        "linear_solves",
        "nonlinear_iterations",
    )
    comparisons = Dict{String, Any}()
    missing_metrics = String[]
    count_regressions = String[]
    for metric_name in count_metrics
        if !haskey(statistics, metric_name)
            push!(missing_metrics, metric_name)
            continue
        end
        reference_value = baseline[metric_name]
        maximum_value = _performance_limit(
            reference_value,
            baseline["count_relative_allowance"],
            baseline["count_absolute_allowance"],
        )
        observed_value = statistics[metric_name]
        comparisons[metric_name] = Dict{String, Any}(
            "reference" => reference_value,
            "maximum" => maximum_value,
            "observed" => observed_value,
        )
        if observed_value > maximum_value
            push!(count_regressions, metric_name)
        end
    end

    minimum_timestep_regression = false
    if !haskey(statistics, "minimum_timestep") || statistics["minimum_timestep"] === nothing
        push!(missing_metrics, "minimum_timestep")
    else
        reference_timestep = baseline["minimum_timestep"]
        minimum_allowed_timestep = reference_timestep / baseline["minimum_timestep_factor"]
        observed_timestep = statistics["minimum_timestep"]
        comparisons["minimum_timestep"] = Dict{String, Any}(
            "reference" => reference_timestep,
            "minimum" => minimum_allowed_timestep,
            "observed" => observed_timestep,
        )
        minimum_timestep_regression = observed_timestep < minimum_allowed_timestep
    end

    runtime_warning = false
    if haskey(statistics, "runtime_seconds")
        reference_runtime = baseline["runtime_seconds"]
        maximum_runtime = max(
            reference_runtime * baseline["runtime_factor"],
            reference_runtime + baseline["runtime_absolute_allowance_seconds"],
        )
        observed_runtime = statistics["runtime_seconds"]
        comparisons["runtime_seconds"] = Dict{String, Any}(
            "reference" => reference_runtime,
            "warning_threshold" => maximum_runtime,
            "observed" => observed_runtime,
            "gating" => false,
        )
        runtime_warning = observed_runtime > maximum_runtime
    else
        push!(missing_metrics, "runtime_seconds")
    end

    passed =
        isempty(missing_metrics) &&
        isempty(count_regressions) &&
        !minimum_timestep_regression
    return _verification_result(
        "performance_$(baseline_name)",
        "performance",
        passed,
        "solver work counters and minimum timestep remain within the persistent baseline";
        metrics = Dict{String, Any}(
            "comparisons" => comparisons,
            "runtime_warning" => runtime_warning,
        ),
        expected = Dict{String, Any}(
            "count_metrics" => "no larger than their baseline allowances",
            "minimum_timestep" => "no smaller than its baseline allowance",
            "runtime" => "non-gating secondary warning",
        ),
        diagnostics = Dict{String, Any}(
            "baseline_path" => PERFORMANCE_BASELINE_PATH,
            "missing_metrics" => missing_metrics,
            "count_regressions" => count_regressions,
            "minimum_timestep_regression" => minimum_timestep_regression,
            "runtime_warning" => runtime_warning,
        ),
    )
end
