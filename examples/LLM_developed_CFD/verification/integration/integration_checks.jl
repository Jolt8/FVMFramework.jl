function _solver_statistics(solution, timed_result, algorithm_name)
    statistics = Dict{String, Any}(
        "algorithm" => algorithm_name,
        "return_code" => string(solution.retcode),
        "runtime_seconds" => timed_result.time,
        "allocated_bytes" => timed_result.bytes,
        "garbage_collection_seconds" => timed_result.gctime,
    )
    statistic_names = Dict(
        :nf => "rhs_evaluations",
        :nf2 => "secondary_rhs_evaluations",
        :njacs => "jacobian_evaluations",
        :nsolve => "linear_solves",
        :nnonliniter => "nonlinear_iterations",
        :nnonlinconvfail => "nonlinear_convergence_failures",
        :naccept => "accepted_steps",
        :nreject => "rejected_steps",
    )
    for (field_name, output_name) in statistic_names
        if hasproperty(solution.destats, field_name)
            statistics[output_name] = getproperty(solution.destats, field_name)
        end
    end

    if length(solution.t) > 1
        saved_steps = diff(solution.t)
        statistics["minimum_timestep"] = minimum(saved_steps)
        statistics["minimum_saved_timestep"] = minimum(saved_steps)
        statistics["maximum_saved_timestep"] = maximum(saved_steps)
    else
        statistics["minimum_timestep"] = nothing
        statistics["minimum_saved_timestep"] = nothing
        statistics["maximum_saved_timestep"] = nothing
    end
    return statistics
end

function _solution_violations(case, solution)
    violations = Dict{String, Any}[]
    for (saved_index, state) in enumerate(solution.u)
        append!(violations, state_violations(
            state,
            case.system.state_axes;
            gamma = case.gamma,
            cv = case.cv,
            time = solution.t[saved_index],
            maximum_violations = max(20 - length(violations), 0),
            species_sum_absolute_tolerance = 5e-8,
            species_sum_relative_tolerance = 5e-8,
        ))
        if length(violations) >= 20
            break
        end
    end
    return violations
end

function _initial_maximum_signal_speed(case)
    state = ComponentVector(case.u0, case.system.state_axes)
    maximum_speed = 0.0
    for cell_id in eachindex(state.density)
        density = state.density[cell_id]
        vel_u = state.momentum_density_u[cell_id] / density
        kinetic_energy_density = 0.5 * (
            state.momentum_density_u[cell_id]^2 +
            state.momentum_density_v[cell_id]^2 +
            state.momentum_density_w[cell_id]^2
        ) / density
        pressure = (case.gamma - 1.0) * (
            state.volumetric_energy[cell_id] - kinetic_energy_density
        )
        sound_speed = sqrt(case.gamma * pressure / density)
        maximum_speed = max(maximum_speed, abs(vel_u) + sound_speed)
    end
    return maximum_speed
end

function run_explicit_smoke(case, final_time)
    problem = ODEProblem(case.rhs!, case.u0, (0.0, final_time), case.p)
    cell_width = minimum(case.geo.cell_volumes)
    stable_timestep = 0.2 * cell_width / _initial_maximum_signal_speed(case)
    timed_result = @timed solve(
        problem,
        SSPRK43();
        adaptive = false,
        dt = min(stable_timestep, final_time),
        save_everystep = true,
        maxiters = 100000,
    )
    solution = timed_result.value
    statistics = _solver_statistics(solution, timed_result, "SSPRK43")
    violations = _solution_violations(case, solution)
    successful = OrdinaryDiffEq.SciMLBase.successful_retcode(solution)
    reached_final_time = isapprox(solution.t[end], final_time; atol = 100.0 * eps(final_time), rtol = 1e-12)
    return (
        solution = solution,
        passed = successful && reached_final_time && isempty(violations),
        statistics = statistics,
        violations = violations,
        stable_timestep = stable_timestep,
    )
end

function run_implicit_smoke(case, final_time)
    sparsity = declared_jacobian_sparsity(case)
    ode_function = ODEFunction(case.rhs!; jac_prototype = float.(sparsity))
    problem = ODEProblem(ode_function, case.u0, (0.0, final_time), case.p)
    invalid_state = function (state, parameters, time)
        return !state_is_valid(
            state,
            case.system.state_axes;
            gamma = case.gamma,
            cv = case.cv,
            time = time,
            species_sum_absolute_tolerance = 1e-8,
            species_sum_relative_tolerance = 1e-6,
        )
    end
    cell_width = minimum(case.geo.cell_volumes)
    explicit_stable_timestep = 0.2 * cell_width / _initial_maximum_signal_speed(case)
    initial_timestep = min(0.05 * explicit_stable_timestep, final_time)
    if case.spatial_method == :muscl
        algorithm = FBDF(autodiff = ADTypes.AutoFiniteDiff())
        algorithm_name = "FBDF(AutoFiniteDiff)"
    else
        algorithm = FBDF()
        algorithm_name = "FBDF(ForwardDiff)"
    end
    timed_result = @timed solve(
        problem,
        algorithm;
        abstol = 1e-8,
        reltol = 1e-6,
        dt = initial_timestep,
        isoutofdomain = invalid_state,
        save_everystep = true,
        maxiters = 100000,
    )
    solution = timed_result.value
    statistics = _solver_statistics(solution, timed_result, algorithm_name)
    violations = _solution_violations(case, solution)
    successful = OrdinaryDiffEq.SciMLBase.successful_retcode(solution)
    reached_final_time = isapprox(solution.t[end], final_time; atol = 100.0 * eps(final_time), rtol = 1e-12)
    return (
        solution = solution,
        passed = successful && reached_final_time && isempty(violations),
        statistics = statistics,
        violations = violations,
        initial_timestep = initial_timestep,
    )
end

function _integration_classification(explicit_run, implicit_run, name)
    if explicit_run.passed && implicit_run.passed
        categories = String[]
        interpretation = "both spatial/time-integration paths completed"
        passed = true
    elseif explicit_run.passed && !implicit_run.passed
        categories = [
            "jacobian",
            "sparsity",
            "nonlinear solver",
            "preconditioner",
            "nonsmooth residual behavior",
        ]
        interpretation = "explicit passed while implicit failed"
        passed = false
    else
        categories = [
            "spatial discretization",
            "state admissibility",
            "boundary conditions",
            "instability",
            "implementation error",
        ]
        interpretation = "explicit failed; spatial or general integration issues are likely"
        passed = false
    end

    if explicit_run.passed
        explicit_status = "pass"
    else
        explicit_status = "fail"
    end
    if implicit_run.passed
        implicit_status = "pass"
    else
        implicit_status = "fail"
    end

    return _verification_result(
        name,
        "integration",
        passed,
        interpretation;
        metrics = Dict{String, Any}(
            "explicit_status" => explicit_status,
            "implicit_status" => implicit_status,
        ),
        expected = Dict{String, Any}(
            "explicit_status" => "pass",
            "implicit_status" => "pass",
        ),
        diagnostics = Dict{String, Any}(
            "likely_categories" => categories,
            "qualification" => "categories are diagnostic leads, not mathematical proof",
        ),
    )
end

function _suffixed_check_name(base_name, suffix)
    if isempty(suffix)
        return base_name
    end
    return "$(base_name)_$(suffix)"
end

function run_integration_checks(case, final_time; name_suffix = "")
    results = VerificationResult[]
    explicit_run = nothing
    implicit_run = nothing
    explicit_name = _suffixed_check_name("explicit_smoke", name_suffix)
    implicit_name = _suffixed_check_name("implicit_smoke", name_suffix)
    classification_name = _suffixed_check_name(
        "explicit_vs_implicit_classification",
        name_suffix,
    )
    benchmark_name = _suffixed_check_name("sod_shock_tube", name_suffix)

    explicit_wrapper = _run_check(explicit_name, "integration") do
        explicit_run = run_explicit_smoke(case, final_time)
        return (
            passed = explicit_run.passed,
            summary = "SSPRK43 completed the diagnostic Sod solve with admissible states",
            metrics = Dict{String, Any}(
                "solver_statistics" => explicit_run.statistics,
                "CFL_timestep" => explicit_run.stable_timestep,
            ),
            expected = Dict{String, Any}(
                "successful_return_code" => true,
                "state_violations" => 0,
                "CFL_number" => 0.2,
            ),
            diagnostics = Dict{String, Any}(
                "state_violations" => explicit_run.violations,
            ),
        )
    end
    push!(results, explicit_wrapper)

    implicit_wrapper = _run_check(implicit_name, "integration") do
        implicit_run = run_implicit_smoke(case, final_time)
        return (
            passed = implicit_run.passed,
            summary = "FBDF completed the diagnostic Sod solve with admissible states",
            metrics = Dict{String, Any}(
                "solver_statistics" => implicit_run.statistics,
                "initial_timestep" => implicit_run.initial_timestep,
            ),
            expected = Dict{String, Any}(
                "successful_return_code" => true,
                "state_violations" => 0,
            ),
            diagnostics = Dict{String, Any}(
                "state_violations" => implicit_run.violations,
            ),
        )
    end
    push!(results, implicit_wrapper)

    if explicit_run !== nothing
        baseline_name = _performance_baseline_name(case, name_suffix, "explicit")
        push!(results, check_performance_baseline(baseline_name, explicit_run.statistics))
    end
    if implicit_run !== nothing
        baseline_name = _performance_baseline_name(case, name_suffix, "implicit")
        push!(results, check_performance_baseline(baseline_name, implicit_run.statistics))
    end

    if explicit_run !== nothing && implicit_run !== nothing
        push!(results, _integration_classification(
            explicit_run,
            implicit_run,
            classification_name,
        ))
    else
        push!(results, VerificationResult(
            classification_name,
            "integration",
            :error,
            "one or both diagnostic integrations raised an exception",
            Dict{String, Any}(),
            Dict{String, Any}(),
            Dict{String, Any}(
                "explicit_available" => explicit_run !== nothing,
                "implicit_available" => implicit_run !== nothing,
            ),
            0.0,
        ))
    end

    if explicit_run !== nothing && explicit_run.passed
        push!(results, _run_check(benchmark_name, "benchmark") do
            check_sod_benchmark(case, explicit_run.solution, final_time)
        end)
    else
        push!(results, VerificationResult(
            benchmark_name,
            "benchmark",
            :error,
            "canonical comparison requires a successful explicit solution",
            Dict{String, Any}(),
            Dict{String, Any}(),
            Dict{String, Any}("explicit_status" => "unavailable_or_failed"),
            0.0,
        ))
    end
    return results
end
