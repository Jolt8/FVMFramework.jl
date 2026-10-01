function _trajectory_sample_indices(n_states; maximum_samples = 64)
    if n_states <= maximum_samples
        return collect(1:n_states)
    end
    return unique(round.(Int, range(1, n_states; length = maximum_samples)))
end

function _solution_health_trajectory(case, solution)
    times = Float64[]
    minimum_densities = Float64[]
    minimum_pressures = Float64[]
    scaled_residuals = Float64[]
    for saved_index in _trajectory_sample_indices(length(solution.u))
        state_vector = solution.u[saved_index]
        state = ComponentVector(state_vector, case.system.state_axes)
        pressure_values = similar(state.density)
        for cell_id in eachindex(state.density)
            density = state.density[cell_id]
            kinetic_energy_density = 0.5 * (
                state.momentum_density_u[cell_id]^2 +
                state.momentum_density_v[cell_id]^2 +
                state.momentum_density_w[cell_id]^2
            ) / density
            pressure_values[cell_id] = (case.gamma - 1.0) * (
                state.volumetric_energy[cell_id] - kinetic_energy_density
            )
        end
        state_residual = residual(case, state_vector, solution.t[saved_index])
        state_scale = max(norm(state_vector, Inf), 1.0)
        push!(times, solution.t[saved_index])
        push!(minimum_densities, minimum(state.density))
        push!(minimum_pressures, minimum(pressure_values))
        push!(scaled_residuals, norm(state_residual, Inf) / state_scale)
    end
    return Dict{String, Any}(
        "times" => times,
        "minimum_density" => minimum_densities,
        "minimum_pressure" => minimum_pressures,
        "scaled_residual_Linf" => scaled_residuals,
    )
end

function check_long_time_known_steady(n_cells, final_time)
    case = build_compressible_case(
        n_cells;
        profile = :stationary_contact,
        boundary_condition = :slip_wall,
        species_diffusion = false,
    )
    explicit_run = run_explicit_smoke(case, final_time)
    implicit_run = run_implicit_smoke(case, final_time)
    explicit_trajectory = _solution_health_trajectory(case, explicit_run.solution)
    implicit_trajectory = _solution_health_trajectory(case, implicit_run.solution)
    explicit_change = norm(explicit_run.solution.u[end] .- case.u0, Inf)
    implicit_change = norm(implicit_run.solution.u[end] .- case.u0, Inf)
    method_difference = norm(
        explicit_run.solution.u[end] .- implicit_run.solution.u[end],
        Inf,
    )
    maximum_scaled_residual = max(
        maximum(explicit_trajectory["scaled_residual_Linf"]),
        maximum(implicit_trajectory["scaled_residual_Linf"]),
    )
    minimum_density = min(
        minimum(explicit_trajectory["minimum_density"]),
        minimum(implicit_trajectory["minimum_density"]),
    )
    minimum_pressure = min(
        minimum(explicit_trajectory["minimum_pressure"]),
        minimum(implicit_trajectory["minimum_pressure"]),
    )
    tolerance = 2e-11
    passed =
        explicit_run.passed &&
        implicit_run.passed &&
        explicit_change <= tolerance &&
        implicit_change <= tolerance &&
        method_difference <= tolerance &&
        maximum_scaled_residual <= tolerance &&
        minimum_density > 0.0 &&
        minimum_pressure > 0.0
    return (
        passed = passed,
        summary = "a nonuniform stationary contact remains an admissible long-time steady solution",
        metrics = Dict{String, Any}(
            "final_time" => final_time,
            "explicit_change_Linf" => explicit_change,
            "implicit_change_Linf" => implicit_change,
            "explicit_implicit_difference_Linf" => method_difference,
            "minimum_density" => minimum_density,
            "minimum_pressure" => minimum_pressure,
            "maximum_scaled_residual_Linf" => maximum_scaled_residual,
            "explicit_trajectory" => explicit_trajectory,
            "implicit_trajectory" => implicit_trajectory,
            "explicit_solver_statistics" => explicit_run.statistics,
            "implicit_solver_statistics" => implicit_run.statistics,
        ),
        expected = Dict{String, Any}(
            "known_steady_state_tolerance" => tolerance,
            "positive_density_and_pressure" => true,
            "explicit_and_implicit_completion" => true,
        ),
        diagnostics = Dict{String, Any}(
            "explicit_violations" => explicit_run.violations,
            "implicit_violations" => implicit_run.violations,
        ),
    )
end

function check_long_time_nontrivial_startup(n_cells, final_time)
    case = build_compressible_case(
        n_cells;
        profile = :smooth_acoustic,
        boundary_condition = :slip_wall,
        species_diffusion = false,
    )
    explicit_run = run_explicit_smoke(case, final_time)
    implicit_run = run_implicit_smoke(case, final_time)
    explicit_trajectory = _solution_health_trajectory(case, explicit_run.solution)
    implicit_trajectory = _solution_health_trajectory(case, implicit_run.solution)
    final_scale = max(norm(explicit_run.solution.u[end], Inf), 1.0)
    method_difference = norm(
        explicit_run.solution.u[end] .- implicit_run.solution.u[end],
        Inf,
    ) / final_scale
    startup_change = norm(explicit_run.solution.u[end] .- case.u0, Inf) /
        max(norm(case.u0, Inf), 1.0)
    minimum_density = min(
        minimum(explicit_trajectory["minimum_density"]),
        minimum(implicit_trajectory["minimum_density"]),
    )
    minimum_pressure = min(
        minimum(explicit_trajectory["minimum_pressure"]),
        minimum(implicit_trajectory["minimum_pressure"]),
    )
    comparison_tolerance = 2e-3
    passed =
        explicit_run.passed &&
        implicit_run.passed &&
        method_difference <= comparison_tolerance &&
        startup_change >= 1e-4 &&
        minimum_density > 0.0 &&
        minimum_pressure > 0.0
    return (
        passed = passed,
        summary = "explicit and implicit paths agree after a nontrivial smooth long-time startup",
        metrics = Dict{String, Any}(
            "final_time" => final_time,
            "relative_explicit_implicit_difference_Linf" => method_difference,
            "relative_startup_change_Linf" => startup_change,
            "minimum_density" => minimum_density,
            "minimum_pressure" => minimum_pressure,
            "explicit_trajectory" => explicit_trajectory,
            "implicit_trajectory" => implicit_trajectory,
            "explicit_solver_statistics" => explicit_run.statistics,
            "implicit_solver_statistics" => implicit_run.statistics,
        ),
        expected = Dict{String, Any}(
            "maximum_relative_method_difference" => comparison_tolerance,
            "minimum_relative_startup_change" => 1e-4,
            "positive_density_and_pressure" => true,
        ),
        diagnostics = Dict{String, Any}(
            "explicit_violations" => explicit_run.violations,
            "implicit_violations" => implicit_run.violations,
        ),
    )
end
