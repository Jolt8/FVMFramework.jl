const _SLAU2_CONSERVATIVE_VARIABLES = (
    :density,
    :momentum_density_u,
    :momentum_density_v,
    :momentum_density_w,
    :volumetric_energy,
    :turbulent_kinetic_energy_density,
    :specific_dissipation_rate_density,
)

function _slau2_relative_error(measured, expected)
    absolute_error = maximum(abs.(collect(measured) .- collect(expected)))
    scale = max(maximum(abs.(collect(expected))), 1.0)
    return absolute_error, absolute_error / scale
end

function check_slau2_identical_state_flux()
    gamma = 1.4
    density = 1.18
    vel_u = 120.0
    vel_v = 11.0
    vel_w = -7.0
    pressure = 101325.0
    normal = Ferrite.Vec{3}((0.6, 0.8, 0.0))
    state = _slau2_conservative_values(
        density,
        vel_u,
        vel_v,
        vel_w,
        pressure,
        gamma,
    )
    expected = _slau2_boundary_physical_flux(
        state...,
        pressure,
        normal,
    )
    measured = slau2_flux(
        state...,
        gamma,
        state...,
        gamma,
        normal,
    )
    absolute_error, relative_error = _slau2_relative_error(measured, expected)
    tolerance = 500.0 * eps(Float64)
    return (
        passed = relative_error <= tolerance,
        summary = "identical-state SLAU2 equals the analytical Euler flux",
        metrics = Dict{String, Any}(
            "absolute_error" => absolute_error,
            "relative_error" => relative_error,
        ),
        expected = Dict{String, Any}(
            "relative_tolerance" => tolerance,
            "oracle" => "independently evaluated physical Euler flux",
        ),
        diagnostics = Dict{String, Any}(
            "measured_flux" => collect(measured),
            "expected_flux" => collect(expected),
        ),
    )
end

function check_slau2_published_pressure_form()
    density_a = 1.1
    density_b = 0.9
    vel_u_a, vel_v_a, vel_w_a = 80.0, 15.0, -4.0
    vel_u_b, vel_v_b, vel_w_b = -35.0, -8.0, 6.0
    pressure_a = 95000.0
    pressure_b = 102000.0
    gamma = 1.4
    normal = Ferrite.Vec{3}((0.6, 0.8, 0.0))
    speed_of_sound_a = sqrt(gamma * pressure_a / density_a)
    speed_of_sound_b = sqrt(gamma * pressure_b / density_b)

    measured_mass_flux, measured_pressure_flux = _slau2_mass_and_pressure_flux(
        density_a,
        vel_u_a,
        vel_v_a,
        vel_w_a,
        pressure_a,
        speed_of_sound_a,
        density_b,
        vel_u_b,
        vel_v_b,
        vel_w_b,
        pressure_b,
        speed_of_sound_b,
        normal,
    )

    # Frozen values independently evaluated from the published SLAU2 mass and
    # pressure formulation. The original SLAU pressure for the same state is included
    # to make accidental implementation of SLAU fail conspicuously.
    expected_mass_flux = 18.7265549202487
    expected_pressure_flux = 102383.988620211
    original_slau_pressure_flux = 103600.574778044
    mass_error = abs(measured_mass_flux - expected_mass_flux)
    pressure_error = abs(measured_pressure_flux - expected_pressure_flux)
    variant_separation = abs(measured_pressure_flux - original_slau_pressure_flux)
    mass_tolerance = 2e-12 * abs(expected_mass_flux)
    pressure_tolerance = 2e-12 * abs(expected_pressure_flux)
    passed = mass_error <= mass_tolerance &&
        pressure_error <= pressure_tolerance &&
        variant_separation >= 1000.0
    return (
        passed = passed,
        summary = "mass split and velocity-scaled pressure dissipation match SLAU2, not SLAU",
        metrics = Dict{String, Any}(
            "mass_flux" => measured_mass_flux,
            "pressure_flux" => measured_pressure_flux,
            "mass_error" => mass_error,
            "pressure_error" => pressure_error,
            "distance_from_original_SLAU_pressure" => variant_separation,
        ),
        expected = Dict{String, Any}(
            "mass_flux" => expected_mass_flux,
            "pressure_flux" => expected_pressure_flux,
            "original_SLAU_pressure_flux" => original_slau_pressure_flux,
            "source" => "Kitamura and Shima, JCP 245 (2013), SLAU2 mass and pressure flux equations",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_slau2_zero_velocity_pressure_jump()
    gamma = 1.4
    density_a, pressure_a = 1.0, 1.0
    density_b, pressure_b = 0.125, 0.1
    speed_of_sound_a = sqrt(gamma * pressure_a / density_a)
    speed_of_sound_b = sqrt(gamma * pressure_b / density_b)
    interface_speed_of_sound = 0.5 * (speed_of_sound_a + speed_of_sound_b)
    normal = Ferrite.Vec{3}((1.0, 0.0, 0.0))
    measured_mass_flux, measured_pressure_flux = _slau2_mass_and_pressure_flux(
        density_a,
        0.0,
        0.0,
        0.0,
        pressure_a,
        speed_of_sound_a,
        density_b,
        0.0,
        0.0,
        0.0,
        pressure_b,
        speed_of_sound_b,
        normal,
    )
    expected_mass_flux = -0.5 * (pressure_b - pressure_a) /
        interface_speed_of_sound
    expected_pressure_flux = 0.5 * (pressure_a + pressure_b)
    error = max(
        abs(measured_mass_flux - expected_mass_flux),
        abs(measured_pressure_flux - expected_pressure_flux),
    )
    tolerance = 100.0 * eps(Float64)
    return (
        passed = error <= tolerance,
        summary = "zero-speed SLAU2 limit retains pressure diffusion and centered pressure",
        metrics = Dict{String, Any}(
            "maximum_absolute_error" => error,
            "mass_flux" => measured_mass_flux,
            "pressure_flux" => measured_pressure_flux,
        ),
        expected = Dict{String, Any}(
            "mass_flux" => expected_mass_flux,
            "pressure_flux" => expected_pressure_flux,
            "absolute_tolerance" => tolerance,
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_slau2_stationary_contact_flux()
    gamma = 1.4
    pressure = 2.5
    normal = Ferrite.Vec{3}((1.0, 0.0, 0.0))
    state_a = _slau2_conservative_values(
        1.0,
        0.0,
        4.0,
        -2.0,
        pressure,
        gamma,
    )
    state_b = _slau2_conservative_values(
        0.2,
        0.0,
        -3.0,
        1.0,
        pressure,
        gamma,
    )
    measured = slau2_flux(
        state_a...,
        gamma,
        state_b...,
        gamma,
        normal,
    )
    expected = (0.0, pressure, 0.0, 0.0, 0.0)
    absolute_error, relative_error = _slau2_relative_error(measured, expected)
    tolerance = 500.0 * eps(Float64)
    return (
        passed = relative_error <= tolerance,
        summary = "SLAU2 preserves a stationary contact and tangential-velocity jump",
        metrics = Dict{String, Any}(
            "absolute_error" => absolute_error,
            "relative_error" => relative_error,
        ),
        expected = Dict{String, Any}(
            "flux" => collect(expected),
            "relative_tolerance" => tolerance,
        ),
        diagnostics = Dict{String, Any}(
            "measured_flux" => collect(measured),
        ),
    )
end

function check_slau2_orientation_symmetry()
    gamma = 1.4
    normal = Ferrite.Vec{3}((0.6, 0.8, 0.0))
    reverse_normal = Ferrite.Vec{3}((-0.6, -0.8, 0.0))
    state_a = _slau2_conservative_values(
        1.1,
        80.0,
        15.0,
        -4.0,
        95000.0,
        gamma,
    )
    state_b = _slau2_conservative_values(
        0.9,
        -35.0,
        -8.0,
        6.0,
        102000.0,
        gamma,
    )
    forward_flux = slau2_flux(
        state_a...,
        gamma,
        state_b...,
        gamma,
        normal,
    )
    reverse_flux = slau2_flux(
        state_b...,
        gamma,
        state_a...,
        gamma,
        reverse_normal,
    )
    imbalance = collect(forward_flux) .+ collect(reverse_flux)
    scale = max(maximum(abs, forward_flux), 1.0)
    normalized_error = maximum(abs, imbalance) / scale
    tolerance = 500.0 * eps(Float64)
    return (
        passed = normalized_error <= tolerance,
        summary = "swapping states and reversing the normal negates the SLAU2 flux",
        metrics = Dict{String, Any}(
            "normalized_error" => normalized_error,
        ),
        expected = Dict{String, Any}(
            "normalized_tolerance" => tolerance,
            "identity" => "F(Ua, Ub, n) = -F(Ub, Ua, -n)",
        ),
        diagnostics = Dict{String, Any}(
            "component_imbalances" => imbalance,
        ),
    )
end

function check_slau2_face_conservation(case)
    derivative = zeros(length(case.u0))
    du, u = unpack_fvm_state(
        derivative,
        case.u0,
        case.p,
        0.0,
        case.system,
    )
    update_region_groups!(du, u, case.p, 0.0, case.system, case.geo)
    connection_group = first(case.system.connection_groups)
    idx_a, neighbor_list = first(connection_group.cell_neighbors)
    idx_b, face_a, face_b = first(neighbor_list)
    connection_group.flux_function!(
        du, u, case.p, 0.0, case.system, case.geo,
        idx_a, face_a,
        idx_b, face_b,
    )

    imbalances = Dict{String, Any}()
    normalized_imbalances = Float64[]
    for variable in _SLAU2_CONSERVATIVE_VARIABLES
        flow = getproperty(du, Symbol(string(variable), "_flow"))
        imbalance = flow[idx_a] + flow[idx_b]
        scale = max(abs(flow[idx_a]), abs(flow[idx_b]), 1.0)
        imbalances[string(variable)] = imbalance
        push!(normalized_imbalances, abs(imbalance) / scale)
    end
    for species_name in propertynames(du.species_density_flow)
        flow = getproperty(du.species_density_flow, species_name)
        imbalance = flow[idx_a] + flow[idx_b]
        scale = max(abs(flow[idx_a]), abs(flow[idx_b]), 1.0)
        imbalances["species_densities.$species_name"] = imbalance
        push!(normalized_imbalances, abs(imbalance) / scale)
    end

    integrated_density_flux = -du.density_flow[idx_a]
    if integrated_density_flux >= 0.0
        upwind_cell = idx_a
    else
        upwind_cell = idx_b
    end
    coupling_errors = Float64[]
    species_flux_sum = 0.0
    for species_name in propertynames(u.mass_fractions)
        species_flow = getproperty(du.species_density_flow, species_name)
        mass_fraction = getproperty(u.mass_fractions, species_name)[upwind_cell]
        expected_flow_a = -integrated_density_flux * mass_fraction
        push!(coupling_errors, abs(species_flow[idx_a] - expected_flow_a))
        species_flux_sum += species_flow[idx_a]
    end
    push!(coupling_errors, abs(species_flux_sum + integrated_density_flux))
    expected_k_flow_a = -integrated_density_flux *
        u.turbulent_kinetic_energy[upwind_cell]
    expected_omega_flow_a = -integrated_density_flux *
        u.specific_dissipation_rate[upwind_cell]
    push!(coupling_errors, abs(
        du.turbulent_kinetic_energy_density_flow[idx_a] - expected_k_flow_a,
    ))
    push!(coupling_errors, abs(
        du.specific_dissipation_rate_density_flow[idx_a] - expected_omega_flow_a,
    ))
    maximum_imbalance = maximum(normalized_imbalances)
    coupling_scale = max(abs(integrated_density_flux), 1.0)
    maximum_coupling_error = maximum(coupling_errors) / coupling_scale
    tolerance = 20.0 * eps(Float64)
    return (
        passed = maximum_imbalance <= tolerance &&
            maximum_coupling_error <= 100.0 * eps(Float64),
        summary = "one SLAU2 face conserves and mass-flux-couples Euler, species, and SST states",
        metrics = Dict{String, Any}(
            "maximum_normalized_imbalance" => maximum_imbalance,
            "maximum_normalized_coupling_error" => maximum_coupling_error,
            "face_cells" => [idx_a, idx_b],
        ),
        expected = Dict{String, Any}(
            "normalized_tolerance" => tolerance,
            "identity" => "flow_into_cell_a + flow_into_cell_b = 0",
        ),
        diagnostics = Dict{String, Any}(
            "imbalances" => imbalances,
        ),
    )
end

function check_slau2_auxiliary_upwinding()
    state = ComponentVector(
        mass_fractions = ComponentVector(
            species_a = [0.2, 0.7],
            species_b = [0.8, 0.3],
        ),
        turbulent_kinetic_energy = [0.4, 1.1],
        specific_dissipation_rate = [3.0, 7.0],
    )
    area = 0.6
    errors = Float64[]
    conservation_errors = Float64[]

    for density_flux in (2.0, -1.5)
        derivative = ComponentVector(
            species_density_flow = ComponentVector(
                species_a = zeros(2),
                species_b = zeros(2),
            ),
            turbulent_kinetic_energy_density_flow = zeros(2),
            specific_dissipation_rate_density_flow = zeros(2),
        )
        add_species_advection_flux!(
            derivative,
            state,
            1,
            2,
            area,
            density_flux,
        )
        add_sst_advection_flux!(
            derivative,
            state,
            1,
            2,
            area,
            density_flux,
        )
        if density_flux >= 0.0
            upwind_cell = 1
        else
            upwind_cell = 2
        end
        integrated_density_flux = area * density_flux
        for species_name in propertynames(state.mass_fractions)
            flow = getproperty(derivative.species_density_flow, species_name)
            expected_flow_a = -integrated_density_flux *
                getproperty(state.mass_fractions, species_name)[upwind_cell]
            push!(errors, abs(flow[1] - expected_flow_a))
            push!(conservation_errors, abs(sum(flow)))
        end
        expected_k_flow_a = -integrated_density_flux *
            state.turbulent_kinetic_energy[upwind_cell]
        expected_omega_flow_a = -integrated_density_flux *
            state.specific_dissipation_rate[upwind_cell]
        push!(errors, abs(
            derivative.turbulent_kinetic_energy_density_flow[1] -
            expected_k_flow_a,
        ))
        push!(errors, abs(
            derivative.specific_dissipation_rate_density_flow[1] -
            expected_omega_flow_a,
        ))
        push!(conservation_errors, abs(sum(
            derivative.turbulent_kinetic_energy_density_flow,
        )))
        push!(conservation_errors, abs(sum(
            derivative.specific_dissipation_rate_density_flow,
        )))
    end

    maximum_error = maximum(errors)
    maximum_conservation_error = maximum(conservation_errors)
    tolerance = 100.0 * eps(Float64)
    return (
        passed = maximum_error <= tolerance &&
            maximum_conservation_error <= tolerance,
        summary = "positive and negative SLAU2 mass fluxes upwind rho*Y, rho*k, and rho*omega",
        metrics = Dict{String, Any}(
            "maximum_upwind_error" => maximum_error,
            "maximum_conservation_error" => maximum_conservation_error,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "transport_rule" => "F_auxiliary = F_density * primitive_upwind",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function _check_slau2_zero_residual(case, summary)
    derivative = slau2_residual(case)
    derivative_state = ComponentVector(derivative, case.system.state_axes)
    state = ComponentVector(case.u0, case.system.state_axes)
    normalized_norms = Dict{String, Any}()
    for variable in _SLAU2_CONSERVATIVE_VARIABLES
        derivative_values = getproperty(derivative_state, variable)
        state_values = getproperty(state, variable)
        scale = max(norm(state_values, Inf), 1.0)
        normalized_norms[string(variable)] = norm(derivative_values, Inf) / scale
    end
    for species_name in propertynames(state.species_densities)
        derivative_values = getproperty(
            derivative_state.species_densities,
            species_name,
        )
        state_values = getproperty(state.species_densities, species_name)
        scale = max(norm(state_values, Inf), 1.0)
        normalized_norms["species_densities.$species_name"] =
            norm(derivative_values, Inf) / scale
    end
    maximum_norm = maximum(values(normalized_norms))
    violations = state_violations(
        case.u0,
        case.system.state_axes;
        gamma = case.gamma,
        cv = case.cv,
    )
    tolerance = 5e-12
    return (
        passed = maximum_norm <= tolerance && isempty(violations),
        summary = summary,
        metrics = Dict{String, Any}(
            "normalized_field_norms" => normalized_norms,
            "maximum_normalized_norm" => maximum_norm,
        ),
        expected = Dict{String, Any}(
            "maximum_normalized_norm" => tolerance,
        ),
        diagnostics = Dict{String, Any}(
            "state_violations" => violations,
        ),
    )
end

function check_slau2_free_stream(case)
    return _check_slau2_zero_residual(
        case,
        "uniform flow with species and SST variables has zero SLAU2 FVM residual",
    )
end

function check_slau2_stationary_auxiliary_contact(case)
    return _check_slau2_zero_residual(
        case,
        "stationary density, composition, k, and omega contacts have zero SLAU2 residual",
    )
end

function run_slau2_sod(case, final_time)
    state = ComponentVector(case.u0, case.system.state_axes)
    maximum_signal_speed = 0.0
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
        speed_of_sound = sqrt(case.gamma * pressure / density)
        maximum_signal_speed = max(
            maximum_signal_speed,
            abs(vel_u) + speed_of_sound,
        )
    end
    stable_timestep = 0.2 * minimum(case.geo.cell_volumes) /
        maximum_signal_speed
    problem = ODEProblem(case.rhs!, case.u0, (0.0, final_time), case.p)
    timed_result = @timed solve(
        problem,
        SSPRK43();
        adaptive = false,
        dt = min(stable_timestep, final_time),
        save_everystep = true,
        maxiters = 100000,
    )
    solution = timed_result.value
    violations = Dict{String, Any}[]
    for (saved_index, saved_state) in enumerate(solution.u)
        append!(violations, state_violations(
            saved_state,
            case.system.state_axes;
            gamma = case.gamma,
            cv = case.cv,
            time = solution.t[saved_index],
            maximum_violations = max(20 - length(violations), 0),
        ))
        if length(violations) >= 20
            break
        end
    end
    successful = OrdinaryDiffEq.SciMLBase.successful_retcode(solution)
    reached_final_time = isapprox(
        solution.t[end],
        final_time;
        atol = 100.0 * eps(final_time),
        rtol = 1e-12,
    )
    return (
        solution = solution,
        passed = successful && reached_final_time && isempty(violations),
        stable_timestep = stable_timestep,
        violations = violations,
        runtime_seconds = timed_result.time,
    )
end

function check_slau2_sod_integration(case, final_time)
    run = run_slau2_sod(case, final_time)
    if !run.passed
        return (
            passed = false,
            summary = "SLAU2 SSPRK43 Sod integration did not complete admissibly",
            metrics = Dict{String, Any}(
                "return_code" => string(run.solution.retcode),
                "final_time" => run.solution.t[end],
                "CFL_timestep" => run.stable_timestep,
            ),
            expected = Dict{String, Any}(
                "return_code" => "Success",
                "state_violations" => 0,
            ),
            diagnostics = Dict{String, Any}(
                "state_violations" => run.violations,
            ),
        )
    end

    benchmark = check_sod_benchmark(case, run.solution, final_time)
    final_state = ComponentVector(run.solution.u[end], case.system.state_axes)
    species_sum_error = norm(
        final_state.species_densities.species_a .+
        final_state.species_densities.species_b .-
        final_state.density,
        Inf,
    )
    species_fraction_error = norm(
        final_state.species_densities.species_a ./ final_state.density .- 0.4,
        Inf,
    )
    turbulent_kinetic_energy_error = norm(
        final_state.turbulent_kinetic_energy_density ./ final_state.density .-
        0.2,
        Inf,
    )
    specific_dissipation_rate_error = norm(
        final_state.specific_dissipation_rate_density ./ final_state.density .-
        3.0,
        Inf,
    )
    auxiliary_tolerance = 2e-10
    metrics = copy(benchmark.metrics)
    metrics["return_code"] = string(run.solution.retcode)
    metrics["CFL_timestep"] = run.stable_timestep
    metrics["runtime_seconds"] = run.runtime_seconds
    metrics["species_density_sum_error_Linf"] = species_sum_error
    metrics["species_fraction_error_Linf"] = species_fraction_error
    metrics["turbulent_kinetic_energy_error_Linf"] =
        turbulent_kinetic_energy_error
    metrics["specific_dissipation_rate_error_Linf"] =
        specific_dissipation_rate_error
    expected = copy(benchmark.expected)
    expected["auxiliary_transport_tolerance"] = auxiliary_tolerance
    return (
        passed = benchmark.passed &&
            species_sum_error <= auxiliary_tolerance &&
            species_fraction_error <= auxiliary_tolerance &&
            turbulent_kinetic_energy_error <= auxiliary_tolerance &&
            specific_dissipation_rate_error <= auxiliary_tolerance,
        summary = "SLAU2 Sod preserves conservative species/SST ratios and matches the exact flow solution",
        metrics = metrics,
        expected = expected,
        diagnostics = benchmark.diagnostics,
    )
end
