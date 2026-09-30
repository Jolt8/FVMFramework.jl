const _CONSERVATIVE_VARIABLES = (
    :density,
    :momentum_density_u,
    :momentum_density_v,
    :momentum_density_w,
    :volumetric_energy,
)

function check_uniform_hllc_flux(case)
    state = ComponentVector(case.u0, case.system.state_axes)
    cell_id = 1
    normal = Ferrite.Vec{3}((1.0, 0.0, 0.0))
    gamma = case.gamma
    expected = physical_flux(
        state.density[cell_id],
        state.momentum_density_u[cell_id],
        state.momentum_density_v[cell_id],
        state.momentum_density_w[cell_id],
        state.volumetric_energy[cell_id],
        (gamma - 1.0) * (
            state.volumetric_energy[cell_id] -
            0.5 * (
                state.momentum_density_u[cell_id]^2 +
                state.momentum_density_v[cell_id]^2 +
                state.momentum_density_w[cell_id]^2
            ) / state.density[cell_id]
        ),
        normal,
    )
    measured = hllc_flux(
        state.density[cell_id],
        state.momentum_density_u[cell_id],
        state.momentum_density_v[cell_id],
        state.momentum_density_w[cell_id],
        state.volumetric_energy[cell_id],
        gamma,
        state.density[cell_id],
        state.momentum_density_u[cell_id],
        state.momentum_density_v[cell_id],
        state.momentum_density_w[cell_id],
        state.volumetric_energy[cell_id],
        gamma,
        normal,
    )

    absolute_error = maximum(abs.(collect(measured) .- collect(expected)))
    flux_scale = max(maximum(abs.(collect(expected))), 1.0)
    relative_error = absolute_error / flux_scale
    tolerance = 200.0 * eps(Float64)
    return (
        passed = relative_error <= tolerance,
        summary = "identical-state HLLC equals the analytical Euler flux",
        metrics = Dict{String, Any}(
            "absolute_error" => absolute_error,
            "relative_error" => relative_error,
        ),
        expected = Dict{String, Any}(
            "relative_tolerance" => tolerance,
            "oracle" => "independently evaluated physical Euler flux",
        ),
        diagnostics = Dict{String, Any}(
            "expected_flux" => collect(expected),
            "measured_flux" => collect(measured),
        ),
    )
end

function check_face_conservation(case)
    state = copy(case.u0)
    derivative = zeros(length(state))
    du, u = unpack_fvm_state(derivative, state, case.p, 0.0, case.system)
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
    for variable in _CONSERVATIVE_VARIABLES
        flow_variable = Symbol(string(variable), "_flow")
        flow = getproperty(du, flow_variable)
        imbalance = flow[idx_a] + flow[idx_b]
        scale = max(abs(flow[idx_a]), abs(flow[idx_b]), 1.0)
        normalized_imbalance = abs(imbalance) / scale
        imbalances[string(variable)] = imbalance
        push!(normalized_imbalances, normalized_imbalance)
    end

    maximum_normalized_imbalance = maximum(normalized_imbalances)
    tolerance = 20.0 * eps(Float64)
    return (
        passed = maximum_normalized_imbalance <= tolerance,
        summary = "one HLLC face contributes equal and opposite integrated flux",
        metrics = Dict{String, Any}(
            "maximum_normalized_imbalance" => maximum_normalized_imbalance,
            "face_cells" => [idx_a, idx_b],
        ),
        expected = Dict{String, Any}(
            "normalized_tolerance" => tolerance,
            "identity" => "flow_into_cell_a + flow_into_cell_b = 0",
        ),
        diagnostics = Dict{String, Any}("imbalances" => imbalances),
    )
end

function _field_residual_norms(case, derivative, state)
    derivative_state = ComponentVector(derivative, case.system.state_axes)
    conservative_state = ComponentVector(state, case.system.state_axes)
    minimum_spacing = minimum(case.geo.cell_volumes)^(1.0 / 3.0)
    speed_scale = 500.0
    norms = Dict{String, Any}()
    normalized_norms = Dict{String, Any}()

    for variable in _CONSERVATIVE_VARIABLES
        derivative_values = getproperty(derivative_state, variable)
        state_values = getproperty(conservative_state, variable)
        absolute_norm = norm(derivative_values, Inf)
        characteristic_rate = max(norm(state_values, Inf) * speed_scale / minimum_spacing, 1.0)
        norms[string(variable)] = absolute_norm
        normalized_norms[string(variable)] = absolute_norm / characteristic_rate
    end
    return norms, normalized_norms
end

function check_free_stream(case)
    derivative = residual(case)
    violations = state_violations(
        case.u0,
        case.system.state_axes;
        gamma = case.gamma,
        cv = case.cv,
    )
    norms, normalized_norms = _field_residual_norms(case, derivative, case.u0)
    maximum_normalized_norm = maximum(values(normalized_norms))
    tolerance = 5e-13

    return (
        passed = maximum_normalized_norm <= tolerance && isempty(violations),
        summary = "uniform three-dimensional flow has a zero semi-discrete residual",
        metrics = Dict{String, Any}(
            "absolute_field_norms" => norms,
            "normalized_field_norms" => normalized_norms,
            "maximum_normalized_norm" => maximum_normalized_norm,
        ),
        expected = Dict{String, Any}(
            "maximum_normalized_norm" => tolerance,
            "boundary_condition" => "transmissive physical flux",
        ),
        diagnostics = Dict{String, Any}(
            "state_violations" => violations,
        ),
    )
end

function check_global_conservation(case)
    derivative = residual(case)
    violations = state_violations(
        case.u0,
        case.system.state_axes;
        gamma = case.gamma,
        cv = case.cv,
    )
    derivative_state = ComponentVector(derivative, case.system.state_axes)
    conservative_state = ComponentVector(case.u0, case.system.state_axes)
    volumes = case.geo.cell_volumes
    net_rates = Dict{String, Any}()
    normalized_rates = Dict{String, Any}()
    normalization_scales = Dict{String, Any}()

    maximum_density = norm(conservative_state.density, Inf)
    maximum_velocity = 0.0
    maximum_pressure = 0.0
    maximum_energy = norm(conservative_state.volumetric_energy, Inf)
    maximum_sound_speed = 0.0
    for cell_id in eachindex(conservative_state.density)
        density = conservative_state.density[cell_id]
        velocity = sqrt(
            conservative_state.momentum_density_u[cell_id]^2 +
            conservative_state.momentum_density_v[cell_id]^2 +
            conservative_state.momentum_density_w[cell_id]^2
        ) / density
        kinetic_energy_density = 0.5 * density * velocity^2
        pressure = (case.gamma - 1.0) * (
            conservative_state.volumetric_energy[cell_id] - kinetic_energy_density
        )
        sound_speed = sqrt(case.gamma * pressure / density)
        maximum_velocity = max(maximum_velocity, velocity)
        maximum_pressure = max(maximum_pressure, pressure)
        maximum_sound_speed = max(maximum_sound_speed, sound_speed)
    end
    signal_speed = maximum_velocity + maximum_sound_speed
    characteristic_scales = Dict(
        :density => maximum_density * signal_speed,
        :momentum_density_u => maximum_density * maximum_velocity * signal_speed + maximum_pressure,
        :momentum_density_v => maximum_density * maximum_velocity * signal_speed + maximum_pressure,
        :momentum_density_w => maximum_density * maximum_velocity * signal_speed + maximum_pressure,
        :volumetric_energy => (maximum_energy + maximum_pressure) * signal_speed,
    )

    for variable in _CONSERVATIVE_VARIABLES
        contributions = volumes .* getproperty(derivative_state, variable)
        net_rate = sum(contributions)
        total_activity = sum(abs, contributions)
        normalization_scale = max(total_activity, characteristic_scales[variable], 1.0)
        normalized_rate = abs(net_rate) / normalization_scale
        net_rates[string(variable)] = net_rate
        normalized_rates[string(variable)] = normalized_rate
        normalization_scales[string(variable)] = normalization_scale
    end

    maximum_normalized_rate = maximum(values(normalized_rates))
    tolerance = 1e-12
    return (
        passed = maximum_normalized_rate <= tolerance && isempty(violations),
        summary = "volume-integrated closed-domain residual conserves all Euler quantities",
        metrics = Dict{String, Any}(
            "net_rates" => net_rates,
            "normalized_rates" => normalized_rates,
            "normalization_scales" => normalization_scales,
            "maximum_normalized_rate" => maximum_normalized_rate,
        ),
        expected = Dict{String, Any}(
            "normalized_tolerance" => tolerance,
            "oracle" => "pairwise internal cancellation plus balanced constant wall pressure",
        ),
        diagnostics = Dict{String, Any}(
            "state_violations" => violations,
        ),
    )
end

function check_state_validation(case)
    valid_violations = state_violations(
        case.u0,
        case.system.state_axes;
        gamma = case.gamma,
        cv = case.cv,
    )

    nonfinite_state = copy(case.u0)
    nonfinite_named = ComponentVector(nonfinite_state, case.system.state_axes)
    nonfinite_named.momentum_density_v[2] = NaN
    nonfinite_violations = state_violations(
        nonfinite_state,
        case.system.state_axes;
        gamma = case.gamma,
        cv = case.cv,
        time = 0.25,
    )

    nonphysical_state = copy(case.u0)
    nonphysical_named = ComponentVector(nonphysical_state, case.system.state_axes)
    nonphysical_named.density[1] = -1.0
    nonphysical_named.volumetric_energy[3] = 0.0
    nonphysical_violations = state_violations(
        nonphysical_state,
        case.system.state_axes;
        gamma = case.gamma,
        cv = case.cv,
        time = 0.5,
    )

    detected_nonfinite = any(
        violation -> violation["variable"] == "momentum_density_v" && violation["cell_index"] == 2,
        nonfinite_violations,
    )
    detected_density = any(
        violation -> violation["variable"] == "density" && violation["cell_index"] == 1,
        nonphysical_violations,
    )
    detected_thermodynamics = any(
        violation -> violation["variable"] in ("pressure", "temperature") && violation["cell_index"] == 3,
        nonphysical_violations,
    )
    passed = isempty(valid_violations) && detected_nonfinite && detected_density && detected_thermodynamics

    return (
        passed = passed,
        summary = "state checker accepts valid flow and localizes non-finite/nonphysical states",
        metrics = Dict{String, Any}(
            "valid_state_violation_count" => length(valid_violations),
            "nonfinite_detection_count" => length(nonfinite_violations),
            "nonphysical_detection_count" => length(nonphysical_violations),
        ),
        expected = Dict{String, Any}(
            "valid_state_violation_count" => 0,
            "required_detections" => ["NaN", "negative density", "nonpositive pressure/temperature"],
        ),
        diagnostics = Dict{String, Any}(
            "valid_state_violations" => valid_violations,
            "nonfinite_state_violations" => nonfinite_violations,
            "nonphysical_state_violations" => nonphysical_violations,
        ),
    )
end
