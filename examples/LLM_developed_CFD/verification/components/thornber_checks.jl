function _thornber_reference(
    velocity_a,
    speed_of_sound_a,
    velocity_b,
    speed_of_sound_b,
    normal,
)
    regularized_mach_a = sqrt(
        sum(abs2, velocity_a) / speed_of_sound_a^2 +
        THORNBER_MACH_REGULARIZATION^2
    )
    regularized_mach_b = sqrt(
        sum(abs2, velocity_b) / speed_of_sound_b^2 +
        THORNBER_MACH_REGULARIZATION^2
    )
    mach_coefficient = min(1.0, max(regularized_mach_a, regularized_mach_b))
    normal_velocity_a = dot(velocity_a, normal)
    normal_velocity_b = dot(velocity_b, normal)
    average_normal_velocity = 0.5 * (normal_velocity_a + normal_velocity_b)
    scaled_half_jump = 0.5 * mach_coefficient * (
        normal_velocity_a - normal_velocity_b
    )
    corrected_normal_a = average_normal_velocity + scaled_half_jump
    corrected_normal_b = average_normal_velocity - scaled_half_jump
    corrected_a = velocity_a + (corrected_normal_a - normal_velocity_a) * normal
    corrected_b = velocity_b + (corrected_normal_b - normal_velocity_b) * normal
    return corrected_a, corrected_b, mach_coefficient
end

function _production_thornber(velocity_a, sound_speed_a, velocity_b, sound_speed_b, normal)
    corrected = thornber_low_mach_correction(
        velocity_a[1],
        velocity_a[2],
        velocity_a[3],
        sound_speed_a,
        velocity_b[1],
        velocity_b[2],
        velocity_b[3],
        sound_speed_b,
        normal,
    )
    return collect(corrected[1:3]), collect(corrected[4:6])
end

function check_thornber_identical_states()
    velocity = [2.0, -0.5, 0.25]
    normal = [1.0, 0.0, 0.0]
    corrected_a, corrected_b = _production_thornber(
        velocity,
        340.0,
        velocity,
        340.0,
        normal,
    )
    maximum_error = max(
        norm(corrected_a - velocity, Inf),
        norm(corrected_b - velocity, Inf),
    )
    tolerance = 20.0 * eps(Float64)
    return (
        passed = maximum_error <= tolerance,
        summary = "Thornber leaves identical velocities unchanged",
        metrics = Dict{String, Any}("maximum_absolute_error" => maximum_error),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "oracle" => "zero left-right velocity jump",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_thornber_high_mach_recovery()
    velocity_a = [420.0, 20.0, -10.0]
    velocity_b = [360.0, -15.0, 5.0]
    normal = [1.0, 0.0, 0.0]
    corrected_a, corrected_b = _production_thornber(
        velocity_a,
        300.0,
        velocity_b,
        300.0,
        normal,
    )
    maximum_error = max(
        norm(corrected_a - velocity_a, Inf),
        norm(corrected_b - velocity_b, Inf),
    )
    tolerance = 100.0 * eps(max(norm(velocity_a, Inf), norm(velocity_b, Inf)))
    return (
        passed = maximum_error <= tolerance,
        summary = "Thornber exactly recovers ordinary reconstruction above Mach one",
        metrics = Dict{String, Any}("maximum_absolute_error" => maximum_error),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "Mach_coefficient" => 1.0,
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_thornber_low_mach_jump()
    velocity_a = [3.0, 4.0, 1.0]
    velocity_b = [-1.0, 2.0, -2.0]
    speed_of_sound_a = 340.0
    speed_of_sound_b = 330.0
    normal = [1.0, 0.0, 0.0]
    expected_a, expected_b, coefficient = _thornber_reference(
        velocity_a,
        speed_of_sound_a,
        velocity_b,
        speed_of_sound_b,
        normal,
    )
    corrected_a, corrected_b = _production_thornber(
        velocity_a,
        speed_of_sound_a,
        velocity_b,
        speed_of_sound_b,
        normal,
    )
    original_jump = dot(velocity_a - velocity_b, normal)
    corrected_jump = dot(corrected_a - corrected_b, normal)
    expected_jump = coefficient * original_jump
    maximum_error = max(
        norm(corrected_a - expected_a, Inf),
        norm(corrected_b - expected_b, Inf),
        abs(corrected_jump - expected_jump),
    )
    tolerance = 200.0 * eps(max(norm(expected_a, Inf), norm(expected_b, Inf), 1.0))
    return (
        passed = maximum_error <= tolerance && abs(corrected_jump) < abs(original_jump),
        summary = "Thornber scales the low-Mach normal jump by the independent Mach coefficient",
        metrics = Dict{String, Any}(
            "Mach_coefficient" => coefficient,
            "original_normal_jump" => original_jump,
            "corrected_normal_jump" => corrected_jump,
            "expected_normal_jump" => expected_jump,
            "maximum_absolute_error" => maximum_error,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "oracle" => "independent Thornber normal-jump decomposition",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_thornber_tangential_preservation()
    velocity_a = [2.0, 5.0, -3.0]
    velocity_b = [-1.0, -4.0, 6.0]
    normal = normalize([1.0, 2.0, -1.0])
    corrected_a, corrected_b = _production_thornber(
        velocity_a,
        330.0,
        velocity_b,
        335.0,
        normal,
    )
    tangential_a = velocity_a - dot(velocity_a, normal) * normal
    tangential_b = velocity_b - dot(velocity_b, normal) * normal
    corrected_tangential_a = corrected_a - dot(corrected_a, normal) * normal
    corrected_tangential_b = corrected_b - dot(corrected_b, normal) * normal
    maximum_error = max(
        norm(corrected_tangential_a - tangential_a, Inf),
        norm(corrected_tangential_b - tangential_b, Inf),
    )
    tolerance = 200.0 * eps(max(norm(velocity_a, Inf), norm(velocity_b, Inf)))
    return (
        passed = maximum_error <= tolerance,
        summary = "Thornber changes only the face-normal velocity components",
        metrics = Dict{String, Any}("maximum_tangential_error" => maximum_error),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "unchanged_quantities" => "both tangential velocity components",
        ),
        diagnostics = Dict{String, Any}(),
    )
end


function check_thornber_zero_velocity_hllc_jacobian()
    gamma = 1.4
    normal = [1.0, 0.0, 0.0]
    state = [
        1.0, 0.0, 0.0, 0.0, 1.0 / (gamma - 1.0),
        0.125, 0.0, 0.0, 0.0, 0.1 / (gamma - 1.0),
    ]
    flux = function (conservative_state)
        return collect(hllc_flux(
            conservative_state[1],
            conservative_state[2],
            conservative_state[3],
            conservative_state[4],
            conservative_state[5],
            gamma,
            conservative_state[6],
            conservative_state[7],
            conservative_state[8],
            conservative_state[9],
            conservative_state[10],
            gamma,
            normal,
            thornber_low_mach_correction,
        ))
    end
    jacobian = ForwardDiff.jacobian(flux, state)
    nonfinite_count = count(!isfinite, jacobian)
    return (
        passed = nonfinite_count == 0,
        summary = "Thornber HLLC has a finite ForwardDiff Jacobian at zero velocity",
        metrics = Dict{String, Any}(
            "Jacobian_rows" => size(jacobian, 1),
            "Jacobian_columns" => size(jacobian, 2),
            "nonfinite_entries" => nonfinite_count,
        ),
        expected = Dict{String, Any}("nonfinite_entries" => 0),
        diagnostics = Dict{String, Any}(),
    )
end


function check_thornber_energy_consistency()
    gamma = 1.4
    density_a = 1.2
    density_b = 0.9
    pressure_a = 101325.0
    pressure_b = 98000.0
    velocity_a = [50.0, 4.0, -2.0]
    velocity_b = [-35.0, -3.0, 1.0]
    normal = [1.0, 0.0, 0.0]
    sound_speed_a = sqrt(gamma * pressure_a / density_a)
    sound_speed_b = sqrt(gamma * pressure_b / density_b)
    energy_a = pressure_a / (gamma - 1.0) + 0.5 * density_a * sum(abs2, velocity_a)
    energy_b = pressure_b / (gamma - 1.0) + 0.5 * density_b * sum(abs2, velocity_b)

    corrected = thornber_low_mach_correction(
        velocity_a[1], velocity_a[2], velocity_a[3], sound_speed_a,
        velocity_b[1], velocity_b[2], velocity_b[3], sound_speed_b,
        normal,
    )
    corrected_energy_a = energy_with_corrected_velocity(
        energy_a,
        density_a,
        velocity_a...,
        corrected[1:3]...,
    )
    corrected_energy_b = energy_with_corrected_velocity(
        energy_b,
        density_b,
        velocity_b...,
        corrected[4:6]...,
    )
    recovered_pressure_a = primitive_from_conservative(
        density_a,
        density_a .* corrected[1:3]...,
        corrected_energy_a,
        gamma,
    )[4]
    recovered_pressure_b = primitive_from_conservative(
        density_b,
        density_b .* corrected[4:6]...,
        corrected_energy_b,
        gamma,
    )[4]
    maximum_pressure_error = max(
        abs(recovered_pressure_a - pressure_a),
        abs(recovered_pressure_b - pressure_b),
    )
    tolerance = 20.0 * eps(max(pressure_a, pressure_b))
    return (
        passed = maximum_pressure_error <= tolerance,
        summary = "Thornber velocity correction preserves reconstructed pressure",
        metrics = Dict{String, Any}(
            "maximum_pressure_error" => maximum_pressure_error,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "preserved_quantity" => "internal energy and pressure",
        ),
        diagnostics = Dict{String, Any}(),
    )
end
