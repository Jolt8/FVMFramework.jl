function _historical_componentwise_gradient_mutation(
    value_a,
    value_b,
    gradient_a,
    gradient_b,
    normal,
    distance,
)
    average_gradient = 0.5 .* (gradient_a .+ gradient_b)
    normal_derivative = (value_b - value_a) / distance
    return [
        average_gradient[dimension] +
        (normal_derivative - average_gradient[dimension] * normal[dimension]) *
        normal[dimension]
        for dimension in 1:3
    ]
end

function check_gradient_mutation_detection()
    normal = normalize([0.5, -0.4, 0.7])
    distance = 0.31
    value_a = 0.2
    normal_derivative = -0.9
    value_b = value_a + distance * normal_derivative
    gradient_a = [0.4, 0.8, -0.3]
    gradient_b = [-0.2, 0.1, 0.9]
    mutated = _historical_componentwise_gradient_mutation(
        value_a,
        value_b,
        gradient_a,
        gradient_b,
        normal,
        distance,
    )
    mutated_normal_error = abs(dot(mutated, normal) - normal_derivative)
    direct_test_tolerance = 100.0 * eps(Float64)
    mutation_detected = mutated_normal_error > direct_test_tolerance
    return (
        passed = mutation_detected,
        summary = "the direct gradient oracle rejects the historical componentwise-projection mutation",
        metrics = Dict{String, Any}(
            "mutated_normal_derivative_error" => mutated_normal_error,
        ),
        expected = Dict{String, Any}(
            "rejection_threshold" => direct_test_tolerance,
            "mutation" => "replace the full dot product with independent scalar products",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_mms_mutation_detection()
    n_cells = 64
    correct_viscous = _viscous_compressible_mms_residuals(n_cells)
    reversed_viscous = _viscous_compressible_mms_residuals(
        n_cells;
        viscous_flux_multiplier = -1.0,
    )
    correct_species = _species_mms_residuals(n_cells)
    reversed_species = _species_mms_residuals(
        n_cells;
        diffusion_flux_multiplier = -1.0,
    )
    interior = 2:(n_cells - 1)
    correct_viscous_error = _mms_l2(correct_viscous.momentum[interior])
    reversed_viscous_error = _mms_l2(reversed_viscous.momentum[interior])
    correct_species_error = _mms_l2(correct_species.species_a[interior])
    reversed_species_error = _mms_l2(reversed_species.species_a[interior])
    viscous_amplification = reversed_viscous_error / correct_viscous_error
    species_amplification = reversed_species_error / correct_species_error
    minimum_amplification = 100.0
    return (
        passed =
            viscous_amplification >= minimum_amplification &&
            species_amplification >= minimum_amplification,
        summary = "MMS rejects reversed viscous and species-diffusion flux signs",
        metrics = Dict{String, Any}(
            "correct_viscous_error" => correct_viscous_error,
            "mutated_viscous_error" => reversed_viscous_error,
            "viscous_error_amplification" => viscous_amplification,
            "correct_species_error" => correct_species_error,
            "mutated_species_error" => reversed_species_error,
            "species_error_amplification" => species_amplification,
        ),
        expected = Dict{String, Any}(
            "minimum_error_amplification" => minimum_amplification,
            "mutations" => ["reverse viscous flux sign", "reverse species-diffusion flux sign"],
        ),
        diagnostics = Dict{String, Any}(),
    )
end
