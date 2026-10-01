function check_stationary_contact(case)
    derivative = residual(case)
    absolute_norms, normalized_norms = _field_residual_norms(
        case,
        derivative,
        case.u0,
    )
    maximum_normalized_norm = maximum(values(normalized_norms))
    violations = state_violations(
        case.u0,
        case.system.state_axes;
        gamma = case.gamma,
        cv = case.cv,
    )
    tolerance = 1e-12
    return (
        passed = maximum_normalized_norm <= tolerance && isempty(violations),
        summary = "HLLC exactly preserves a stationary ideal-gas contact discontinuity",
        metrics = Dict{String, Any}(
            "absolute_field_norms" => absolute_norms,
            "normalized_field_norms" => normalized_norms,
            "maximum_normalized_norm" => maximum_normalized_norm,
            "spatial_method" => string(case.spatial_method),
        ),
        expected = Dict{String, Any}(
            "maximum_normalized_norm" => tolerance,
            "oracle" => "exact stationary Euler contact with equal pressure and velocity",
        ),
        diagnostics = Dict{String, Any}(
            "state_violations" => violations,
        ),
    )
end
