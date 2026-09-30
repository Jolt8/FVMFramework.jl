function _verification_environment(target, level)
    benchmark_cells = if level == :full
        target.full_cells
    elseif level == :integration
        target.integration_cells
    else
        0
    end
    return Dict{String, Any}(
        "julia_version" => string(VERSION),
        "FVMFramework_version" => string(Base.pkgversion(FVMFramework)),
        "OrdinaryDiffEq_version" => string(Base.pkgversion(OrdinaryDiffEq)),
        "ForwardDiff_version" => string(Base.pkgversion(ForwardDiff)),
        "random_seed" => target.random_seed,
        "unit_mesh_cells" => target.unit_cells,
        "benchmark_mesh_cells" => benchmark_cells,
        "spatial_method" => "first-order finite volume with HLLC",
        "explicit_algorithm" => "SSPRK43, fixed CFL=0.2",
        "implicit_algorithm" => "FBDF, reltol=1e-6, abstol=1e-8",
    )
end

function validate(
    target::CompressibleFlowVerification = CompressibleFlowVerification();
    level = :unit,
    output_directory = nothing,
    io = stdout,
    print_report = true,
)
    if !(level in (:unit, :integration, :full))
        throw(ArgumentError("verification level must be :unit, :integration, or :full"))
    end

    started_at = now()
    results = VerificationResult[]
    uniform_case = build_compressible_case(
        target.unit_cells;
        profile = :uniform,
        boundary_condition = :transmissive,
    )
    closed_case = build_compressible_case(
        target.unit_cells;
        profile = :stationary_smooth,
        boundary_condition = :slip_wall,
    )
    jacobian_case = build_compressible_case(
        target.unit_cells;
        profile = :jacobian,
        boundary_condition = :transmissive,
    )

    push!(results, _run_check("uniform_hllc_flux", "component") do
        check_uniform_hllc_flux(uniform_case)
    end)
    push!(results, _run_check("internal_face_conservation", "component") do
        check_face_conservation(jacobian_case)
    end)
    push!(results, _run_check("free_stream_residual", "residual") do
        check_free_stream(uniform_case)
    end)
    push!(results, _run_check("global_conservation", "residual") do
        check_global_conservation(closed_case)
    end)
    push!(results, _run_check("state_admissibility", "state") do
        check_state_validation(uniform_case)
    end)
    push!(results, _run_check("empirical_dependency", "jacobian") do
        check_empirical_sparsity(jacobian_case)
    end)
    push!(results, _run_check("ad_vs_finite_difference_jvp", "jacobian") do
        check_ad_jvp(jacobian_case, target.random_seed)
    end)

    if level in (:integration, :full)
        if level == :full
            benchmark_cells = target.full_cells
        else
            benchmark_cells = target.integration_cells
        end
        benchmark_case = build_compressible_case(
            benchmark_cells;
            profile = :sod,
            boundary_condition = :slip_wall,
        )
        append!(results, run_integration_checks(benchmark_case, 0.12))
    end

    report = VerificationReport(
        "CompressibleFlow(HLLC, first_order)",
        level,
        started_at,
        now(),
        _verification_environment(target, level),
        results,
    )
    if print_report
        println(io, report_string(report))
    end
    if output_directory !== nothing
        write_report(report, output_directory)
    end
    return report
end
