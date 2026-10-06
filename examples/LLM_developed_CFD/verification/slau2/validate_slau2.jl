struct SLAU2VerificationTarget
    unit_cells::Int
    integration_cells::Int
    full_cells::Int
end

function SLAU2VerificationTarget(;
    unit_cells = 5,
    integration_cells = 48,
    full_cells = 96,
)
    return SLAU2VerificationTarget(
        unit_cells,
        integration_cells,
        full_cells,
    )
end

function _slau2_verification_environment(target, level)
    benchmark_cells = if level == :unit
        0
    elseif level == :integration
        target.integration_cells
    else
        target.full_cells
    end
    return Dict{String, Any}(
        "julia_version" => string(VERSION),
        "FVMFramework_version" => string(Base.pkgversion(FVMFramework)),
        "OrdinaryDiffEq_version" => string(Base.pkgversion(OrdinaryDiffEq)),
        "unit_mesh_cells" => target.unit_cells,
        "benchmark_mesh_cells" => benchmark_cells,
        "inviscid_flux" => "SLAU2 (not SLAU)",
        "conservative_state" => "Euler variables, rho*Y_k, rho*k, and rho*omega",
        "auxiliary_advection" => "signed SLAU2 mass flux with upwind primitive ratios",
        "reference" => "Kitamura and Shima, Journal of Computational Physics 245 (2013), 62-83",
        "reference_doi" => "10.1016/j.jcp.2013.02.046",
        "cross_check" => "SU2 official SLAU2 implementation, ausm_slau.cpp",
        "integration_algorithm" => "SSPRK43 with fixed acoustic CFL=0.2",
    )
end

function validate_slau2(
    target::SLAU2VerificationTarget = SLAU2VerificationTarget();
    level = :unit,
    output_directory = nothing,
    io = stdout,
    print_report = true,
)
    if !(level in (:unit, :integration, :full))
        throw(ArgumentError("SLAU2 verification level must be :unit, :integration, or :full"))
    end

    started_at = now()
    results = VerificationResult[]
    smooth_case = build_slau2_case(
        target.unit_cells;
        profile = :smooth,
        boundary_condition = :transmissive,
    )
    uniform_case = build_slau2_case(
        target.unit_cells;
        profile = :uniform,
        boundary_condition = :transmissive,
    )
    stationary_contact_case = build_slau2_case(
        target.unit_cells + 1;
        profile = :stationary_contact,
        boundary_condition = :slip_wall,
    )

    push!(results, _run_check("slau2_identical_state_flux", "component") do
        check_slau2_identical_state_flux()
    end)
    push!(results, _run_check("slau2_published_pressure_form", "component") do
        check_slau2_published_pressure_form()
    end)
    push!(results, _run_check("slau2_zero_velocity_pressure_jump", "component") do
        check_slau2_zero_velocity_pressure_jump()
    end)
    push!(results, _run_check("slau2_stationary_contact_flux", "component") do
        check_slau2_stationary_contact_flux()
    end)
    push!(results, _run_check("slau2_orientation_symmetry", "component") do
        check_slau2_orientation_symmetry()
    end)
    push!(results, _run_check("slau2_internal_face_conservation", "component") do
        check_slau2_face_conservation(smooth_case)
    end)
    push!(results, _run_check("slau2_auxiliary_upwinding", "component") do
        check_slau2_auxiliary_upwinding()
    end)
    push!(results, _run_check("slau2_free_stream_residual", "residual") do
        check_slau2_free_stream(uniform_case)
    end)
    push!(results, _run_check("slau2_stationary_auxiliary_contact", "residual") do
        check_slau2_stationary_auxiliary_contact(stationary_contact_case)
    end)

    if level in (:integration, :full)
        benchmark_cells = if level == :full
            target.full_cells
        else
            target.integration_cells
        end
        sod_case = build_slau2_case(
            benchmark_cells;
            profile = :sod,
            boundary_condition = :transmissive,
        )
        push!(results, _run_check("slau2_sod_shock_tube", "integration") do
            check_slau2_sod_integration(sod_case, 0.2)
        end)
    end

    finished_at = now()
    report = VerificationReport(
        "SLAU2 compressible flux with conservative species and SST advection",
        level,
        started_at,
        finished_at,
        _slau2_verification_environment(target, level),
        results,
    )
    if print_report
        print(io, report_string(report))
    end
    if output_directory !== nothing
        write_report(report, output_directory)
    end
    return report
end
