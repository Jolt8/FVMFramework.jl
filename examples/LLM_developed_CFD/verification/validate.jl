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
        "spatial_methods" => ["first-order HLLC", "MUSCL HLLC"],
        "low_mach_method" => "Thornber velocity reconstruction correction",
        "species_model" => "conservative rho*Y state with HLLC upwinding and mixture-corrected Fickian diffusion",
        "explicit_algorithm" => "SSPRK43, fixed CFL=0.2",
        "implicit_algorithm" => "FBDF; ForwardDiff for first order and finite differences for shock-limited MUSCL",
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
    muscl_uniform_case = build_compressible_case(
        target.unit_cells;
        profile = :uniform,
        boundary_condition = :transmissive,
        spatial_method = :muscl,
    )
    muscl_jacobian_case = build_compressible_case(
        target.unit_cells + 2;
        profile = :jacobian,
        boundary_condition = :transmissive,
        spatial_method = :muscl,
    )
    contact_case = build_compressible_case(
        target.unit_cells + 1;
        profile = :stationary_contact,
        boundary_condition = :slip_wall,
        species_diffusion = false,
    )
    muscl_contact_case = build_compressible_case(
        target.unit_cells + 1;
        profile = :stationary_contact,
        boundary_condition = :slip_wall,
        spatial_method = :muscl,
        species_diffusion = false,
    )
    thornber_uniform_case = build_compressible_case(
        target.unit_cells;
        profile = :uniform,
        boundary_condition = :transmissive,
        low_mach_method = :thornber,
    )
    muscl_thornber_uniform_case = build_compressible_case(
        target.unit_cells;
        profile = :uniform,
        boundary_condition = :transmissive,
        spatial_method = :muscl,
        low_mach_method = :thornber,
    )
    thornber_zero_velocity_case = build_compressible_case(
        target.unit_cells;
        profile = :stationary_smooth,
        boundary_condition = :slip_wall,
        low_mach_method = :thornber,
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
    push!(results, _run_check("muscl_constant_reconstruction", "muscl") do
        check_muscl_constant_reconstruction(16)
    end)
    push!(results, _run_check("muscl_linear_reconstruction", "muscl") do
        check_muscl_linear_reconstruction(32)
    end)
    convergence_grids = if level == :full
        [24, 48, 96, 192]
    else
        [16, 32, 64, 128]
    end
    push!(results, _run_check("muscl_smooth_convergence", "muscl") do
        check_muscl_smooth_convergence(convergence_grids)
    end)
    push!(results, _run_check("muscl_discontinuity", "muscl") do
        check_muscl_discontinuity(32)
    end)
    push!(results, _run_check("muscl_free_stream", "muscl") do
        check_free_stream(muscl_uniform_case)
    end)
    push!(results, _run_check("muscl_stencil_expansion", "jacobian") do
        check_muscl_stencil_expansion(target.unit_cells + 2)
    end)
    push!(results, _run_check("muscl_empirical_dependency", "jacobian") do
        check_empirical_sparsity(muscl_jacobian_case)
    end)
    push!(results, _run_check("muscl_ad_vs_finite_difference_jvp", "jacobian") do
        check_ad_jvp(muscl_jacobian_case, target.random_seed + 1)
    end)
    push!(results, _run_check("thornber_identical_states", "low_mach") do
        check_thornber_identical_states()
    end)
    push!(results, _run_check("thornber_high_mach_recovery", "low_mach") do
        check_thornber_high_mach_recovery()
    end)
    push!(results, _run_check("thornber_low_mach_jump", "low_mach") do
        check_thornber_low_mach_jump()
    end)
    push!(results, _run_check("thornber_tangential_preservation", "low_mach") do
        check_thornber_tangential_preservation()
    end)
    push!(results, _run_check("thornber_zero_velocity_hllc_jacobian", "low_mach") do
        check_thornber_zero_velocity_hllc_jacobian()
    end)
    push!(results, _run_check("thornber_energy_consistency", "low_mach") do
        check_thornber_energy_consistency()
    end)
    push!(results, _run_check("thornber_zero_velocity_ad_vs_finite_difference_jvp", "jacobian") do
        check_ad_jvp(thornber_zero_velocity_case, target.random_seed + 2)
    end)
    push!(results, _run_check("thornber_free_stream", "low_mach") do
        check_free_stream(thornber_uniform_case)
    end)
    push!(results, _run_check("muscl_thornber_free_stream", "low_mach") do
        check_free_stream(muscl_thornber_uniform_case)
    end)
    push!(results, _run_check("species_zero_diffusion", "species") do
        check_species_zero_diffusion()
    end)
    push!(results, _run_check("species_linear_diffusion", "species") do
        check_species_linear_diffusion()
    end)
    push!(results, _run_check("species_diffusive_conservation", "species") do
        check_species_diffusive_conservation()
    end)
    push!(results, _run_check("species_advection_sum_constraint", "species") do
        check_species_advection_and_sum_constraint()
    end)
    push!(results, _run_check("conservative_species_state", "species") do
        check_conservative_species_state(jacobian_case)
    end)
    push!(results, _run_check("conservative_species_face_advection", "species") do
        check_conservative_species_face_advection()
    end)
    push!(results, _run_check("mixture_corrected_species_diffusion", "species") do
        check_mixture_corrected_species_diffusion()
    end)
    push!(results, _run_check("stationary_contact", "benchmark") do
        check_stationary_contact(contact_case)
    end)
    push!(results, _run_check("stationary_contact_muscl", "benchmark") do
        check_stationary_contact(muscl_contact_case)
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
        muscl_benchmark_case = build_compressible_case(
            benchmark_cells;
            profile = :sod,
            boundary_condition = :slip_wall,
            spatial_method = :muscl,
        )
        append!(results, run_integration_checks(
            muscl_benchmark_case,
            0.12;
            name_suffix = "muscl",
        ))
        thornber_benchmark_case = build_compressible_case(
            benchmark_cells;
            profile = :sod,
            boundary_condition = :slip_wall,
            low_mach_method = :thornber,
        )
        append!(results, run_integration_checks(
            thornber_benchmark_case,
            0.12;
            name_suffix = "thornber",
        ))
        species_advection_case = build_compressible_case(
            benchmark_cells;
            profile = :species_advection,
            boundary_condition = :transmissive,
            species_diffusion = false,
        )
        push!(results, _run_check("coupled_species_advection", "integration") do
            check_coupled_species_advection_integration(
                species_advection_case,
                0.1,
            )
        end)
    end

    report = VerificationReport(
        "CompressibleFlow(HLLC, first_order + MUSCL + Thornber + species)",
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
