include(joinpath(
    @__DIR__,
    "..", "..",
    "examples", "LLM_developed_CFD", "verification", "FVMVerification.jl",
))
using .FVMVerification

@testset "CFD verification framework" begin
    report = validate(; level = :unit, print_report = false)
    @test !has_failures(report)
    result_names = Set(result.name for result in report.results)
    @test "muscl_smooth_convergence" in result_names
    @test "thornber_low_mach_jump" in result_names
    @test "thornber_zero_velocity_hllc_jacobian" in result_names
    @test "thornber_energy_consistency" in result_names
    @test "thornber_zero_velocity_ad_vs_finite_difference_jvp" in result_names
    @test "species_diffusive_conservation" in result_names
    @test "conservative_species_state" in result_names
    @test "conservative_species_face_advection" in result_names
    @test "mixture_corrected_species_diffusion" in result_names
    @test "stationary_contact_muscl" in result_names
    @test "muscl_thornber_free_stream" in result_names
    @test all(result -> !isempty(result.diagnostics) || result.status == :pass, report.results)

    @test error_norm([1.0, -1.0], [1.0, 1.0], :L1) == 1.0
    @test error_norm([3.0, 4.0], [1.0, 1.0], :Linf) == 4.0
    orders = observed_convergence_orders([0.25, 0.0625, 0.015625], [0.5, 0.25, 0.125])
    @test orders == Union{Nothing, Float64}[nothing, 2.0, 2.0]

    machine_report = report_json(report)
    @test occursin("\"schema_version\": \"1.0\"", machine_report)
    @test occursin("\"status\": \"pass\"", machine_report)
    @test occursin("FVMFramework Verification Report", report_string(report))
end
