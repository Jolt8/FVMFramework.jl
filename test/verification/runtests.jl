include(joinpath(
    @__DIR__,
    "..", "..",
    "examples", "LLM_developed_CFD", "verification", "FVMVerification.jl",
))
using .FVMVerification

@testset "CFD verification framework" begin
    report = validate(; level = :unit, print_report = false)
    @test !has_failures(report)
    @test length(report.results) == 7
    @test all(result -> !isempty(result.diagnostics) || result.status == :pass, report.results)

    machine_report = report_json(report)
    @test occursin("\"schema_version\": \"1.0\"", machine_report)
    @test occursin("\"status\": \"pass\"", machine_report)
    @test occursin("FVMFramework Verification Report", report_string(report))
end
