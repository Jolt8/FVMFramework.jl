module SLAU2Verification

using ComponentArrays
using Dates
using Ferrite
using FVMFramework
using LinearAlgebra
using OrdinaryDiffEq
using Printf
using Unitful

include(joinpath(@__DIR__, "..", "utilities", "reporting.jl"))
include(joinpath(@__DIR__, "..", "utilities", "state_validation.jl"))
include(joinpath(
    @__DIR__,
    "..", "..", "navier_stokes_testing_grounds",
    "face_reconstructors", "first_order_face_reconstruction.jl",
))
include(joinpath(
    @__DIR__,
    "..", "..", "navier_stokes_testing_grounds",
    "navier_stokes_fluid_property_update_functions",
    "fluid_property_update_functions.jl",
))
include(joinpath(
    @__DIR__,
    "..", "..", "navier_stokes_testing_grounds",
    "sum_and_cap_functions", "cap_functions.jl",
))
include(joinpath(
    @__DIR__,
    "..", "..", "navier_stokes_testing_grounds",
    "species_transport", "conservative_species_transport.jl",
))
include(joinpath(
    @__DIR__,
    "..", "..", "navier_stokes_testing_grounds",
    "turbulence", "sst_k_omega.jl",
))
include(joinpath(
    @__DIR__,
    "..", "..", "navier_stokes_testing_grounds",
    "riemann_solvers", "SLAU2.jl",
))
include(joinpath(@__DIR__, "..", "benchmarks", "sod_shock_tube.jl"))
include("slau2_case.jl")
include("slau2_checks.jl")
include("validate_slau2.jl")

export SLAU2VerificationTarget
export VerificationReport, VerificationResult
export has_failures, report_dict, report_json, report_string
export validate_slau2, write_report

end
