module FVMVerification

using ComponentArrays
using Dates
using Ferrite
using ForwardDiff
using FVMFramework
using LinearAlgebra
using OrdinaryDiffEq
using Printf
using Random
using SparseArrays
using TOML
using Unitful

import ADTypes
import SparseConnectivityTracer

include("utilities/reporting.jl")
include("utilities/state_validation.jl")
include("utilities/convergence.jl")
include("components/compressible_case.jl")
include("components/viscous_checks.jl")
include("components/muscl_checks.jl")
include("components/thornber_checks.jl")
include("components/species_checks.jl")
include("turbulence/sst_case.jl")
include("turbulence/sst_checks.jl")
include("residual/residual_checks.jl")
include("jacobian/jacobian_checks.jl")
include("benchmarks/contact_discontinuity.jl")
include("benchmarks/viscous_analytical.jl")
include("benchmarks/sod_shock_tube.jl")
include("mms/manufactured_solutions.jl")
include("mutation/mutation_checks.jl")
include("performance/performance_checks.jl")
include("integration/integration_checks.jl")
include("integration/long_time_checks.jl")
include("validate.jl")

export CompressibleFlowVerification
export VerificationReport, VerificationResult
export ConvergenceStudy, convergence_table, error_norm
export observed_convergence_orders, run_convergence_study
export has_failures, report_dict, report_json, report_string
export validate, write_report

end
