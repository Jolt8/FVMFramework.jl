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
using Unitful

import ADTypes
import SparseConnectivityTracer

include("utilities/reporting.jl")
include("utilities/state_validation.jl")
include("utilities/convergence.jl")
include("components/compressible_case.jl")
include("components/muscl_checks.jl")
include("components/thornber_checks.jl")
include("components/species_checks.jl")
include("residual/residual_checks.jl")
include("jacobian/jacobian_checks.jl")
include("benchmarks/contact_discontinuity.jl")
include("benchmarks/sod_shock_tube.jl")
include("integration/integration_checks.jl")
include("validate.jl")

export CompressibleFlowVerification
export VerificationReport, VerificationResult
export ConvergenceStudy, convergence_table, error_norm
export observed_convergence_orders, run_convergence_study
export has_failures, report_dict, report_json, report_string
export validate, write_report

end
