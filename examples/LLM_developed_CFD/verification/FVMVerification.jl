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
include("components/compressible_case.jl")
include("residual/residual_checks.jl")
include("jacobian/jacobian_checks.jl")
include("benchmarks/sod_shock_tube.jl")
include("integration/integration_checks.jl")
include("validate.jl")

export CompressibleFlowVerification
export VerificationReport, VerificationResult
export has_failures, report_dict, report_json, report_string
export validate, write_report

end
