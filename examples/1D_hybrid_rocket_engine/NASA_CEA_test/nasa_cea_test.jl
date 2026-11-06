using Revise

using FVMFramework

using Logging
using Unitful
using OrdinaryDiffEq
using NonlinearSolve
using Sparspak
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using ILUZero
using LinearAlgebra
using StaticArrays
using Clapeyron
using PreallocationTools
using DataInterpolations
using JSON

Revise.includet(joinpath(@__DIR__, "main_nasa_cea_test_setup.jl"))

t0 = 0.0
tMax = 10.0

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity)) #old version, sparsity pattern from t = 0 only

startup_prob = ODEProblem(ode_func, u0_vec, (t0, tMax), p_guess)

@time sol = solve(
    startup_prob,
    FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff())
    #callback = approximate_time_to_finish_cb,
)

@time sol = solve(
    startup_prob,
    #Tsit5(),
    #FBDF(linsolve = SparSpakFactorization())
    #FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff())
    FBDF(linsolve = KLUFactorization()),
    #callback = approximate_time_to_finish_cb,
)

@time sol = solve(
    startup_prob,
    Tsit5(),
    #callback = approximate_time_to_finish_cb,
)



