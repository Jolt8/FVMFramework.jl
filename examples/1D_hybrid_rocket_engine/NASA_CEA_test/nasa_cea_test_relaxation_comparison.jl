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
using BenchmarkTools
using Plots

Revise.includet(joinpath(@__DIR__, "main_nasa_cea_test_setup.jl"))

t0 = 0.0
tMax = 10.0

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity)) #old version, sparsity pattern from t = 0 only

startup_prob = ODEProblem(ode_func, u0_vec, (t0, tMax), p_guess)

sparspak_times = []
klu_times = []
tsit_5_times = []

relaxation_factors_to_test = collect(exp10.(range(log10(1e-4), log10(100), 20)))

for relaxation_factor in relaxation_factors_to_test
    original_p_named = ComponentVector(startup_prob.p, system.p_axes)

    original_p_named.relaxation_factor = relaxation_factor

    new_prob = remake(startup_prob; p = original_p_named)

    @show relaxation_factor
    
    sparspak_time = @belapsed sol = solve(
        $new_prob,
        FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff()),
        #callback = approximate_time_to_finish_cb,
    )
    push!(sparspak_times, sparspak_time)

    klu_time = @belapsed sol = solve(
        $new_prob,
        FBDF(linsolve = KLUFactorization()),
        #callback = approximate_time_to_finish_cb,
    )
    push!(klu_times, klu_time)

    #=
    tsit_5_time = @belapsed sol = solve(
        $new_prob,
        Tsit5(),
        callback = approximate_time_to_finish_cb,
    )
    push!(tsit_5_times, tsit_5_time)
    =#
end

plot(relaxation_factors_to_test, klu_times, label = "KLU")
plot!(relaxation_factors_to_test, sparspak_times, label = "Sparspak")
#plot!(relaxation_factors_to_test, tsit_5_times, xaxis=:log10, label = "Tsit5")