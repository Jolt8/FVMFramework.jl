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

startup_prob.p.relaxation_factor = 0.1

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

desired_shown_steps = 250
save_fraction = desired_shown_steps / length(sol.t)

reduced_sol = (
    u = [[sol.u[i] for i in 1:(Int(ceil(1/save_fraction))):length(sol.u)]..., sol.u[end]],
    t = [sol.t[1:(Int(ceil(1/save_fraction))):length(sol.t)]..., sol.t[end]]
)

du_named, u_named = regenerate_fvm_state(reduced_sol, system, solve_system!, config.geo, p_guess, save_properties_and_caches = true, track_progress = true);

using Plots

plot(reduced_sol.t, [u_named[i].mass_fractions.nitrous_oxide[1] for i in eachindex(reduced_sol.t)], label = "nitrous_oxide")
plot!(reduced_sol.t, [u_named[i].mass_fractions.ethylene[1] for i in eachindex(reduced_sol.t)], label = "ethylene")
plot!(reduced_sol.t, [u_named[i].mass_fractions.nitrogen[1] for i in eachindex(reduced_sol.t)], label = "nitrogen")
plot!(reduced_sol.t, [u_named[i].mass_fractions.oxygen[1] for i in eachindex(reduced_sol.t)], label = "oxygen")
plot!(reduced_sol.t, [u_named[i].mass_fractions.carbon_dioxide[1] for i in eachindex(reduced_sol.t)], label = "carrbon_dioxide")
plot!(reduced_sol.t, [u_named[i].mass_fractions.carbon_monoxide[1] for i in eachindex(reduced_sol.t)], label = "carbon_monoxide")
plot!(reduced_sol.t, [u_named[i].mass_fractions.water[1] for i in eachindex(reduced_sol.t)], label = "water")
plot!(reduced_sol.t, [u_named[i].mass_fractions.hydrogen[1] for i in eachindex(reduced_sol.t)], label = "hydrogen")