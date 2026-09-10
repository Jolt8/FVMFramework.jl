using DynamicQuantities
using OrdinaryDiffEq
using Ferrite
#using FerriteGmsh
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using NonlinearSolve
using Sparspak
using Dates
using DataInterpolations
using Random
using Lux

using XLSX
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using ForwardDiff
using DataFrames
using CSV
using JLD2
using Plots

using SymbolicRegression

using FVMFramework

model_path = joinpath(@__DIR__, "learned_model.jld2")

model_data = load(model_path)
loaded_nn = model_data["nn"]
loaded_p = model_data["p_trained"]
loaded_st = model_data["st_nn"]

T_ref = 300.0
k_ref = 50.0
d_ref = 0.2
A_ref = 1.0

x = ComponentVector(
    temp_a = collect(290:1:310),
    temp_b = collect(290:1:310),
    k_a = collect(1.0:1.0:50.0),
    k_b = collect(1.0:1.0:50.0),
    distance = collect(0.05:0.01:0.5),
    area = collect(1.0:0.5:5.0)
)

test_data_1 = [loaded_nn([temp_a/T_ref, x.temp_b[1]/T_ref, x.k_a[1]/k_ref, x.k_b[1]/k_ref, x.distance[1]/d_ref, x.area[1]/A_ref], loaded_p, loaded_st)[1][1] for temp_a in x.temp_a]
test_data_2 = [loaded_nn([x.temp_a[1]/T_ref, temp_b/T_ref, x.k_a[1]/k_ref, x.k_b[1]/k_ref, x.distance[1]/d_ref, x.area[1]/A_ref], loaded_p, loaded_st)[1][1] for temp_b in x.temp_b]

plot(x.temp_a, test_data_1, label = "varied_temp_a")
plot!(x.temp_b, test_data_2, label = "varied_temp_b")

n_datapoints = 1000
x_data = Vector{Float64}[]
y_data = Float64[]

for i in 1:n_datapoints
    rand_temp_a = rand(x.temp_a)
    rand_temp_b = rand(x.temp_b)
    rand_k_a = rand(x.k_a)
    rand_k_b = rand(x.k_b)
    rand_distance = rand(x.distance)
    rand_area = rand(x.area)

    x_input = [rand_temp_a, rand_temp_b, rand_k_a, rand_k_b, rand_distance, rand_area]
    x_reference = x_input ./ [T_ref, T_ref, k_ref, k_ref, d_ref, A_ref]
    
    nn_val, _ = loaded_nn(x_reference, loaded_p, loaded_st)

    push!(x_data, x_input)
    push!(y_data, nn_val[1])
end

x_units = [u"K", u"K", u"W/(m*K)", u"W/(m*K)", u"m", u"m^2"]
y_unit = u"W"

x_matrix = reduce(hcat, x_data)

options = Options(
    binary_operators = [+, *, /, -],
    unary_operators = [],
    input_stream = devnull,
    early_stop_condition = 1e-20,
    maxsize = 20,
    timeout_in_seconds = 10,
)

hall_of_fame = equation_search(
    x_matrix, y_data, niterations = 40, options = options,
    parallelism = :serial,
    X_units = x_units, y_units = y_unit
)
