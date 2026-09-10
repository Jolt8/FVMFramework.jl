using DynamicQuantities
using DataDrivenDiffEq
using DataDrivenSR
#using SymbolicRegression
using Random

# Function to calculate pressure drop per unit length using the Ergun equation
function ergun_pressure_drop(epsilon, d_p, mu, v_s, rho)
    term1 = 150.0 * ((1.0 - epsilon) ^ 2) / (epsilon ^ 3) * mu * v_s / (d_p ^ 2)
    term2 = 1.75 * (1.0 - epsilon) / (epsilon ^ 3) * rho * (v_s ^ 2) / d_p
    return term1 + term2
end

# Generate synthetic dataset
n_datapoints = 1000
x_data = Vector{Float64}[]
y_data = Float64[]

# Seed for reproducibility
Random.seed!(42)

for i in 1:n_datapoints
    # void fraction (dimensionless)
    epsilon = rand() * (0.5 - 0.3) + 0.3
    
    # particle diameter (m)
    d_p = rand() * (0.01 - 0.001) + 0.001
    
    # dynamic viscosity (kg/(m*s))
    mu = 10 ^ (rand() * (log10(1e-2) - log10(1e-5)) + log10(1e-5))
    
    # superficial velocity (m/s)
    v_s = rand() * (5.0 - 0.1) + 0.1
    
    # fluid density (kg/m^3)
    rho = 10 ^ (rand() * (log10(1000.0) - log10(1.0)) + log10(1.0))
    
    dp_L = ergun_pressure_drop(epsilon, d_p, mu, v_s, rho)
    
    push!(x_data, [epsilon, d_p, mu, v_s, rho])
    push!(y_data, dp_L)
end

# Format input data as matrices (features x samples)
x_matrix = reduce(hcat, x_data)
y_matrix = reshape(y_data, 1, :)

# Define physical units using DynamicQuantities
# Column 1: epsilon (dimensionless)
# Column 2: d_p (m)
# Column 3: mu (kg/(m*s))
# Column 4: v_s (m/s)
# Column 5: rho (kg/m^3)
x_units = [u"1", u"m", u"kg/(m*s)", u"m/s", u"kg/m^3"]
y_unit = u"kg/(m^2*s^2)"

# Define the DirectDataDrivenProblem
problem = DirectDataDrivenProblem(x_matrix, y_matrix, name = :ErgunProblem)

# Define search options for SymbolicRegression.jl
options = Options(
    binary_operators = [+, *, /, -],
    unary_operators = [],
    input_stream = devnull,
    early_stop_condition = 1e-20,
    maxsize = 25,
    timeout_in_seconds = 300,
)

# Initialize the EQSearch solver algorithm
alg = EQSearch(eq_options = options)

# Solve the data-driven problem with dimensional constraints passed as kwargs
res = solve(problem, alg, X_units = x_units, y_units = y_unit)

# Print the discovered equation / basis
println(res)
println(get_basis(res))
