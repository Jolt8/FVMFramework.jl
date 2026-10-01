struct ConvergenceStudy
    grid_sizes::Vector{Int}
    grid_spacings::Vector{Float64}
    errors::Vector{Float64}
    observed_orders::Vector{Union{Nothing, Float64}}
    norm_name::Symbol
end

function error_norm(errors, weights, norm_name)
    if length(errors) != length(weights)
        throw(DimensionMismatch("errors and weights must have equal lengths"))
    end
    if any(weight -> weight < 0.0, weights)
        throw(ArgumentError("error-norm weights must be non-negative"))
    end
    weight_sum = sum(weights)
    if weight_sum <= 0.0
        throw(ArgumentError("error-norm weights must have a positive sum"))
    end

    if norm_name == :L1
        return sum(weights .* abs.(errors)) / weight_sum
    elseif norm_name == :L2
        return sqrt(sum(weights .* abs2.(errors)) / weight_sum)
    elseif norm_name == :Linf
        return norm(errors, Inf)
    end
    throw(ArgumentError("norm_name must be :L1, :L2, or :Linf"))
end

function observed_convergence_orders(errors, grid_spacings)
    if length(errors) != length(grid_spacings)
        throw(DimensionMismatch("errors and grid spacings must have equal lengths"))
    end
    orders = Union{Nothing, Float64}[nothing]
    for level in 2:length(errors)
        if errors[level] <= 0.0 || errors[level - 1] <= 0.0
            push!(orders, nothing)
            continue
        end
        refinement_ratio = grid_spacings[level - 1] / grid_spacings[level]
        if refinement_ratio <= 1.0
            throw(ArgumentError("grid spacings must decrease between convergence levels"))
        end
        push!(orders, log(errors[level - 1] / errors[level]) / log(refinement_ratio))
    end
    return orders
end

function run_convergence_study(grid_sizes, evaluator; norm_name = :L2)
    if length(grid_sizes) < 2
        throw(ArgumentError("a convergence study needs at least two grid levels"))
    end
    if any(grid_size -> grid_size <= 0, grid_sizes)
        throw(ArgumentError("grid sizes must be positive"))
    end

    spacings = 1.0 ./ Float64.(grid_sizes)
    errors = Float64[evaluator(grid_size, norm_name) for grid_size in grid_sizes]
    orders = observed_convergence_orders(errors, spacings)
    return ConvergenceStudy(
        collect(grid_sizes),
        spacings,
        errors,
        orders,
        norm_name,
    )
end

function convergence_table(study::ConvergenceStudy)
    rows = Dict{String, Any}[]
    for level in eachindex(study.grid_sizes)
        push!(rows, Dict{String, Any}(
            "grid" => study.grid_sizes[level],
            "spacing" => study.grid_spacings[level],
            "error" => study.errors[level],
            "observed_order" => study.observed_orders[level],
        ))
    end
    return rows
end
