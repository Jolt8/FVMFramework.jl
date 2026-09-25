# Overload power_by_squaring to support tracer types inside Clapeyron's complex root solver
#Base.power_by_squaring(x::Any, p::SparseConnectivityTracer.Dual) = p
Base.power_by_squaring(x::SparseConnectivityTracer.Dual, p::SparseConnectivityTracer.Dual) = x+p
#Base.power_by_squaring(x::Any, p::SparseConnectivityTracer.Dual) = x^p #this is the only thing left that fails, the above works though

# Overload ldexp for SparseConnectivityTracer.Dual
Base.ldexp(x::SparseConnectivityTracer.Dual, i::Int) = x

# Define explicit promotion rules to resolve ambiguity
Base.promote_rule(::Type{SparseConnectivityTracer.Dual{P, T}}, ::Type{ForwardDiff.Dual{Tx, Vx, Nx}}) where {P, T, Tx, Vx, Nx} = 
    ForwardDiff.Dual{Tx, promote_type(SparseConnectivityTracer.Dual{P, T}, Vx), Nx}

Base.promote_rule(::Type{ForwardDiff.Dual{Tx, Vx, Nx}}, ::Type{SparseConnectivityTracer.Dual{P, T}}) where {P, T, Tx, Vx, Nx} = 
    ForwardDiff.Dual{Tx, promote_type(SparseConnectivityTracer.Dual{P, T}, Vx), Nx}

# Define convert overload to resolve conversion ambiguity
Base.convert(::Type{ForwardDiff.Dual{Tx, Vx, Nx}}, y::D) where {Tx, Vx, Nx, P, T, D <: SparseConnectivityTracer.Dual{P, T}} = 
    ForwardDiff.Dual{Tx, Vx, Nx}(convert(Vx, y), ForwardDiff.Partials{Nx, Vx}(ntuple(i -> zero(Vx), Nx)))

# Resolve ambiguities for common binary operators and comparisons (+, -, *, /, ^, ==, <, <=)
for TracerType in (SparseConnectivityTracer.GradientTracer, SparseConnectivityTracer.HessianTracer)
    # Define ^(::D, ::D) using the exp(y * log(x)) identity
    @eval begin
        function Base.:^(x::D, y::D) where {P, T <: $TracerType, D <: SparseConnectivityTracer.Dual{P, T}}
            return exp(y * log(x))
        end
    end
    for op in (:+, :-, :*, :/, :^, Symbol("=="), Symbol("<"), Symbol("<="))
        @eval begin
            function Base.$op(x::ForwardDiff.Dual{Tx, Vx, Nx}, y::D) where {Tx, Vx, Nx, P, T <: $TracerType, D <: SparseConnectivityTracer.Dual{P, T}}
                x_prom, y_prom = promote(x, y)
                return Base.$op(x_prom, y_prom)
            end
            function Base.$op(x::D, y::ForwardDiff.Dual{Tx, Vx, Nx}) where {Tx, Vx, Nx, P, T <: $TracerType, D <: SparseConnectivityTracer.Dual{P, T}}
                x_prom, y_prom = promote(x, y)
                return Base.$op(x_prom, y_prom)
            end
        end
    end
end