
#=
How these files flow into each other (in order)
    - tracer.jl traces a flux_function or region_function
    - tracer_merging.jl merges the results from the tracer and records the times read and written to du and u as well as access history
    - tracer_classifying.jl classifies the variables as state, cache, or fixed
    - tracer_to_CompArray.jl converts the classified variables to a ComponentArray
=#

#Structs
struct VarPath
    root::Symbol #:du or :u
    path::Vector{Symbol} #[:reactions, :MSR_rxn, :kf_A]
end

struct VarAccess
    path::VarPath
    access_type::Symbol #:read or :write
end

struct TracerContext
    access_logs::Vector{VarAccess}
    pending_reads::Vector{VarPath}
    write_dependencies::Dict{VarPath,Set{VarPath}}
    loop_context::Dict{Symbol,Vector{Symbol}}
end

struct TracerVariable
    path::VarPath
    tracer_ref::Ref{TracerContext}
    value::Float64
end

#Setup 
function create_tracer_context(ctx_loop_context)
    ctx = TracerContext(VarAccess[], VarPath[], Dict(), ctx_loop_context)
    du_tracer = TracerVariable(
        VarPath(:du, Symbol[]),
        ctx,
        1.0
    )
    u_tracer = TracerVariable(
        VarPath(:u, Symbol[]),
        ctx,
        1.0
    )
    return du_tracer, u_tracer, ctx
end

#Functions and Overloads

function Base.getproperty(t::TracerVariable, sym::Symbol)
    parent_path = getfield(t, :path)
    new_path = VarPath(parent_path.root, vcat(parent_path.path, sym))

    return TracerVariable(new_path, getfield(t, :tracer_ref)[], rand())
end

function Base.setproperty!(t::TracerVariable, sym::Symbol, val)
    ctx = getfield(t, :tracer_ref)[]
    parent_path = getfield(t, :path)
    write_path = VarPath(parent_path.root, vcat(parent_path.path, sym))

    push!(ctx.access_logs, VarAccess(write_path, :write))

    if !haskey(ctx.write_dependencies, write_path)
        ctx.write_dependencies[write_path] = Set{VarPath}()
    end
    union!(ctx.write_dependencies[write_path], ctx.pending_reads)

    empty!(ctx.pending_reads)
end

#hopefully it's safe to assume that whenever it's indexed by an integer, it's being indexed by face
function Base.getindex(t::TracerVariable, idx::Int)
    ctx = getfield(t, :tracer_ref)[]
    context = ctx.loop_context[:face_idxs][idx]
    parent_path = getfield(t, :path)
    return TracerVariable(
        VarPath(parent_path.root, vcat(parent_path.path, context)),
        getfield(t, :tracer_ref),
        getfield(t, :value)
    )
end


function Base.getindex(t::TracerVariable, sym::Symbol)
    parent_path = getfield(t, :path)
    return TracerVariable(
        VarPath(parent_path.root, vcat(parent_path.path, sym)),
        getfield(t, :tracer_ref),
        getfield(t, :value)
    )
end

function Base.getindex(t::TracerVariable, s::String)
    parent_path = getfield(t, :path)
    sym = Symbol(s)
    return TracerVariable(
        VarPath(parent_path.root, vcat(parent_path.path, sym)),
        getfield(t, :tracer_ref),
        getfield(t, :value)
    )
end
#FIXME: this needs to be fixed 
#errors with no method matching ndims(::Type{Nothing})
#this is unused for now, we'll come back to this later
function Base.getindex(t::TracerVariable, ::Colon, i::Int)
    ctx = getfield(t, :tracer_ref)[]
    context = ctx.loop_context[getfield(t, :path).path[end]]

    for c in context
        getindex(getindex(t, i), c)
    end
end

function Base.setindex!(t::TracerVariable, val, i::Int)
    ctx = getfield(t, :tracer_ref)[]
    context = ctx.loop_context[:face_idxs][i]
    parent_path = getfield(t, :path)
    write_path = VarPath(parent_path.root, vcat(parent_path.path, context))

    push!(ctx.access_logs, VarAccess(write_path, :write))

    if !haskey(ctx.write_dependencies, write_path)
        ctx.write_dependencies[write_path] = Set{VarPath}()
    end
    union!(ctx.write_dependencies[write_path], ctx.pending_reads)
    empty!(ctx.pending_reads)
end


function Base.setindex!(t::TracerVariable, val, sym::Symbol)
    ctx = getfield(t, :tracer_ref)[]
    parent_path = getfield(t, :path)
    write_path = VarPath(parent_path.root, vcat(parent_path.path, sym))

    push!(ctx.access_logs, VarAccess(write_path, :write))

    if !haskey(ctx.write_dependencies, write_path)
        ctx.write_dependencies[write_path] = Set{VarPath}()
    end
    union!(ctx.write_dependencies[write_path], ctx.pending_reads)
    empty!(ctx.pending_reads)
end

function Base.setindex!(t::TracerVariable, val, s::String)
    ctx = getfield(t, :tracer_ref)[]
    parent_path = getfield(t, :path)
    write_path = VarPath(parent_path.root, vcat(parent_path.path, Symbol(s)))

    push!(ctx.access_logs, VarAccess(write_path, :write))

    if !haskey(ctx.write_dependencies, write_path)
        ctx.write_dependencies[write_path] = Set{VarPath}()
    end
    union!(ctx.write_dependencies[write_path], ctx.pending_reads)
    empty!(ctx.pending_reads)
end


function logread(t::TracerVariable)
    ctx = getfield(t, :tracer_ref)[]
    vp = getfield(t, :path)
    push!(ctx.access_logs, VarAccess(vp, :read))
    push!(ctx.pending_reads, vp)
end

function Base.propertynames(t::TracerVariable) #this is just to make sure it doesn't yield path, tracer_ref, and value
    ctx = getfield(t, :tracer_ref)[]
    context = ctx.loop_context[getfield(t, :path).path[end]]
    return context
end

#this is for when properties are summed by face
function Base.eachindex(t::TracerVariable)
    ctx = getfield(t, :tracer_ref)[]
    context = ctx.loop_context[:face_idxs]
    return context
end

#MATH
Base.:+(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) + b
)
Base.:-(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) - b
)
Base.:*(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) * b
)
Base.:/(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value) / b
)
Base.:^(a::TracerVariable, b::Number) = TracerVariable(
    getfield(a, :path), getfield(a, :tracer_ref), getfield(a, :value)^b
)
Base.:isless(x::TracerVariable, y::Number) = (logread(x); return getfield(x, :value) < y)

Base.:+(a::Number, b::TracerVariable) = (logread(b); return a + getfield(b, :value))
Base.:-(a::Number, b::TracerVariable) = (logread(b); return a - getfield(b, :value))
Base.:*(a::Number, b::TracerVariable) = (logread(b); return a * getfield(b, :value))
Base.:/(a::Number, b::TracerVariable) = (logread(b); return a / getfield(b, :value))
Base.:^(a::Number, b::TracerVariable) = (logread(b); return a^getfield(b, :value))
Base.:isless(x::Number, y::TracerVariable) = (logread(y); return x < getfield(y, :value))

Base.:+(a::TracerVariable, b::TracerVariable) = (logread(a);
logread(b);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value) + getfield(b, :value)
)
)
Base.:-(a::TracerVariable, b::TracerVariable) = (logread(a);
logread(b);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value) - getfield(b, :value)
)
)
Base.:*(a::TracerVariable, b::TracerVariable) = (logread(a);
logread(b);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value) * getfield(b, :value)
)
)
Base.:/(a::TracerVariable, b::TracerVariable) = (logread(a);
logread(b);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value) / getfield(b, :value)
)
)
Base.:^(a::TracerVariable, b::TracerVariable) = (logread(a);
logread(b);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    getfield(a, :value)^getfield(b, :value)
)
)
Base.:isless(x::TracerVariable, y::TracerVariable) = (logread(x);
logread(y);
return getfield(x, :value) < getfield(y, :value)
)
Base.:sqrt(a::TracerVariable) = (logread(a);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    sqrt(getfield(a, :value))
)
)
Base.:exp(a::TracerVariable) = (logread(a);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    exp(getfield(a, :value))
)
)
Base.:-(a::TracerVariable) = (logread(a);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    -getfield(a, :value)
)
)


#MISC:
#I really wish there was a way to shorten this
Base.:clamp(a::TracerVariable, b::TracerVariable, c::TracerVariable) = (logread(a); logread(b); logread(c);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    clamp(getfield(a, :value), getfield(b, :value), getfield(c, :value))
)
)
Base.:clamp(a::TracerVariable, b::TracerVariable, c::Number) = (logread(a); logread(b);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    clamp(getfield(a, :value), getfield(b, :value), c)
)
)
Base.:clamp(a::TracerVariable, b::Number, c::TracerVariable) = (logread(a); logread(c);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    clamp(getfield(a, :value), b, getfield(c, :value))
)
)
Base.:clamp(a::Number, b::TracerVariable, c::TracerVariable) = (logread(b); logread(c);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(b, :tracer_ref),
    clamp(a, getfield(b, :value), getfield(c, :value))
)
)
Base.:clamp(a::TracerVariable, b::Number, c::Number) = (logread(a);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(a, :tracer_ref),
    clamp(getfield(a, :value), b, c)
)
)
Base.:clamp(a::Number, b::TracerVariable, c::Number) = (logread(b);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(b, :tracer_ref),
    clamp(a, getfield(b, :value), c)
)
)
Base.:clamp(a::Number, b::Number, c::TracerVariable) = (logread(c);
return TracerVariable(
    VarPath(:_expr, Symbol[]),
    getfield(c, :tracer_ref),
    clamp(a, b, getfield(c, :value))
)
)

#Example required for input:
#=
ctx_loop_context = Dict(
    :mass_fractions => [:methanol, :water, :carbon_monoxide, :hydrogen, :carbon_dioxide],
    :molar_concentrations => [:methanol, :water, :carbon_monoxide, :hydrogen, :carbon_dioxide],
    :reforming_reactions => [:MSR_rxn, :MD_rxn, :WGS_rxn]
)

ctx = TracerContext(VarAccess[], VarPath[], Dict(), ctx_loop_context)

u_tracer = TracerVariable(
    VarPath(:u, Symbol[]),
    ctx,
    1.0
)

du_tracer = TracerVariable(
    VarPath(:du, Symbol[]),
    ctx,
    1.0
)
=#