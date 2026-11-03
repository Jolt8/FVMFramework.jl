using PreallocationTools: GeneralLazyBufferCache, get_tmp

#=
Per-cell species vectors (moles and mole fractions) that Clapeyron needs as plain AbstractVectors.

Why this exists:
Clapeyron takes the cell composition as a vector (vt_flash(model, V, T, moles, method)), but the state lives
in a VirtualFVMArray, so each cell's composition has to be gathered into a scratch vector first. That scratch
vector must have the same number type as the current RHS evaluation:
    - Float64 during a normal explicit step
    - ForwardDiff.Dual (any tag / chunk size) during Jacobian or parameter-sensitivity evaluations
    - SparseConnectivityTracer tracers during sparsity detection

How it works (same idea as ADGeometryCache in fvm/ad_compatible_mesh/new_ad_compatible_mesh.jl):
    1. SpeciesVectorCache owns one flat buffer per number type, created lazily the first time that type is seen.
    2. Once per RHS call, fvm_operator_with_species_vectors! picks the buffer matching the current number type
       and wraps the base system in a SystemWithSpeciesVectors view.
    3. Physics code keeps writing system.moles_vec[cell_id][species_idx] and passing system.moles_vec[cell_id]
       to Clapeyron. Inside the RHS that now resolves to a concretely typed view of the flat buffer, so the
       per-cell loop does no allocation and no dynamic dispatch.
=#


# ---------------------------------------------------------------------------------------------------------
# CellSpeciesVectors: a flat buffer laid out cell by cell, indexed as [cell_id] -> species vector of that cell
# ---------------------------------------------------------------------------------------------------------

struct CellSpeciesVectors{T, V <: AbstractVector{T}}
    flat_values::V
    n_cells::Int
    n_species::Int
end

function CellSpeciesVectors(flat_values::AbstractVector, n_species::Int)
    length(flat_values) % n_species == 0 || throw(DimensionMismatch("the flat species buffer length must be a multiple of n_species"))

    n_cells = length(flat_values) ÷ n_species

    return CellSpeciesVectors{eltype(flat_values), typeof(flat_values)}(flat_values, n_cells, n_species)
end

Base.length(species_vectors::CellSpeciesVectors) = species_vectors.n_cells
Base.eachindex(species_vectors::CellSpeciesVectors) = Base.OneTo(species_vectors.n_cells)
Base.eltype(::Type{<:CellSpeciesVectors{T}}) where {T} = T

# Species values for one cell are stored contiguously:
# [cell 1: species 1..n_species | cell 2: species 1..n_species | ...]
@inline function cell_species_range(species_vectors::CellSpeciesVectors, cell_id::Int)
    first_species_index = (cell_id - 1) * species_vectors.n_species + 1
    last_species_index = cell_id * species_vectors.n_species

    return first_species_index:last_species_index
end

# Returns a contiguous SubArray (no copy, no heap allocation). Writes through it go straight into the buffer,
# so system.moles_vec[cell_id][species_idx] = value updates the cache, and the same view can be handed to Clapeyron.
@inline function Base.getindex(species_vectors::CellSpeciesVectors, cell_id::Int)
    @boundscheck 1 <= cell_id <= species_vectors.n_cells || throw(BoundsError(species_vectors, cell_id))

    species_range = cell_species_range(species_vectors, cell_id)

    return @inbounds view(species_vectors.flat_values, species_range)
end

# Allows overwriting a whole cell at once, e.g. system.moles_vec[cell_id] = some_vector
@inline function Base.setindex!(species_vectors::CellSpeciesVectors, cell_values, cell_id::Int)
    @boundscheck 1 <= cell_id <= species_vectors.n_cells || throw(BoundsError(species_vectors, cell_id))

    species_range = cell_species_range(species_vectors, cell_id)

    copyto!(view(species_vectors.flat_values, species_range), cell_values)

    return cell_values
end

function Base.show(io::IO, species_vectors::CellSpeciesVectors{T}) where {T}
    print(io, "CellSpeciesVectors{", T, "}(", species_vectors.n_cells, " cells × ", species_vectors.n_species, " species)")
end


# ---------------------------------------------------------------------------------------------------------
# SpeciesVectorCache: one CellSpeciesVectors buffer per number type (Float64, Dual, tracer, ...)
# ---------------------------------------------------------------------------------------------------------

struct SpeciesVectorCache{B <: GeneralLazyBufferCache}
    n_cells::Int
    n_species::Int
    typed_buffers::B
end

function create_species_vector_cache(n_cells::Int, n_species::Int)
    # GeneralLazyBufferCache keys on typeof(scalar_seed), so every ForwardDiff tag / chunk size and every
    # tracer type gets its own buffer the first time it appears and reuses it afterwards.
    # (DiffCache is fixed to one chunk size and allocates a fresh array whenever the chunk size differs.)
    typed_buffers = GeneralLazyBufferCache() do scalar_seed
        flat_values = zeros(typeof(scalar_seed), n_cells * n_species)
        CellSpeciesVectors(flat_values, n_species)
    end

    return SpeciesVectorCache(n_cells, n_species, typed_buffers)
end

function get_species_vectors(species_vector_cache::SpeciesVectorCache, scalar_seed)
    return get_tmp(species_vector_cache.typed_buffers, scalar_seed)
end

# Outside the RHS (e.g. in callbacks) the base FVMSystem is used, so system.moles_vec is the cache itself.
# Indexing it then returns the Float64 buffer, which holds the values from the most recent Float64 RHS call.
@inline function Base.getindex(species_vector_cache::SpeciesVectorCache, cell_id::Int)
    float64_species_vectors = get_species_vectors(species_vector_cache, 0.0)::CellSpeciesVectors{Float64, Vector{Float64}}

    return float64_species_vectors[cell_id]
end

Base.length(species_vector_cache::SpeciesVectorCache) = species_vector_cache.n_cells
Base.eachindex(species_vector_cache::SpeciesVectorCache) = Base.OneTo(species_vector_cache.n_cells)


# ---------------------------------------------------------------------------------------------------------
# SystemWithSpeciesVectors: the FVMSystem plus concretely typed species vectors for the current RHS call
# ---------------------------------------------------------------------------------------------------------

struct SystemWithSpeciesVectors{S, M <: CellSpeciesVectors, X <: CellSpeciesVectors}
    base_system::S
    moles_vec::M
    mole_fractions_vec::X
end

# system.moles_vec and system.mole_fractions_vec resolve to the typed views.
# Everything else (region_groups, model, clapeyron_method, cad_data, ...) is forwarded to the base system.
@inline function Base.getproperty(system::SystemWithSpeciesVectors, name::Symbol)
    if hasfield(SystemWithSpeciesVectors, name)
        return getfield(system, name)
    end

    return getproperty(getfield(system, :base_system), name)
end

function Base.propertynames(system::SystemWithSpeciesVectors, private::Bool = false)
    return (fieldnames(SystemWithSpeciesVectors)..., propertynames(getfield(system, :base_system), private)...)
end

"""
    fvm_operator_with_species_vectors!(du_vec, u_vec, p, t, system, geo, solve_groups!)

Drop-in replacement for `fvm_operator!` that selects the species vectors matching the number type of the
current evaluation before calling `solve_groups!`. Requires `moles_vec` and `mole_fractions_vec`
(both `SpeciesVectorCache`s) in the system's additional data.
"""
function fvm_operator_with_species_vectors!(du_vec, u_vec, p, t, system, geo, solve_groups!)
    du, u = unpack_fvm_state(du_vec, u_vec, p, t, system)

    # Cell moles = density(u) * mass fraction(u) * cell volume(geo) / molecular weight, and p can carry duals
    # during parameter fitting, so the scratch number type has to cover the state, the parameters and the geometry.
    scalar_seed = first(u_vec) + first(p) + first(geo.cell_volumes)

    moles_vec = get_species_vectors(system.moles_vec, scalar_seed)
    mole_fractions_vec = get_species_vectors(system.mole_fractions_vec, scalar_seed)

    system_with_species_vectors = SystemWithSpeciesVectors(system, moles_vec, mole_fractions_vec)

    # Function barrier: inside solve_groups! the species vectors have a concrete type
    solve_groups!(du, u, p, t, system_with_species_vectors, geo)

    return nothing
end


# ---------------------------------------------------------------------------------------------------------
# Per-cell updates (called from fluid_property_update_ideal! / fluid_property_update_clapeyron!)
# ---------------------------------------------------------------------------------------------------------

function update_cell_moles_vec!(du, u, p, t, system, geo, cell_id)
    # species_idx is a Ref because a plain Int that is reassigned inside the do-block closure gets wrapped in
    # an untyped Core.Box, which allocates and dispatches dynamically on every species
    species_idx = Ref(1)

    cell_mass = u.density[cell_id] * geo.cell_volumes[cell_id]

    cell_moles = system.moles_vec[cell_id]

    for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
        cell_moles[species_idx[]] = (u_mass_fractions[species[cell_id]] * cell_mass) / u_molecular_weights[species[cell_id]]

        species_idx[] += 1
    end
end

function update_cell_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    species_idx = Ref(1)

    cell_mole_fractions = system.mole_fractions_vec[cell_id]

    for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
        cell_mole_fractions[species_idx[]] = u_mass_fractions[species[cell_id]] * (u.mw_avg[cell_id] / u_molecular_weights[species[cell_id]])

        species_idx[] += 1
    end
end

function update_cell_moles_and_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_vec!(du, u, p, t, system, geo, cell_id)
    update_cell_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
end


# ---------------------------------------------------------------------------------------------------------
# Setup: attach the caches to config.data so they end up in system.additional_data
# ---------------------------------------------------------------------------------------------------------

n_species = length(propertynames(config.u_proto.mass_fractions))

system_moles_vec_cache = create_species_vector_cache(config.geo.n_cells, n_species)
system_mole_fractions_vec_cache = create_species_vector_cache(config.geo.n_cells, n_species)

config.data = merge(config.data, (
    n_species = n_species,
    moles_vec = system_moles_vec_cache,
    mole_fractions_vec = system_mole_fractions_vec_cache,
))
