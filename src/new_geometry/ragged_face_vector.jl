"""
    CellFaceLayout(n_faces_per_cell)

Map a cell-local face index to a packed face index. Each cell owns one entry for
each of its faces, so an interior mesh face appears once for each adjacent cell.
"""
struct CellFaceLayout
    offsets::Vector{Int}

    function CellFaceLayout(offsets::Vector{Int}, ::Val{:offsets})
        isempty(offsets) && throw(ArgumentError("cell-face offsets cannot be empty"))
        first(offsets) == 1 || throw(ArgumentError("cell-face offsets must start at 1"))
        issorted(offsets) || throw(ArgumentError("cell-face offsets must be sorted"))
        return new(offsets)
    end
end

function CellFaceLayout(n_faces_per_cell::AbstractVector{<:Integer})
    offsets = Vector{Int}(undef, length(n_faces_per_cell) + 1)
    offsets[1] = 1

    for cell_id in eachindex(n_faces_per_cell)
        n_faces = n_faces_per_cell[cell_id]
        n_faces >= 0 || throw(ArgumentError("the number of faces cannot be negative"))
        offsets[cell_id + 1] = offsets[cell_id] + n_faces
    end

    return CellFaceLayout(offsets, Val(:offsets))
end

Base.length(layout::CellFaceLayout) = length(layout.offsets) - 1
Base.eachindex(layout::CellFaceLayout) = Base.OneTo(length(layout))

"""Return the packed indices belonging to `cell_id`."""
@inline function Base.getindex(layout::CellFaceLayout, cell_id::Int)
    @boundscheck checkbounds(layout.offsets, cell_id:(cell_id + 1))
    return layout.offsets[cell_id]:(layout.offsets[cell_id + 1] - 1)
end

"""Return the packed index of `face_idx` on `cell_id`."""
@inline function Base.getindex(layout::CellFaceLayout, cell_id::Int, face_idx::Int)
    @boundscheck begin
        checkbounds(layout.offsets, cell_id:(cell_id + 1))
        1 <= face_idx <= nfaces(layout, cell_id) ||
            throw(BoundsError(layout, (cell_id, face_idx)))
    end
    return layout.offsets[cell_id] + face_idx - 1
end

@inline nfaces(layout::CellFaceLayout, cell_id::Int) =
    layout.offsets[cell_id + 1] - layout.offsets[cell_id]

@inline total_faces(layout::CellFaceLayout) = last(layout.offsets) - 1

"""
    CellFaceVector(layout, data)

A flat vector with an additional `(cell_id, face_idx)` indexing mode. The
one-index form is deliberately ordinary flat-vector indexing; cell-local access
is always written as `values[cell_id, face_idx]`.
"""
struct CellFaceVector{T, V <: AbstractVector{T}} <: AbstractVector{T}
    layout::CellFaceLayout
    data::V

    function CellFaceVector(layout::CellFaceLayout, data::V) where {T, V <: AbstractVector{T}}
        length(data) == total_faces(layout) || throw(DimensionMismatch(
            "the face data has $(length(data)) entries but the layout requires $(total_faces(layout))",
        ))
        return new{T, V}(layout, data)
    end
end

Base.size(values::CellFaceVector) = size(values.data)
Base.axes(values::CellFaceVector) = axes(values.data)
Base.IndexStyle(::Type{<:CellFaceVector}) = IndexLinear()
Base.parent(values::CellFaceVector) = values.data

@inline Base.getindex(values::CellFaceVector, index::Int) = values.data[index]
@inline Base.setindex!(values::CellFaceVector, value, index::Int) =
    (values.data[index] = value)

@inline function Base.getindex(values::CellFaceVector, cell_id::Int, face_idx::Int)
    packed_face_id = values.layout[cell_id, face_idx]
    return @inbounds values.data[packed_face_id]
end

@inline function Base.setindex!(values::CellFaceVector, value, cell_id::Int, face_idx::Int)
    packed_face_id = values.layout[cell_id, face_idx]
    return @inbounds values.data[packed_face_id] = value
end

@inline Base.getindex(values::CellFaceVector, cell_id::Int, ::Colon) =
    @views values.data[values.layout[cell_id]]

@inline Base.view(values::CellFaceVector, cell_id::Int, ::Colon) =
    view(values.data, values.layout[cell_id])

function Base.similar(values::CellFaceVector, ::Type{S}) where {S}
    return CellFaceVector(values.layout, similar(values.data, S))
end

Base.copy(values::CellFaceVector) = CellFaceVector(values.layout, copy(values.data))

cell_face_zeros(::Type{T}, layout::CellFaceLayout) where {T} =
    CellFaceVector(layout, zeros(T, total_faces(layout)))

cell_face_fill(value, layout::CellFaceLayout) =
    CellFaceVector(layout, fill(value, total_faces(layout)))
