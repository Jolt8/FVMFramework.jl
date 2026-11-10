using FVMFramework
using LinearAlgebra: dot, norm
using PreallocationTools: GeneralLazyBufferCache, get_tmp
using StaticArrays: SVector

struct ADGeometryBufferLayout
    cell_volumes::UnitRange{Int}
    cell_centroids::UnitRange{Int}
    cell_neighbor_normals::UnitRange{Int}
    cell_neighbor_distances::UnitRange{Int}
    cell_face_areas::UnitRange{Int}
    cell_face_normals::UnitRange{Int}
    cell_face_distances::UnitRange{Int}
    cell_face_centroids::UnitRange{Int}
    node_coordinates::UnitRange{Int}
end

struct ADGeometryCache{G, V, C}
    base_geo::G
    initial_values::V
    buffers::C
    layout::ADGeometryBufferLayout
end

struct ADGeometryView{G, V, CV, CC, CNN, CND, CFA, CFN, CFD, CFC, NC}
    base_geo::G
    values::V
    cell_volumes::CV
    cell_centroids::CC
    cell_neighbor_normals::CNN
    cell_neighbor_distances::CND
    cell_face_areas::CFA
    cell_face_normals::CFN
    cell_face_distances::CFD
    cell_face_centroids::CFC
    node_coordinates::NC
end

struct SVector3View{T, V <: AbstractVector{T}} <: AbstractVector{SVector{3, T}}
    data::V

    function SVector3View(data::V) where {T, V <: AbstractVector{T}}
        length(data) % 3 == 0 || throw(DimensionMismatch(
            "the coordinate buffer must contain a multiple of three values",
        ))
        return new{T, V}(data)
    end
end

Base.size(values::SVector3View) = (length(values.data) ÷ 3,)
Base.axes(values::SVector3View) = (Base.OneTo(length(values)),)
Base.IndexStyle(::Type{<:SVector3View}) = IndexLinear()
Base.parent(values::SVector3View) = values.data

@inline function Base.getindex(values::SVector3View, index::Int)
    @boundscheck checkbounds(values, index)
    offset = 3 * (index - 1)
    return @inbounds SVector(
        values.data[offset + 1],
        values.data[offset + 2],
        values.data[offset + 3],
    )
end

@inline function Base.setindex!(values::SVector3View, value, index::Int)
    @boundscheck checkbounds(values, index)
    offset = 3 * (index - 1)
    @inbounds begin
        values.data[offset + 1] = value[1]
        values.data[offset + 2] = value[2]
        values.data[offset + 3] = value[3]
    end
    return value
end

function Base.getproperty(geo::ADGeometryView, name::Symbol)
    if hasfield(typeof(geo), name)
        return getfield(geo, name)
    end

    base_geo = getfield(geo, :base_geo)
    if hasfield(typeof(base_geo), name)
        return getfield(base_geo, name)
    end

    return getfield(geo, name)
end

function Base.propertynames(geo::ADGeometryView, private::Bool = false)
    return (fieldnames(typeof(geo))..., fieldnames(typeof(getfield(geo, :base_geo)))...)
end

@inline function _take_geometry_range(cursor, count)
    range = cursor:(cursor + count - 1)
    return range, last(range) + 1
end

function _create_geometry_buffer_layout(geo)
    n_cells = length(geo.cell_volumes)
    n_cell_faces = length(geo.cell_face_areas)
    n_nodes = length(geo.node_coordinates)
    cursor = 1

    cell_volumes, cursor = _take_geometry_range(cursor, n_cells)
    cell_centroids, cursor = _take_geometry_range(cursor, 3 * n_cells)
    cell_neighbor_normals, cursor = _take_geometry_range(cursor, 3 * n_cell_faces)
    cell_neighbor_distances, cursor = _take_geometry_range(cursor, n_cell_faces)
    cell_face_areas, cursor = _take_geometry_range(cursor, n_cell_faces)
    cell_face_normals, cursor = _take_geometry_range(cursor, 3 * n_cell_faces)
    cell_face_distances, cursor = _take_geometry_range(cursor, n_cell_faces)
    cell_face_centroids, cursor = _take_geometry_range(cursor, 3 * n_cell_faces)
    node_coordinates, cursor = _take_geometry_range(cursor, 3 * n_nodes)

    return ADGeometryBufferLayout(
        cell_volumes,
        cell_centroids,
        cell_neighbor_normals,
        cell_neighbor_distances,
        cell_face_areas,
        cell_face_normals,
        cell_face_distances,
        cell_face_centroids,
        node_coordinates,
    ), cursor - 1
end

function _copy_scalar_values!(destination, range, source)
    copyto!(@view(destination[range]), source)
    return nothing
end

function _copy_coordinate_values!(destination, range, source)
    destination_values = @view destination[range]
    destination_idx = firstindex(destination_values)
    for coordinate in source
        @inbounds begin
            destination_values[destination_idx] = coordinate[1]
            destination_values[destination_idx + 1] = coordinate[2]
            destination_values[destination_idx + 2] = coordinate[3]
        end
        destination_idx += 3
    end
    return nothing
end

function _pack_geometry(geo, layout, buffer_length)
    scalar_type = eltype(geo.cell_volumes)
    values = Vector{scalar_type}(undef, buffer_length)

    _copy_scalar_values!(values, layout.cell_volumes, geo.cell_volumes)
    _copy_coordinate_values!(values, layout.cell_centroids, geo.cell_centroids)
    _copy_coordinate_values!(values, layout.cell_neighbor_normals, parent(geo.cell_neighbor_normals))
    _copy_scalar_values!(values, layout.cell_neighbor_distances, parent(geo.cell_neighbor_distances))
    _copy_scalar_values!(values, layout.cell_face_areas, parent(geo.cell_face_areas))
    _copy_coordinate_values!(values, layout.cell_face_normals, parent(geo.cell_face_normals))
    _copy_scalar_values!(values, layout.cell_face_distances, parent(geo.cell_face_distances))
    _copy_coordinate_values!(values, layout.cell_face_centroids, parent(geo.cell_face_centroids))
    _copy_coordinate_values!(values, layout.node_coordinates, geo.node_coordinates)

    return values
end

@inline _scalar_view(values, range) = @view values[range]

@inline function _coordinate_view(values, range)
    return SVector3View(@view values[range])
end

@inline function _cell_face_scalar_view(values, range, cell_face_layout)
    return FVMFramework.CellFaceVector(
        cell_face_layout,
        _scalar_view(values, range),
    )
end

@inline function _cell_face_coordinate_view(values, range, cell_face_layout)
    return FVMFramework.CellFaceVector(
        cell_face_layout,
        _coordinate_view(values, range),
    )
end

function _geometry_view(base_geo, values, layout)
    cell_face_layout = base_geo.cell_face_layout
    return ADGeometryView(
        base_geo,
        values,
        _scalar_view(values, layout.cell_volumes),
        _coordinate_view(values, layout.cell_centroids),
        _cell_face_coordinate_view(
            values, layout.cell_neighbor_normals, cell_face_layout,
        ),
        _cell_face_scalar_view(
            values, layout.cell_neighbor_distances, cell_face_layout,
        ),
        _cell_face_scalar_view(values, layout.cell_face_areas, cell_face_layout),
        _cell_face_coordinate_view(values, layout.cell_face_normals, cell_face_layout),
        _cell_face_scalar_view(values, layout.cell_face_distances, cell_face_layout),
        _cell_face_coordinate_view(
            values, layout.cell_face_centroids, cell_face_layout,
        ),
        _coordinate_view(values, layout.node_coordinates),
    )
end

"""
    create_ad_geometry_cache(geo)

Create a geometry cache whose numeric storage is selected from the scalar type of
the current state. A new buffer is allocated only the first time a scalar type is
encountered, so different `ForwardDiff` tags and chunk sizes do not require a new
geometry struct definition.
"""
function create_ad_geometry_cache(geo)
    layout, buffer_length = _create_geometry_buffer_layout(geo)
    initial_values = _pack_geometry(geo, layout, buffer_length)
    buffers = GeneralLazyBufferCache() do seed
        values = similar(initial_values, typeof(seed))
        _geometry_view(geo, values, layout)
    end
    return ADGeometryCache(geo, initial_values, buffers, layout)
end

"""
    get_ad_geometry(cache, scalar_seed)

Return a typed geometry view backed by the reusable buffer associated with
`typeof(scalar_seed)`. The buffer is reset to the original geometry each call.
The returned view aliases the cache and should only be used for the current RHS
evaluation.
"""
function get_ad_geometry(cache::ADGeometryCache, scalar_seed)
    geo = get_tmp(cache.buffers, scalar_seed)
    copyto!(geo.values, cache.initial_values)
    return geo
end

get_ad_geometry(cache::ADGeometryCache, state::AbstractArray) =
    get_ad_geometry(cache, first(state))

function get_ad_geometry(cache::ADGeometryCache, state::VirtualFVMArray)
    state_values = first(getfield(state, :data))
    return get_ad_geometry(cache, first(state_values))
end

create_diff_cache_geo(geo, _grid) = create_ad_geometry_cache(geo)
create_diff_cache_geo(geo, _grid, _chunk_size) = create_ad_geometry_cache(geo)
regenerate_ad_geo!(state, cache::ADGeometryCache) = get_ad_geometry(cache, state)

function update_fvm_geometry!(geo::ADGeometryView, node_coordinates = geo.node_coordinates)
    for cell_id in eachindex(geo.cell_node_ids)
        volume, centroid = FVMFramework._cell_volume_and_centroid(
            node_coordinates,
            geo.cell_node_ids[cell_id],
            geo.cell_face_map[cell_id],
            geo.cell_faces_node_ids,
        )
        geo.cell_volumes[cell_id] = volume
        geo.cell_centroids[cell_id] = centroid
    end

    for cell_id in eachindex(geo.cell_node_ids)
        for face_idx in 1:geo.cell_n_faces[cell_id]
            packed_face_id = geo.cell_face_layout[cell_id, face_idx]
            face_node_ids = geo.cell_faces_node_ids[packed_face_id]
            face_area, face_normal, face_distance, face_centroid =
                FVMFramework._face_geometry(
                    node_coordinates,
                    face_node_ids,
                    geo.cell_centroids[cell_id],
                )

            geo.cell_face_areas[cell_id, face_idx] = face_area
            geo.cell_face_normals[cell_id, face_idx] = face_normal
            geo.cell_face_distances[cell_id, face_idx] = face_distance
            geo.cell_face_centroids[cell_id, face_idx] = face_centroid
        end
    end

    fill!(parent(geo.cell_neighbor_normals), zero(eltype(geo.cell_centroids)))
    fill!(parent(geo.cell_neighbor_distances), zero(eltype(geo.cell_volumes)))

    for (cell_id, neighbors) in geo.cell_neighbors
        for (neighbor_id, face_idx, _) in neighbors
            neighbor_direction =
                geo.cell_centroids[neighbor_id] - geo.cell_centroids[cell_id]
            neighbor_normal = geo.cell_face_normals[cell_id, face_idx]

            if dot(neighbor_normal, neighbor_direction) < 0
                neighbor_normal = -neighbor_normal
            end

            geo.cell_neighbor_normals[cell_id, face_idx] = neighbor_normal
            geo.cell_neighbor_distances[cell_id, face_idx] = norm(neighbor_direction)
        end
    end
    return nothing
end

"""
    fvm_operator!(du, u, p, t, system, geo_cache, solve_groups!)

Operator overload that materializes the correctly typed geometry view before
calling the normal system RHS. Pass an `ADGeometryCache` in place of the ordinary
`FVMGeometry` when constructing the ODE closure.
"""
function fvm_operator!(
        du_vec,
        u_vec,
        p,
        t,
        system::FVMSystem,
        geo_cache::ADGeometryCache,
        solve_groups!,
    )
    du, u = unpack_fvm_state(du_vec, u_vec, p, t, system)
    geo = get_ad_geometry(geo_cache, first(u_vec))
    solve_groups!(du, u, p, t, system, geo)
    return nothing
end
