cell_face_layout::CellFaceLayout

    #Geometry information
    #   Cell-centered information
    cell_volumes::Vector{T}
    cell_centroids::Vector{CoordType}

    #   Neighbor information, stored in the same half-face layout
    cell_neighbor_normals::CellFaceVector{CoordType, Vector{CoordType}}
    cell_neighbor_distances::CellFaceVector{T, Vector{T}}

    #   Face information
    cell_face_areas::CellFaceVector{T, Vector{T}}
    cell_face_normals::CellFaceVector{CoordType, Vector{CoordType}}
    cell_face_distances::CellFaceVector{T, Vector{T}}
    cell_face_centroids::CellFaceVector{CoordType, Vector{CoordType}}

    # Node information
    node_coordinates::Vector{CoordType}
end

struct FVMGeometryAD{T, CoordType, DualType <: Real, DualCoordType <: SVector{3, DualType}}
    base_geo::FVMGeometry{T, CoordType}

    cell_face_layout::CellFaceLayout

    cell_volumes_diff_cache::DiffCache
    cell_volumes_values::Vector{T}
    cell_volumes::Vector{DualType}

    cell_centroids_diff_cache::DiffCache
    cell_centroids_values::Vector{CoordType}
    cell_centroids::Vector{DualCoordType}

    cell_neighbor_normals_diff_cache::DiffCache
    cell_neighbor_normals_values::CellFaceVector{CoordType, Vector{CoordType}}
    cell_neighbor_normals::CellFaceVector{DualCoordType, Vector{DualCoordType}}

    cell_neighbor_distances_diff_cache::DiffCache
    cell_neighbor_distances_values::CellFaceVector{T, Vector{T}}
    cell_neighbor_distances::CellFaceVector{DualType, Vector{DualType}}

    cell_face_areas_diff_cache::DiffCache
    cell_face_areas_values::CellFaceVector{T, Vector{T}}
    cell_face_areas::CellFaceVector{DualType, Vector{DualType}}

    cell_face_normals_diff_cache::DiffCache
    cell_face_normals_values::CellFaceVector{CoordType, Vector{CoordType}}
    cell_face_normals::CellFaceVector{CoordType, Vector{DualCoordType}}

    cell_face_distances_diff_cache::DiffCache
    cell_face_distances_values::CellFaceVector{T, Vector{T}}
    cell_face_distances::CellFaceVector{DualType, Vector{DualType}}

    cell_face_centroids_diff_cache::DiffCache
    cell_face_centroids_values::CellFaceVector{CoordType, Vector{CoordType}}
    cell_face_centroids::CellFaceVector{DualCoordType, Vector{DualCoordType}}

    node_coordinates_diff_cache::DiffCache
    node_coordinates_values::Vector{CoordType}
    node_coordinates::Vector{DualCoordType}

    field_syms::Vector{Symbol}
end

function Base.getproperty(geo::FVMGeometryAD, sym::Symbol)
    # 1. Fall back to FVMSystem's own fields
    if sym in fieldnames(FVMGeometryAD)
        return getfield(geo, sym)
    end
    
    # 2. Otherwise look inside the additional_data NamedTuple
    base_geo = getfield(geo, :base_geo)
    if hasfield(typeof(base_geo), sym)
        return getfield(base_geo, sym)
    end
    
    return getfield(geo, sym) # triggers standard error
end

function create_diff_cache_geo(geo, grid, N)
    recorded_vectors = []
    field_syms = Symbol[]

    geometry_field_names = [
        :cell_volumes,
        :cell_centroids,
        :cell_neighbor_normals,
        :cell_neighbor_distances,
        :cell_face_areas,
        :cell_face_normals,
        :cell_face_distances,
        :cell_face_centroids,
        :node_coordinates
    ]

    for field_sym in geometry_field_names
        @show field_sym
        @show typeof(getproperty(geo, field_sym))

        diff_cache = DiffCache(getproperty(geo, field_sym), N)
        values = getproperty(geo, field_sym)

        T = eltype(geo.cell_volumes)
        CoordType = eltype(geo.cell_centroids)
        DualType = ForwardDiff.Dual{Nothing, T, N}
        DualCoordType = SVector{3, DualType}

        n_cells = length(geo.cell_volumes)

        if values isa Vector{T}
            dual_values = zeros(DualType, length(values))
        elseif values isa Vector{CoordType}
            dual_values = zeros(DualCoordType, length(values))
        elseif values isa CellFaceVector{CoordType, Vector{CoordType}}
            dual_values = cell_face_layout_zeros(DualCoordType, geo.cell_face_layout)
        elseif values isa CellFaceVector{T, Vector{T}}
            dual_values = cell_face_layout_zeros(DualType, geo.cell_face_layout)
        else 
            error("Unsupported type")
        end

        push!(recorded_vectors, diff_cache)
        push!(recorded_vectors, values)
        push!(recorded_vectors, dual_values)
        push!(field_syms, field_sym)
    end

    @show length([geo, geo.cell_face_layout, recorded_vectors..., field_syms])
    @show fieldnames(FVMGeometryAD)

    return FVMGeometryAD(
        geo, 
        geo.cell_face_layout,
        recorded_vectors...,
        field_syms,
    )
end

function regenerate_ad_geo!(u, ad_geo::FVMGeometryAD)
    for field_sym in ad_geo.field_syms
        diff_cache_sym = Symbol(string(field_sym) * "_diff_cache")
        values_sym = Symbol(string(field_sym) * "_values")

        ad_cache = get_tmp(getproperty(ad_geo, diff_cache_sym), first(u))
        ad_cache .= getproperty(ad_geo, values_sym)
        
        getproperty(ad_geo, field_sym) .= ad_cache
    end
end

        




