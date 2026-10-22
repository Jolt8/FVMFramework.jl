struct FVMGeometry{T, CoordType}
    #Topology information
    #   Cell Node Information
    cell_node_ids::Vector{Vector{Int}}

    #   Each neighbor entry is (neighbor cell, this local face, neighbor local face).
    cell_neighbors::Vector{Tuple{Int, Vector{NTuple{3, Int}}}}
    cell_neighbors_node_ids::Vector{Vector{Vector{Int}}}

    #   `cell_face_map[cell][local_face]` is the packed half-face index.
    cell_face_map::Vector{Vector{Int}}
    cell_faces_node_ids::Vector{Vector{Int}}
    cell_n_faces::Vector{Int}
    unconnected_cell_face_map::Vector{Tuple{Int, Int}}

    #convenience variables
    n_cells::Int

    #   Accessing faces
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
end

"""
    generate_fvm_topology_and_empty_geometry(grid, top; node_coordinates)

Build the mesh topology and allocate zeroed geometry arrays. Faces are stored as
cell-local half-faces: an interface shared by two cells therefore has two packed
entries, preserving Ferrite's local face numbering on both sides.

Returns `geo`.
"""
function generate_fvm_geometry(
        grid,
        top;
        node_coordinates = get_node_coordinates(grid),
    )
    n_cells = length(grid.cells)
    cell_node_ids = [collect(cell.nodes) for cell in grid.cells]
    cell_n_faces = [nfacets(cell) for cell in grid.cells]
    cell_face_layout = CellFaceLayout(cell_n_faces)

    cell_neighbors = [
        (cell_id, NTuple{3, Int}[])
        for cell_id in 1:n_cells
    ]
    cell_neighbors_node_ids = [
        [Int[] for _ in 1:cell_n_faces[cell_id]]
        for cell_id in 1:n_cells
    ]

    cell_face_map = [collect(cell_face_layout[cell_id]) for cell_id in 1:n_cells]
    cell_faces_node_ids = Vector{Vector{Int}}(undef, total_faces(cell_face_layout))
    unconnected_cell_face_map = Tuple{Int, Int}[]

    for cell_id in 1:n_cells
        cell_faces = Ferrite.faces(grid.cells[cell_id])

        for face_idx in 1:cell_n_faces[cell_id]
            face_node_ids = collect(cell_faces[face_idx])
            packed_face_id = cell_face_layout[cell_id, face_idx]
            cell_faces_node_ids[packed_face_id] = face_node_ids

            neighbor_info = top.face_face_neighbor[cell_id, face_idx]
            if isempty(neighbor_info)
                push!(unconnected_cell_face_map, (cell_id, face_idx))
                continue
            end

            neighbor_id, neighbor_face_idx = first(neighbor_info).idx
            push!(
                cell_neighbors[cell_id][2],
                (neighbor_id, face_idx, neighbor_face_idx),
            )
            cell_neighbors_node_ids[cell_id][face_idx] = copy(face_node_ids)
        end
    end

    CoordType = eltype(node_coordinates)
    T = eltype(CoordType)
    geo = FVMGeometry(
        #Topology information
        cell_node_ids,
        cell_neighbors,
        cell_neighbors_node_ids,
        cell_face_map,
        cell_faces_node_ids,
        cell_n_faces,
        unconnected_cell_face_map,
        n_cells,

        #Geometry Information
        cell_face_layout,
        zeros(T, n_cells),
        zeros(CoordType, n_cells),
        cell_face_layout_zeros(CoordType, cell_face_layout),
        cell_face_layout_zeros(T, cell_face_layout),
        cell_face_layout_zeros(T, cell_face_layout),
        cell_face_layout_zeros(CoordType, cell_face_layout),
        cell_face_layout_zeros(T, cell_face_layout),
        cell_face_layout_zeros(CoordType, cell_face_layout),
    )

    return geo
end

@inline function _cross_product(a, b)
    return SVector(
        a[2] * b[3] - a[3] * b[2],
        a[3] * b[1] - a[1] * b[3],
        a[1] * b[2] - a[2] * b[1],
    )
end

function _cell_volume_and_centroid(
        node_coordinates,
        cell_node_ids,
        packed_face_ids,
        cell_faces_node_ids,
    )
    CoordType = eltype(node_coordinates)
    T = eltype(CoordType)

    reference_point = zero(CoordType)
    for node_id in cell_node_ids
        reference_point += node_coordinates[node_id]
    end
    reference_point /= length(cell_node_ids)

    volume = zero(T)
    centroid_moment = zero(CoordType)

    # Each polygonal face is triangulated as a fan. Together with an interior
    # reference point, those triangles partition every supported convex cell.
    for packed_face_id in packed_face_ids
        face_node_ids = cell_faces_node_ids[packed_face_id]
        point_a = node_coordinates[first(face_node_ids)]

        for face_node_idx in 2:(length(face_node_ids) - 1)
            point_b = node_coordinates[face_node_ids[face_node_idx]]
            point_c = node_coordinates[face_node_ids[face_node_idx + 1]]

            tetra_volume = abs(dot(
                point_a - reference_point,
                _cross_product(point_b - reference_point, point_c - reference_point),
            )) / 6
            tetra_centroid = (reference_point + point_a + point_b + point_c) / 4

            volume += tetra_volume
            centroid_moment += tetra_volume * tetra_centroid
        end
    end

    return volume, centroid_moment / volume
end

function _face_geometry(node_coordinates, face_node_ids, cell_centroid)
    CoordType = eltype(node_coordinates)
    T = eltype(CoordType)

    point_a = node_coordinates[first(face_node_ids)]
    area_vector = zero(CoordType)
    surface_area = zero(T)
    centroid_moment = zero(CoordType)

    for face_node_idx in 2:(length(face_node_ids) - 1)
        point_b = node_coordinates[face_node_ids[face_node_idx]]
        point_c = node_coordinates[face_node_ids[face_node_idx + 1]]
        triangle_area_vector = _cross_product(point_b - point_a, point_c - point_a) / 2
        triangle_area = norm(triangle_area_vector)

        area_vector += triangle_area_vector
        surface_area += triangle_area
        centroid_moment += triangle_area * ((point_a + point_b + point_c) / 3)
    end

    face_centroid = centroid_moment / surface_area
    distance_to_face = face_centroid - cell_centroid

    if dot(area_vector, distance_to_face) < 0
        area_vector = -area_vector
    end

    face_area = norm(area_vector)
    raw_normal = area_vector / face_area
    clean_normal = map(component -> abs(component) < 1e-12 ? zero(component) : component, raw_normal)
    face_normal = normalize(clean_normal)
    face_distance = abs(dot(distance_to_face, face_normal))

    return face_area, face_normal, face_distance, face_centroid
end

"""
    update_fvm_geometry!(geo, fvm_topology, node_coordinates)

Recompute all cell and face geometry in place. This is the entry point intended
for future geometry optimization after `node_coordinates` changes.
"""
function update_fvm_geometry!(geo::FVMGeometry, node_coordinates)
    for cell_id in eachindex(geo.cell_node_ids)
        volume, centroid = _cell_volume_and_centroid(
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
            face_area, face_normal, face_distance, face_centroid = _face_geometry(
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

    return geo
end

"""
    build_fvm_geo_into_struct(grid, top)

Construct and populate both halves of the new geometry representation. Returns
`geo`.
"""
function build_fvm_geo_into_struct(grid, top)
    node_coordinates = get_node_coordinates(grid)
    geo = generate_fvm_geometry(grid, top; node_coordinates)
    update_fvm_geometry!(geo, node_coordinates)
    return geo
end
