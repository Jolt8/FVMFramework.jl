abstract type FVMGeometry end

struct FVMGeometryTetra{T, CoordType} <: FVMGeometry
    cell_volumes::Vector{T}
    cell_centroids::Vector{CoordType}
    cell_neighbors::Vector{Tuple{Int64, Vector{Tuple{Int64, Int64, Int64}}}}
    cell_neighbor_areas::Vector{MVector{4, T}}
    cell_neighbor_normals::Vector{MVector{4, CoordType}}
    cell_neighbor_distances::Vector{MVector{4, T}}
    unconnected_cell_face_map::Vector{Tuple{Int, Int}}
    cell_face_areas::Vector{MVector{4, T}}
    cell_face_normals::Vector{MVector{4, CoordType}}
    cell_face_distances::Vector{MVector{4, T}}
end

struct FVMGeometryHexa{T, CoordType} <: FVMGeometry
    cell_volumes::Vector{T}
    cell_centroids::Vector{CoordType}
    cell_neighbors::Vector{Tuple{Int64, Vector{Tuple{Int64, Int64, Int64}}}}
    cell_neighbor_areas::Vector{MVector{6, T}}
    cell_neighbor_normals::Vector{MVector{6, CoordType}}
    cell_neighbor_distances::Vector{MVector{6, T}}
    unconnected_cell_face_map::Vector{Tuple{Int, Int}}
    cell_face_areas::Vector{MVector{6, T}}
    cell_face_normals::Vector{MVector{6, CoordType}}
    cell_face_distances::Vector{MVector{6, T}}
end

function compress_geo_to_struct(
        grid::Grid{3, Tetrahedron, Float64}, top,
        cell_volumes, cell_centroids,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_face_areas, cell_face_normals, cell_face_distances,
        initial_node_coordinates, nodes_of_cells,
        cell_neighbors, cell_neighbors_node_ids,
        cell_face_map, map_respective_node_ids
    )

    rebuild_fvm_geometry_tetra!(
        #mutated vars
        cell_volumes, cell_centroids,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_face_areas, cell_face_normals, cell_face_distances,

        #non-mutated vars
        initial_node_coordinates, nodes_of_cells,
        cell_neighbors, cell_neighbors_node_ids,
        cell_face_map, map_respective_node_ids
    )

    unconnected_cell_face_map = get_unconnected_map(grid, top)

    return FVMGeometryTetra(
        cell_volumes, cell_centroids,
        cell_neighbors,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        unconnected_cell_face_map, cell_face_areas, cell_face_normals, cell_face_distances
    )
end

function compress_geo_to_struct(
        grid::Grid{3, Hexahedron, Float64}, top,
        cell_volumes, cell_centroids,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_face_areas, cell_face_normals, cell_face_distances, 
        initial_node_coordinates, nodes_of_cells,
        cell_neighbors, cell_neighbors_node_ids,
        cell_face_map, map_respective_node_ids
    )

    rebuild_fvm_geometry_hexa!(
        #mutated vars
        cell_volumes, cell_centroids,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_face_areas, cell_face_normals, cell_face_distances,

        #non-mutated vars
        initial_node_coordinates, nodes_of_cells,
        cell_neighbors, cell_neighbors_node_ids,
        cell_face_map, map_respective_node_ids
    )

    unconnected_cell_face_map = get_unconnected_map(grid, top)

    return FVMGeometryHexa(
        cell_volumes, cell_centroids,
        cell_neighbors,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        unconnected_cell_face_map, cell_face_areas, cell_face_normals, cell_face_distances
    )
end

function build_fvm_geo_into_struct(grid, top)
    initial_node_coordinates = get_node_coordinates(grid)
    nodes_of_cells = get_nodes_of_cells(grid)

    cell_neighbors, cell_neighbors_node_ids = get_cell_neighbors(grid, top)

    cell_face_map, map_respective_node_ids = get_cell_face_map(grid)

    CoordType = eltype(initial_node_coordinates)
    T = eltype(CoordType)

    n_cells = length(nodes_of_cells)
    cell_volumes = zeros(T, n_cells)
    cell_centroids = zeros(CoordType, n_cells)

    n_facets = nfacets(grid.cells[1])

    n_cell_neighbors = length(cell_neighbors)
    cell_neighbor_areas = [zero(MVector{n_facets, T}) for _ in 1:n_cell_neighbors]
    cell_neighbor_normals = [zero(MVector{n_facets, CoordType}) for _ in 1:n_cell_neighbors]
    cell_neighbor_distances = [zero(MVector{n_facets, T}) for _ in 1:n_cell_neighbors]

    cell_face_areas = [zero(MVector{n_facets, T}) for _ in 1:n_cells]
    cell_face_normals = [zero(MVector{n_facets, CoordType}) for _ in 1:n_cells]
    cell_face_distances = [zero(MVector{n_facets, T}) for _ in 1:n_cells]

    #we use the dynamic dispatch of compress_geo_to_struct to prevent having to check the grid's type with if statements
    #the reason this mutates is that it will be later used for geometry optimization which is one of the biggest applications of this framework
    return compress_geo_to_struct(
        grid, top,
        cell_volumes, cell_centroids,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_face_areas, cell_face_normals, cell_face_distances,
        initial_node_coordinates, nodes_of_cells,
        cell_neighbors, cell_neighbors_node_ids,
        cell_face_map, map_respective_node_ids
    )
end
