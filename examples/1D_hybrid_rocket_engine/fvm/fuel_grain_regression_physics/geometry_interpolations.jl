function get_min_max_x_fuel_grain_bounds(cells, node_coordinates)
    cells_node_coordinates = []
    for cell_id in cells 
        append!(cells_node_coordinates, [node_coordinates[node_id][3] for node_id in grid.cells[cell_id].nodes])
    end

    sort!(cells_node_coordinates)

    min_x_fuel_grain_bound = cells_node_coordinates[1]
    max_x_fuel_grian_bound = cells_node_coordinates[end]

    return min_x_fuel_grain_bound, max_x_fuel_grian_bound
end

function precompute_fuel_grain_regression_geometry_interpolations(
    grid, geo, 
    dominant_fuel_grain_cells, 
    fuel_grain_subsidiary_cell_ids_to_dominant_cell_id_map,
    fuel_grain_dominant_cell_id_to_subsidiaries_map,
    fluid_domain_subsidiary_cell_ids_to_dominant_cell_id_map,
    fluid_domain_dominant_cell_id_to_subsidiaries_map,
    dominant_fuel_grain_closest_fluid_domain_cell_map, #access with [fuel_grain_dominant_cell_id] --> closest_fluid_domain_cell_id
)
    node_coordinates = get_node_coordinates(grid)

    min_x_fuel_grain_node = 0
    max_x_fuel_grian_node = 0

    first_dominant_cell_id = dominant_fuel_grian_cells[1]
    first_subsidiary_cells = fuel_grain_dominant_cell_id_to_subsidiaries_map[first_dominant_cell_id]
    first_fluid_cells = fluid_domain_dominant_cell_id_to_subsidiaries_map[first_dominant_cell_id]

    min_x_fuel_grain_bound, max_x_fuel_grian_bound = get_min_max_x_fuel_grain_bounds(vcat(first_dominant_cell_id, first_subsidiary_cells...), node_coordinates)
    min_x_fluid_bound, max_x_fluid_bound = get_min_max_x_fuel_grain_bounds(first_fluid_cells, node_coordinates)
    #max_x_fuel_grain_bound and min_x_fluid_bound should be the same

    #since every cell is a rectangular prism with 90* corners, we can use some simple logic here to get the interpolations we need
    #However, the geometry of these cells have been modified previously with geometry_editing so this might be more challenging that I think
    #A lazy way would just be to rerun update_fvm_geometry!(geo, node_coordinates)
    #then rerun the geometry_editing_script
    #this might actually be the most efficient way to do this since it's really easy and probably quite fast
    #yep, we're going with this

    #update_fvm_geometry!(geo, node_coordinates)
    #correct_hybrid_geo(geo, grid)

    #interpolations we need (for all cells)
    cell_volumes = Vector{DataInterpolation.LinearInterpolation}(undef, n_cells)
    cell_centroids = Vector{SVector{3, DataInterpolation.LinearInterpolation}}(undef, n_cells)

    #   Neighbor information, stored in the same half-face layout
    cell_neighbor_normals = Vector{SVector{6, SVector{3, DataInterpolation.LinearInterpolation}}}(undef, n_cells) #although we have wedges, we'll use 6 faces per cell for simplicity
    cell_neighbor_distances = Vector{CellFaceVector{DataInterpolation.LinearInterpolation, Vector{DataInterpolation.LinearInterpolation}}}(undef, n_cells)
    
    #   Face information
    cell_face_areas = Vector{CellFaceVector{DataInterpolation.LinearInterpolation, Vector{DataInterpolation.LinearInterpolation}}}(undef, n_cells)
    cell_face_normals = Vector{CellFaceVector{SVector{3, DataInterpolation.LinearInterpolation}, Vector{SVector{3, DataInterpolation.LinearInterpolation}}}}(undef, n_cells)
    cell_face_distances = Vector{CellFaceVector{DataInterpolation.LinearInterpolation, Vector{DataInterpolation.LinearInterpolation}}}(undef, n_cells)
    cell_face_centroids = Vector{CellFaceVector{SVector{3, DataInterpolation.LinearInterpolation}, Vector{SVector{3, DataInterpolation.LinearInterpolation}}}}(undef, n_cells)

    cell_volumes_over_steps = []
    cell_centroids_over_steps = []

    cell_neighbor_normals_over_steps = []
    cell_neighbor_distances_over_steps = []

    cell_face_areas_over_steps = []
    cell_face_normals_over_steps = []
    cell_face_distances_over_steps = []
    cell_face_centroids_over_steps = []

    n_steps = 5

    starting_node_coordinates = deepcopy(node_coordinates)
    node_coordinates_over_steps = []

    for i in 1:n_steps
        for dominant_cell_id in dominant_fuel_grain_cells #we're iterating through each 'row' of fuel grain and fluid_domain cells
            sorted_closest_to_wall_fluid_domain_cells = []
            sorted_closest_to_wall_fluid_domain_node_ids = []

            last_fluid_domain_cell_nodes_x = min_x_fluid_bound

            for fluid_domain_cell_id in fluid_domain_dominant_cell_id_to_subsidiaries_map[dominant_cell_id]
                fluid_domain_cell_nodes_x = [node_coordinates[node_id][1] for node_id in grid.cells[fluid_domain_cell_id].nodes]

                
                push!(sorted_closest_to_wall_fluid_domain_node_ids, [node_coordinates[node_id][1] for node_id in grid.cells[fluid_domain_cell_id].nodes])
                
                
                push!(sorted_closest_to_wall_fluid_domain_cells, fluid_domain_cell_id)
            end
        end
    end
end

    

