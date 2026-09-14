#COMMITING BEFORE SWITCHING TO A CELL CENTRIC CONNECTION MAP

function calculate_tetra_volume(p)
    volume = (1 / 6) * abs(dot(p[1] - p[4], cross_product((p[2] - p[4]), (p[3] - p[4]))))

    return volume
end

function calculate_tri_face_area(face_node_coordiantes)
    node_1_coords = face_node_coordiantes[1]
    node_2_coords = face_node_coordiantes[2]
    node_3_coords = face_node_coordiantes[3]

    total_area_vec = 0.5 * cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)

    return norm(total_area_vec)
end

function rebuild_fvm_geometry_tetra!(
        #mutated vars
        cell_volumes, cell_centroids, 
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances, 
        cell_face_areas, cell_face_normals, cell_face_distances,

        #unmutated vars
        node_coordinates, nodes_of_cells,
        cell_neighbors, cell_neighbors_node_ids,
        all_cell_face_map, map_respective_node_ids::Vector{NTuple{3, Int}}
    )
    CoordType = eltype(node_coordinates)
    T = eltype(CoordType)

    for cell_id in eachindex(nodes_of_cells)
        cell_nodes = nodes_of_cells[cell_id]

        n_nodes = length(cell_nodes)

        p = ntuple(4) do i
            @inbounds node_coordinates[cell_nodes[i]]
        end

        vol = calculate_tetra_volume(p)

        cent = sum(p) / T(n_nodes)

        cell_volumes[cell_id] = vol
        cell_centroids[cell_id] = cent
    end

    #for geometry optimization in the future, we might want to store a separate version of cell_neighbors that doesn't contain duplicates
    #this would only require making normals negative for the neighboring cell and would reduce the math needed
    for (cell_id, this_cell_neighbors) in cell_neighbors
        for (neighbor_id, face_idx) in this_cell_neighbors
            face_node_indices = cell_neighbors_node_ids[cell_id][face_idx] #cell_neighbors_node_ids[face_idx] looks like (1, 4, 7) 
            node_1_coords = node_coordinates[face_node_indices[1]]
            node_2_coords = node_coordinates[face_node_indices[2]]
            node_3_coords = node_coordinates[face_node_indices[3]]

            #get area
            total_area_vec = 0.5 * cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)

            total_area = norm(total_area_vec)

            #get normal
            #make sure the cell's normal is actually pointing away
            vec_AB = cell_centroids[neighbor_id] - cell_centroids[cell_id]
            if dot(total_area_vec, vec_AB) < 0
                total_area_vec = -total_area_vec
            end

            cell_normal = normalize(total_area_vec)

            #get distance
            dist = norm(cell_centroids[cell_id] - cell_centroids[neighbor_id])

            cell_neighbor_areas[cell_id][face_idx] = total_area
            cell_neighbor_normals[cell_id][face_idx] = cell_normal
            cell_neighbor_distances[cell_id][face_idx] = dist
        end
    end

    for (i, (cell_id, face_idx)) in enumerate(all_cell_face_map)
        face_node_indices = map_respective_node_ids[i] #unconnected_map_respective_node_ids[i] looks like (1, 4, 7) 
        node_1_coords = node_coordinates[face_node_indices[1]]
        node_2_coords = node_coordinates[face_node_indices[2]]
        node_3_coords = node_coordinates[face_node_indices[3]]

        total_area_vec = 0.5 * cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)

        total_area = norm(total_area_vec)

        cell_face_areas[cell_id][face_idx] = total_area

        dist_to_face_vec = (node_1_coords + node_2_coords + node_3_coords) / 3 - cell_centroids[cell_id]

        cell_face_normals[cell_id][face_idx] = normalize(dist_to_face_vec)
        cell_face_distances[cell_id][face_idx] = norm(dist_to_face_vec)
    end
end