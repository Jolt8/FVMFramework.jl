function calculate_hex_volume(p)
    c = sum(p) / 8.0
    
    # Faces of Hex8 (node indices for each face)
    faces = (
        (1, 4, 3, 2),  # bottom
        (1, 2, 6, 5),  # front
        (2, 3, 7, 6),  # right
        (3, 4, 8, 7),  # back
        (1, 5, 8, 4),  # left 
        (5, 6, 7, 8)   # top
    )
    
    total_vol = 0.0
    for face in faces
        node_1, node_2, node_3, node_4 = p[face[1]], p[face[2]], p[face[3]], p[face[4]]
        
        # Split quad face into 2 triangles, form tetrahedra with centroid
        total_vol += dot(node_1 - c, cross_product(node_2 - c, node_3 - c))
        total_vol += dot(node_1 - c, cross_product(node_3 - c, node_4 - c))
    end
    
    return abs(total_vol) / 6.0
end

function calculate_quad_face_area(face_node_coordiantes)
    node_1_coords = face_node_coordiantes[1]
    node_2_coords = face_node_coordiantes[2]
    node_3_coords = face_node_coordiantes[3]
    node_4_coords = face_node_coordiantes[4]

    cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
    cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
    
    area_vec_1 = 0.5 * cross_a
    area_vec_2 = 0.5 * cross_b

    total_area_vec = area_vec_1 + area_vec_2
    return norm(total_area_vec)
end

function rebuild_fvm_geometry_hexa!(
        #mutated vars
        cell_volumes, cell_centroids, 
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances, 
        cell_face_areas, cell_face_normals, cell_face_distances,

        #unmutated vars
        node_coordinates, nodes_of_cells,
        cell_neighbors, cell_neighbors_node_ids, 
        all_cell_face_map, map_respective_node_ids::Vector{NTuple{4, Int}}
    )
    CoordType = eltype(node_coordinates)
    T = eltype(CoordType)

    for cell_id in eachindex(nodes_of_cells)
        cell_nodes = nodes_of_cells[cell_id]

        n_nodes = length(cell_nodes)

        p = ntuple(8) do i
            @inbounds node_coordinates[cell_nodes[i]]
        end

        vol = calculate_hex_volume(p)
        
        cent = sum(p) / T(n_nodes)

        cell_volumes[cell_id] = vol
        cell_centroids[cell_id] = cent
    end

    for (cell_id, this_cell_neighbors) in cell_neighbors
        for (neighbor_id, face_idx) in this_cell_neighbors
            face_node_indices = cell_neighbors_node_ids[cell_id][face_idx] #cell_neighbors_node_ids[face_idx] looks like (1, 4, 7, 21) 
            node_1_coords = node_coordinates[face_node_indices[1]]
            node_2_coords = node_coordinates[face_node_indices[2]]
            node_3_coords = node_coordinates[face_node_indices[3]]
            node_4_coords = node_coordinates[face_node_indices[4]]

            #get_area
            cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
            cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
            
            area_vec_1 = 0.5 * cross_a
            area_vec_2 = 0.5 * cross_b

            total_area_vec = area_vec_1 + area_vec_2
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
        face_node_indices = map_respective_node_ids[i] #unconnected_map_respective_node_ids[i] looks like (1, 4, 7, 21) 
        node_1_coords = node_coordinates[face_node_indices[1]]
        node_2_coords = node_coordinates[face_node_indices[2]]
        node_3_coords = node_coordinates[face_node_indices[3]]
        node_4_coords = node_coordinates[face_node_indices[4]]

        #get_area
        cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
        cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
        
        area_vec_1 = 0.5 * cross_a
        area_vec_2 = 0.5 * cross_b

        total_area_vec = area_vec_1 + area_vec_2
        total_area = norm(total_area_vec)

        cell_face_areas[cell_id][face_idx] = total_area 

        dist_to_face_vec = (node_1_coords + node_2_coords + node_3_coords + node_4_coords) / 4 - cell_centroids[cell_id]

        cell_face_normals[cell_id][face_idx] = normalize(dist_to_face_vec)
        cell_face_distances[cell_id][face_idx] = norm(dist_to_face_vec)
    end
end