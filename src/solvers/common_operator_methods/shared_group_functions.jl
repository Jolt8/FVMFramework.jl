function solve_connection_group!(
    du, u, p, t,
    flux!::F, cell_neighbors,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
) where {F}

    for (idx_a, neighbor_list) in cell_neighbors
        for (idx_b, face_idx) in neighbor_list
            flux!(
                du, u, p, t, 
                idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                cell_volumes
            )
        end
    end
end

function update_region_group!(
    du, u, p, t,
    property_update_function!::F, region_cells,
    cell_volumes, system
) where {F}
    for cell_id in region_cells
        property_update_function!(
            du, u, p, t, 
            cell_id,
            cell_volumes[cell_id],
            system
        )
    end
end

function solve_region_group!(
    du, u, p, t,
    region_function!::G, region_cells,
    cell_volumes
) where {G}
    for cell_id in region_cells
        region_function!(
            du, u, p, t, 
            cell_id,
            cell_volumes[cell_id]
        )
    end
end

function solve_patch_group!(
    du, u, p, t,
    patch_physics!::F, cell_neighbors,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
) where {F} 
    for (idx_a, neighbor_list) in cell_neighbors
        for (idx_b, face_idx) in neighbor_list
            patch_physics!(
                du, u, p, t,
                idx_a, idx_b, face_idx,
                cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
                cell_volumes
            )
        end
    end
end