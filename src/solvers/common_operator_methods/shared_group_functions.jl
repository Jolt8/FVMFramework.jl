function update_region_group!(
    du, u, p, t, system, geo,
    property_update_function!::F, region_cells
) where {F}
    for cell_id in region_cells
        property_update_function!(
            du, u, p, t, system, geo,
            cell_id
        )
    end
end

function solve_connection_group!(
    du, u, p, t, system, geo,
    flux!::F, cell_neighbors
) where {F}
    for (idx_a, neighbor_list) in cell_neighbors
        for (idx_b, face_idx_a, face_idx_b) in neighbor_list
            flux!(
                du, u, p, t, system, geo,
                idx_a, face_idx_a,
                idx_b, face_idx_b
            )
        end
    end
end

function solve_patch_group!(
    du, u, p, t, system, geo,
    patch_physics!::F, cell_neighbors
) where {F} 
    for (idx_a, neighbor_list) in cell_neighbors
        for (idx_b, face_idx_a, face_idx_b) in neighbor_list
            patch_physics!(
                du, u, p, t, system, geo,
                idx_a, face_idx_a,
                idx_b, face_idx_b
            )
        end
    end
end

function solve_region_group!(
    du, u, p, t, system, geo,
    region_function!::G, region_cells
) where {G}
    for cell_id in region_cells
        region_function!(
            du, u, p, t, system, geo,
            cell_id
        )
    end
end