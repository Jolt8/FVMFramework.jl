
function update_region_groups!(du, u, p, t, geo, system)
    for reg in system.region_groups
        update_region_group!(
            du, u, p, t,
            reg.property_update_function!, reg.region_cells,
            geo.cell_volumes, system
        )
    end
end

function solve_connection_groups!(du, u, p, t, geo, system)
    for conn in system.connection_groups
        solve_connection_group!(
            du, u, p, t,
            conn.flux_function!, conn.cell_neighbors,
            geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
            geo.cell_volumes
        )
    end
end

function solve_patch_groups!(du, u, p, t, geo, system)
    for patch in system.patch_groups
        solve_patch_group!(
            du, u, p, t, 
            patch.patch_function!, patch.cell_neighbors,
            geo.cell_face_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
            geo.cell_volumes
        )
    end
end

function solve_region_groups!(du, u, p, t, geo, system)
    for reg in system.region_groups
        solve_region_group!(
            du, u, p, t,
            reg.region_function!, reg.region_cells,
            geo.cell_volumes
        )
    end
end