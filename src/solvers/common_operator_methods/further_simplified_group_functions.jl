function update_region_groups!(du, u, p, t, system, geo)
    for reg in system.region_groups
        update_region_group!(
            du, u, p, t, system, geo,
            reg.property_update_function!, reg.region_cells
        )
    end
end

function solve_connection_groups!(du, u, p, t, system, geo)
    for conn in system.connection_groups
        solve_connection_group!(
            du, u, p, t, system, geo,
            conn.flux_function!, conn.cell_neighbors
        )
    end
end

function solve_patch_groups!(du, u, p, t, system, geo)
    for patch in system.patch_groups
        solve_patch_group!(
            du, u, p, t, system, geo,
            patch.patch_function!, patch.cell_neighbors
        )
    end
end

function solve_region_groups!(du, u, p, t, system, geo)
    for reg in system.region_groups
        solve_region_group!(
            du, u, p, t, system, geo,
            reg.region_function!, reg.region_cells
        )
    end
end