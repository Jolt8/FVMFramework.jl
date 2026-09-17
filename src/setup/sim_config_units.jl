
function run_and_check_units(du0_vec_units, u0_vec_units, system, geo, du_unitful_cache_vec, u_unitful_cache_vec, properties_vec_units, p_vec_units)
    du_unitful_cache_vec = upreferred.(du_unitful_cache_vec)
    du_unitful_cache_vec .*= 0.0
    du_unitful_cache_vec = du_unitful_cache_vec ./ 1.0u"s"
    u_unitful_cache_vec = upreferred.(u_unitful_cache_vec)
    u_unitful_cache_vec .*= 0.0

    du0_vec_units .*= 0.0
    du0_vec_units = du0_vec_units ./ 1.0u"s"

    u = VirtualFVMArray((u0_vec_units, u_unitful_cache_vec, properties_vec_units), system.u_virtual_axes)
    du = VirtualFVMArray((du0_vec_units, du_unitful_cache_vec), system.du_virtual_axes)

    t = 0.0u"s"
    p = p_vec_units

    #applying units 
    cell_volumes = geo.cell_volumes .* u"m^3"
    cell_centroids = geo.cell_centroids .* u"m"
    
    cell_neighbor_areas = geo.cell_neighbor_areas .* u"m^2"
    cell_neighbor_normals = geo.cell_neighbor_normals .* u"m"
    cell_neighbor_distances = geo.cell_neighbor_distances .* u"m"

    cell_face_areas = geo.cell_face_areas .* u"m^2"
    cell_face_normals = geo.cell_face_normals .* u"m"

    for reg in system.region_groups
        update_region_group!(
            du, u, p, t,
            reg.property_update_function!, 
            reg.region_cells,
            cell_volumes
        )
    end

    for conn in system.connection_groups
        solve_connection_group!(
            du, u, p, t,
            conn.flux_function!, conn.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
            cell_volumes
        )
    end

    for patch in system.patch_groups
        solve_patch_group!(
            du, u, p, t,
            patch.patch_function!, patch.cell_neighbors,
            cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
            cell_volumes
        )
    end
    
    for reg in system.region_groups
        solve_region_group!(
            du, u, p, t,
            reg.region_function!, reg.region_cells,
            cell_volumes
        )
    end

    return du, u
end

#just an idea, but we could also run a separate run-through of the solver that doesn't contain units 
#and then at the end strip the units from the unitful run-through and compare the two
#while I'm pretty sure that upreferred would take care of most of this, it would make you extra
#confident that you're not missing any unit conversions

