
function create_sst_wall_distances(grid, geo, wall_patch_names)
    wall_distances = []

    for wall_patch_name in wall_patch_names
        face_centers = typeof(geo.cell_centroids[1])[]
        for (cell_id, face_id) in getfacetset(grid, wall_patch_name)
            face_node_ids = Ferrite.faces(grid.cells[cell_id])[face_id]
            face_center = sum([grid.nodes[node_id].x for node_id in face_node_ids]) /
            length(face_node_ids)
            push!(face_centers, face_center)
        end
    end

    isempty(airfoil_face_centers) && error("No fluid-airfoil interface was found in the mesh")
    sst_wall_distances = fill(1.0, n_cells)
    for cell_id in eachindex(geo.cell_centroids)
        sst_wall_distances[cell_id] = minimum(
            norm(geo.cell_centroids[cell_id] - face_center)
            for face_center in airfoil_face_centers
        )
    end

    return sst_wall_distances
end