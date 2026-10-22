function create_sst_wall_distances(grid, geo, wall_patch_names)
    n_cells = length(geo.cell_volumes)
    face_centers = typeof(geo.cell_centroids[1])[]

    for wall_patch_name in wall_patch_names
        for (cell_id, face_id) in getfacetset(grid, wall_patch_name)
            face_node_ids = Ferrite.faces(grid.cells[cell_id])[face_id]
            face_center = sum(grid.nodes[node_id].x for node_id in face_node_ids) / length(face_node_ids)
            push!(face_centers, face_center)
        end
    end

    isempty(face_centers) && error("No wall faces were found in the mesh for the specified wall patches.")

    sst_wall_distances = zeros(Float64, n_cells)
    for cell_id in 1:n_cells
        sst_wall_distances[cell_id] = minimum(
            norm(geo.cell_centroids[cell_id] - fc)
            for fc in face_centers
        )
    end

    if any(isnan, sst_wall_distances) || any(isinf, sst_wall_distances) || any(x -> x < 0.0, sst_wall_distances)
        error("An unphysical wall distance was calculated.")
    end

    return sst_wall_distances
end
