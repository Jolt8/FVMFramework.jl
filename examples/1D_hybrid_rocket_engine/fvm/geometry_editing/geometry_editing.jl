# Cell Volume:    V_i = 2π * r_cell * Δx * Δz = V_cart * (2π * r_cell / Δy)
# Axial Face Area: A_z = 2π * r_face * Δx     = A_cart * (2π * r_face / Δy)
# Radial Face Area: A_x = 2π * r_face * Δz     = A_cart * (2π * r_face / Δy)

function correct_hybrid_geo!(geo, grid)
    Δy = 0.001  # Extrusion depth = 1.0 mm = 0.001 m

    # 1. Scale cell volumes
    for cell_id in 1:length(geo.cell_volumes)
        r_cell = geo.cell_centroids[cell_id][1]  # X-coordinate is radius
        geo.cell_volumes[cell_id] *= (2π * r_cell / Δy)
    end

    # 2. Scale internal neighbor face areas
    for (cell_id, this_cell_neighbors) in geo.cell_neighbors
        for (neighbor_id, face_idx, _) in this_cell_neighbors
            normal = geo.cell_neighbor_normals[cell_id, face_idx]
            # Only scale in-plane faces (radial or axial, normal not along Y)
            if abs(normal[2]) < 0.5
                r_face = geo.cell_face_centroids[cell_id, face_idx][1]
                geo.cell_face_areas[cell_id, face_idx] *= (2π * r_face / Δy)
            end
        end
    end

    # 3. Scale boundary face areas (inlet, outlet, wall, symmetry)
    for (cell_id, face_idx) in geo.unconnected_cell_face_map
        normal = geo.cell_face_normals[cell_id, face_idx]
        if abs(normal[2]) < 0.5
            r_face = geo.cell_face_centroids[cell_id, face_idx][1]

            r_face = abs(r_face) < 1e-12 ? zero(r_face) : r_face
            @assert r_face >= 0
            
            geo.cell_face_areas[cell_id, face_idx] *= 2π * r_face / Δy
        end
    end
end
