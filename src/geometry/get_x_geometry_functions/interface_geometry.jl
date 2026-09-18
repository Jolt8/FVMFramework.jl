
"""
Typical usage:
(
    dist,
    face_area_a, face_normal_a, face_distance_a, vol_a,
    face_area_b, face_normal_b, face_distance_b, vol_b
) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)
"""
function interface_geometry(geo, idx_a, face_a, idx_b, face_b)
    return (
        geo.cell_neighbor_distances[idx_a][face_a],
        geo.cell_face_areas[idx_a][face_a], geo.cell_face_normals[idx_a][face_a], geo.cell_face_distances[idx_a][face_a], geo.cell_volumes[idx_a],
        geo.cell_face_areas[idx_b][face_b], geo.cell_face_normals[idx_b][face_b], geo.cell_face_distances[idx_b][face_b], geo.cell_volumes[idx_b]
    )
end