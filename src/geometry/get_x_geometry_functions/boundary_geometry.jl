
"""
Typical usage:
face_area, face_normal, face_distance, vol = boundary_geometry(geo, cell_id, face_idx)
"""
function boundary_geometry(geo, cell_id, face_idx)
    return (
        geo.cell_face_areas[cell_id][face_idx], geo.cell_face_normals[cell_id][face_idx], geo.cell_face_distances[cell_id][face_idx], geo.cell_volumes[cell_id]
    )
end