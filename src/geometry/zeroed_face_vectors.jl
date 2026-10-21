function zeroed_face_vectors(geo::FVMGeometry)
    return fill!(deepcopy(geo.cell_face_areas), 0.0)
end

function fill_face_vectors(geo::FVMGeometry, value)
    return fill!(deepcopy(geo.cell_face_areas), value)
end
    