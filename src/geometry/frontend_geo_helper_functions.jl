face_fill(value, geo::FVMGeometry) = [fill(value, N) for N in geo.cell_n_faces]

face_zeros(::Type{T}, geo::FVMGeometry) where {T} = [zeros(T, N) for N in geo.cell_n_faces]

face_zeros(geo::FVMGeometry) = [zeros(Float64, N) for N in geo.cell_n_faces]