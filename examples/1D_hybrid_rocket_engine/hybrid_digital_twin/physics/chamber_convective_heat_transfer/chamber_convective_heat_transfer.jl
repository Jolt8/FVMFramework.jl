function chamber_convective_heat_flux_custom_geometry!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b, face_area_a)    
    du.heat[idx_a] -= p.chamber_convective_heat_transfer_coefficient * face_area_a * (u.temperature[idx_a] - u.temperature[idx_b])
    du.heat[idx_b] += p.chamber_convective_heat_transfer_coefficient * face_area_a * (u.temperature[idx_a] - u.temperature[idx_b])
end

function chamber_convective_heat_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    (
        dist,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    chamber_convective_heat_flux_custom_geometry!(
        du, u, p, t, system, geo, 
        idx_a, face_a, idx_b, face_b, 
        face_area_a
    )
end