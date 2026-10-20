function heat_diffusion!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b
)
    (
        distance,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / distance

    du.heat[idx_a] -= -k_effective * grad_T * face_area_a
    du.heat[idx_b] += -k_effective * grad_T * face_area_b
end