function heat_diffusion!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / dist

    du.heat[idx_a] -= -k_effective * grad_T * area
end