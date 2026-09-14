function first_order_face_reconstruction!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    area, cell_face_normal, dist_to_face,
    cell_neighbor_normal, cell_neighbor_dist,
    cell_volume_a,
)
    return (
        u.density[idx_a],
        u.momentum_density_u[idx_a],
        u.momentum_density_v[idx_a],
        u.momentum_density_w[idx_a],
        u.volumetric_energy[idx_a]
    )
end