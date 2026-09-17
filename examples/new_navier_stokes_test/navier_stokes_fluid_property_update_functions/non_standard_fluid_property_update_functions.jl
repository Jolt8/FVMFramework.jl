function update_k_from_prandtl!(du, u, p, t, cell_id, vol)
    u.k[cell_id] = (u.mu[cell_id] * u.cp[cell_id]) / u.prandtl_number[cell_id]
end