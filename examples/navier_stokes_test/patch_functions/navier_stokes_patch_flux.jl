
# Boundary patch flux
function wall_patch_flux_generic!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_volume_a,
    area, cell_face_normal, dist_to_face,
    cell_neighbor_normal, cell_neighbor_dist
)
    # Viscous flux at boundary (drag force)
    F_visc_u = u.nu[idx_a] * area * (u.u_bc[idx_a] - u.u_vel[idx_a]) / dist_to_face
    F_visc_v = u.nu[idx_a] * area * (u.v_bc[idx_a] - u.v_vel[idx_a]) / dist_to_face
    F_visc_w = u.nu[idx_a] * area * (u.w_bc[idx_a] - u.w_vel[idx_a]) / dist_to_face

    # Pressure force from boundary
    F_press_u = (1.0 / u.rho[idx_a]) * u.pressure[idx_a] * area * cell_face_normal[1]
    F_press_v = (1.0 / u.rho[idx_a]) * u.pressure[idx_a] * area * cell_face_normal[2]
    F_press_w = (1.0 / u.rho[idx_a]) * u.pressure[idx_a] * area * cell_face_normal[3]

    # Accumulate
    du.u_flow[idx_a] += F_visc_u - F_press_u
    du.v_flow[idx_a] += F_visc_v - F_press_v
    du.w_flow[idx_a] += F_visc_w - F_press_w
end