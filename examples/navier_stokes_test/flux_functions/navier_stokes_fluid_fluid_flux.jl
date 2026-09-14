# Connection flux (internal faces)
function fluid_fluid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    area = cell_neighbor_areas[idx_a][face_idx]
    norm = cell_neighbor_normals[idx_a][face_idx]
    dist = cell_neighbor_distances[idx_a][face_idx]

    # Velocity at cells a and b
    v_a = Ferrite.Vec{3}((u.u_vel[idx_a], u.v_vel[idx_a], u.w_vel[idx_a]))
    v_b = Ferrite.Vec{3}((u.u_vel[idx_b], u.v_vel[idx_b], u.w_vel[idx_b]))

    # Average velocity at face
    v_f = 0.5 * (v_a + v_b)

    rho_eff = 0.5 * (u.rho[idx_a] + u.rho[idx_b])

    # Volumetric flow rate across face from a to b with Rhie-Chow pressure stabilization
    # We subtract a pressure difference term to couple pressure and velocity and prevent checkerboarding
    chi = 0.05 / rho_eff  # stabilization coefficient (can be tuned, typical range is 0.01 - 0.1)
    dp_dn = (u.pressure[idx_b] - u.pressure[idx_a]) / dist
    Q_f = (dot(v_f, norm) - chi * dp_dn) * area

    # Convective flux of velocity using upwind scheme
    if Q_f > 0.0
        F_conv_u = Q_f * u.u_vel[idx_a]
        F_conv_v = Q_f * u.v_vel[idx_a]
        F_conv_w = Q_f * u.w_vel[idx_a]
    else
        F_conv_u = Q_f * u.u_vel[idx_b]
        F_conv_v = Q_f * u.v_vel[idx_b]
        F_conv_w = Q_f * u.w_vel[idx_b]
    end

    # Pressure force on cell a (using physical density)
    p_f = 0.5 * (u.pressure[idx_a] + u.pressure[idx_b])

    F_press_u = (1.0 / rho_eff) * p_f * area * norm[1]
    F_press_v = (1.0 / rho_eff) * p_f * area * norm[2]
    F_press_w = (1.0 / rho_eff) * p_f * area * norm[3]

    # Viscous flux from cell b to cell a
    nu_eff = 0.5 * (u.nu[idx_a] + u.nu[idx_b])

    F_visc_u = nu_eff * area * (u.u_vel[idx_b] - u.u_vel[idx_a]) / dist
    F_visc_v = nu_eff * area * (u.v_vel[idx_b] - u.v_vel[idx_a]) / dist
    F_visc_w = nu_eff * area * (u.w_vel[idx_b] - u.w_vel[idx_a]) / dist

    # Continuity flux (compressibility)
    beta_eff = 0.5 * (u.beta[idx_a] + u.beta[idx_b])
    F_mass = beta_eff * Q_f

    # Accumulate into du for cell a
    du.u_flow[idx_a] -= F_conv_u + F_press_u - F_visc_u
    du.v_flow[idx_a] -= F_conv_v + F_press_v - F_visc_v
    du.w_flow[idx_a] -= F_conv_w + F_press_w - F_visc_w
    du.p_flow[idx_a] -= F_mass
end