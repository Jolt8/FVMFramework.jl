function viscous_fluid_and_diffusive_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    area, cell_face_normal, dist_to_face,
    cell_neighbor_normal, cell_neighbor_dist,
    vol
)
    # ------------------------------------------------------------
    # 1. Interpolate primitive quantities/properties to the face
    # ------------------------------------------------------------

    # Velocity should NOT use a harmonic mean.
    # A simple arithmetic interpolation is appropriate here.
    #TODO: use a better interpolation here in the future if necessary
    vel_u_face = 0.5 * (u.vel_u[idx_a] + u.vel_u[idx_b])
    vel_v_face = 0.5 * (u.vel_v[idx_a] + u.vel_v[idx_b])
    vel_w_face = 0.5 * (u.vel_w[idx_a] + u.vel_w[idx_b])

    # Harmonic averaging is reasonable for diffusive transport
    # coefficients, especially if properties vary spatially.
    mu_face = harmonic_mean(u.mu[idx_a], u.mu[idx_b])
    k_face  = harmonic_mean(u.k[idx_a],  u.k[idx_b])

    # ------------------------------------------------------------
    # 2. Interpolate cell-centered gradients to the face
    # ------------------------------------------------------------

    grad_u_face = 0.5 * (
        u.grad_vel_u[idx_a] +
        u.grad_vel_u[idx_b]
    )

    grad_v_face = 0.5 * (
        u.grad_vel_v[idx_a] +
        u.grad_vel_v[idx_b]
    )

    grad_w_face = 0.5 * (
        u.grad_vel_w[idx_a] +
        u.grad_vel_w[idx_b]
    )

    grad_T_face = 0.5 * (
        u.grad_temperature[idx_a] +
        u.grad_temperature[idx_b]
    )

    # ------------------------------------------------------------
    # 3. Velocity divergence at the face
    #
    # ∇⋅v = ∂u/∂x + ∂v/∂y + ∂w/∂z
    # ------------------------------------------------------------

    div_v =
        grad_u_face[1] +
        grad_v_face[2] +
        grad_w_face[3]

    # ------------------------------------------------------------
    # 4. Newtonian viscous stress tensor
    #
    # τ = μ[∇v + (∇v)' - (2/3)(∇⋅v)I]
    # ------------------------------------------------------------

    tau_xx =
        2 * mu_face * grad_u_face[1] -
        (2 / 3) * mu_face * div_v

    tau_yy =
        2 * mu_face * grad_v_face[2] -
        (2 / 3) * mu_face * div_v

    tau_zz =
        2 * mu_face * grad_w_face[3] -
        (2 / 3) * mu_face * div_v

    tau_xy =
        mu_face * (
            grad_u_face[2] +
            grad_v_face[1]
        )

    tau_xz =
        mu_face * (
            grad_u_face[3] +
            grad_w_face[1]
        )

    tau_yz =
        mu_face * (
            grad_v_face[3] +
            grad_w_face[2]
        )

    # ------------------------------------------------------------
    # 5. Viscous traction τ⋅n
    # ------------------------------------------------------------

    n_x = cell_face_normal[1]
    n_y = cell_face_normal[2]
    n_z = cell_face_normal[3]

    traction_x =
        tau_xx * n_x +
        tau_xy * n_y +
        tau_xz * n_z

    traction_y =
        tau_xy * n_x +
        tau_yy * n_y +
        tau_yz * n_z

    traction_z =
        tau_xz * n_x +
        tau_yz * n_y +
        tau_zz * n_z

    # ------------------------------------------------------------
    # 6. Fourier heat conduction
    #
    # q = -k∇T
    #
    # In the Navier-Stokes viscous flux,
    # the energy term contains:
    #
    #   v⋅(τ⋅n) + k∇T⋅n
    #
    # because the minus sign is already contained in q = -k∇T.
    # ------------------------------------------------------------

    conductive_energy_flux =
        k_face * (
            grad_T_face[1] * n_x +
            grad_T_face[2] * n_y +
            grad_T_face[3] * n_z
        )

    viscous_work_flux =
        vel_u_face * traction_x +
        vel_v_face * traction_y +
        vel_w_face * traction_z

    viscous_energy_flux =
        viscous_work_flux +
        conductive_energy_flux

    # ------------------------------------------------------------
    # 7. Add integrated viscous flux to cell A
    #
    # Governing equation:
    #
    #   ∂U/∂t = -∇⋅F_conv + ∇⋅F_visc
    #
    # so viscous flux enters with a PLUS sign here.
    # ------------------------------------------------------------

    du.momentum_density_u_flow[idx_a] += traction_x * area
    du.momentum_density_v_flow[idx_a] += traction_y * area
    du.momentum_density_w_flow[idx_a] += traction_z * area

    du.volumetric_energy_flow[idx_a] +=
        viscous_energy_flux * area

    return nothing
end