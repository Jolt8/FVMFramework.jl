
function _correct_wall_gradient(
    cell_gradients,
    cell_id,
    cell_value,
    wall_value,
    normal_x,
    normal_y,
    normal_z,
    dist_to_face,
)
    gradient_x = cell_gradients[cell_id, 1]
    gradient_y = cell_gradients[cell_id, 2]
    gradient_z = cell_gradients[cell_id, 3]

    projected_gradient =
        gradient_x * normal_x +
        gradient_y * normal_y +
        gradient_z * normal_z
    required_normal_gradient = (wall_value - cell_value) / dist_to_face
    correction = required_normal_gradient - projected_gradient

    return (
        gradient_x + correction * normal_x,
        gradient_y + correction * normal_y,
        gradient_z + correction * normal_z,
    )
end

function non_moving_wall_viscous_and_diffusive_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    return _wall_viscous_and_diffusive_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, 
        idx_b, face_b,
        0.0, 0.0, 0.0
    )
end


function moving_wall_viscous_and_diffusive_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    return _wall_viscous_and_diffusive_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, 
        idx_b, face_b,
        u.wall_vel_x[idx_b], u.wall_vel_y[idx_b], u.wall_vel_z[idx_b]
    )
end


function _wall_viscous_and_diffusive_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
    wall_vel_x, wall_vel_y, wall_vel_z,
)
    face_area, face_normal, face_distance, vol = boundary_geometry(geo, idx_a, face_a)

    normal_magnitude = sqrt(
        face_normal[1]^2 +
        face_normal[2]^2 +
        face_normal[3]^2
    )
    iszero(normal_magnitude) && throw(ArgumentError("the wall normal must be nonzero"))
    face_distance > zero(face_distance) || throw(ArgumentError(
        "the cell-centroid-to-wall distance must be positive",
    ))

    normal_x = face_normal[1] / normal_magnitude
    normal_y = face_normal[2] / normal_magnitude
    normal_z = face_normal[3] / normal_magnitude

    # Retain the tangential part of each cell-centred gradient and replace its
    # wall-normal component with the one-sided Dirichlet gradient required by
    # the no-slip condition.
    grad_vel_u_x, grad_vel_u_y, grad_vel_u_z = _correct_wall_gradient(
        u.grad_vel_u, idx_a, u.vel_u[idx_a], wall_vel_x,
        normal_x, normal_y, normal_z, face_distance,
    )
    grad_vel_v_x, grad_vel_v_y, grad_vel_v_z = _correct_wall_gradient(
        u.grad_vel_v, idx_a, u.vel_v[idx_a], wall_vel_y,
        normal_x, normal_y, normal_z, face_distance,
    )
    grad_vel_w_x, grad_vel_w_y, grad_vel_w_z = _correct_wall_gradient(
        u.grad_vel_w, idx_a, u.vel_w[idx_a], wall_vel_z,
        normal_x, normal_y, normal_z, face_distance,
    )

    dynamic_viscosity = u.mu[idx_a]
    velocity_divergence = grad_vel_u_x + grad_vel_v_y + grad_vel_w_z
    isotropic_stress = (2 / 3) * dynamic_viscosity * velocity_divergence

    tau_xx = 2 * dynamic_viscosity * grad_vel_u_x - isotropic_stress
    tau_yy = 2 * dynamic_viscosity * grad_vel_v_y - isotropic_stress
    tau_zz = 2 * dynamic_viscosity * grad_vel_w_z - isotropic_stress
    tau_xy = dynamic_viscosity * (grad_vel_u_y + grad_vel_v_x)
    tau_xz = dynamic_viscosity * (grad_vel_u_z + grad_vel_w_x)
    tau_yz = dynamic_viscosity * (grad_vel_v_z + grad_vel_w_y)

    traction_x = tau_xx * normal_x + tau_xy * normal_y + tau_xz * normal_z
    traction_y = tau_xy * normal_x + tau_yy * normal_y + tau_yz * normal_z
    traction_z = tau_xz * normal_x + tau_yz * normal_y + tau_zz * normal_z

    # No temperature is prescribed, so the wall is adiabatic and its
    # conductive heat flux is zero.
    viscous_energy_flux =
        wall_vel_x * traction_x +
        wall_vel_y * traction_y +
        wall_vel_z * traction_z

    du.momentum_density_u_flow[idx_a] += face_area * traction_x
    du.momentum_density_v_flow[idx_a] += face_area * traction_y
    du.momentum_density_w_flow[idx_a] += face_area * traction_z
    du.volumetric_energy_flow[idx_a] += face_area * viscous_energy_flux

    return nothing
end
