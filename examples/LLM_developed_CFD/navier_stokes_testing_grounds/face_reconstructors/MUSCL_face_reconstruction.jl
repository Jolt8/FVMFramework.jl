"""
    update_MUSCL_gradients!(u, stencil)

Populate the density gradient required by the primitive-variable MUSCL
reconstruction. Velocity and temperature gradients are populated by
`update_weighted_least_squares_gradients!`.
"""
function update_MUSCL_gradients!(u, stencil)
    populate_weighted_least_squares_gradient!(u.grad_density, u.density, stencil)
    return nothing
end


"""
    MUSCL_face_reconstruction!(du, u, p, t, system, geo,
                               idx_a, face_a, idx_b, face_b)

Reconstruct the state on the `idx_a` side of an internal face using a
cell-centred, piecewise-linear MUSCL reconstruction. Density, velocity, and
temperature are reconstructed with weighted-least-squares gradients and a
Venkatakrishnan limiter, then converted back to the conservative variables
expected by HLLC.

HLLC calls this function once from each side of the face, so the returned
state is always oriented from `idx_a` toward the shared face.
"""
function MUSCL_face_reconstruction!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b,
)
    (
        dist,
        face_area_a, face_normal_a, face_distance_a, vol_a,
        face_area_b, face_normal_b, face_distance_b, vol_b
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    density_face = _MUSCL_limited_face_value(
        u.density, u.grad_density, idx_a,
        face_normal_a, face_distance_a, geo,
    )
    vel_u_face = _MUSCL_limited_face_value(
        u.vel_u, u.grad_vel_u, idx_a,
        face_normal_a, face_distance_a, geo,
    )
    vel_v_face = _MUSCL_limited_face_value(
        u.vel_v, u.grad_vel_v, idx_a,
        face_normal_a, face_distance_a, geo,
    )
    vel_w_face = _MUSCL_limited_face_value(
        u.vel_w, u.grad_vel_w, idx_a,
        face_normal_a, face_distance_a, geo,
    )
    temperature_face = _MUSCL_limited_face_value(
        u.temperature, u.grad_temperature, idx_a,
        face_normal_a, face_distance_a, geo,
    )

    momentum_density_u_face = density_face * vel_u_face
    momentum_density_v_face = density_face * vel_v_face
    momentum_density_w_face = density_face * vel_w_face

    specific_internal_energy_face = u.cv[idx_a] * temperature_face
    specific_kinetic_energy_face = 0.5 * (
        vel_u_face^2 + vel_v_face^2 + vel_w_face^2
    )
    volumetric_energy_face = density_face * (
        specific_internal_energy_face + specific_kinetic_energy_face
    )

    return (
        density_face,
        momentum_density_u_face,
        momentum_density_v_face,
        momentum_density_w_face,
        volumetric_energy_face,
    )
end


function _MUSCL_limited_face_value(
    cell_values,
    cell_gradients,
    cell_id,
    face_normal,
    face_distance,
    geo,
)
    cell_value = cell_values[cell_id]
    unlimited_increment = face_distance * (
        cell_gradients[cell_id, 1] * face_normal[1] +
        cell_gradients[cell_id, 2] * face_normal[2] +
        cell_gradients[cell_id, 3] * face_normal[3]
    )

    local_minimum = cell_value
    local_maximum = cell_value
    for neighbor_data in geo.cell_neighbors[cell_id][2]
        neighbor_id = neighbor_data[1]
        neighbor_id > 0 || continue
        neighbor_value = cell_values[neighbor_id]
        local_minimum = min(local_minimum, neighbor_value)
        local_maximum = max(local_maximum, neighbor_value)
    end

    limiter = _MUSCL_venkatakrishnan_limiter(
        cell_value,
        unlimited_increment,
        local_minimum,
        local_maximum,
        geo.cell_volumes[cell_id],
    )

    return cell_value + limiter * unlimited_increment
end

"""
    _MUSCL_venkatakrishnan_limiter(
        cell_value, unlimited_increment, local_minimum, local_maximum,
        cell_volume; K=3.0,
    )

Return the smooth Venkatakrishnan limiter for one face extrapolation. The
regularization is `epsilon^2 = (K * h)^3`, where `h` is the cube root of the
cell volume. The state variables used by this example are unit-stripped before
the spatial operator is evaluated.

The result is capped at one so the limiter never steepens the
weighted-least-squares gradient.
"""
function _MUSCL_venkatakrishnan_limiter(
    cell_value,
    unlimited_increment,
    local_minimum,
    local_maximum,
    cell_volume;
    K=3.0,
)
    if unlimited_increment > zero(unlimited_increment)
        bound_increment = local_maximum - cell_value
    elseif unlimited_increment < zero(unlimited_increment)
        bound_increment = local_minimum - cell_value
    else
        return 1.0
    end

    K >= zero(K) || throw(ArgumentError("K must be non-negative"))
    cell_volume > zero(cell_volume) || throw(ArgumentError(
        "cell_volume must be positive",
    ))

    characteristic_length = cbrt(cell_volume)
    epsilon_squared = (K * characteristic_length)^3
    bound_increment_squared = bound_increment^2
    unlimited_increment_squared = unlimited_increment^2
    cross_increment = bound_increment * unlimited_increment

    limiter = (
        bound_increment_squared + 2 * cross_increment + epsilon_squared
    ) / (
        bound_increment_squared + cross_increment +
        2 * unlimited_increment_squared + epsilon_squared
    )

    return max(zero(limiter), min(one(limiter), limiter))
end

function _MUSCL_bound_limiter(
    cell_value,
    unlimited_increment,
    local_minimum,
    local_maximum,
)
    if unlimited_increment > zero(unlimited_increment)
        ratio = (local_maximum - cell_value) / unlimited_increment
    elseif unlimited_increment < zero(unlimited_increment)
        ratio = (local_minimum - cell_value) / unlimited_increment
    else
        return 1.0
    end

    return max(zero(ratio), min(one(ratio), ratio))
end

