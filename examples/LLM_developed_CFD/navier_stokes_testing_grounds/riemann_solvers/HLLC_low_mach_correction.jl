"""
    no_low_mach_correction(
        vel_u_a, vel_v_a, vel_w_a, speed_of_sound_a,
        vel_u_b, vel_v_b, vel_w_b, speed_of_sound_b,
        cell_face_normal,
    )

Return the left and right velocities unchanged. Pass this function to `HLLC!`
when the standard, uncorrected HLLC flux is desired.
"""
function no_low_mach_correction(
    vel_u_a,
    vel_v_a,
    vel_w_a,
    speed_of_sound_a,
    vel_u_b,
    vel_v_b,
    vel_w_b,
    speed_of_sound_b,
    cell_face_normal,
)
    return (
        vel_u_a,
        vel_v_a,
        vel_w_a,
        vel_u_b,
        vel_v_b,
        vel_w_b,
    )
end


const THORNBER_MACH_REGULARIZATION = 1.0e-3


"""
    regularized_velocity_magnitude(vel_u, vel_v, vel_w, speed_of_sound)

Return a smoothly regularized velocity magnitude. The small acoustic-speed
term removes the undefined derivative of the Euclidean norm at zero velocity
and provides a small amount of velocity-jump damping in stagnant cells.
"""
function regularized_velocity_magnitude(
    vel_u,
    vel_v,
    vel_w,
    speed_of_sound,
)
    regularization_speed = THORNBER_MACH_REGULARIZATION * speed_of_sound
    return sqrt(
        vel_u^2 + vel_v^2 + vel_w^2 + regularization_speed^2
    )
end


function c1_absolute_value(value, transition_width)
    absolute_value = abs(value)
    if absolute_value >= transition_width
        return absolute_value
    end
    return 0.5 * (value^2 / transition_width + transition_width)
end


function c1_maximum(value_a, value_b; transition_width = 1.0e-8)
    width = one(value_a + value_b) * transition_width
    return 0.5 * (
        value_a + value_b + c1_absolute_value(value_a - value_b, width)
    )
end


function c1_unit_cap(value; transition_width = 1.0e-6)
    one_value = one(value)
    width = one_value * transition_width
    return 0.5 * (
        one_value + value - c1_absolute_value(one_value - value, width)
    )
end


"""
    thornber_low_mach_correction(
        vel_u_a, vel_v_a, vel_w_a, speed_of_sound_a,
        vel_u_b, vel_v_b, vel_w_b, speed_of_sound_b,
        cell_face_normal,
    )

Apply the velocity-reconstruction correction of Thornber et al. (2008) to
the two states supplied to HLLC.

Only the face-normal velocity jump is reduced. The average normal velocity
and both tangential velocity components are preserved. The jump is multiplied
by the local Mach-number coefficient

    z = min(1, max(norm(velocity_a) / speed_of_sound_a,
                   norm(velocity_b) / speed_of_sound_b)).

The velocity magnitudes have a `1e-3` Mach regularization, and the maximum and
unit cap use compact C1 transition bands around their kinks. Away from the
regularization and transition bands the expression recovers the Thornber
coefficient. This gives automatic differentiation a continuous derivative at
zero velocity, when the two Mach numbers are equal, and when the local Mach
number crosses one.

Consequently, the correction strongly damps the normal velocity jump as Mach
number tends to zero and exactly recovers standard HLLC when `z == 1`.
"""
function thornber_low_mach_correction(
    vel_u_a,
    vel_v_a,
    vel_w_a,
    speed_of_sound_a,
    vel_u_b,
    vel_v_b,
    vel_w_b,
    speed_of_sound_b,
    cell_face_normal,
)
    velocity_magnitude_a = regularized_velocity_magnitude(
        vel_u_a,
        vel_v_a,
        vel_w_a,
        speed_of_sound_a,
    )
    velocity_magnitude_b = regularized_velocity_magnitude(
        vel_u_b,
        vel_v_b,
        vel_w_b,
        speed_of_sound_b,
    )

    local_mach_number = c1_maximum(
        velocity_magnitude_a / speed_of_sound_a,
        velocity_magnitude_b / speed_of_sound_b,
    )
    velocity_jump_scaling = c1_unit_cap(local_mach_number)

    normal_velocity_a =
        vel_u_a * cell_face_normal[1] +
        vel_v_a * cell_face_normal[2] +
        vel_w_a * cell_face_normal[3]

    normal_velocity_b =
        vel_u_b * cell_face_normal[1] +
        vel_v_b * cell_face_normal[2] +
        vel_w_b * cell_face_normal[3]

    average_normal_velocity = 0.5 * (normal_velocity_a + normal_velocity_b)
    half_normal_velocity_jump = 0.5 * (normal_velocity_a - normal_velocity_b)

    corrected_normal_velocity_a =
        average_normal_velocity + velocity_jump_scaling * half_normal_velocity_jump
    corrected_normal_velocity_b =
        average_normal_velocity - velocity_jump_scaling * half_normal_velocity_jump

    normal_velocity_change_a = corrected_normal_velocity_a - normal_velocity_a
    normal_velocity_change_b = corrected_normal_velocity_b - normal_velocity_b

    return (
        vel_u_a + normal_velocity_change_a * cell_face_normal[1],
        vel_v_a + normal_velocity_change_a * cell_face_normal[2],
        vel_w_a + normal_velocity_change_a * cell_face_normal[3],
        vel_u_b + normal_velocity_change_b * cell_face_normal[1],
        vel_v_b + normal_velocity_change_b * cell_face_normal[2],
        vel_w_b + normal_velocity_change_b * cell_face_normal[3],
    )
end
