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

Consequently, the correction approaches a common normal velocity as Mach
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
    velocity_magnitude_a = sqrt(vel_u_a^2 + vel_v_a^2 + vel_w_a^2)
    velocity_magnitude_b = sqrt(vel_u_b^2 + vel_v_b^2 + vel_w_b^2)

    local_mach_number = max(
        velocity_magnitude_a / speed_of_sound_a,
        velocity_magnitude_b / speed_of_sound_b,
    )
    velocity_jump_scaling = min(one(local_mach_number), local_mach_number)

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
