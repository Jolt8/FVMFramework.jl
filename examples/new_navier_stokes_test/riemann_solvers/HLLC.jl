function physical_flux(
    density,
    momentum_density_u,
    momentum_density_v,
    momentum_density_w,
    volumetric_energy,
    pressure,
    normal
)
    vel_u = momentum_density_u / density
    vel_v = momentum_density_v / density
    vel_w = momentum_density_w / density

    normal_velocity =
        vel_u * normal[1] +
        vel_v * normal[2] +
        vel_w * normal[3]

    return (
        density * normal_velocity,

        momentum_density_u * normal_velocity +
            pressure * normal[1],

        momentum_density_v * normal_velocity +
            pressure * normal[2],

        momentum_density_w * normal_velocity +
            pressure * normal[3],

        (volumetric_energy + pressure) * normal_velocity
    )
end

function primitive_from_conservative(
    density,
    momentum_density_u,
    momentum_density_v,
    momentum_density_w,
    volumetric_energy,
    specific_heat_ratio
)
    vel_u = momentum_density_u / density
    vel_v = momentum_density_v / density
    vel_w = momentum_density_w / density

    kinetic_energy_density = 0.5 * (
        momentum_density_u^2 +
        momentum_density_v^2 +
        momentum_density_w^2
    ) / density

    pressure = (specific_heat_ratio - 1) * (
        volumetric_energy - kinetic_energy_density
    )

    speed_of_sound = sqrt(specific_heat_ratio * pressure / density)

    return vel_u, vel_v, vel_w, pressure, speed_of_sound
end

function hllc_flux(
    density_a,
    momentum_density_u_a,
    momentum_density_v_a,
    momentum_density_w_a,
    volumetric_energy_a,
    gamma_a,

    density_b,
    momentum_density_u_b,
    momentum_density_v_b,
    momentum_density_w_b,
    volumetric_energy_b,
    gamma_b,

    cell_face_normal
)
    vel_u_a, vel_v_a, vel_w_a, pressure_a, speed_of_sound_a = primitive_from_conservative(
        density_a,
        momentum_density_u_a,
        momentum_density_v_a,
        momentum_density_w_a,
        volumetric_energy_a,
        gamma_a
    )

    vel_u_b, vel_v_b, vel_w_b, pressure_b, speed_of_sound_b = primitive_from_conservative(
        density_b,
        momentum_density_u_b,
        momentum_density_v_b,
        momentum_density_w_b,
        volumetric_energy_b,
        gamma_b
    )

    normal_velocity_a = 
        vel_u_a * cell_face_normal[1] + 
        vel_v_a * cell_face_normal[2] + 
        vel_w_a * cell_face_normal[3]

    normal_velocity_b = 
        vel_u_b * cell_face_normal[1] + 
        vel_v_b * cell_face_normal[2] + 
        vel_w_b * cell_face_normal[3]

    S_a = min(normal_velocity_a - speed_of_sound_a, normal_velocity_b - speed_of_sound_b)
    S_b = max(normal_velocity_a + speed_of_sound_a, normal_velocity_b + speed_of_sound_b)

    S_M = (
        (pressure_b - pressure_a) + 
        (density_a * normal_velocity_a * (S_a - normal_velocity_a)) - 
        (density_b * normal_velocity_b * (S_b - normal_velocity_b))
    ) / (density_a * (S_a - normal_velocity_a) - density_b * (S_b - normal_velocity_b))

    pressure_star = pressure_a + density_a * (S_a - normal_velocity_a) * (S_M - normal_velocity_a)

    density_star_a = density_a * (S_a - normal_velocity_a) / (S_a - S_M)
    density_star_b = density_b * (S_b - normal_velocity_b) / (S_b - S_M)
    
    vel_u_star_a = vel_u_a + (S_M - normal_velocity_a) * cell_face_normal[1]
    vel_v_star_a = vel_v_a + (S_M - normal_velocity_a) * cell_face_normal[2]
    vel_w_star_a = vel_w_a + (S_M - normal_velocity_a) * cell_face_normal[3]

    vel_u_star_b = vel_u_b + (S_M - normal_velocity_b) * cell_face_normal[1]
    vel_v_star_b = vel_v_b + (S_M - normal_velocity_b) * cell_face_normal[2]
    vel_w_star_b = vel_w_b + (S_M - normal_velocity_b) * cell_face_normal[3]

    momentum_density_u_star_a = density_star_a * vel_u_star_a
    momentum_density_v_star_a = density_star_a * vel_v_star_a
    momentum_density_w_star_a = density_star_a * vel_w_star_a

    momentum_density_u_star_b = density_star_b * vel_u_star_b
    momentum_density_v_star_b = density_star_b * vel_v_star_b
    momentum_density_w_star_b = density_star_b * vel_w_star_b

    volumetric_energy_star_a = (
        volumetric_energy_a * (S_a - normal_velocity_a) - 
        pressure_a * normal_velocity_a +
        pressure_star * S_M
    ) / (S_a - S_M)

    volumetric_energy_star_b = (
        volumetric_energy_b * (S_b - normal_velocity_b) - 
        pressure_b * normal_velocity_b +
        pressure_star * S_M
    ) / (S_b - S_M)


    F_density_a, 
    F_momentum_density_u_a,
    F_momentum_density_v_a,
    F_momentum_density_w_a,
    F_volumetric_energy_a = 
    physical_flux(
        density_a,
        momentum_density_u_a,
        momentum_density_v_a,
        momentum_density_w_a,
        volumetric_energy_a,
        pressure_a,
        cell_face_normal
    )

    F_density_b, 
    F_momentum_density_u_b,
    F_momentum_density_v_b,
    F_momentum_density_w_b,
    F_volumetric_energy_b = 
    physical_flux(
        density_b,
        momentum_density_u_b,
        momentum_density_v_b,
        momentum_density_w_b,
        volumetric_energy_b,
        pressure_b,
        cell_face_normal
    )

    F_density_hllc = 0.0
    F_momentum_density_u_hllc = 0.0
    F_momentum_density_v_hllc = 0.0
    F_momentum_density_w_hllc = 0.0
    F_volumetric_energy_hllc = 0.0

    if 0.0 <= S_a
        F_density_hllc = F_density_a
        F_momentum_density_u_hllc = F_momentum_density_u_a
        F_momentum_density_v_hllc = F_momentum_density_v_a
        F_momentum_density_w_hllc = F_momentum_density_w_a
        F_volumetric_energy_hllc = F_volumetric_energy_a
    elseif S_a < 0.0 && 0.0 <= S_M
        F_density_hllc = F_density_a + S_a * (density_star_a - density_a)
        F_momentum_density_u_hllc = F_momentum_density_u_a + S_a * (momentum_density_u_star_a - momentum_density_u_a)
        F_momentum_density_v_hllc = F_momentum_density_v_a + S_a * (momentum_density_v_star_a - momentum_density_v_a)
        F_momentum_density_w_hllc = F_momentum_density_w_a + S_a * (momentum_density_w_star_a - momentum_density_w_a)
        F_volumetric_energy_hllc = F_volumetric_energy_a + S_a * (volumetric_energy_star_a - volumetric_energy_a)
    elseif S_M <= 0.0 && 0.0 <= S_b
        F_density_hllc = F_density_b + S_b * (density_star_b - density_b)
        F_momentum_density_u_hllc = F_momentum_density_u_b + S_b * (momentum_density_u_star_b - momentum_density_u_b)
        F_momentum_density_v_hllc = F_momentum_density_v_b + S_b * (momentum_density_v_star_b - momentum_density_v_b)
        F_momentum_density_w_hllc = F_momentum_density_w_b + S_b * (momentum_density_w_star_b - momentum_density_w_b)
        F_volumetric_energy_hllc = F_volumetric_energy_b + S_b * (volumetric_energy_star_b - volumetric_energy_b)
    else
        F_density_hllc = F_density_b
        F_momentum_density_u_hllc = F_momentum_density_u_b
        F_momentum_density_v_hllc = F_momentum_density_v_b
        F_momentum_density_w_hllc = F_momentum_density_w_b
        F_volumetric_energy_hllc = F_volumetric_energy_b
    end

    return (
        F_density_hllc,
        F_momentum_density_u_hllc,
        F_momentum_density_v_hllc,
        F_momentum_density_w_hllc,
        F_volumetric_energy_hllc
    )
end

function HLLC!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    area, cell_face_normal, dist_to_face,
    cell_neighbor_normal, cell_neighbor_dist,
    vol,
    face_reconstructor!
)
    density_a, momentum_density_u_a, momentum_density_v_a, momentum_density_w_a, volumetric_energy_a = 
    face_reconstructor!(
        du, u, p, t, 
        idx_a, idx_b, face_idx, 
        area, cell_face_normal, dist_to_face, 
        cell_neighbor_normal, cell_neighbor_dist,
        vol
    )

    density_b, momentum_density_u_b, momentum_density_v_b, momentum_density_w_b, volumetric_energy_b = 
    face_reconstructor!(
        du, u, p, t, 
        idx_b, idx_a, face_idx, 
        area, cell_face_normal, dist_to_face, 
        cell_neighbor_normal, cell_neighbor_dist,
        vol
    )

    F_density, 
    F_momentum_density_u,
    F_momentum_density_v,
    F_momentum_density_w,
    F_volumetric_energy = 
    hllc_flux(
        density_a,
        momentum_density_u_a,
        momentum_density_v_a,
        momentum_density_w_a,
        volumetric_energy_a,
        (u.cp[idx_a] / u.cv[idx_b]),

        density_b,
        momentum_density_u_b,
        momentum_density_v_b,
        momentum_density_w_b,
        volumetric_energy_b,
        (u.cp[idx_b] / u.cv[idx_b]),

        cell_face_normal
    )

    du.density_flow[idx_a] -= area * F_density
    du.momentum_density_u_flow[idx_a] -= area * F_momentum_density_u
    du.momentum_density_v_flow[idx_a] -= area * F_momentum_density_v
    du.momentum_density_w_flow[idx_a] -= area * F_momentum_density_w
    du.volumetric_energy_flow[idx_a] -= area * F_volumetric_energy
end