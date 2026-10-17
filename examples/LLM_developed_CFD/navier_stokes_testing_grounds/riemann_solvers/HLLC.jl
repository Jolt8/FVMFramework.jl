include(joinpath(
    @__DIR__, "..", "species_transport", "conservative_species_transport.jl",
))

function energy_with_corrected_velocity(
    volumetric_energy,
    density,
    original_vel_u,
    original_vel_v,
    original_vel_w,
    corrected_vel_u,
    corrected_vel_v,
    corrected_vel_w,
)
    original_squared_speed =
        original_vel_u^2 + original_vel_v^2 + original_vel_w^2
    corrected_squared_speed =
        corrected_vel_u^2 + corrected_vel_v^2 + corrected_vel_w^2
    return volumetric_energy +
           0.5 * density * (corrected_squared_speed - original_squared_speed)
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

    cell_face_normal,
    low_mach_correction = no_low_mach_correction,
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

    original_vel_u_a = vel_u_a
    original_vel_v_a = vel_v_a
    original_vel_w_a = vel_w_a
    original_vel_u_b = vel_u_b
    original_vel_v_b = vel_v_b
    original_vel_w_b = vel_w_b

    (
        vel_u_a,
        vel_v_a,
        vel_w_a,
        vel_u_b,
        vel_v_b,
        vel_w_b,
    ) = low_mach_correction(
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

    volumetric_energy_a = energy_with_corrected_velocity(
        volumetric_energy_a,
        density_a,
        original_vel_u_a,
        original_vel_v_a,
        original_vel_w_a,
        vel_u_a,
        vel_v_a,
        vel_w_a,
    )
    volumetric_energy_b = energy_with_corrected_velocity(
        volumetric_energy_b,
        density_b,
        original_vel_u_b,
        original_vel_v_b,
        original_vel_w_b,
        vel_u_b,
        vel_v_b,
        vel_w_b,
    )

    momentum_density_u_a = density_a * vel_u_a
    momentum_density_v_a = density_a * vel_v_a
    momentum_density_w_a = density_a * vel_w_a

    momentum_density_u_b = density_b * vel_u_b
    momentum_density_v_b = density_b * vel_v_b
    momentum_density_w_b = density_b * vel_w_b

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

    if zero(S_a) <= S_a
        F_density_hllc = F_density_a
        F_momentum_density_u_hllc = F_momentum_density_u_a
        F_momentum_density_v_hllc = F_momentum_density_v_a
        F_momentum_density_w_hllc = F_momentum_density_w_a
        F_volumetric_energy_hllc = F_volumetric_energy_a
    elseif S_a < zero(S_a) && zero(S_M) <= S_M
        F_density_hllc = F_density_a + S_a * (density_star_a - density_a)
        F_momentum_density_u_hllc = F_momentum_density_u_a + S_a * (momentum_density_u_star_a - momentum_density_u_a)
        F_momentum_density_v_hllc = F_momentum_density_v_a + S_a * (momentum_density_v_star_a - momentum_density_v_a)
        F_momentum_density_w_hllc = F_momentum_density_w_a + S_a * (momentum_density_w_star_a - momentum_density_w_a)
        F_volumetric_energy_hllc = F_volumetric_energy_a + S_a * (volumetric_energy_star_a - volumetric_energy_a)
    elseif S_M <= zero(S_M) && zero(S_b) <= S_b
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
#=
hllc_flux(
    1.18u"kg/m^3",
    600.0u"m/s" * 1.18u"kg/m^3",
    0.0u"kg/(m^2*s)",
    0.0u"kg/(m^2*s)",
    466572.0u"J/m^3",
    1.4,

    1.18u"kg/m^3",
    600.0u"m/s" * 1.18u"kg/m^3",
    0.0u"kg/(m^2*s)",
    0.0u"kg/(m^2*s)",
    466572.0u"J/m^3",
    1.4,

    Ferrite.Vec{3}((1.0, 0.0, 0.0))
) 

hllc_flux(
    1.18u"kg/m^3",
    600.0u"m/s" * 1.18u"kg/m^3",
    0.0u"kg/(m^2*s)",
    0.0u"kg/(m^2*s)",
    466572.0u"J/m^3",
    1.4,

    1.18u"kg/m^3",
    600.0u"m/s" * 1.18u"kg/m^3",
    0.0u"kg/(m^2*s)",
    0.0u"kg/(m^2*s)",
    466572.0u"J/m^3",
    1.4,

    Ferrite.Vec{3}((-1.0, 0.0, 0.0))
)=#

function HLLC!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
    face_reconstructor!,
    low_mach_correction = no_low_mach_correction,
)
    (
        dist,
        face_area_a, face_normal_a, face_distance_a, vol_a,
        face_area_b, face_normal_b, face_distance_b, vol_b
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    density_a, momentum_density_u_a, momentum_density_v_a, momentum_density_w_a, volumetric_energy_a = 
    face_reconstructor!(
        du, u, p, t, system, geo,
        idx_a, face_a, 
        idx_b, face_b,
    )

    density_b, momentum_density_u_b, momentum_density_v_b, momentum_density_w_b, volumetric_energy_b = 
    face_reconstructor!(
        du, u, p, t, system, geo,
        idx_b, face_b,
        idx_a, face_a,
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
        (u.cp[idx_a] / u.cv[idx_a]),

        density_b,
        momentum_density_u_b,
        momentum_density_v_b,
        momentum_density_w_b,
        volumetric_energy_b,
        (u.cp[idx_b] / u.cv[idx_b]),

        face_normal_a,
        low_mach_correction,
    )

    du.density_flow[idx_a] -= face_area_a * F_density
    du.momentum_density_u_flow[idx_a] -= face_area_a * F_momentum_density_u
    du.momentum_density_v_flow[idx_a] -= face_area_a * F_momentum_density_v
    du.momentum_density_w_flow[idx_a] -= face_area_a * F_momentum_density_w
    du.volumetric_energy_flow[idx_a] -= face_area_a * F_volumetric_energy

    du.density_flow[idx_b] += face_area_a * F_density
    du.momentum_density_u_flow[idx_b] += face_area_a * F_momentum_density_u
    du.momentum_density_v_flow[idx_b] += face_area_a * F_momentum_density_v
    du.momentum_density_w_flow[idx_b] += face_area_a * F_momentum_density_w
    du.volumetric_energy_flow[idx_b] += face_area_a * F_volumetric_energy

    if hasproperty(u, :species_densities)
        add_hllc_species_advection_flux!(
            du,
            u,
            idx_a,
            idx_b,
            face_area_a,
            F_density,
        )
    end
    if hasproperty(u, :turbulent_kinetic_energy_density)
        add_hllc_sst_advection_flux!(
            du,
            u,
            idx_a,
            idx_b,
            face_area_a,
            F_density,
        )
    end
end
