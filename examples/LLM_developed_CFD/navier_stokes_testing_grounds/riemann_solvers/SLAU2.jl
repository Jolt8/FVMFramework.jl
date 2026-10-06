"""
    slau2_flux(...)

Evaluate the five-equation Euler flux with the SLAU2 all-speed AUSM-family
scheme of Kitamura and Shima (JCP 245, 2013, pp. 62--83,
https://doi.org/10.1016/j.jcp.2013.02.046).

This is SLAU2, not the original SLAU scheme. In particular, the last pressure
dissipation term is scaled by the RMS velocity, interface sound speed, and
mean density. Conservative species and SST variables use the signed SLAU2 mass
flux for upwind advection.
"""
function slau2_flux(
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
)
    (
        vel_u_a,
        vel_v_a,
        vel_w_a,
        pressure_a,
        speed_of_sound_a,
    ) = _slau2_primitive_from_conservative(
        density_a,
        momentum_density_u_a,
        momentum_density_v_a,
        momentum_density_w_a,
        volumetric_energy_a,
        gamma_a,
    )
    (
        vel_u_b,
        vel_v_b,
        vel_w_b,
        pressure_b,
        speed_of_sound_b,
    ) = _slau2_primitive_from_conservative(
        density_b,
        momentum_density_u_b,
        momentum_density_v_b,
        momentum_density_w_b,
        volumetric_energy_b,
        gamma_b,
    )

    mass_flux, pressure_flux = _slau2_mass_and_pressure_flux(
        density_a,
        vel_u_a,
        vel_v_a,
        vel_w_a,
        pressure_a,
        speed_of_sound_a,
        density_b,
        vel_u_b,
        vel_v_b,
        vel_w_b,
        pressure_b,
        speed_of_sound_b,
        cell_face_normal,
    )

    if mass_flux >= 0.0
        face_vel_u = vel_u_a
        face_vel_v = vel_v_a
        face_vel_w = vel_w_a
        face_total_enthalpy = (volumetric_energy_a + pressure_a) / density_a
    else
        face_vel_u = vel_u_b
        face_vel_v = vel_v_b
        face_vel_w = vel_w_b
        face_total_enthalpy = (volumetric_energy_b + pressure_b) / density_b
    end

    return (
        mass_flux,
        mass_flux * face_vel_u + pressure_flux * cell_face_normal[1],
        mass_flux * face_vel_v + pressure_flux * cell_face_normal[2],
        mass_flux * face_vel_w + pressure_flux * cell_face_normal[3],
        mass_flux * face_total_enthalpy,
    )
end

function _slau2_primitive_from_conservative(
    density,
    momentum_density_u,
    momentum_density_v,
    momentum_density_w,
    volumetric_energy,
    gamma,
)
    density = max(density, 1e-10)
    vel_u = momentum_density_u / density
    vel_v = momentum_density_v / density
    vel_w = momentum_density_w / density
    kinetic_energy_density = 0.5 * (
        momentum_density_u^2 +
        momentum_density_v^2 +
        momentum_density_w^2
    ) / density
    internal_energy_density = max(
        volumetric_energy - kinetic_energy_density,
        1e-10,
    )
    pressure = (gamma - 1.0) * internal_energy_density
    speed_of_sound = sqrt(gamma * pressure / density)
    return vel_u, vel_v, vel_w, pressure, speed_of_sound
end

function _slau2_pressure_weight_left(mach_number)
    if abs(mach_number) < 1.0
        return 0.25 * (2.0 - mach_number) * (mach_number + 1.0)^2
    elseif mach_number >= 0.0
        return 1.0
    else
        return 0.0
    end
end

function _slau2_pressure_weight_right(mach_number)
    if abs(mach_number) < 1.0
        return 0.25 * (2.0 + mach_number) * (mach_number - 1.0)^2
    elseif mach_number >= 0.0
        return 0.0
    else
        return 1.0
    end
end

function _slau2_mass_and_pressure_flux(
    density_a,
    vel_u_a,
    vel_v_a,
    vel_w_a,
    pressure_a,
    speed_of_sound_a,
    density_b,
    vel_u_b,
    vel_v_b,
    vel_w_b,
    pressure_b,
    speed_of_sound_b,
    cell_face_normal,
)
    normal_velocity_a =
        vel_u_a * cell_face_normal[1] +
        vel_v_a * cell_face_normal[2] +
        vel_w_a * cell_face_normal[3]
    normal_velocity_b =
        vel_u_b * cell_face_normal[1] +
        vel_v_b * cell_face_normal[2] +
        vel_w_b * cell_face_normal[3]

    interface_speed_of_sound = 0.5 * (speed_of_sound_a + speed_of_sound_b)
    mach_number_a = normal_velocity_a / interface_speed_of_sound
    mach_number_b = normal_velocity_b / interface_speed_of_sound

    squared_speed_a = vel_u_a^2 + vel_v_a^2 + vel_w_a^2
    squared_speed_b = vel_u_b^2 + vel_v_b^2 + vel_w_b^2
    rms_speed = sqrt(0.5 * (squared_speed_a + squared_speed_b))
    limited_mach_number = min(1.0, rms_speed / interface_speed_of_sound)
    pressure_diffusion_scaling = (1.0 - limited_mach_number)^2

    density_weighting = -max(min(mach_number_a, 0.0), -1.0) *
        min(max(mach_number_b, 0.0), 1.0)
    mean_normal_speed = (
        density_a * abs(normal_velocity_a) +
        density_b * abs(normal_velocity_b)
    ) / (density_a + density_b)
    normal_speed_a = (1.0 - density_weighting) * mean_normal_speed +
        density_weighting * abs(normal_velocity_a)
    normal_speed_b = (1.0 - density_weighting) * mean_normal_speed +
        density_weighting * abs(normal_velocity_b)

    mass_flux = 0.5 * (
        density_a * (normal_velocity_a + normal_speed_a) +
        density_b * (normal_velocity_b - normal_speed_b) -
        pressure_diffusion_scaling / interface_speed_of_sound *
            (pressure_b - pressure_a)
    )

    pressure_weight_a = _slau2_pressure_weight_left(mach_number_a)
    pressure_weight_b = _slau2_pressure_weight_right(mach_number_b)
    mean_density = 0.5 * (density_a + density_b)

    # SLAU2 differs from SLAU in this final term: SLAU uses a pressure-based
    # factor, whereas SLAU2 uses rms_speed * rho_bar * c_bar.
    pressure_flux = 0.5 * (pressure_a + pressure_b) +
        0.5 * (pressure_weight_a - pressure_weight_b) *
            (pressure_a - pressure_b) +
        rms_speed * (pressure_weight_a + pressure_weight_b - 1.0) *
            mean_density * interface_speed_of_sound

    return mass_flux, pressure_flux
end

function SLAU2!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b,
    face_reconstructor!,
)
    (
        distance,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    (
        density_a,
        momentum_density_u_a,
        momentum_density_v_a,
        momentum_density_w_a,
        volumetric_energy_a,
    ) = face_reconstructor!(
        du, u, p, t, system, geo,
        idx_a, face_a,
        idx_b, face_b,
    )
    (
        density_b,
        momentum_density_u_b,
        momentum_density_v_b,
        momentum_density_w_b,
        volumetric_energy_b,
    ) = face_reconstructor!(
        du, u, p, t, system, geo,
        idx_b, face_b,
        idx_a, face_a,
    )

    (
        density_flux,
        momentum_density_u_flux,
        momentum_density_v_flux,
        momentum_density_w_flux,
        volumetric_energy_flux,
    ) = slau2_flux(
        density_a,
        momentum_density_u_a,
        momentum_density_v_a,
        momentum_density_w_a,
        volumetric_energy_a,
        u.cp[idx_a] / u.cv[idx_a],
        density_b,
        momentum_density_u_b,
        momentum_density_v_b,
        momentum_density_w_b,
        volumetric_energy_b,
        u.cp[idx_b] / u.cv[idx_b],
        face_normal_a,
    )

    du.density_flow[idx_a] -= face_area_a * density_flux
    du.momentum_density_u_flow[idx_a] -= face_area_a * momentum_density_u_flux
    du.momentum_density_v_flow[idx_a] -= face_area_a * momentum_density_v_flux
    du.momentum_density_w_flow[idx_a] -= face_area_a * momentum_density_w_flux
    du.volumetric_energy_flow[idx_a] -= face_area_a * volumetric_energy_flux

    du.density_flow[idx_b] += face_area_a * density_flux
    du.momentum_density_u_flow[idx_b] += face_area_a * momentum_density_u_flux
    du.momentum_density_v_flow[idx_b] += face_area_a * momentum_density_v_flux
    du.momentum_density_w_flow[idx_b] += face_area_a * momentum_density_w_flux
    du.volumetric_energy_flow[idx_b] += face_area_a * volumetric_energy_flux

    if hasproperty(u, :species_densities)
        add_species_advection_flux!(
            du,
            u,
            idx_a,
            idx_b,
            face_area_a,
            density_flux,
        )
    end
    if hasproperty(u, :turbulent_kinetic_energy_density)
        add_sst_advection_flux!(
            du,
            u,
            idx_a,
            idx_b,
            face_area_a,
            density_flux,
        )
    end
    return nothing
end
