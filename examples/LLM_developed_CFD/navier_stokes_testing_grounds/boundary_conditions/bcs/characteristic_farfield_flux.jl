
"""
    characteristic_farfield_flux!(...)

Open far-field boundary for subsonic or supersonic normal flow.  At subsonic
faces the outgoing acoustic invariant is taken from the interior and the
incoming invariant from the prescribed freestream.  Entropy and tangential
velocity are taken from the upwind side.  Supersonic outflow extrapolates the
interior state, while supersonic inflow prescribes the complete freestream
state.
"""
function characteristic_farfield_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b;
    farfield_conditions = farfield_conditions
)
    face_area, face_normal, face_distance, vol = boundary_geometry(geo, idx_a, face_a)

    vel_u, vel_v, vel_w, pressure, sound = primitive_from_conservative(
        u.density[idx_a],
        u.momentum_density_u[idx_a],
        u.momentum_density_v[idx_a],
        u.momentum_density_w[idx_a],
        u.volumetric_energy[idx_a],
        u.gamma[idx_a],
    )
    
    normal_vel =
        vel_u * face_normal[1] +
        vel_v * face_normal[2] +
        vel_w * face_normal[3]

    density_inf = farfield_conditions.density
    momentum_u_inf = farfield_conditions.momentum_density_u
    momentum_v_inf = farfield_conditions.momentum_density_v
    momentum_w_inf = farfield_conditions.momentum_density_w
    energy_inf = farfield_conditions.volumetric_energy
    
    vel_u_inf, vel_v_inf, vel_w_inf, pressure_inf, sound_inf = 
    primitive_from_conservative(
        density_inf,
        momentum_u_inf,
        momentum_v_inf,
        momentum_w_inf,
        energy_inf,
        u.gamma[idx_a],
    )

    normal_vel_inf =
        vel_u_inf * face_normal[1] +
        vel_v_inf * face_normal[2] +
        vel_w_inf * face_normal[3]

    if normal_vel >= sound
        # Every characteristic leaves the domain.
        density_b = u.density[idx_a]
        momentum_u_b = u.momentum_density_u[idx_a]
        momentum_v_b = u.momentum_density_v[idx_a]
        momentum_w_b = u.momentum_density_w[idx_a]
        energy_b = u.volumetric_energy[idx_a]
    elseif normal_vel <= -sound
        # Every characteristic enters the domain.
        density_b = density_inf
        momentum_u_b = momentum_u_inf
        momentum_v_b = momentum_v_inf
        momentum_w_b = momentum_w_inf
        energy_b = energy_inf
    else
        # R+ leaves the domain and R- enters it.
        invariant_plus = normal_vel + 2.0 * sound / (u.gamma[idx_a] - 1.0)
        invariant_minus = normal_vel_inf - 2.0 * sound_inf / (u.gamma[idx_a] - 1.0)
        normal_vel_b = 0.5 * (invariant_plus + invariant_minus)
        sound_b = max(
            0.25 * (u.gamma[idx_a] - 1.0) * (invariant_plus - invariant_minus),
            1e-8,
        )

        if normal_vel_b >= 0.0
            # Subsonic outflow: entropy and tangential velocity leave the domain.
            entropy_b = pressure / u.density[idx_a]^u.gamma[idx_a]
            tangential_u = vel_u - normal_vel * face_normal[1]
            tangential_v = vel_v - normal_vel * face_normal[2]
            tangential_w = vel_w - normal_vel * face_normal[3]
        else
            # Subsonic inflow: entropy and tangential velocity enter from freestream.
            entropy_b = pressure_inf / density_inf^u.gamma[idx_a]
            tangential_u = vel_u_inf - normal_vel_inf * face_normal[1]
            tangential_v = vel_v_inf - normal_vel_inf * face_normal[2]
            tangential_w = vel_w_inf - normal_vel_inf * face_normal[3]
        end

        density_b = (sound_b^2 / (u.gamma[idx_a] * entropy_b))^(1.0 / (u.gamma[idx_a] - 1.0))
        pressure_b = entropy_b * density_b^u.gamma[idx_a]
        vel_u_b = tangential_u + normal_vel_b * face_normal[1]
        vel_v_b = tangential_v + normal_vel_b * face_normal[2]
        vel_w_b = tangential_w + normal_vel_b * face_normal[3]
        momentum_u_b = density_b * vel_u_b
        momentum_v_b = density_b * vel_v_b
        momentum_w_b = density_b * vel_w_b
        energy_b = pressure_b / (u.gamma[idx_a] - 1.0) + 0.5 * density_b * (vel_u_b^2 + vel_v_b^2 + vel_w_b^2)
    end

    (
        F_density,
        F_momentum_density_u,
        F_momentum_density_v,
        F_momentum_density_w,
        F_volumetric_energy,
    ) = slau2_flux(
        u.density[idx_a],
        u.momentum_density_u[idx_a],
        u.momentum_density_v[idx_a],
        u.momentum_density_w[idx_a],
        u.volumetric_energy[idx_a],
        u.gamma[idx_a],
        density_b,
        momentum_u_b,
        momentum_v_b,
        momentum_w_b,
        energy_b,
        u.gamma[idx_a],
        face_normal,
    )

    du.density_flow[idx_a] -= face_area * F_density
    du.momentum_density_u_flow[idx_a] -= face_area * F_momentum_density_u
    du.momentum_density_v_flow[idx_a] -= face_area * F_momentum_density_v
    du.momentum_density_w_flow[idx_a] -= face_area * F_momentum_density_w
    du.volumetric_energy_flow[idx_a] -= face_area * F_volumetric_energy

    if compressible_solver_options.use_sst_turbulence == true
        if F_density >= 0.0
            add_boundary_sst_advection_flux!(du, u, idx_a, face_area, F_density)
        else
            add_prescribed_boundary_sst_advection_flux!(
                du, idx_a, face_area, F_density,
                farfield_conditions.turbulent_kinetic_energy_density / farfield_conditions.density,
                farfield_conditions.specific_dissipation_rate_density / farfield_conditions.density,
            )
        end
    end
end