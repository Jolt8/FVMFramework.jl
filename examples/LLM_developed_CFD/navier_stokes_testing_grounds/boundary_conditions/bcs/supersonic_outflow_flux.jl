function supersonic_outlet_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    face_area, face_normal, face_distance, vol = boundary_geometry(geo, idx_a, face_a)

    vel_u, vel_v, vel_w, pressure, speed_of_sound =
        primitive_from_conservative(
            u.density[idx_a],
            u.momentum_density_u[idx_a],
            u.momentum_density_v[idx_a],
            u.momentum_density_w[idx_a],
            u.volumetric_energy[idx_a],
            u.gamma[idx_a],
        )

    (
        F_density,
        F_momentum_density_u,
        F_momentum_density_v,
        F_momentum_density_w,
        F_volumetric_energy,
    ) = physical_flux(
        u.density[idx_a],
        u.momentum_density_u[idx_a],
        u.momentum_density_v[idx_a],
        u.momentum_density_w[idx_a],
        u.volumetric_energy[idx_a],
        pressure,
        face_normal
    )

    du.density_flow[idx_a] -= face_area * F_density

    du.momentum_density_u_flow[idx_a] -= face_area * F_momentum_density_u
    du.momentum_density_v_flow[idx_a] -= face_area * F_momentum_density_v
    du.momentum_density_w_flow[idx_a] -= face_area * F_momentum_density_w

    du.volumetric_energy_flow[idx_a] -= face_area * F_volumetric_energy

    if system.compressible_solver_options.use_sst_turbulence == true
        add_boundary_sst_advection_flux!(du, u, idx_a, face_area, F_density)
    end
end