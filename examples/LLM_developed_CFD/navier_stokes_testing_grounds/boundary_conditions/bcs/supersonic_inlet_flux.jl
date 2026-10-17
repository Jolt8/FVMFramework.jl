function supersonic_inlet_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b;
    inlet_conditions = inlet_conditions
)
    face_area, face_normal, face_distance, vol = boundary_geometry(geo, idx_a, face_a)

    density_bc = inlet_conditions.density
    momentum_density_u_bc = inlet_conditions.momentum_density_u
    momentum_density_v_bc = inlet_conditions.momentum_density_v
    momentum_density_w_bc = inlet_conditions.momentum_density_w
    volumetric_energy_bc = inlet_conditions.volumetric_energy

    vel_u_bc, vel_v_bc, vel_w_bc, pressure_bc, speed_of_sound_bc =
        primitive_from_conservative(
            density_bc,
            momentum_density_u_bc,
            momentum_density_v_bc,
            momentum_density_w_bc,
            volumetric_energy_bc,
            u.gamma[idx_a],
        )

    (
        F_density,
        F_momentum_density_u,
        F_momentum_density_v,
        F_momentum_density_w,
        F_volumetric_energy,
    ) = physical_flux(
        density_bc,
        momentum_density_u_bc,
        momentum_density_v_bc,
        momentum_density_w_bc,
        volumetric_energy_bc,
        pressure_bc,
        face_normal,
    )

    du.density_flow[idx_a] -= face_area * F_density

    du.momentum_density_u_flow[idx_a] -= face_area * F_momentum_density_u
    du.momentum_density_v_flow[idx_a] -= face_area * F_momentum_density_v
    du.momentum_density_w_flow[idx_a] -= face_area * F_momentum_density_w

    du.volumetric_energy_flow[idx_a] -= face_area * F_volumetric_energy

    if system.compressible_solver_options.use_sst_turbulence == true
        add_prescribed_boundary_sst_advection_flux!(
            du, idx_a, face_area,
            F_density,
            inlet_conditions.turbulent_kinetic_energy_density / inlet_conditions.density,
            inlet_conditions.specific_dissipation_rate_density / inlet_conditions.density,
        )
    end
end