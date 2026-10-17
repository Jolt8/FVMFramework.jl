
function no_slip_adiabatic_wall_flux!(du, u, p, t, system, geo, cell_id, face_id)
    face_area, face_normal, face_distance, vol = boundary_geometry(geo, cell_id, face_id)

    vel_u, vel_v, vel_w, pressure, speed_of_sound = primitive_from_conservative(
        u.density[cell_id],
        u.momentum_density_u[cell_id],
        u.momentum_density_v[cell_id],
        u.momentum_density_w[cell_id],
        u.volumetric_energy[cell_id],
        u.gamma[cell_id],
    )

    du.momentum_density_u_flow[cell_id] -= face_area * pressure * face_normal[1]
    du.momentum_density_v_flow[cell_id] -= face_area * pressure * face_normal[2]
    du.momentum_density_w_flow[cell_id] -= face_area * pressure * face_normal[3]
    
    non_moving_wall_viscous_and_diffusive_flux!(du, u, p, t, system, geo, cell_id, face_id, 0, 0)
    
    if system.compressible_solver_options.use_sst_turbulence
        add_sst_wall_flux!(du, u, p, t, system, geo, cell_id, face_id)
    end
end