function slip_wall_flux!(
    du, u, p, t, system, geo, 
    idx_a, face_a, idx_b, face_b
)
    face_area, face_normal, _, _ = boundary_geometry(geo, idx_a, face_a)
    
    vel_u, vel_v, vel_w, pressure, speed_of_sound = primitive_from_conservative(
        u.density[idx_a],
        u.momentum_density_u[idx_a],
        u.momentum_density_v[idx_a],
        u.momentum_density_w[idx_a],
        u.volumetric_energy[idx_a],
        u.gamma[idx_a],
    )
    du.momentum_density_u_flow[idx_a] -= face_area * pressure * face_normal[1]
    du.momentum_density_v_flow[idx_a] -= face_area * pressure * face_normal[2]
    du.momentum_density_w_flow[idx_a] -= face_area * pressure * face_normal[3]
end

function symmetry_flux!(
    du, u, p, t, system, geo, 
    idx_a, face_a, idx_b, face_b
)
    #identical to a slip wall flux 
    slip_wall_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
end

