function cap_navier_stokes_flow!(du, u, p, t, system, geo, cell_id)
    cell_volume = cell_geometry(geo, cell_id)
    du.density[cell_id] += du.density_flow[cell_id] / cell_volume

    du.momentum_density_u[cell_id] += du.momentum_density_u_flow[cell_id] / cell_volume
    du.momentum_density_v[cell_id] += du.momentum_density_v_flow[cell_id] / cell_volume
    du.momentum_density_w[cell_id] += du.momentum_density_w_flow[cell_id] / cell_volume
    
    du.volumetric_energy[cell_id] += du.volumetric_energy_flow[cell_id] / cell_volume
end