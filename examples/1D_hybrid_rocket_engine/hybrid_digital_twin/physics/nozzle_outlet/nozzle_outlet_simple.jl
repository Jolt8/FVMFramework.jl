function nozzle_outlet_simple!(
    du, u, p, t, system, geo, 
    idx_a, face_a, idx_b, face_b,
)
    (
        dist,
        face_area_a, face_normal_a, face_distance_a, vol_a,
        face_area_b, face_normal_b, face_distance_b, vol_b
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    nozzle_throat_area = system.cad_data.nozzle_throat_area

    gamma = u.cp[idx_a] / u.cv[idx_a]

    mach_number = sqrt(
        (2 / (gamma - 1)) * 
        ((u.pressure[idx_a] / p.ambient_pressure)^((gamma - 1) / gamma) - 1)
    )

    mass_flow = p.nozzle_throat_area * 
        u.pressure[idx_a] / sqrt(u.temperature[idx_a]) * 
        sqrt(gamma / R_gas[idx_a]) *
        mach_number * (1 + ((gamma - 1) / 2) * mach_number^2)^(-(gamma + 1) / (2 * (gamma - 1)))


    #=
    #Old equation which always assumed supersonic flow:
    mass_flow = p.nozzle_discharge_coefficient * nozzle_throat_area * u.pressure[idx_a] *
        sqrt(
            (k / u.R_gas[idx_a] * u.temperature[idx_a]) * 
            (2 / (gamma + 1)) ^ ((k + 1) / (k-1))
    )
    =#

    du.density[idx_a] -= mass_flow / geo.cell_volumes[idx_a]

    # Energy, on a sensible basis (e = cv * T, h = cp * T). Heats of reaction are added
    # separately by the combustion functions.
    # Upwind cell:   d(M e)/dt = -ṁ h,  dM/dt = -ṁ   ->   M de/dt = -ṁ (h_up - e_up)
    #   The leaving gas does flow work pushing out, which cools the donor (tank blowdown cooling).
    # Downwind cell: d(M e)/dt = +ṁ h_up, dM/dt = +ṁ  ->   M de/dt = ṁ (h_up - e_down)
    #   Incoming gas carries its enthalpy and compresses the cell (tank-filling heating).
    # Both are added to du.heat (W) and converted to temperature by cap_heat_flux_to_temp_change!
    upwind_specific_enthalpy = u.cp[idx_a] * u.temperature[idx_a]

    upwind_specific_internal_energy = u.cv[idx_a] * u.temperature[idx_a]

    upwind_flow_work_power = mass_flow * (upwind_specific_enthalpy - upwind_specific_internal_energy)

    du.heat[idx_a] -= upwind_flow_work_power
end