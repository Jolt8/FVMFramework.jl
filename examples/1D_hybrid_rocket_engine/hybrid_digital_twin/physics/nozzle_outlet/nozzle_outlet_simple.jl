function nozzle_outlet_simple!(
    du, u, p, t, system, geo, 
    idx_a, face_a, idx_b, face_b,
)
    # This is a boundary face: idx_b and face_b are both 0, so there is no neighbour to look up geometry for
    # (interface_geometry would index cell 0). The throat area comes from the CAD data instead.
    nozzle_throat_area = system.cad_data.nozzle_throat_area

    gamma = u.cp[idx_a] / u.cv[idx_a]

    # Specific gas constant of the mixture in the cell: R = R_u / MW (J/(kg*K)), with mw_avg in kg/mol
    universal_gas_constant = 8.314462618

    specific_gas_constant = universal_gas_constant / u.mw_avg[idx_a]

    chamber_pressure = u.pressure[idx_a]

    # The port cell is treated as the stagnation state ahead of the throat. Below ambient pressure nothing
    # flows out (backflow into the motor isn't modelled), and the square root below would be undefined.
    if chamber_pressure <= p.ambient_pressure
        return nothing
    end

    # Isentropic throat Mach number for expansion from chamber pressure down to ambient. It can't exceed 1 at
    # the throat: once the pressure ratio passes the critical ratio, ((γ + 1) / 2)^(γ / (γ - 1)) ≈ 1.8, the throat
    # chokes and the Mach number there stays at 1. Without the cap, the formula below would give a falling mass
    # flow as the chamber pressure rises.
    unchoked_throat_mach_number_squared = (2 / (gamma - 1)) * ((chamber_pressure / p.ambient_pressure)^((gamma - 1) / gamma) - 1)

    throat_mach_number = min(sqrt(unchoked_throat_mach_number_squared), 1.0)

    mass_flow = p.nozzle_discharge_coefficient * nozzle_throat_area *
        chamber_pressure / sqrt(u.temperature[idx_a]) *
        sqrt(gamma / specific_gas_constant) *
        throat_mach_number * (1 + ((gamma - 1) / 2) * throat_mach_number^2)^(-(gamma + 1) / (2 * (gamma - 1)))


    #=
    #Old equation which always assumed supersonic flow:
    mass_flow = p.nozzle_discharge_coefficient * nozzle_throat_area * u.pressure[idx_a] *
        sqrt(
            (k / u.R_gas[idx_a] * u.temperature[idx_a]) * 
            (2 / (gamma + 1)) ^ ((k + 1) / (k-1))
    )
    =#

    # idx_a is the port cell at the nozzle end, whose volume grows as the port regresses
    du.density[idx_a] -= mass_flow / u.cell_volume[idx_a]

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