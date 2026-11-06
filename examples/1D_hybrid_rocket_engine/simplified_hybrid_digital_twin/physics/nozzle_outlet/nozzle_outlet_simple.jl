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

    # A Newton iterate the solver will reject can hand this cell a negative temperature or gas constant, and the
    # square roots below would throw and end the run. Marking the outflow NaN makes the solver reject the iterate and
    # shorten the step instead. Written as !(x > 0) so NaN is caught too.
    if !(u.temperature[idx_a] > 0) || !(gamma / specific_gas_constant > 0)
        du.density[idx_a] += oftype(chamber_pressure, NaN)
        return nothing
    end

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

    # sqrt(M^2) has an infinite slope at M^2 = 0, which is exactly where the chamber sits before the valve opens
    # (chamber pressure = ambient). ForwardDiff turns that into NaN across the whole Jacobian row and the implicit
    # solver stops before its first step. The same regularization as the orifice equations is used:
    #     M = M^2 / ((M^2)^2 + ε^2)^(1/4)
    # which equals sqrt(M^2) to within 0.5% once M^2 > 7ε, and is linear below ε. ε = 0.01 corresponds to a chamber
    # only about 0.7 kPa above ambient (M^2 ≈ (2/γ) ΔP / P_ambient), so it never touches a real operating point.
    mach_number_squared_regularization = 0.01

    regularized_throat_mach_number = unchoked_throat_mach_number_squared / (unchoked_throat_mach_number_squared^2 + mach_number_squared_regularization^2)^0.25

    throat_mach_number = min(regularized_throat_mach_number, 1.0)

    #throat_mach_number = min(sqrt(unchoked_throat_mach_number_squared), 1.0) #old version, infinite slope at chamber pressure = ambient

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