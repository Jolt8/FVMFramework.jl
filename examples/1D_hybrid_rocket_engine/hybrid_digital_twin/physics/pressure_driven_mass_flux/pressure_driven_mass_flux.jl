"""
Starter quasi-steady mass flux between two low-speed fluid cells (pre-chamber,
port cells, post-combustion chamber). Not for the injector, valve, or nozzle:
those faces need choked or two-phase closures.

The face is treated as a short circular duct. The pressure difference is split
into a viscous (Hagen-Poiseuille) part and an inertial minor-loss part:

    |ΔP| = R * U + 0.5 * ρ_up * K_L * U^2        R = 8 * μ * Δx / r_avg^2

ΔP: pressure difference between the two cell centroids (Pa); U: face velocity
in the smaller of the two face areas (m/s); R: viscous resistance (Pa*s/m);
μ: average dynamic viscosity (Pa*s); Δx: centroid-to-centroid distance (m);
r_avg: average face radius (m); ρ_up: upwind gas density (kg/m^3);
K_L: minor-loss coefficient for a sudden expansion or contraction (-).

The positive root is solved in closed form. Density, composition and enthalpy
of the gas crossing the face are taken from the upwind (donor) cell.
"""
function pressure_driven_mass_flux_custom_geometry!(
    du, u, p, t, system, geo,
    idx_a, face_a, idx_b, face_b,
    dist, 
    face_area_a, vol_a, 
    face_area_b, vol_b
)
    radius_a = sqrt(face_area_a / pi)
    radius_b = sqrt(face_area_b / pi)

    # Gas flows from high to low pressure. The upwind cell is the one the gas leaves.
    pressure_difference = u.pressure[idx_a] - u.pressure[idx_b]

    if pressure_difference >= 0
        upwind_cell = idx_a
        downwind_cell = idx_b
        upwind_radius = radius_a
        downwind_radius = radius_b
    else
        upwind_cell = idx_b
        downwind_cell = idx_a
        upwind_radius = radius_b
        downwind_radius = radius_a
    end

    pressure_difference_magnitude = abs(pressure_difference)

    # Minor loss from the area change, referenced to the velocity in the smaller area
    small_radius = min(radius_a, radius_b)
    large_radius = max(radius_a, radius_b)

    small_face_area = pi * small_radius^2

    area_ratio = (small_radius / large_radius)^2

    if abs(radius_a - radius_b) < 1e-6
        minor_loss_coefficient = 0.0
    elseif downwind_radius > upwind_radius
        # Sudden expansion (Borda-Carnot)
        minor_loss_coefficient = (1.0 - area_ratio)^2
    else
        # Sudden contraction (common approximation for a sharp-edged step)
        minor_loss_coefficient = 0.5 * (1.0 - area_ratio)
    end

    # Viscous resistance of a duct with the average face radius (Hagen-Poiseuille)
    average_radius = 0.5 * (radius_a + radius_b)

    average_viscosity = 0.5 * (u.mu[idx_a] + u.mu[idx_b])

    viscous_resistance = 8.0 * average_viscosity * dist / average_radius^2

    # Solve |ΔP| = R * U + 0.5 * ρ_up * K_L * U^2 for the positive root.
    # Written as U = 2|ΔP| / (R + sqrt(R^2 + 2 ρ_up K_L |ΔP|)) so it stays accurate when the
    # inertial term is small, and reduces to U = |ΔP| / R when K_L = 0.
    upwind_density = u.density[upwind_cell]

    # A trial stage the solver will reject can hand this flux a cell with no mass. sqrt() would throw on the
    # negative inertial term and end the run, so the flux is marked NaN instead: the step's error estimate becomes
    # NaN and the solver rejects it and retries with a smaller step. Nothing is clamped.
    if !(upwind_density > 0)
        not_a_number = oftype(upwind_density, NaN)

        du.density[idx_a] += not_a_number
        du.density[idx_b] += not_a_number

        return nothing
    end

    inertial_term = 2.0 * upwind_density * minor_loss_coefficient * pressure_difference_magnitude

    face_velocity = 2.0 * pressure_difference_magnitude / (viscous_resistance + sqrt(viscous_resistance^2 + inertial_term))

    # Always >= 0, directed from the upwind cell to the downwind cell
    face_mass_flow = upwind_density * face_velocity * small_face_area

    # Mass conservation. Cell volumes are treated as fixed here; the -ρ dV/dt term for
    # regressing port cells belongs with the fuel regression geometry update.
    if upwind_cell == idx_a
        upwind_cell_vol = vol_a
        downwind_cell_vol = vol_b
    else
        upwind_cell_vol = vol_b
        downwind_cell_vol = vol_a
    end

    du.density[upwind_cell] -= face_mass_flow / upwind_cell_vol
    du.density[downwind_cell] += face_mass_flow / downwind_cell_vol

    # Species. Outflow leaves the upwind cell's composition unchanged. The downwind cell
    # mixes toward the upwind composition:
    #     dY_i/dt = (Y_i,upwind - Y_i,downwind) * ṁ / M_downwind
    downwind_gas_mass = u.density[downwind_cell] * downwind_cell_vol

    for_fields!(du.mass_fractions, u.mass_fractions) do species, du_mass_fractions, u_mass_fractions
        upwind_mass_fraction = u_mass_fractions[species[upwind_cell]]
        downwind_mass_fraction = u_mass_fractions[species[downwind_cell]]

        du_mass_fractions[species[downwind_cell]] += (upwind_mass_fraction - downwind_mass_fraction) * face_mass_flow / downwind_gas_mass
    end

    # Energy, on a sensible basis (e = cv * T, h = cp * T). Heats of reaction are added
    # separately by the combustion functions.
    # Upwind cell:   d(M e)/dt = -ṁ h,  dM/dt = -ṁ   ->   M de/dt = -ṁ (h_up - e_up)
    #   The leaving gas does flow work pushing out, which cools the donor (tank blowdown cooling).
    # Downwind cell: d(M e)/dt = +ṁ h_up, dM/dt = +ṁ  ->   M de/dt = ṁ (h_up - e_down)
    #   Incoming gas carries its enthalpy and compresses the cell (tank-filling heating).
    # Both are added to du.heat (W) and converted to temperature by cap_heat_flux_to_temp_change!
    upwind_specific_enthalpy = u.cp[upwind_cell] * u.temperature[upwind_cell]

    upwind_specific_internal_energy = u.cv[upwind_cell] * u.temperature[upwind_cell]

    downwind_specific_internal_energy = u.cv[downwind_cell] * u.temperature[downwind_cell]

    upwind_flow_work_power = face_mass_flow * (upwind_specific_enthalpy - upwind_specific_internal_energy)

    downwind_enthalpy_inflow_power = face_mass_flow * (upwind_specific_enthalpy - downwind_specific_internal_energy)

    du.heat[upwind_cell] -= upwind_flow_work_power
    du.heat[downwind_cell] += downwind_enthalpy_inflow_power
end


function pressure_driven_mass_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, idx_b, face_b
)
    (
        dist,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    # Same argument order as the custom-geometry signature: each face area is followed by its cell's volume.
    # The current volume comes from u.cell_volume, so this stays correct for cells whose volume changes.
    pressure_driven_mass_flux_custom_geometry!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
        dist,
        face_area_a, u.cell_volume[idx_a],
        face_area_b, u.cell_volume[idx_b]
    )
end

