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

    # Positive when gas flows from cell a to cell b
    pressure_difference = u.pressure[idx_a] - u.pressure[idx_b]

    pressure_difference_magnitude = abs(pressure_difference)


    # Smooth upwinding. The gas crossing the face carries the donor cell's density, composition and enthalpy.
    # A hard switch of donor at ΔP = 0 makes the right-hand side jump there, and the laminar resistance below keeps
    # neighbouring port cells within a few Pa of each other, so faces sit near ΔP = 0 almost all the time. Every
    # flicker of direction then made Newton fail. Instead the donor is blended over a narrow band around ΔP = 0:
    #     φ_face = w φ_a + (1 - w) φ_b,     w = smoothstep(0.5 + ΔP / (2 ΔP_blend))
    # φ: density, mass fractions or enthalpy; w: weight of cell a (-), 1 for clear a -> b flow, 0 for clear b -> a;
    # ΔP_blend: blending width (Pa). Outside ±ΔP_blend the flux is exactly upwind, as before.
    upwind_blending_pressure_width = 1.0

    blend_position = clamp(0.5 + pressure_difference / (2 * upwind_blending_pressure_width), 0.0, 1.0)

    weight_a = blend_position^2 * (3.0 - 2.0 * blend_position)
    weight_b = 1.0 - weight_a


    # Minor loss from the area change, referenced to the velocity in the smaller area. Whether the step is an
    # expansion or a contraction depends on the flow direction, so the two coefficients are blended the same way.
    small_radius = min(radius_a, radius_b)
    large_radius = max(radius_a, radius_b)

    small_face_area = pi * small_radius^2

    area_ratio = (small_radius / large_radius)^2

    sudden_expansion_coefficient = (1.0 - area_ratio)^2      # Borda-Carnot
    sudden_contraction_coefficient = 0.5 * (1.0 - area_ratio) # common approximation for a sharp-edged step

    if abs(radius_a - radius_b) < 1e-6
        minor_loss_coefficient_a_to_b = 0.0
        minor_loss_coefficient_b_to_a = 0.0
    elseif radius_b > radius_a
        minor_loss_coefficient_a_to_b = sudden_expansion_coefficient
        minor_loss_coefficient_b_to_a = sudden_contraction_coefficient
    else
        minor_loss_coefficient_a_to_b = sudden_contraction_coefficient
        minor_loss_coefficient_b_to_a = sudden_expansion_coefficient
    end

    minor_loss_coefficient = weight_a * minor_loss_coefficient_a_to_b + weight_b * minor_loss_coefficient_b_to_a


    # Viscous resistance of a duct with the average face radius (Hagen-Poiseuille)
    average_radius = 0.5 * (radius_a + radius_b)

    # Turbulent wall friction between the two cell centroids, added to the minor loss:
    #     ΔP_friction = f (L / D) 0.5 ρ U^2
    # f: port_darcy_friction_factor (-), about 0.02 for a smooth tube at Re ~ 1e5; L: centroid-to-centroid distance (m);
    # D: average face diameter (m). At these Reynolds numbers this term carries most of the pressure drop (about 30 Pa per
    # port cell instead of ~1 Pa from the laminar term alone). Without it the laminar resistance was ~50x too weak, so
    # neighbouring port cells equalized pressure in nanoseconds and the port rang step to step after the valve opened.
    average_diameter = 2 * average_radius

    friction_loss_coefficient = p.port_darcy_friction_factor * dist / average_diameter

    minor_loss_coefficient = minor_loss_coefficient + friction_loss_coefficient

    average_viscosity = 0.5 * (u.mu[idx_a] + u.mu[idx_b])

    viscous_resistance = 8.0 * average_viscosity * dist / average_radius^2

    face_density = weight_a * u.density[idx_a] + weight_b * u.density[idx_b]

    # A trial stage the solver will reject can hand this flux a cell with no mass. sqrt() would throw on the
    # negative inertial term and end the run, so the flux is marked NaN instead: the step's error estimate becomes
    # NaN and the solver rejects it and retries with a smaller step. Nothing is clamped.
    if !(face_density > 0)
        not_a_number = oftype(face_density, NaN)

        du.density[idx_a] += not_a_number
        du.density[idx_b] += not_a_number

        return nothing
    end

    # Solve |ΔP| = R * U + 0.5 * ρ_face * K_L * U^2 for the positive root, written as
    #     U = 2 ΔP / (R + sqrt(R^2 + 2 ρ_face K_L |ΔP|))
    # so it stays accurate when the inertial term is small, carries the sign of ΔP (positive for a -> b), and
    # reduces to U = ΔP / R when K_L = 0.
    inertial_term = 2.0 * face_density * minor_loss_coefficient * pressure_difference_magnitude

    face_velocity = 2.0 * pressure_difference / (viscous_resistance + sqrt(viscous_resistance^2 + inertial_term))

    # Signed: positive from a to b
    face_mass_flow = face_density * face_velocity * small_face_area


    # Mass conservation (geo.cell_volumes already holds the current volumes of the moving port cells)
    du.density[idx_a] -= face_mass_flow / vol_a
    du.density[idx_b] += face_mass_flow / vol_b


    # Species. Each cell mixes toward the face composition at the flow crossing its face:
    #     cell a:  M_a dY_i,a/dt = -ṁ (Y_i,face - Y_i,a)
    #     cell b:  M_b dY_i,b/dt = +ṁ (Y_i,face - Y_i,b)
    # With clear a -> b flow, Y_face = Y_a, so cell a's composition is unchanged by its outflow and cell b mixes
    # toward Y_a, the same as the previous upwind form.
    gas_mass_a = u.density[idx_a] * vol_a
    gas_mass_b = u.density[idx_b] * vol_b

    for_fields!(du.mass_fractions, u.mass_fractions) do species, du_mass_fractions, u_mass_fractions
        mass_fraction_a = u_mass_fractions[species[idx_a]]
        mass_fraction_b = u_mass_fractions[species[idx_b]]

        face_mass_fraction = weight_a * mass_fraction_a + weight_b * mass_fraction_b

        du_mass_fractions[species[idx_a]] -= face_mass_flow * (face_mass_fraction - mass_fraction_a) / gas_mass_a
        du_mass_fractions[species[idx_b]] += face_mass_flow * (face_mass_fraction - mass_fraction_b) / gas_mass_b
    end


    # Energy, on a sensible basis (e = cv * T, h = cp * T). Heats of reaction are added separately by the combustion functions.
    #     cell a:  M_a de_a/dt = -ṁ (h_face - e_a)
    #     cell b:  M_b de_b/dt = +ṁ (h_face - e_b)
    # With clear a -> b flow, h_face = h_a: cell a loses its flow work, ṁ (h_a - e_a) (blowdown cooling), and cell b
    # gains ṁ (h_a - e_b) (filling heating), the same as the previous upwind form.
    # Both are added to du.heat (W) and converted to temperature by cap_heat_flux_to_temp_change!
    specific_enthalpy_a = u.cp[idx_a] * u.temperature[idx_a]
    specific_enthalpy_b = u.cp[idx_b] * u.temperature[idx_b]

    face_specific_enthalpy = weight_a * specific_enthalpy_a + weight_b * specific_enthalpy_b

    specific_internal_energy_a = u.cv[idx_a] * u.temperature[idx_a]
    specific_internal_energy_b = u.cv[idx_b] * u.temperature[idx_b]

    du.heat[idx_a] -= face_mass_flow * (face_specific_enthalpy - specific_internal_energy_a)
    du.heat[idx_b] += face_mass_flow * (face_specific_enthalpy - specific_internal_energy_b)

    #Old version: a hard upwind switch at ΔP = 0 (density, composition, enthalpy and the expansion/contraction
    #coefficient all jumped when the flow direction flipped). See git history before this change for the full code.
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
    # The current volume comes from geo.cell_volumes, so this stays correct for cells whose volume changes.
    pressure_driven_mass_flux_custom_geometry!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
        dist,
        face_area_a, geo.cell_volumes[idx_a],
        face_area_b, geo.cell_volumes[idx_b]
    )
end

