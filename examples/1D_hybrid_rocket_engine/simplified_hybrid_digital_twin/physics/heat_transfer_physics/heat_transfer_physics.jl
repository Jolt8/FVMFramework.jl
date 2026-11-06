function heat_diffusion!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b
)
    (
        dist,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

    grad_T = (u.temperature[idx_b] - u.temperature[idx_a]) / dist

    du.heat[idx_a] -= -k_effective * grad_T * face_area_a
    du.heat[idx_b] += -k_effective * grad_T * face_area_b
end

"""
Conduction across an interface as two thermal resistances in series, one per
half-cell:

    R_a = d_a / (k_a * A_a)        R_b = d_b / (k_b * A_b)        Q = (T_b - T_a) / (R_a + R_b)

R_a, R_b: half-cell thermal resistances (K/W); d_a, d_b: half-cell distances
from each cell node to the face (m, `cell_face_distances`); k_a, k_b: thermal
conductivities (W/m/K); A_a, A_b: face areas (m^2); Q: heat flow from b into a (W).

Unlike `heat_diffusion!` (harmonic mean of k over the total distance) this
is exact when the two halves differ, e.g. a 0.8 mm liner against a 3.5 mm
housing wall. With d_a = d_b and A_a = A_b both functions agree. For a
coaxial cylindrical face of a ROM export, d is d_eq = r_f |ln(r_f / r_c)|,
which makes this flat-wall form exact for a cylindrical shell.
"""

function heat_diffusion_series_resistance_custom_geometry!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b,
    face_area_a, face_distance_a,
    face_area_b, face_distance_b,
)
    half_cell_resistance_a = face_distance_a / (u.k[idx_a] * face_area_a)
    half_cell_resistance_b = face_distance_b / (u.k[idx_b] * face_area_b)

    conduction_heat_flow = (u.temperature[idx_b] - u.temperature[idx_a]) / (half_cell_resistance_a + half_cell_resistance_b)

    du.heat[idx_a] += conduction_heat_flow
    du.heat[idx_b] -= conduction_heat_flow
end

function heat_diffusion_series_resistance!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b
)
    (
        dist,
        face_area_a, face_normal_a, face_distance_a, volume_a,
        face_area_b, face_normal_b, face_distance_b, volume_b,
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    heat_diffusion_series_resistance_custom_geometry!(
        du, u, p, t, system, geo,
        idx_a, face_a,
        idx_b, face_b,
        face_area_a, face_distance_a,
        face_area_b, face_distance_b,
    )
end