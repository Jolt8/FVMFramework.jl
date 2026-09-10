"""
    pressure_driven_mass_flux!(du, u, idx_a, idx_b, face_idx, area, norm, dist)

Calculates the mass flow rate across a face in a porous medium using Darcy's Law.

Parameters:
    - du: the derivative of the state variable vector
    - u: the state variable vector
    - idx_a: the index of the first cell
    - idx_b: the index of the second cell
    - face_idx: the index of the face
    - area: the area of the face
    - norm: the normal of the face
    - dist: the distance between the two cells
    - rho_avg: average density of the two cells sharing the face in [kg/m^3]
    - permeability: permeability of the porous medium in [m^2]
    - viscosity: viscosity of the fluid in [Pa*s]
    - pressure_a: pressure in the first cell in [Pa]
    - pressure_b: pressure in the second cell in [Pa]
    - area: area of the face in [m^2]
    - dist: distance between the two cells in [m]

Required variables in u:
    - pressure: pressure in [Pa]
    - rho: density in [kg/m^3]
    - mu: dynamic viscosity in [Pa*s]
    - permeability: permeability of the porous medium in [m^2]

Required in du:
    - mass_face: mass flow rate across faces in [kg/s]

Returns:
    None
end
"""
function pressure_driven_mass_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist
)
    rho_avg = 0.5 * (u.rho[idx_a] + u.rho[idx_b])
    mu_avg = 0.5 * (u.mu[idx_a] + u.mu[idx_b])
    permeability_avg = 0.5 * (u.permeability[idx_a] + u.permeability[idx_b])

    p_grad = (u.pressure[idx_b] - u.pressure[idx_a]) / dist
    face_m_dot = -rho_avg * (permeability_avg / mu_avg) * p_grad * area

    du.mass_face[idx_a, face_idx] -= face_m_dot
end

