"""
    mass_fraction_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        area, norm, dist,
    )

Computes the mass fraction diffusion flux between two adjacent cells.

Parameters:
    - du: the derivative of the state variable vector
    - u: the state variable vector
    - idx_a: the index of the first cell
    - idx_b: the index of the second cell
    - face_idx: the index of the face
    - area: the area of the face
    - norm: the normal of the face pointing away from cell a and towards cell b
    - dist: the distance between the two cells

Required variables in u:
    - rho: density in [kg/m^3]
    - mass_fractions: mass fractions of each species in [kg/kg]
    - diffusion_coefficients: diffusion coefficients of each species in [m^2/s]

Required in du:
    - mass_face: mass flow rate across faces in [kg/s]

Returns:
    - nothing

"""
function mass_fraction_diffusion!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
)
    rho_avg = 0.5 * (u.rho[idx_a] + u.rho[idx_b])

    for_fields!(u.mass_fractions, du.mass_fractions, u.diffusion_coefficients) do species, mass_fractions, du_mass_fractions, diffusion_coefficients

        diffusion_coeff_effective = 
        diffusion_coeff_effective = harmonic_mean(diffusion_coefficients[species[idx_a]], diffusion_coefficients[species[idx_b]])
        
        concentration_gradient = (mass_fractions[species[idx_b]] - mass_fractions[species[idx_a]]) / dist
        diffusive_flux = -rho_avg * diffusion_coeff_effective * concentration_gradient * area
        
        du_mass_fractions[species[idx_a]] -= diffusive_flux
    end
end