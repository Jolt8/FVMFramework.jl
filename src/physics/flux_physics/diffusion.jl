"""
    species_numerical_flux(
        density_a, density_b,
        mass_fraction_a, mass_fraction_b,
        diffusion_coefficient_a, diffusion_coefficient_b,
        area, distance,
    )

Return the integrated Fickian species flux in the normal direction from cell
`a` toward cell `b`. Density is arithmetically averaged and the diffusion
coefficient is harmonically averaged. A positive result is transport from
`a` to `b`.
"""
function species_numerical_flux(
    density_a,
    density_b,
    mass_fraction_a,
    mass_fraction_b,
    diffusion_coefficient_a,
    diffusion_coefficient_b,
    area,
    distance,
)
    if distance <= 0.0
        throw(ArgumentError("distance must be positive"))
    end
    if area < 0.0
        throw(ArgumentError("area must be non-negative"))
    end
    if density_a < 0.0 || density_b < 0.0
        throw(ArgumentError("densities must be non-negative"))
    end
    if diffusion_coefficient_a < 0.0 || diffusion_coefficient_b < 0.0
        throw(ArgumentError("diffusion coefficients must be non-negative"))
    end

    density_average = 0.5 * (density_a + density_b)
    if diffusion_coefficient_a == 0.0 || diffusion_coefficient_b == 0.0
        diffusion_coefficient_effective = 0.0
    else
        diffusion_coefficient_effective = harmonic_mean(
            diffusion_coefficient_a,
            diffusion_coefficient_b,
        )
    end
    mass_fraction_gradient = (mass_fraction_b - mass_fraction_a) / distance
    return -density_average * diffusion_coefficient_effective * mass_fraction_gradient * area
end

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
    - mass_fractions: oriented integrated species-flow accumulator

This is an oriented-cell update: call it once for each side of an internal
face, with the cell indices reversed on the second call, when the surrounding
operator requires equal-and-opposite cell contributions.

Returns:
    - nothing

"""
function mass_fraction_diffusion!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
)
    for_fields!(u.mass_fractions, du.mass_fractions, u.diffusion_coefficients) do species, mass_fractions, du_mass_fractions, diffusion_coefficients
        diffusive_flux = species_numerical_flux(
            u.rho[idx_a],
            u.rho[idx_b],
            mass_fractions[species[idx_a]],
            mass_fractions[species[idx_b]],
            diffusion_coefficients[species[idx_a]],
            diffusion_coefficients[species[idx_b]],
            area,
            dist,
        )
        du_mass_fractions[species[idx_a]] -= diffusive_flux
    end
end
