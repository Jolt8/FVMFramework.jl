"""
    update_species_mass_fractions!(du, u, p, t, system, geo, cell_id)

Recover the primitive species mass fractions from the conservative species
densities `rho * Y_k`. Species thermodynamics are currently passive: the
mixture heat capacities are still supplied by the fluid properties.
"""
function update_species_mass_fractions!(du, u, p, t, system, geo, cell_id)
    density = max(u.density[cell_id], 1e-10)
    foreach_field_at!(cell_id, u.species_densities, u.mass_fractions) do local_cell_id, species_densities, mass_fractions
        mass_fractions[local_cell_id] = species_densities[local_cell_id] / density
    end
    return nothing
end

"""
    add_species_advection_flux!(du, u, idx_a, idx_b, area, density_flux)

Advect every conservative species density with a numerical mixture mass flux.
The mass fraction is selected from the upwind side of the face. Consequently,
when the cell fractions sum to one, the species fluxes sum exactly to the
mixture density flux.
"""
function add_species_advection_flux!(
    du,
    u,
    idx_a,
    idx_b,
    area,
    density_flux,
)
    foreach_field_at!(idx_a, u.mass_fractions, du.species_density_flow) do local_cell_id, mass_fractions, species_density_flow
        if density_flux >= 0.0
            face_mass_fraction = mass_fractions[idx_a]
        else
            face_mass_fraction = mass_fractions[idx_b]
        end
        integrated_species_flux = area * density_flux * face_mass_fraction
        species_density_flow[idx_a] -= integrated_species_flux
        species_density_flow[idx_b] += integrated_species_flux
    end
    return nothing
end

function add_hllc_species_advection_flux!(
    du,
    u,
    idx_a,
    idx_b,
    area,
    density_flux,
)
    add_species_advection_flux!(
        du,
        u,
        idx_a,
        idx_b,
        area,
        density_flux,
    )
    return nothing
end

"""
    add_boundary_species_advection_flux!(du, u, cell_id, area, density_flux)

Apply a one-sided convective species flux at a boundary. This is appropriate
for transmissive/outflow boundaries, where the boundary composition is the
adjacent cell composition.
"""
function add_boundary_species_advection_flux!(
    du,
    u,
    cell_id,
    area,
    density_flux,
)
    foreach_field_at!(cell_id, u.mass_fractions, du.species_density_flow) do local_cell_id, mass_fractions, species_density_flow
        species_density_flow[local_cell_id] -=
            area * density_flux * mass_fractions[local_cell_id]
    end
    return nothing
end

"""
    add_prescribed_boundary_species_advection_flux!(
        du, cell_id, area, density_flux, prescribed_mass_fractions,
    )

Apply a convective species flux using a prescribed boundary composition. The
fields of `prescribed_mass_fractions` must match `du.species_density_flow`.
"""
function add_prescribed_boundary_species_advection_flux!(
    du,
    cell_id,
    area,
    density_flux,
    prescribed_mass_fractions,
)
    foreach_field_at!(cell_id, prescribed_mass_fractions, du.species_density_flow) do local_cell_id, mass_fractions, species_density_flow
        species_density_flow[local_cell_id] -=
            area * density_flux * mass_fractions[local_cell_id]
    end
    return nothing
end

"""
    add_conservative_species_diffusion_flux!(
        du, u, idx_a, idx_b, area, distance,
    )

Add an equal-and-opposite mixture-corrected Fickian flux to two cells. The
correction `J_k <- J_k - Y_k * sum(J)` enforces `sum(J_k) = 0`, so diffusion
cannot change total mixture density when the face mass fractions sum to one.
"""
function add_conservative_species_diffusion_flux!(
    du,
    u,
    idx_a,
    idx_b,
    area,
    distance,
)
    uncorrected_flux_sum = 0.0
    foreach_field_at!(idx_a, u.mass_fractions, u.diffusion_coefficients) do local_cell_id, mass_fractions, diffusion_coefficients
        uncorrected_flux_sum += species_numerical_flux(
            u.density[idx_a],
            u.density[idx_b],
            mass_fractions[idx_a],
            mass_fractions[idx_b],
            diffusion_coefficients[idx_a],
            diffusion_coefficients[idx_b],
            area,
            distance,
        )
    end

    foreach_field_at!(idx_a, u.mass_fractions, u.diffusion_coefficients, du.species_density_flow) do local_cell_id, mass_fractions, diffusion_coefficients, species_density_flow
        uncorrected_flux = species_numerical_flux(
            u.density[idx_a],
            u.density[idx_b],
            mass_fractions[idx_a],
            mass_fractions[idx_b],
            diffusion_coefficients[idx_a],
            diffusion_coefficients[idx_b],
            area,
            distance,
        )
        face_mass_fraction = 0.5 * (
            mass_fractions[idx_a] + mass_fractions[idx_b]
        )
        corrected_flux = uncorrected_flux -
            face_mass_fraction * uncorrected_flux_sum
        species_density_flow[idx_a] -= corrected_flux
        species_density_flow[idx_b] += corrected_flux
    end
    return nothing
end

"""
    cap_species_density_flow!(du, u, p, t, system, geo, cell_id)

Convert integrated conservative species flows to the volume-normalized time
derivatives evolved by the ODE integrator.
"""
function cap_species_density_flow!(du, u, p, t, system, geo, cell_id)
    cell_volume = cell_geometry(geo, cell_id)
    foreach_field_at!(cell_id, du.species_densities, du.species_density_flow) do local_cell_id, species_densities, species_density_flow
        species_densities[local_cell_id] +=
            species_density_flow[local_cell_id] / cell_volume
    end
    return nothing
end
