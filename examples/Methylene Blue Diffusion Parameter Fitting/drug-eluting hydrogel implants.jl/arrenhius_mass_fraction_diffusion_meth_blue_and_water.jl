function diffusion_arrenhius_equation(du, u, cell_id, temp_avg)
    return u.diffusion_pre_exponential_factor[1] * exp(-u.diffusion_activation_energy[1] / (R_gas * temp_avg))
end

#the reason these look so weird is that optimizing a the diffusion_pre_exponential_factor and diffusion_activation_energy
#for both water and methylene blue would give 4 parameters to optimize, which is too many for the data we have.
#so instead, we are only optimizing the diffusion_pre_exponential_factor and diffusion_activation_energy
#for methylene blue, and assuming that each m3 of methylene_blue diffused causes an equal volume of water to diffuse
#in the opposite direction.

function arrenhius_mass_fraction_diffusion_meth_blue_and_water!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    rho_avg = 0.5 * (u.rho[idx_a] + u.rho[idx_b])

    temp_avg = 0.5 * (u.temp[idx_a] + u.temp[idx_b])

    meth_blue_species_diffusion_coefficient = diffusion_arrenhius_equation(du, u, idx_a, temp_avg)

    meth_blue_concentration_gradient = (u.mass_fractions.methylene_blue[idx_b] - u.mass_fractions.methylene_blue[idx_a]) / dist
    meth_blue_diffusion = -rho_avg * meth_blue_species_diffusion_coefficient * meth_blue_concentration_gradient * area

    du.species_masses.methylene_blue[idx_a] -= meth_blue_diffusion
    du.species_masses.water[idx_a] += meth_blue_diffusion
end