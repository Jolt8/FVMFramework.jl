
function construct_fluid_initial_conditions_from_intuitive_inputs(compressible_solver_options, u)
    u = upreferred.(u)
    momentum_density_u = u.density * u.vel_u
    momentum_density_v = u.density * u.vel_v
    momentum_density_w = u.density * u.vel_w

    internal_energy = u.cv * u.temperature
    kinetic_energy = 0.5 * (u.vel_u^2 + u.vel_v^2 + u.vel_w^2)
    
    volumetric_energy = u.density * (internal_energy + kinetic_energy)
    turbulent_kinetic_energy_density =
        u.density * u.turbulent_kinetic_energy
    specific_dissipation_rate_density =
        u.density * u.specific_dissipation_rate
    species_densities = ComponentVector(
        nitrous_oxide = u.density * u.mass_fractions.nitrous_oxide,
        ethylene = u.density * u.mass_fractions.ethylene,
        nitrogen = u.density * u.mass_fractions.nitrogen,
        oxygen = u.density * u.mass_fractions.oxygen,
        carbon_dioxide = u.density * u.mass_fractions.carbon_dioxide,
        carbon_monoxide = u.density * u.mass_fractions.carbon_monoxide,
        hydrogen = u.density * u.mass_fractions.hydrogen,
        water = u.density * u.mass_fractions.water,
    )
    
    properties_comp_vec = ComponentVector(
        cp = u.cp,
        cv = u.cv,
        R_gas = u.R_gas,
        mw = u.mw,
        mu = u.mu, 
        molecular_viscosity = u.mu,
        prandtl_number = u.prandtl_number,
        diffusion_coefficients_ambient = u.diffusion_coefficients_ambient,
        diffusion_coefficients_ambient_temperature = u.diffusion_coefficients_ambient_temperature,
        diffusion_coefficients_ambient_pressure = u.diffusion_coefficients_ambient_pressure
    )

    if compressible_solver_options.use_sst_turbulence == true
        return ustrip.(upreferred.(ComponentVector(
            density = u.density,
            momentum_density_u = momentum_density_u,
            momentum_density_v = momentum_density_v,
            momentum_density_w = momentum_density_w,
            volumetric_energy = volumetric_energy,
            turbulent_kinetic_energy_density = turbulent_kinetic_energy_density,
            specific_dissipation_rate_density = specific_dissipation_rate_density,
            species_densities = species_densities,
        ))), ustrip.(upreferred.(properties_comp_vec))
    else
        return ustrip.(upreferred.(ComponentVector(
            density = u.density,
            momentum_density_u = momentum_density_u,
            momentum_density_v = momentum_density_v,
            momentum_density_w = momentum_density_w,
            volumetric_energy = volumetric_energy,
            species_densities = species_densities,
        ))), ustrip.(upreferred.(properties_comp_vec))
    end
end

function construct_solid_initial_conditions_from_intuitive_inputs(u)
    """
    Requires:
    - density
    - temperature
    - cp
    - cv
    - k 
    """
    u = upreferred.(u)
    
    state_comp_vec = ComponentVector(
        density = u.density,
        momentum_density_u = 0.0,
        momentum_density_v = 0.0,
        momentum_density_w = 0.0,
        volumetric_energy = u.density * u.cv * u.temperature,
    )

    properties_comp_vec = ComponentVector(
        cp = u.cp,
        k = u.k
    )

    return ustrip.(upreferred.(state_comp_vec)), ustrip.(upreferred.(properties_comp_vec))
end