Base.@kwdef struct SSTKOmegaConstants
    beta_star::Float64 = 0.09
    a1::Float64 = 0.31
    beta_1::Float64 = 0.075
    beta_2::Float64 = 0.0828
    gamma_1::Float64 = 5.0 / 9.0
    gamma_2::Float64 = 0.44
    sigma_k_1::Float64 = 0.85
    sigma_k_2::Float64 = 1.0
    sigma_omega_1::Float64 = 0.5
    sigma_omega_2::Float64 = 0.856
    turbulent_prandtl::Float64 = 0.9
    production_limiter::Float64 = 10.0
end

const DEFAULT_SST_CONSTANTS = SSTKOmegaConstants()

function _sst_blend(first_value, second_value, blending_function)
    return blending_function * first_value + (1.0 - blending_function) * second_value
end

function sst_strain_rate_squared(gradient_u, gradient_v, gradient_w)
    strain_xx = gradient_u[1]
    strain_yy = gradient_v[2]
    strain_zz = gradient_w[3]
    strain_xy = 0.5 * (gradient_u[2] + gradient_v[1])
    strain_xz = 0.5 * (gradient_u[3] + gradient_w[1])
    strain_yz = 0.5 * (gradient_v[3] + gradient_w[2])
    return 2.0 * (
        strain_xx^2 + strain_yy^2 + strain_zz^2 +
        2.0 * (strain_xy^2 + strain_xz^2 + strain_yz^2)
    )
end

function sst_blending_functions(
    density,
    turbulent_kinetic_energy,
    specific_dissipation_rate,
    molecular_viscosity,
    wall_distance,
    gradient_k,
    gradient_omega;
    constants = DEFAULT_SST_CONSTANTS,
)
    tiny = 1e-14
    positive_density = max(density, tiny)
    positive_k = max(turbulent_kinetic_energy, tiny)
    positive_omega = max(specific_dissipation_rate, tiny)
    positive_distance = max(wall_distance, tiny)
    kinematic_viscosity = molecular_viscosity / positive_density
    cross_diffusion = max(
        2.0 * positive_density * constants.sigma_omega_2 *
        dot(gradient_k, gradient_omega) / positive_omega,
        1e-20,
    )
    argument_1 = min(
        max(
            sqrt(positive_k) /
            (constants.beta_star * positive_omega * positive_distance),
            500.0 * kinematic_viscosity /
            (positive_distance^2 * positive_omega),
        ),
        4.0 * positive_density * constants.sigma_omega_2 * positive_k /
        (cross_diffusion * positive_distance^2),
    )
    argument_2 = max(
        2.0 * sqrt(positive_k) /
        (constants.beta_star * positive_omega * positive_distance),
        500.0 * kinematic_viscosity /
        (positive_distance^2 * positive_omega),
    )
    return (
        F1 = tanh(argument_1^4),
        F2 = tanh(argument_2^2),
        cross_diffusion_coefficient = cross_diffusion,
    )
end

function sst_eddy_viscosity(
    density,
    turbulent_kinetic_energy,
    specific_dissipation_rate,
    strain_rate_squared,
    F2;
    constants = DEFAULT_SST_CONSTANTS,
)
    positive_k = max(turbulent_kinetic_energy, 0.0)
    positive_omega = max(specific_dissipation_rate, 1e-14)
    strain_rate = sqrt(max(strain_rate_squared, 0.0))
    denominator = max(constants.a1 * positive_omega, strain_rate * F2, 1e-14)
    return density * constants.a1 * positive_k / denominator
end

function sst_source_terms(
    density,
    turbulent_kinetic_energy,
    specific_dissipation_rate,
    turbulent_viscosity,
    strain_rate_squared,
    F1,
    gradient_k,
    gradient_omega;
    constants = DEFAULT_SST_CONSTANTS,
)
    beta = _sst_blend(constants.beta_1, constants.beta_2, F1)
    gamma = _sst_blend(constants.gamma_1, constants.gamma_2, F1)
    positive_k = max(turbulent_kinetic_energy, 0.0)
    positive_omega = max(specific_dissipation_rate, 1e-14)
    raw_production = turbulent_viscosity * strain_rate_squared
    production_limit = constants.production_limiter * constants.beta_star *
        density * positive_k * positive_omega
    production = min(raw_production, production_limit)
    k_destruction = constants.beta_star * density *
        positive_k * positive_omega
    if turbulent_viscosity > 1e-20
        omega_production = gamma * density * production / turbulent_viscosity
    else
        omega_production = 0.0
    end
    omega_destruction = beta * density * positive_omega^2
    cross_diffusion = 2.0 * (1.0 - F1) * density *
        constants.sigma_omega_2 * dot(gradient_k, gradient_omega) /
        positive_omega
    return (
        k_source = production - k_destruction,
        omega_source = omega_production - omega_destruction + cross_diffusion,
        production = production,
        k_destruction = k_destruction,
        omega_production = omega_production,
        omega_destruction = omega_destruction,
        cross_diffusion = cross_diffusion,
    )
end

function update_sst_primitives!(du, u, p, t, system, geo, cell_id)
    density = max(u.density[cell_id], 1e-14)
    u.turbulent_kinetic_energy[cell_id] =
        u.turbulent_kinetic_energy_density[cell_id] / density
    u.specific_dissipation_rate[cell_id] =
        u.specific_dissipation_rate_density[cell_id] / density
    return nothing
end

function update_sst_closure!(u, wall_distances; constants = DEFAULT_SST_CONSTANTS)
    for cell_id in eachindex(u.density)
        gradient_u = @view u.grad_vel_u[cell_id, :]
        gradient_v = @view u.grad_vel_v[cell_id, :]
        gradient_w = @view u.grad_vel_w[cell_id, :]
        gradient_k = @view u.grad_turbulent_kinetic_energy[cell_id, :]
        gradient_omega = @view u.grad_specific_dissipation_rate[cell_id, :]
        strain_squared = sst_strain_rate_squared(gradient_u, gradient_v, gradient_w)
        blending = sst_blending_functions(
            u.density[cell_id],
            u.turbulent_kinetic_energy[cell_id],
            u.specific_dissipation_rate[cell_id],
            u.molecular_viscosity[cell_id],
            wall_distances[cell_id],
            gradient_k,
            gradient_omega;
            constants = constants,
        )
        turbulent_viscosity = sst_eddy_viscosity(
            u.density[cell_id],
            u.turbulent_kinetic_energy[cell_id],
            u.specific_dissipation_rate[cell_id],
            strain_squared,
            blending.F2;
            constants = constants,
        )
        u.wall_distance[cell_id] = wall_distances[cell_id]
        u.strain_rate_squared[cell_id] = strain_squared
        u.sst_F1[cell_id] = blending.F1
        u.sst_F2[cell_id] = blending.F2
        u.turbulent_viscosity[cell_id] = turbulent_viscosity
        u.effective_dynamic_viscosity[cell_id] =
            u.molecular_viscosity[cell_id] + turbulent_viscosity
        u.effective_thermal_conductivity[cell_id] = u.k[cell_id] +
            turbulent_viscosity * u.cp[cell_id] / constants.turbulent_prandtl
    end
    return nothing
end

function add_hllc_sst_advection_flux!(du, u, idx_a, idx_b, area, density_flux)
    if density_flux >= 0.0
        face_k = u.turbulent_kinetic_energy[idx_a]
        face_omega = u.specific_dissipation_rate[idx_a]
    else
        face_k = u.turbulent_kinetic_energy[idx_b]
        face_omega = u.specific_dissipation_rate[idx_b]
    end
    integrated_k_flux = area * density_flux * face_k
    integrated_omega_flux = area * density_flux * face_omega
    du.turbulent_kinetic_energy_density_flow[idx_a] -= integrated_k_flux
    du.turbulent_kinetic_energy_density_flow[idx_b] += integrated_k_flux
    du.specific_dissipation_rate_density_flow[idx_a] -= integrated_omega_flux
    du.specific_dissipation_rate_density_flow[idx_b] += integrated_omega_flux
    return nothing
end

function add_boundary_sst_advection_flux!(du, u, cell_id, area, density_flux)
    du.turbulent_kinetic_energy_density_flow[cell_id] -=
        area * density_flux * u.turbulent_kinetic_energy[cell_id]
    du.specific_dissipation_rate_density_flow[cell_id] -=
        area * density_flux * u.specific_dissipation_rate[cell_id]
    return nothing
end

function add_prescribed_boundary_sst_advection_flux!(
    du,
    cell_id,
    area,
    density_flux,
    prescribed_k,
    prescribed_omega,
)
    du.turbulent_kinetic_energy_density_flow[cell_id] -=
        area * density_flux * prescribed_k
    du.specific_dissipation_rate_density_flow[cell_id] -=
        area * density_flux * prescribed_omega
    return nothing
end

function _sst_face_gradient(value_a, value_b, gradient_a, gradient_b, normal, distance)
    average_gradient = 0.5 .* (gradient_a .+ gradient_b)
    normal_derivative = (value_b - value_a) / distance
    return average_gradient .+ (
        normal_derivative - dot(average_gradient, normal)
    ) .* normal
end

function add_sst_diffusion_flux!(du, u, idx_a, idx_b, area, normal, distance;
    constants = DEFAULT_SST_CONSTANTS,
)
    F1_face = 0.5 * (u.sst_F1[idx_a] + u.sst_F1[idx_b])
    sigma_k = _sst_blend(constants.sigma_k_1, constants.sigma_k_2, F1_face)
    sigma_omega = _sst_blend(
        constants.sigma_omega_1,
        constants.sigma_omega_2,
        F1_face,
    )
    molecular_viscosity = harmonic_mean(
        u.molecular_viscosity[idx_a],
        u.molecular_viscosity[idx_b],
    )
    turbulent_viscosity = harmonic_mean(
        max(u.turbulent_viscosity[idx_a], 1e-30),
        max(u.turbulent_viscosity[idx_b], 1e-30),
    )
    gradient_k = _sst_face_gradient(
        u.turbulent_kinetic_energy[idx_a],
        u.turbulent_kinetic_energy[idx_b],
        @view(u.grad_turbulent_kinetic_energy[idx_a, :]),
        @view(u.grad_turbulent_kinetic_energy[idx_b, :]),
        normal,
        distance,
    )
    gradient_omega = _sst_face_gradient(
        u.specific_dissipation_rate[idx_a],
        u.specific_dissipation_rate[idx_b],
        @view(u.grad_specific_dissipation_rate[idx_a, :]),
        @view(u.grad_specific_dissipation_rate[idx_b, :]),
        normal,
        distance,
    )
    k_flux = (molecular_viscosity + sigma_k * turbulent_viscosity) *
        dot(gradient_k, normal) * area
    omega_flux = (molecular_viscosity + sigma_omega * turbulent_viscosity) *
        dot(gradient_omega, normal) * area
    du.turbulent_kinetic_energy_density_flow[idx_a] += k_flux
    du.turbulent_kinetic_energy_density_flow[idx_b] -= k_flux
    du.specific_dissipation_rate_density_flow[idx_a] += omega_flux
    du.specific_dissipation_rate_density_flow[idx_b] -= omega_flux
    return nothing
end

function add_sst_wall_flux!(du, u, geo, cell_id, face_id;
    constants = DEFAULT_SST_CONSTANTS,
)
    area, normal, face_distance, volume = boundary_geometry(geo, cell_id, face_id)
    density = u.density[cell_id]
    molecular_viscosity = u.molecular_viscosity[cell_id]
    kinematic_viscosity = molecular_viscosity / density
    wall_omega = 60.0 * kinematic_viscosity /
        (constants.beta_1 * face_distance^2)
    k_normal_gradient = -u.turbulent_kinetic_energy[cell_id] / face_distance
    omega_normal_gradient =
        (wall_omega - u.specific_dissipation_rate[cell_id]) / face_distance
    du.turbulent_kinetic_energy_density_flow[cell_id] +=
        area * molecular_viscosity * k_normal_gradient
    du.specific_dissipation_rate_density_flow[cell_id] +=
        area * molecular_viscosity * omega_normal_gradient
    return nothing
end

function cap_sst_transport!(du, u, geo, cell_id; constants = DEFAULT_SST_CONSTANTS)
    volume = cell_geometry(geo, cell_id)
    sources = sst_source_terms(
        u.density[cell_id],
        u.turbulent_kinetic_energy[cell_id],
        u.specific_dissipation_rate[cell_id],
        u.turbulent_viscosity[cell_id],
        u.strain_rate_squared[cell_id],
        u.sst_F1[cell_id],
        @view(u.grad_turbulent_kinetic_energy[cell_id, :]),
        @view(u.grad_specific_dissipation_rate[cell_id, :]);
        constants = constants,
    )
    du.turbulent_kinetic_energy_density[cell_id] +=
        du.turbulent_kinetic_energy_density_flow[cell_id] / volume +
        sources.k_source
    du.specific_dissipation_rate_density[cell_id] +=
        du.specific_dissipation_rate_density_flow[cell_id] / volume +
        sources.omega_source
    return nothing
end
