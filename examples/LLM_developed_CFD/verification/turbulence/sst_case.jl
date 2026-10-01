struct SSTVerificationFluid <: AbstractPhysics end

struct SSTVerificationCase{S, G, F}
    system::S
    geo::G
    rhs!::F
    u0::Vector{Float64}
    p::Float64
    gamma::Float64
    cv::Float64
end

function _sst_verification_internal_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, idx_b, face_b,
)
    HLLC!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
        first_order_face_reconstruction!,
        no_low_mach_correction,
    )
    fluid_viscous_and_diffusive_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
    )
    distance, area, normal, _, _, _, _, _, _ = interface_geometry(
        geo, idx_a, face_a, idx_b, face_b,
    )
    add_sst_diffusion_flux!(
        du, u, idx_a, idx_b, area, normal, distance,
    )
    return nothing
end

function _sst_verification_boundary_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, idx_b, face_b,
)
    area, normal, face_distance, volume = boundary_geometry(geo, idx_a, face_a)
    density_flux, momentum_u_flux, momentum_v_flux, momentum_w_flux, energy_flux =
        physical_flux(
            u.density[idx_a],
            u.momentum_density_u[idx_a],
            u.momentum_density_v[idx_a],
            u.momentum_density_w[idx_a],
            u.volumetric_energy[idx_a],
            u.pressure[idx_a],
            normal,
        )
    du.density_flow[idx_a] -= area * density_flux
    du.momentum_density_u_flow[idx_a] -= area * momentum_u_flux
    du.momentum_density_v_flow[idx_a] -= area * momentum_v_flux
    du.momentum_density_w_flow[idx_a] -= area * momentum_w_flux
    du.volumetric_energy_flow[idx_a] -= area * energy_flux
    add_boundary_sst_advection_flux!(du, u, idx_a, area, density_flux)

    turbulent_normal_stress = -(2.0 / 3.0) * u.density[idx_a] *
        u.turbulent_kinetic_energy[idx_a]
    traction_u = turbulent_normal_stress * normal[1]
    traction_v = turbulent_normal_stress * normal[2]
    traction_w = turbulent_normal_stress * normal[3]
    du.momentum_density_u_flow[idx_a] += area * traction_u
    du.momentum_density_v_flow[idx_a] += area * traction_v
    du.momentum_density_w_flow[idx_a] += area * traction_w
    du.volumetric_energy_flow[idx_a] += area * (
        u.vel_u[idx_a] * traction_u +
        u.vel_v[idx_a] * traction_v +
        u.vel_w[idx_a] * traction_w
    )
    return nothing
end

function _sst_verification_connection_map(physics_a, physics_b)
    if physics_a isa SSTVerificationFluid && physics_b isa SSTVerificationFluid
        return _sst_verification_internal_flux!
    end
    throw(ArgumentError("unsupported SST verification region pair"))
end

function build_sst_verification_case(n_cells; wall_distance = 1e-5)
    gamma = 1.4
    cv = 717.5
    cp = gamma * cv
    density = 1.0
    pressure = 1.0
    turbulent_kinetic_energy = 0.2
    specific_dissipation_rate = 3.0
    grid = generate_grid(
        Hexahedron,
        (n_cells, 1, 1),
        Ferrite.Vec{3}((0.0, 0.0, 0.0)),
        Ferrite.Vec{3}((1.0, 1.0, 1.0)),
    )
    addcellset!(grid, "fluid", xyz -> true)
    _add_boundary_sets!(grid)
    state_prototype = ComponentVector(
        density = zeros(n_cells)u"kg/m^3",
        momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
        momentum_density_v = zeros(n_cells)u"kg/(m^2*s)",
        momentum_density_w = zeros(n_cells)u"kg/(m^2*s)",
        volumetric_energy = zeros(n_cells)u"J/m^3",
        turbulent_kinetic_energy_density = zeros(n_cells)u"J/m^3",
        specific_dissipation_rate_density = zeros(n_cells)u"kg/(m^3*s)",
    )
    config = create_fvm_config(grid, state_prototype)
    add_setup_syms!(
        config;
        cache_syms_and_units = (
            vel_u = u"m/s",
            vel_v = u"m/s",
            vel_w = u"m/s",
            pressure = u"Pa",
            temperature = u"K",
            speed_of_sound = u"m/s",
            density_flow = u"kg/s",
            momentum_density_u_flow = u"N",
            momentum_density_v_flow = u"N",
            momentum_density_w_flow = u"N",
            volumetric_energy_flow = u"W",
            turbulent_kinetic_energy = u"m^2/s^2",
            specific_dissipation_rate = u"1/s",
            turbulent_viscosity = u"Pa*s",
            effective_dynamic_viscosity = u"Pa*s",
            effective_thermal_conductivity = u"W/(m*K)",
            wall_distance = u"m",
            strain_rate_squared = u"1/s^2",
            sst_F1 = u"1",
            sst_F2 = u"1",
            turbulent_kinetic_energy_density_flow = u"W",
            specific_dissipation_rate_density_flow = u"kg/s^2",
        ),
        special_caches = ComponentVector(
            grad_density = zeros(n_cells, 3)u"kg/m^4",
            grad_pressure = zeros(n_cells, 3)u"Pa/m",
            grad_vel_u = zeros(n_cells, 3)u"1/s",
            grad_vel_v = zeros(n_cells, 3)u"1/s",
            grad_vel_w = zeros(n_cells, 3)u"1/s",
            grad_temperature = zeros(n_cells, 3)u"K/m",
            grad_turbulent_kinetic_energy = zeros(n_cells, 3)u"m/s^2",
            grad_specific_dissipation_rate = zeros(n_cells, 3)u"1/(m*s)",
        ),
        second_order_syms = [],
        optimized_parameters = ComponentVector(),
    )
    add_region!(
        config,
        "fluid";
        type = SSTVerificationFluid(),
        initial_conditions = ComponentVector(
            density = density * u"kg/m^3",
            momentum_density_u = 0.0u"kg/(m^2*s)",
            momentum_density_v = 0.0u"kg/(m^2*s)",
            momentum_density_w = 0.0u"kg/(m^2*s)",
            volumetric_energy = pressure / (gamma - 1.0) * u"J/m^3",
            turbulent_kinetic_energy_density =
                density * turbulent_kinetic_energy * u"J/m^3",
            specific_dissipation_rate_density =
                density * specific_dissipation_rate * u"kg/(m^3*s)",
        ),
        properties = ComponentVector(
            cp = cp * u"J/(kg*K)",
            cv = cv * u"J/(kg*K)",
            mu = 1e-5u"Pa*s",
            molecular_viscosity = 1e-5u"Pa*s",
            k = 0.02u"W/(m*K)",
        ),
        property_update_function = function (du, u, p, t, system, geo, cell_id)
            overall_navier_stokes_property_update!(du, u, p, t, system, geo, cell_id)
            update_sst_primitives!(du, u, p, t, system, geo, cell_id)
        end,
        region_function = function (du, u, p, t, system, geo, cell_id)
            cap_navier_stokes_flow!(du, u, p, t, system, geo, cell_id)
            cap_sst_transport!(du, u, geo, cell_id)
        end,
    )
    for patch_name in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
        add_patch!(
            config,
            patch_name;
            properties = ComponentVector(),
            patch_function = _sst_verification_boundary_flux!,
        )
    end
    wall_distances = fill(wall_distance, n_cells)
    additional_data = (wall_distances = wall_distances,)
    du0, u0, system, geo = finish_fvm_config(
        config,
        _sst_verification_connection_map,
        additional_data;
        check_units = false,
    )
    stencil = build_weighted_least_squares_stencil(geo)
    function solve_groups!(du, u, p, t, system, geo)
        update_region_groups!(du, u, p, t, system, geo)
        update_weighted_least_squares_gradients!(u, stencil)
        update_sst_closure!(u, wall_distances)
        solve_connection_groups!(du, u, p, t, system, geo)
        solve_patch_groups!(du, u, p, t, system, geo)
        solve_region_groups!(du, u, p, t, system, geo)
        return nothing
    end
    rhs! = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_groups!)
    return SSTVerificationCase(system, geo, rhs!, u0, 0.0, gamma, cv)
end

function check_sst_conservative_operator()
    case = build_sst_verification_case(5)
    derivative = similar(case.u0)
    case.rhs!(derivative, case.u0, case.p, 0.0)
    state = ComponentVector(case.u0, case.system.state_axes)
    derivative_state = ComponentVector(derivative, case.system.state_axes)
    mean_flow_residual = maximum(abs.([
        derivative_state.density;
        derivative_state.momentum_density_u;
        derivative_state.momentum_density_v;
        derivative_state.momentum_density_w;
        derivative_state.volumetric_energy;
    ]))
    expected_k_rate = -DEFAULT_SST_CONSTANTS.beta_star .* state.density .*
        (state.turbulent_kinetic_energy_density ./ state.density) .*
        (state.specific_dissipation_rate_density ./ state.density)
    expected_omega_rate = -DEFAULT_SST_CONSTANTS.beta_1 .* state.density .*
        (state.specific_dissipation_rate_density ./ state.density).^2
    k_error = norm(
        derivative_state.turbulent_kinetic_energy_density .- expected_k_rate,
        Inf,
    )
    omega_error = norm(
        derivative_state.specific_dissipation_rate_density .- expected_omega_rate,
        Inf,
    )
    tolerance = 2e-12
    return (
        passed =
            mean_flow_residual <= tolerance &&
            k_error <= tolerance &&
            omega_error <= tolerance,
        summary = "the full FVM operator evolves conservative rho*k and rho*omega with the SST sources",
        metrics = Dict{String, Any}(
            "mean_flow_residual_Linf" => mean_flow_residual,
            "rho_k_source_error_Linf" => k_error,
            "rho_omega_source_error_Linf" => omega_error,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "state" => "conservative rho*k and rho*omega",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_sst_homogeneous_decay_integration()
    case = build_sst_verification_case(3)
    initial_state = ComponentVector(case.u0, case.system.state_axes)
    initial_k = initial_state.turbulent_kinetic_energy_density[1] /
        initial_state.density[1]
    initial_omega = initial_state.specific_dissipation_rate_density[1] /
        initial_state.density[1]
    final_time = 0.2
    problem = ODEProblem(case.rhs!, case.u0, (0.0, final_time), case.p)
    solution = solve(
        problem,
        Tsit5();
        abstol = 1e-11,
        reltol = 1e-10,
        isoutofdomain = (state, parameters, time) -> begin
            named = ComponentVector(state, case.system.state_axes)
            return minimum(named.turbulent_kinetic_energy_density) <= 0.0 ||
                minimum(named.specific_dissipation_rate_density) <= 0.0
        end,
    )
    final_state = ComponentVector(solution.u[end], case.system.state_axes)
    beta = DEFAULT_SST_CONSTANTS.beta_1
    beta_star = DEFAULT_SST_CONSTANTS.beta_star
    denominator = 1.0 + beta * initial_omega * final_time
    expected_omega = initial_omega / denominator
    expected_k = initial_k * denominator^(-beta_star / beta)
    calculated_k = final_state.turbulent_kinetic_energy_density[1] /
        final_state.density[1]
    calculated_omega = final_state.specific_dissipation_rate_density[1] /
        final_state.density[1]
    relative_k_error = abs(calculated_k - expected_k) / expected_k
    relative_omega_error = abs(calculated_omega - expected_omega) / expected_omega
    tolerance = 2e-8
    passed =
        OrdinaryDiffEq.SciMLBase.successful_retcode(solution) &&
        relative_k_error <= tolerance &&
        relative_omega_error <= tolerance &&
        minimum(final_state.turbulent_kinetic_energy_density) > 0.0 &&
        minimum(final_state.specific_dissipation_rate_density) > 0.0
    return (
        passed = passed,
        summary = "conservative SST turbulence follows the analytical homogeneous-decay solution",
        metrics = Dict{String, Any}(
            "relative_k_error" => relative_k_error,
            "relative_omega_error" => relative_omega_error,
            "final_k" => calculated_k,
            "final_omega" => calculated_omega,
            "accepted_steps" => solution.destats.naccept,
            "rejected_steps" => solution.destats.nreject,
            "rhs_evaluations" => solution.destats.nf,
        ),
        expected = Dict{String, Any}(
            "relative_tolerance" => tolerance,
            "positive_k_and_omega" => true,
        ),
        diagnostics = Dict{String, Any}("return_code" => string(solution.retcode)),
    )
end
