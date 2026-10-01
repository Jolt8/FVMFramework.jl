function check_species_zero_diffusion()
    uniform_flux = species_numerical_flux(
        1.0,
        1.4,
        0.35,
        0.35,
        1e-5,
        2e-5,
        0.8,
        0.2,
    )
    zero_coefficient_flux = species_numerical_flux(
        1.0,
        1.4,
        0.2,
        0.8,
        0.0,
        2e-5,
        0.8,
        0.2,
    )
    maximum_flux = max(abs(uniform_flux), abs(zero_coefficient_flux))
    return (
        passed = maximum_flux == 0.0,
        summary = "uniform composition and an impermeable side give zero species diffusion",
        metrics = Dict{String, Any}(
            "uniform_composition_flux" => uniform_flux,
            "zero_coefficient_flux" => zero_coefficient_flux,
        ),
        expected = Dict{String, Any}(
            "exact_flux" => 0.0,
            "oracle" => "Fick's law with zero composition gradient or zero effective diffusivity",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_species_linear_diffusion()
    density_a = 1.0
    density_b = 1.2
    mass_fraction_a = 0.2
    mass_fraction_b = 0.5
    diffusion_coefficient_a = 2e-5
    diffusion_coefficient_b = 4e-5
    area = 0.7
    distance = 0.1
    measured_flux = species_numerical_flux(
        density_a,
        density_b,
        mass_fraction_a,
        mass_fraction_b,
        diffusion_coefficient_a,
        diffusion_coefficient_b,
        area,
        distance,
    )
    density_average = 0.5 * (density_a + density_b)
    diffusion_coefficient_effective = 2.0 * diffusion_coefficient_a *
        diffusion_coefficient_b / (diffusion_coefficient_a + diffusion_coefficient_b)
    expected_flux = -density_average * diffusion_coefficient_effective *
        ((mass_fraction_b - mass_fraction_a) / distance) * area
    absolute_error = abs(measured_flux - expected_flux)
    tolerance = 50.0 * eps(abs(expected_flux))
    correct_direction = measured_flux < 0.0
    return (
        passed = absolute_error <= tolerance && correct_direction,
        summary = "linear composition gradient has the independent Fickian magnitude and direction",
        metrics = Dict{String, Any}(
            "measured_flux" => measured_flux,
            "expected_flux" => expected_flux,
            "absolute_error" => absolute_error,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "direction" => "from high toward low mass fraction",
            "oracle" => "independent arithmetic-density/harmonic-diffusivity Fick law",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function _two_species_state()
    return ComponentVector(
        rho = [1.0, 1.0],
        mass_fractions = ComponentVector(
            species_a = [0.2, 0.8],
            species_b = [0.8, 0.2],
        ),
        diffusion_coefficients = ComponentVector(
            species_a = [1e-5, 1e-5],
            species_b = [1e-5, 1e-5],
        ),
    )
end

function check_species_diffusive_conservation()
    state = _two_species_state()
    derivative = ComponentVector(
        mass_fractions = ComponentVector(
            species_a = zeros(2),
            species_b = zeros(2),
        ),
    )
    area = 1.0
    distance = 0.5
    normal = [1.0, 0.0, 0.0]
    mass_fraction_diffusion!(
        derivative,
        state,
        1,
        2,
        1,
        area,
        normal,
        distance,
    )
    mass_fraction_diffusion!(
        derivative,
        state,
        2,
        1,
        1,
        area,
        -normal,
        distance,
    )

    species_a_global_rate = sum(derivative.mass_fractions.species_a)
    species_b_global_rate = sum(derivative.mass_fractions.species_b)
    cell_one_sum_rate = derivative.mass_fractions.species_a[1] +
        derivative.mass_fractions.species_b[1]
    cell_two_sum_rate = derivative.mass_fractions.species_a[2] +
        derivative.mass_fractions.species_b[2]
    maximum_imbalance = maximum(abs.((
        species_a_global_rate,
        species_b_global_rate,
        cell_one_sum_rate,
        cell_two_sum_rate,
    )))
    tolerance = 100.0 * eps(Float64)
    return (
        passed = maximum_imbalance <= tolerance,
        summary = "oriented Fickian face updates conserve each species and the equal-D fraction sum",
        metrics = Dict{String, Any}(
            "species_a_global_rate" => species_a_global_rate,
            "species_b_global_rate" => species_b_global_rate,
            "cell_one_fraction_sum_rate" => cell_one_sum_rate,
            "cell_two_fraction_sum_rate" => cell_two_sum_rate,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "constraint" => "equal diffusion coefficients and sum(Y)=1 imply sum(J)=0",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_species_advection_and_sum_constraint()
    state = _two_species_state()
    derivative = ComponentVector(
        mass_face = zeros(2, 1),
        mass = [-2.0, 2.0],
        mass_fractions = ComponentVector(
            species_a = zeros(2),
            species_b = zeros(2),
        ),
        species_masses = ComponentVector(
            species_a = zeros(2),
            species_b = zeros(2),
        ),
    )
    derivative.mass_face[1, 1] = -2.0
    derivative.mass_face[2, 1] = 2.0
    normal = [1.0, 0.0, 0.0]
    all_species_advection!(
        derivative,
        state,
        1,
        2,
        1,
        1.0,
        normal,
        0.5,
        1.0,
        1.0,
    )
    all_species_advection!(
        derivative,
        state,
        2,
        1,
        1,
        1.0,
        -normal,
        0.5,
        1.0,
        1.0,
    )

    species_a_global_rate = sum(derivative.species_masses.species_a)
    species_b_global_rate = sum(derivative.species_masses.species_b)
    cell_one_species_rate = derivative.species_masses.species_a[1] +
        derivative.species_masses.species_b[1]
    cell_two_species_rate = derivative.species_masses.species_a[2] +
        derivative.species_masses.species_b[2]

    cap_species_mass_flux_to_mass_fraction_change!(derivative, state, 1, 1.0)
    cap_species_mass_flux_to_mass_fraction_change!(derivative, state, 2, 1.0)
    cell_one_fraction_sum_rate = derivative.mass_fractions.species_a[1] +
        derivative.mass_fractions.species_b[1]
    cell_two_fraction_sum_rate = derivative.mass_fractions.species_a[2] +
        derivative.mass_fractions.species_b[2]
    imbalances = (
        species_a_global_rate,
        species_b_global_rate,
        cell_one_species_rate - derivative.mass[1],
        cell_two_species_rate - derivative.mass[2],
        cell_one_fraction_sum_rate,
        cell_two_fraction_sum_rate,
    )
    maximum_imbalance = maximum(abs.(imbalances))
    tolerance = 100.0 * eps(Float64)
    return (
        passed = maximum_imbalance <= tolerance,
        summary = "upwind species advection conserves species and preserves sum(Y)=1 after capping",
        metrics = Dict{String, Any}(
            "species_a_global_rate" => species_a_global_rate,
            "species_b_global_rate" => species_b_global_rate,
            "cell_one_fraction_sum_rate" => cell_one_fraction_sum_rate,
            "cell_two_fraction_sum_rate" => cell_two_fraction_sum_rate,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "upwind_source_cell" => 1,
            "constraint" => "sum of species mass fluxes equals total mass flux",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_conservative_species_state(case)
    state = ComponentVector(case.u0, case.system.state_axes)
    derivative = residual(case)
    derivative_state = ComponentVector(derivative, case.system.state_axes)

    cache_derivative = zeros(length(case.u0))
    du, u = unpack_fvm_state(
        cache_derivative,
        case.u0,
        case.p,
        0.0,
        case.system,
    )
    update_region_groups!(du, u, case.p, 0.0, case.system, case.geo)

    maximum_state_sum_error = 0.0
    maximum_cache_error = 0.0
    maximum_residual_sum_error = 0.0
    for cell_id in eachindex(state.density)
        species_density_sum = 0.0
        species_derivative_sum = 0.0
        for species_name in propertynames(state.species_densities)
            species_density = getproperty(
                state.species_densities,
                species_name,
            )[cell_id]
            species_density_sum += species_density
            species_derivative_sum += getproperty(
                derivative_state.species_densities,
                species_name,
            )[cell_id]
            cached_mass_fraction = getproperty(
                u.mass_fractions,
                species_name,
            )[cell_id]
            expected_mass_fraction = species_density / state.density[cell_id]
            maximum_cache_error = max(
                maximum_cache_error,
                abs(cached_mass_fraction - expected_mass_fraction),
            )
        end
        maximum_state_sum_error = max(
            maximum_state_sum_error,
            abs(species_density_sum - state.density[cell_id]),
        )
        maximum_residual_sum_error = max(
            maximum_residual_sum_error,
            abs(species_derivative_sum - derivative_state.density[cell_id]),
        )
    end

    scale = max(norm(derivative_state.density, Inf), 1.0)
    normalized_residual_error = maximum_residual_sum_error / scale
    tolerance = 2e-12
    return (
        passed = maximum_state_sum_error <= tolerance &&
            maximum_cache_error <= tolerance &&
            normalized_residual_error <= tolerance,
        summary = "rho*Y is evolved in the ODE state and its summed residual equals the density residual",
        metrics = Dict{String, Any}(
            "state_sum_error" => maximum_state_sum_error,
            "mass_fraction_cache_error" => maximum_cache_error,
            "normalized_residual_sum_error" => normalized_residual_error,
            "species_degrees_of_freedom" => length(state.species_densities),
        ),
        expected = Dict{String, Any}(
            "absolute_or_normalized_tolerance" => tolerance,
            "state_identity" => "sum(rho*Y_k) = rho",
            "residual_identity" => "sum(d(rho*Y_k)/dt) = d(rho)/dt",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_conservative_species_face_advection()
    state = ComponentVector(
        mass_fractions = ComponentVector(
            species_a = [0.2, 0.7],
            species_b = [0.8, 0.3],
        ),
    )
    derivative = ComponentVector(
        species_density_flow = ComponentVector(
            species_a = zeros(2),
            species_b = zeros(2),
        ),
    )
    area = 0.7
    density_flux = 2.0
    add_hllc_species_advection_flux!(
        derivative,
        state,
        1,
        2,
        area,
        density_flux,
    )

    expected_integrated_density_flux = area * density_flux
    cell_a_species_rate =
        derivative.species_density_flow.species_a[1] +
        derivative.species_density_flow.species_b[1]
    cell_b_species_rate =
        derivative.species_density_flow.species_a[2] +
        derivative.species_density_flow.species_b[2]
    species_a_global_rate = sum(derivative.species_density_flow.species_a)
    species_b_global_rate = sum(derivative.species_density_flow.species_b)
    maximum_error = maximum(abs.((
        cell_a_species_rate + expected_integrated_density_flux,
        cell_b_species_rate - expected_integrated_density_flux,
        species_a_global_rate,
        species_b_global_rate,
    )))
    tolerance = 100.0 * eps(Float64)
    return (
        passed = maximum_error <= tolerance,
        summary = "HLLC mass flux advects conservative species densities with exact face conservation",
        metrics = Dict{String, Any}(
            "maximum_error" => maximum_error,
            "cell_a_species_rate" => cell_a_species_rate,
            "cell_b_species_rate" => cell_b_species_rate,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "species_flux_sum" => expected_integrated_density_flux,
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_mixture_corrected_species_diffusion()
    state = ComponentVector(
        density = [1.0, 1.3],
        mass_fractions = ComponentVector(
            species_a = [0.2, 0.7],
            species_b = [0.8, 0.3],
        ),
        diffusion_coefficients = ComponentVector(
            species_a = [1e-5, 2e-5],
            species_b = [4e-5, 3e-5],
        ),
    )
    derivative = ComponentVector(
        species_density_flow = ComponentVector(
            species_a = zeros(2),
            species_b = zeros(2),
        ),
    )
    add_conservative_species_diffusion_flux!(
        derivative,
        state,
        1,
        2,
        0.8,
        0.25,
    )

    cell_a_sum = derivative.species_density_flow.species_a[1] +
        derivative.species_density_flow.species_b[1]
    cell_b_sum = derivative.species_density_flow.species_a[2] +
        derivative.species_density_flow.species_b[2]
    species_a_global_rate = sum(derivative.species_density_flow.species_a)
    species_b_global_rate = sum(derivative.species_density_flow.species_b)
    activity = maximum(abs.(derivative.species_density_flow))
    maximum_imbalance = maximum(abs.((
        cell_a_sum,
        cell_b_sum,
        species_a_global_rate,
        species_b_global_rate,
    )))
    tolerance = 1e-14
    return (
        passed = maximum_imbalance <= tolerance && activity > tolerance,
        summary = "mixture-corrected unequal-D diffusion conserves every species and total density",
        metrics = Dict{String, Any}(
            "maximum_imbalance" => maximum_imbalance,
            "maximum_species_flux_activity" => activity,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "constraint" => "sum(J_k)=0 with unequal diffusion coefficients",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_coupled_species_advection_integration(case, final_time)
    initial_state = ComponentVector(case.u0, case.system.state_axes)
    initial_species_a_mass = sum(
        case.geo.cell_volumes .*
        initial_state.species_densities.species_a,
    )
    integration = run_explicit_smoke(case, final_time)
    final_state = ComponentVector(
        integration.solution.u[end],
        case.system.state_axes,
    )
    final_species_a_mass = sum(
        case.geo.cell_volumes .*
        final_state.species_densities.species_a,
    )

    density = initial_state.density[1]
    velocity = initial_state.momentum_density_u[1] / density
    inlet_mass_fraction = 0.8
    outlet_mass_fraction = 0.2
    cross_sectional_area = 1.0
    expected_mass_change = density * velocity * cross_sectional_area *
        (inlet_mass_fraction - outlet_mass_fraction) * final_time
    measured_mass_change = final_species_a_mass - initial_species_a_mass
    conservation_error = abs(measured_mass_change - expected_mass_change)

    maximum_sum_error = 0.0
    minimum_mass_fraction = Inf
    maximum_mass_fraction = -Inf
    for cell_id in eachindex(final_state.density)
        species_density_sum =
            final_state.species_densities.species_a[cell_id] +
            final_state.species_densities.species_b[cell_id]
        maximum_sum_error = max(
            maximum_sum_error,
            abs(species_density_sum - final_state.density[cell_id]),
        )
        species_a_mass_fraction =
            final_state.species_densities.species_a[cell_id] /
            final_state.density[cell_id]
        minimum_mass_fraction = min(
            minimum_mass_fraction,
            species_a_mass_fraction,
        )
        maximum_mass_fraction = max(
            maximum_mass_fraction,
            species_a_mass_fraction,
        )
    end

    tolerance = 2e-11
    bounded = minimum_mass_fraction >= -tolerance &&
        maximum_mass_fraction <= 1.0 + tolerance
    return (
        passed = integration.passed && conservation_error <= tolerance &&
            maximum_sum_error <= tolerance && bounded,
        summary = "the integrated compressible ODE advects rho*Y conservatively and preserves sum(Y)=1",
        metrics = Dict{String, Any}(
            "solver_statistics" => integration.statistics,
            "measured_species_a_mass_change" => measured_mass_change,
            "expected_species_a_mass_change" => expected_mass_change,
            "conservation_error" => conservation_error,
            "maximum_species_density_sum_error" => maximum_sum_error,
            "minimum_species_a_mass_fraction" => minimum_mass_fraction,
            "maximum_species_a_mass_fraction" => maximum_mass_fraction,
        ),
        expected = Dict{String, Any}(
            "absolute_tolerance" => tolerance,
            "oracle" => "boundary flux balance for a passively advected two-species contact",
        ),
        diagnostics = Dict{String, Any}(
            "state_violations" => integration.violations,
        ),
    )
end
