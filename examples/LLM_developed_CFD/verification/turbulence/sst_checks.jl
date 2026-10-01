function check_sst_blending_limits()
    density = 1.0
    molecular_viscosity = 1e-5
    zero_gradient = [0.0, 0.0, 0.0]
    near_wall = sst_blending_functions(
        density,
        0.1,
        100.0,
        molecular_viscosity,
        1e-5,
        zero_gradient,
        zero_gradient,
    )
    free_stream = sst_blending_functions(
        density,
        1e-6,
        100.0,
        molecular_viscosity,
        10.0,
        zero_gradient,
        zero_gradient,
    )
    tolerance = 1e-8
    return (
        passed =
            near_wall.F1 >= 1.0 - tolerance &&
            near_wall.F2 >= 1.0 - tolerance &&
            free_stream.F1 <= tolerance &&
            free_stream.F2 <= tolerance,
        summary = "SST blending selects the k-omega wall limit and outer-flow limit",
        metrics = Dict{String, Any}(
            "near_wall_F1" => near_wall.F1,
            "near_wall_F2" => near_wall.F2,
            "free_stream_F1" => free_stream.F1,
            "free_stream_F2" => free_stream.F2,
        ),
        expected = Dict{String, Any}(
            "near_wall_limit" => 1.0,
            "free_stream_limit" => 0.0,
            "absolute_tolerance" => tolerance,
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_sst_homogeneous_decay_sources()
    sources = sst_source_terms(
        1.2,
        0.4,
        3.0,
        0.0,
        0.0,
        1.0,
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    )
    expected_k_source = -DEFAULT_SST_CONSTANTS.beta_star * 1.2 * 0.4 * 3.0
    expected_omega_source = -DEFAULT_SST_CONSTANTS.beta_1 * 1.2 * 3.0^2
    maximum_error = max(
        abs(sources.k_source - expected_k_source),
        abs(sources.omega_source - expected_omega_source),
    )
    tolerance = 1e-14
    return (
        passed =
            maximum_error <= tolerance &&
            sources.k_source < 0.0 &&
            sources.omega_source < 0.0,
        summary = "homogeneous zero-strain SST turbulence decays at the analytical rates",
        metrics = Dict{String, Any}(
            "k_source" => sources.k_source,
            "omega_source" => sources.omega_source,
            "maximum_error" => maximum_error,
        ),
        expected = Dict{String, Any}(
            "k_source" => expected_k_source,
            "omega_source" => expected_omega_source,
            "absolute_tolerance" => tolerance,
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_sst_channel_log_layer()
    constants = DEFAULT_SST_CONSTANTS
    density = 1.0
    molecular_viscosity = 1e-5
    friction_velocity = 0.05
    von_karman = 0.41
    y_plus_values = [50.0, 100.0, 200.0]
    viscosity_ratios = Float64[]
    production_destruction_ratios = Float64[]
    F2_values = Float64[]
    for y_plus in y_plus_values
        wall_distance = y_plus * molecular_viscosity /
            (density * friction_velocity)
        turbulent_kinetic_energy =
            friction_velocity^2 / sqrt(constants.beta_star)
        specific_dissipation_rate = friction_velocity /
            (sqrt(constants.beta_star) * von_karman * wall_distance)
        shear_rate = friction_velocity / (von_karman * wall_distance)
        blending = sst_blending_functions(
            density,
            turbulent_kinetic_energy,
            specific_dissipation_rate,
            molecular_viscosity,
            wall_distance,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        )
        turbulent_viscosity = sst_eddy_viscosity(
            density,
            turbulent_kinetic_energy,
            specific_dissipation_rate,
            shear_rate^2,
            blending.F2,
        )
        sources = sst_source_terms(
            density,
            turbulent_kinetic_energy,
            specific_dissipation_rate,
            turbulent_viscosity,
            shear_rate^2,
            blending.F1,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        )
        expected_viscosity = density * von_karman * wall_distance * friction_velocity
        push!(viscosity_ratios, turbulent_viscosity / expected_viscosity)
        push!(production_destruction_ratios, sources.production / sources.k_destruction)
        push!(F2_values, blending.F2)
    end
    viscosity_error = maximum(abs.(viscosity_ratios .- 1.0))
    equilibrium_error = maximum(abs.(production_destruction_ratios .- 1.0))
    tolerance = 2e-12
    return (
        passed =
            viscosity_error <= tolerance &&
            equilibrium_error <= tolerance &&
            minimum(F2_values) >= 1.0 - tolerance,
        summary = "SST recovers log-layer channel eddy viscosity and k-equilibrium",
        metrics = Dict{String, Any}(
            "y_plus" => y_plus_values,
            "eddy_viscosity_ratios" => viscosity_ratios,
            "production_destruction_ratios" => production_destruction_ratios,
            "F2" => F2_values,
            "maximum_eddy_viscosity_error" => viscosity_error,
            "maximum_k_equilibrium_error" => equilibrium_error,
        ),
        expected = Dict{String, Any}(
            "eddy_viscosity_ratio" => 1.0,
            "production_destruction_ratio" => 1.0,
            "absolute_tolerance" => tolerance,
            "regime" => "equilibrium turbulent channel log layer",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function _read_flat_plate_reference()
    reference_path = joinpath(
        @__DIR__,
        "..",
        "reference_data",
        "turbulent_flat_plate.csv",
    )
    rows = NamedTuple[]
    for (line_number, line) in enumerate(eachline(reference_path))
        if line_number == 1
            continue
        end
        fields = split(line, ',')
        push!(rows, (
            reynolds_x = parse(Float64, fields[1]),
            skin_friction_coefficient = parse(Float64, fields[2]),
            source = fields[3],
        ))
    end
    return rows, reference_path
end

function check_sst_flat_plate_equilibrium()
    rows, reference_path = _read_flat_plate_reference()
    density = 1.0
    free_stream_velocity = 1.0
    streamwise_location = 1.0
    von_karman = 0.41
    y_plus = 100.0
    relative_errors = Float64[]
    predicted_coefficients = Float64[]
    for row in rows
        molecular_viscosity = density * free_stream_velocity *
            streamwise_location / row.reynolds_x
        friction_velocity = free_stream_velocity *
            sqrt(0.5 * row.skin_friction_coefficient)
        wall_distance = y_plus * molecular_viscosity /
            (density * friction_velocity)
        turbulent_kinetic_energy = friction_velocity^2 /
            sqrt(DEFAULT_SST_CONSTANTS.beta_star)
        specific_dissipation_rate = friction_velocity /
            (sqrt(DEFAULT_SST_CONSTANTS.beta_star) * von_karman * wall_distance)
        shear_rate = friction_velocity / (von_karman * wall_distance)
        blending = sst_blending_functions(
            density,
            turbulent_kinetic_energy,
            specific_dissipation_rate,
            molecular_viscosity,
            wall_distance,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        )
        turbulent_viscosity = sst_eddy_viscosity(
            density,
            turbulent_kinetic_energy,
            specific_dissipation_rate,
            shear_rate^2,
            blending.F2,
        )
        modeled_wall_stress = turbulent_viscosity * shear_rate
        predicted_coefficient = 2.0 * modeled_wall_stress /
            (density * free_stream_velocity^2)
        push!(predicted_coefficients, predicted_coefficient)
        push!(relative_errors, abs(
            predicted_coefficient - row.skin_friction_coefficient
        ) / row.skin_friction_coefficient)
    end
    maximum_relative_error = maximum(relative_errors)
    tolerance = 2e-12
    return (
        passed = maximum_relative_error <= tolerance,
        summary = "SST log-layer stress matches the turbulent flat-plate skin-friction reference",
        metrics = Dict{String, Any}(
            "reynolds_x" => [row.reynolds_x for row in rows],
            "reference_skin_friction" => [row.skin_friction_coefficient for row in rows],
            "predicted_skin_friction" => predicted_coefficients,
            "relative_errors" => relative_errors,
            "maximum_relative_error" => maximum_relative_error,
        ),
        expected = Dict{String, Any}(
            "maximum_relative_error" => tolerance,
            "reference" => "one-fifth-power-law smooth turbulent flat plate",
            "reference_data_path" => reference_path,
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function check_sst_face_conservation()
    state = ComponentVector(
        density = [1.0, 1.1],
        turbulent_kinetic_energy = [0.2, 0.35],
        specific_dissipation_rate = [3.0, 2.0],
        molecular_viscosity = [1e-5, 1.2e-5],
        turbulent_viscosity = [0.01, 0.015],
        sst_F1 = [0.8, 0.6],
        grad_turbulent_kinetic_energy = [0.1 0.2 0.0; -0.2 0.3 0.0],
        grad_specific_dissipation_rate = [0.5 -0.4 0.0; 0.2 0.7 0.0],
    )
    derivative = ComponentVector(
        turbulent_kinetic_energy_density_flow = zeros(2),
        specific_dissipation_rate_density_flow = zeros(2),
    )
    add_hllc_sst_advection_flux!(derivative, state, 1, 2, 0.7, 0.4)
    add_sst_diffusion_flux!(
        derivative,
        state,
        1,
        2,
        0.7,
        [1.0, 0.0, 0.0],
        0.25,
    )
    k_imbalance = sum(derivative.turbulent_kinetic_energy_density_flow)
    omega_imbalance = sum(derivative.specific_dissipation_rate_density_flow)
    activity = norm(derivative, Inf)
    tolerance = 1e-14
    return (
        passed =
            abs(k_imbalance) <= tolerance &&
            abs(omega_imbalance) <= tolerance &&
            activity > tolerance,
        summary = "SST conservative advection and diffusion are equal and opposite across a face",
        metrics = Dict{String, Any}(
            "k_global_imbalance" => k_imbalance,
            "omega_global_imbalance" => omega_imbalance,
            "flux_activity" => activity,
        ),
        expected = Dict{String, Any}("absolute_tolerance" => tolerance),
        diagnostics = Dict{String, Any}(),
    )
end
