function _viscous_compressible_mms_residuals(n_cells; viscous_flux_multiplier = 1.0)
    cell_width = 1.0 / n_cells
    centers = [(cell_id - 0.5) * cell_width for cell_id in 1:n_cells]
    density = 1.0
    dynamic_viscosity = 0.03
    mean_velocity = 0.2
    velocity_amplitude = 0.1
    wave_number = 2.0 * pi
    velocity(y) = mean_velocity + velocity_amplitude * sin(wave_number * y)
    velocity_gradient(y) = velocity_amplitude * wave_number * cos(wave_number * y)
    velocity_laplacian(y) = -velocity_amplitude * wave_number^2 * sin(wave_number * y)
    cell_velocity = velocity.(centers)

    density_residual = zeros(n_cells)
    momentum_residual = [
        -dynamic_viscosity * velocity_laplacian(y) for y in centers
    ]
    energy_residual = [
        -dynamic_viscosity * (
            velocity_gradient(y)^2 + velocity(y) * velocity_laplacian(y)
        ) for y in centers
    ]
    normal = [0.0, 1.0, 0.0]
    zero_gradient = (0.0, 0.0, 0.0)

    left_momentum_flux = viscous_flux_multiplier * dynamic_viscosity * velocity_gradient(0.0)
    left_energy_flux = left_momentum_flux * velocity(0.0)
    momentum_residual[1] -= left_momentum_flux / cell_width
    energy_residual[1] -= left_energy_flux / cell_width

    for left_cell in 1:(n_cells - 1)
        right_cell = left_cell + 1
        face_gradient = corrected_face_gradient(
            cell_velocity[left_cell],
            cell_velocity[right_cell],
            zero_gradient...,
            zero_gradient...,
            normal,
            cell_width,
        )[2]
        face_velocity = 0.5 * (
            cell_velocity[left_cell] + cell_velocity[right_cell]
        )
        momentum_flux = viscous_flux_multiplier * dynamic_viscosity * face_gradient
        energy_flux = face_velocity * momentum_flux
        momentum_residual[left_cell] += momentum_flux / cell_width
        momentum_residual[right_cell] -= momentum_flux / cell_width
        energy_residual[left_cell] += energy_flux / cell_width
        energy_residual[right_cell] -= energy_flux / cell_width
    end

    right_momentum_flux = viscous_flux_multiplier * dynamic_viscosity * velocity_gradient(1.0)
    right_energy_flux = right_momentum_flux * velocity(1.0)
    momentum_residual[end] += right_momentum_flux / cell_width
    energy_residual[end] += right_energy_flux / cell_width
    return (
        density = density_residual,
        momentum = momentum_residual ./ density,
        energy = energy_residual,
    )
end

function _mms_l2(values)
    return sqrt(sum(abs2, values) / length(values))
end

function _mms_orders(errors, grid_sizes)
    return [
        log(errors[index - 1] / errors[index]) /
        log(grid_sizes[index] / grid_sizes[index - 1])
        for index in 2:length(errors)
    ]
end

function check_viscous_compressible_mms(grid_sizes)
    density_errors = Float64[]
    momentum_errors = Float64[]
    energy_errors = Float64[]
    for n_cells in grid_sizes
        residuals = _viscous_compressible_mms_residuals(n_cells)
        push!(density_errors, _mms_l2(residuals.density))
        push!(momentum_errors, _mms_l2(residuals.momentum[2:(end - 1)]))
        push!(energy_errors, _mms_l2(residuals.energy[2:(end - 1)]))
    end
    momentum_orders = _mms_orders(momentum_errors, grid_sizes)
    energy_orders = _mms_orders(energy_errors, grid_sizes)
    asymptotic_order = min(
        minimum(last(momentum_orders, min(2, length(momentum_orders)))),
        minimum(last(energy_orders, min(2, length(energy_orders)))),
    )
    minimum_order = 1.8
    passed =
        maximum(density_errors) == 0.0 &&
        asymptotic_order >= minimum_order
    return (
        passed = passed,
        summary = "viscous compressible shear MMS converges in momentum and total energy",
        metrics = Dict{String, Any}(
            "grid_sizes" => grid_sizes,
            "density_L2_errors" => density_errors,
            "interior_momentum_L2_errors" => momentum_errors,
            "interior_energy_L2_errors" => energy_errors,
            "momentum_orders" => momentum_orders,
            "energy_orders" => energy_orders,
            "asymptotic_minimum_order" => asymptotic_order,
        ),
        expected = Dict{String, Any}(
            "minimum_asymptotic_order" => minimum_order,
            "density_residual" => 0.0,
            "manufactured_field" => "rho,p,T constant; u(y)=u0+A*sin(2*pi*y)",
            "norm_domain" => "interior cells; analytical boundary fluxes are tested separately",
        ),
        diagnostics = Dict{String, Any}(),
    )
end

function _species_mms_residuals(n_cells; diffusion_flux_multiplier = 1.0)
    cell_width = 1.0 / n_cells
    centers = [(cell_id - 0.5) * cell_width for cell_id in 1:n_cells]
    density = 1.3
    diffusion_coefficient = 0.04
    mean_fraction = 0.4
    fraction_amplitude = 0.1
    wave_number = 2.0 * pi
    species_a_fraction(y) = mean_fraction + fraction_amplitude * sin(wave_number * y)
    species_a_gradient(y) = fraction_amplitude * wave_number * cos(wave_number * y)
    species_a_laplacian(y) = -fraction_amplitude * wave_number^2 * sin(wave_number * y)
    fractions_a = species_a_fraction.(centers)
    fractions_b = 1.0 .- fractions_a
    state = ComponentVector(
        density = fill(density, n_cells),
        mass_fractions = ComponentVector(
            species_a = fractions_a,
            species_b = fractions_b,
        ),
        diffusion_coefficients = ComponentVector(
            species_a = fill(diffusion_coefficient, n_cells),
            species_b = fill(diffusion_coefficient, n_cells),
        ),
    )
    derivative = ComponentVector(
        species_density_flow = ComponentVector(
            species_a = zeros(n_cells),
            species_b = zeros(n_cells),
        ),
    )
    for left_cell in 1:(n_cells - 1)
        right_cell = left_cell + 1
        add_conservative_species_diffusion_flux!(
            derivative,
            state,
            left_cell,
            right_cell,
            1.0,
            cell_width,
        )
    end
    derivative.species_density_flow.species_a .*= diffusion_flux_multiplier
    derivative.species_density_flow.species_b .*= diffusion_flux_multiplier

    left_flux_a = density * diffusion_coefficient * species_a_gradient(0.0)
    right_flux_a = density * diffusion_coefficient * species_a_gradient(1.0)
    derivative.species_density_flow.species_a[1] -=
        diffusion_flux_multiplier * left_flux_a
    derivative.species_density_flow.species_a[end] +=
        diffusion_flux_multiplier * right_flux_a
    derivative.species_density_flow.species_b[1] +=
        diffusion_flux_multiplier * left_flux_a
    derivative.species_density_flow.species_b[end] -=
        diffusion_flux_multiplier * right_flux_a

    source_a = -density * diffusion_coefficient .* species_a_laplacian.(centers)
    residual_a = derivative.species_density_flow.species_a ./ cell_width .+ source_a
    residual_b = derivative.species_density_flow.species_b ./ cell_width .- source_a
    return (species_a = residual_a, species_b = residual_b)
end

function check_species_mms(grid_sizes)
    species_a_errors = Float64[]
    species_b_errors = Float64[]
    sum_errors = Float64[]
    for n_cells in grid_sizes
        residuals = _species_mms_residuals(n_cells)
        push!(species_a_errors, _mms_l2(residuals.species_a[2:(end - 1)]))
        push!(species_b_errors, _mms_l2(residuals.species_b[2:(end - 1)]))
        push!(sum_errors, norm(residuals.species_a .+ residuals.species_b, Inf))
    end
    species_a_orders = _mms_orders(species_a_errors, grid_sizes)
    species_b_orders = _mms_orders(species_b_errors, grid_sizes)
    asymptotic_order = min(
        minimum(last(species_a_orders, min(2, length(species_a_orders)))),
        minimum(last(species_b_orders, min(2, length(species_b_orders)))),
    )
    minimum_order = 1.8
    sum_tolerance = 1e-12
    return (
        passed = asymptotic_order >= minimum_order && maximum(sum_errors) <= sum_tolerance,
        summary = "conservative two-species diffusion MMS converges while preserving the mixture sum",
        metrics = Dict{String, Any}(
            "grid_sizes" => grid_sizes,
            "interior_species_a_L2_errors" => species_a_errors,
            "interior_species_b_L2_errors" => species_b_errors,
            "species_a_orders" => species_a_orders,
            "species_b_orders" => species_b_orders,
            "species_sum_Linf_errors" => sum_errors,
            "asymptotic_minimum_order" => asymptotic_order,
        ),
        expected = Dict{String, Any}(
            "minimum_asymptotic_order" => minimum_order,
            "species_sum_tolerance" => sum_tolerance,
            "manufactured_field" => "Y_a=0.4+0.1*sin(2*pi*y), Y_b=1-Y_a",
            "norm_domain" => "interior cells; analytical boundary fluxes are tested separately",
        ),
        diagnostics = Dict{String, Any}(),
    )
end
