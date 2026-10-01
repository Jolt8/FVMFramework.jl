function _state_snapshot(state, cell_id)
    snapshot = Dict{String, Any}()
    for variable in (
        :density,
        :momentum_density_u,
        :momentum_density_v,
        :momentum_density_w,
        :volumetric_energy,
        :turbulent_kinetic_energy_density,
        :specific_dissipation_rate_density,
    )
        if hasproperty(state, variable)
            values = getproperty(state, variable)
            snapshot[string(variable)] = values[cell_id]
        end
    end
    if hasproperty(state, :species_densities)
        species_snapshot = Dict{String, Any}()
        for species_name in propertynames(state.species_densities)
            species_snapshot[string(species_name)] =
                getproperty(state.species_densities, species_name)[cell_id]
        end
        snapshot["species_densities"] = species_snapshot
    end
    return snapshot
end

function state_violations(
    state_vector,
    state_axes;
    gamma,
    cv,
    time = 0.0,
    maximum_violations = 20,
    species_sum_absolute_tolerance = 1e-12,
    species_sum_relative_tolerance = 1e-12,
)
    state = ComponentVector(state_vector, state_axes)
    violations = Dict{String, Any}[]
    n_cells = length(state.density)

    function record_violation(variable, cell_id, value, reason)
        if length(violations) >= maximum_violations
            return
        end
        push!(violations, Dict{String, Any}(
            "variable" => string(variable),
            "cell_index" => cell_id,
            "time" => time,
            "value" => value,
            "reason" => reason,
            "nearby_state" => _state_snapshot(state, cell_id),
        ))
    end

    for cell_id in 1:n_cells
        conservative_values = (
            state.density[cell_id],
            state.momentum_density_u[cell_id],
            state.momentum_density_v[cell_id],
            state.momentum_density_w[cell_id],
            state.volumetric_energy[cell_id],
        )
        conservative_names = (
            :density,
            :momentum_density_u,
            :momentum_density_v,
            :momentum_density_w,
            :volumetric_energy,
        )

        for (variable, value) in zip(conservative_names, conservative_values)
            if !isfinite(value)
                record_violation(variable, cell_id, value, "non-finite conservative state")
            end
        end

        density = state.density[cell_id]
        if !isfinite(density) || density <= 0.0
            record_violation(:density, cell_id, density, "density must be finite and positive")
            continue
        end

        kinetic_energy_density = 0.5 * (
            state.momentum_density_u[cell_id]^2 +
            state.momentum_density_v[cell_id]^2 +
            state.momentum_density_w[cell_id]^2
        ) / density
        internal_energy_density = state.volumetric_energy[cell_id] - kinetic_energy_density
        pressure = (gamma - 1.0) * internal_energy_density
        temperature = internal_energy_density / (density * cv)

        if !isfinite(pressure) || pressure <= 0.0
            record_violation(:pressure, cell_id, pressure, "pressure derived from the conservative state must be finite and positive")
        end
        if !isfinite(temperature) || temperature <= 0.0
            record_violation(:temperature, cell_id, temperature, "temperature derived from the conservative state must be finite and positive")
        end

        for turbulence_variable in (
            :turbulent_kinetic_energy,
            :density_turbulent_kinetic_energy,
            :turbulent_kinetic_energy_density,
            :specific_dissipation_rate,
            :density_specific_dissipation_rate,
            :specific_dissipation_rate_density,
        )
            if hasproperty(state, turbulence_variable)
                turbulence_value = getproperty(state, turbulence_variable)[cell_id]
                if !isfinite(turbulence_value) || turbulence_value <= 0.0
                    record_violation(
                        turbulence_variable,
                        cell_id,
                        turbulence_value,
                        "SST k-omega state variables must be finite and strictly positive",
                    )
                end
            end
        end

        if hasproperty(state, :species_densities)
            species_density_sum = 0.0
            for species_name in propertynames(state.species_densities)
                species_density = getproperty(
                    state.species_densities,
                    species_name,
                )[cell_id]
                species_density_sum += species_density
                if !isfinite(species_density) || species_density < 0.0
                    record_violation(
                        species_name,
                        cell_id,
                        species_density,
                        "conservative species density must be finite and non-negative",
                    )
                end
            end
            if !isapprox(
                species_density_sum,
                density;
                atol = species_sum_absolute_tolerance,
                rtol = species_sum_relative_tolerance,
            )
                record_violation(
                    :species_density_sum,
                    cell_id,
                    species_density_sum,
                    "conservative species densities must sum to mixture density",
                )
            end
        end
    end

    if hasproperty(state, :mass_fractions)
        mass_fractions = state.mass_fractions
        for cell_id in 1:n_cells
            fraction_sum = 0.0
            for species_name in propertynames(mass_fractions)
                fraction = getproperty(mass_fractions, species_name)[cell_id]
                fraction_sum += fraction
                if !isfinite(fraction) || fraction < 0.0 || fraction > 1.0
                    record_violation(species_name, cell_id, fraction, "mass fraction must be finite and lie in [0, 1]")
                end
            end
            if !isapprox(fraction_sum, 1.0; atol = 1e-12, rtol = 1e-12)
                record_violation(:mass_fraction_sum, cell_id, fraction_sum, "mass fractions must sum to one")
            end
        end
    end

    return violations
end

function state_is_valid(
    state_vector,
    state_axes;
    gamma,
    cv,
    time = 0.0,
    species_sum_absolute_tolerance = 1e-12,
    species_sum_relative_tolerance = 1e-12,
)
    return isempty(state_violations(
        state_vector,
        state_axes;
        gamma = gamma,
        cv = cv,
        time = time,
        maximum_violations = 1,
        species_sum_absolute_tolerance = species_sum_absolute_tolerance,
        species_sum_relative_tolerance = species_sum_relative_tolerance,
    ))
end
