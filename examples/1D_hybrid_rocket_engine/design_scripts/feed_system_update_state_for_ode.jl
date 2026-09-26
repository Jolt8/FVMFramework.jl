include(joinpath(@__DIR__, "CEA_lookup_table.jl"))
#includes isp_interpolator_Pa(pressure_Pa, oxidizer_to_fuel_ratio) and cstar_interpolator_Pa(pressure_Pa, oxidizer_to_fuel_ratio)

function update_state!(du, u, p, t, oxidizer_model, chamber_model)
    du .= 0.0 
    
    tank_n_moles = u.tank_oxidizer_mass / p.nitrous_oxide_molecular_weight

    show_debug = false

    if show_debug
        @show "tank"
        @show u.tank_oxidizer_mass
        @show tank_n_moles
        @show u.tank_oxidizer_internal_energy
        @show p.tank_volume
    end

    result = uv_flash_via_vt(
        oxidizer_model,
        u.tank_oxidizer_internal_energy,
        p.tank_volume,
        [tank_n_moles],
    )

    p.tank_temperature = result.data.T
    p.tank_pressure = pressure(oxidizer_model, result)
    vapor_phase = argmax(result.volumes)
    p.tank_vapor_fraction = result.fractions[vapor_phase] / sum(result.fractions)

    #this increases stiffness significantly especially near the end
    if p.tank_vapor_fraction <= 0.999
        p.tank_density = mass_density(oxidizer_model, result, 1) #get liquid density because we drawing from the bottom of the tank
        p.tank_specific_enthalpy = mass_enthalpy(oxidizer_model, result, 1)
    else
        p.tank_density = mass_density(oxidizer_model, result) #otherwise, we will be drawing from the remaining vapor in the tank
        p.tank_specific_enthalpy = mass_enthalpy(oxidizer_model, result)
    end

    #p.tank_density = mass_density(oxidizer_model, result)
    #p.tank_specific_enthalpy = mass_enthalpy(oxidizer_model, result)

    mid_section_n_moles = u.mid_section_mass / p.nitrous_oxide_molecular_weight

    if show_debug
        @show "midsection"
        @show u.mid_section_mass
        @show mid_section_n_moles
        @show u.mid_section_internal_energy
        @show p.mid_section_volume
    end
    result = uv_flash_via_vt(
        oxidizer_model,
        u.mid_section_internal_energy,
        p.mid_section_volume,
        [mid_section_n_moles],
    )

    p.mid_section_temperature = result.data.T
    p.mid_section_pressure = pressure(oxidizer_model, result)
    vapor_phase = argmax(result.volumes)
    p.mid_section_vapor_fraction = result.fractions[vapor_phase] / sum(result.fractions)
    p.mid_section_density = mass_density(oxidizer_model, result)
    p.mid_section_specific_enthalpy = mass_enthalpy(oxidizer_model, result)

    if show_debug
        @show "chamber"
        @show u.chamber_gas_mass
        @show u.chamber_nitrous_oxide_mass_fraction
        @show u.chamber_hdpe_mass_fraction
        @show u.chamber_gas_internal_energy
        @show p.chamber_volume
    end

    oxidizer_moles = (u.chamber_gas_mass * u.chamber_nitrous_oxide_mass_fraction) / p.nitrous_oxide_molecular_weight
    fuel_moles = (u.chamber_gas_mass * u.chamber_hdpe_mass_fraction) / p.ethylene_molecular_weight
    chamber_moles = [oxidizer_moles, fuel_moles]

    # The combustion chamber is assumed to be a homogeneous gas, so recover its
    # temperature directly from U(V, T, n) instead of performing a phase flash.
    
    chamber_energy_residual = (temperature, _) -> Clapeyron.VT0.internal_energy(chamber_model, p.chamber_volume, temperature, chamber_moles) - u.chamber_gas_internal_energy

    
    temperature_problem = NonlinearProblem(chamber_energy_residual, p.chamber_temperature)
    temperature_solution = solve(temperature_problem, NewtonRaphson(); abstol = 1e-8, reltol = 1e-8)
    p.chamber_temperature = temperature_solution.u
    

    #=
    energy_residual(T) =
        Clapeyron.VT0.internal_energy(
            chamber_model,
            p.chamber_volume,
            T,
            chamber_moles,
        ) - u.chamber_gas_internal_energy

    p.chamber_temperature = find_zero(
        energy_residual,
        (200.0, 5000.0),
        Roots.Brent(),
    )
        =#
    

    #=
    energy_residual(T) =
        Clapeyron.VT0.internal_energy(
            chamber_model,
            p.chamber_volume,
            T,
            chamber_moles
        ) - u.chamber_gas_internal_energy

    p.chamber_temperature = find_zero(
        energy_residual,
        (300.0, 5000.0),
        Roots.Brent()
    )
        =#

    p.chamber_pressure = pressure(chamber_model, p.chamber_volume, p.chamber_temperature, chamber_moles)
    p.chamber_vapor_fraction = one(p.chamber_temperature)
    p.chamber_density = Clapeyron.VT0.mass_density(chamber_model, p.chamber_volume, p.chamber_temperature, chamber_moles)
    p.chamber_specific_enthalpy = Clapeyron.VT0.mass_enthalpy(chamber_model, p.chamber_volume, p.chamber_temperature, chamber_moles)

    #Other state updates:
    #Overall rocket
    p.desired_impulse = p.desired_average_thrust * p.desired_burn_time

    #Propellant properties
    #we don't need ISP yet, we only need that for the objective function, we do need propellant_characteristic_velocity for the mass flow out of the nozzle however
    #p.propellant_characteristic_velocity = cstar_interpolator_Pa(p.chamber_pressure, p.oxidizer_to_fuel_ratio)
    
    #Adjustable Valve
    p.valve_opening = valve_opening_at_t(t)
    p.valve_flow_capacity_factor = valve_flow_capacity_factor(p.valve_opening, p)

    #Injector
    #p.injector_orifice_area = p.injector_number_of_orifices * (pi / 4) * (p.injector_orifice_diameter^2)

    #Fuel Grain
    p.final_fuel_grain_void_diameter = p.u0_fuel_grain_void_diameter + p.additional_fuel_grain_void_diameter
    
    p.fuel_grain_average_cross_sectional_area = pi * (u.port_diameter / 2)^2
    
    p.fuel_grain_burning_surface_area = pi * u.port_diameter * p.fuel_grain_length
    
    p.fuel_mass = p.fuel_density * (pi * (p.final_fuel_grain_void_diameter / 2)^2 - pi * (u.port_diameter / 2)^2) * p.fuel_grain_length

    #Nozzle
    p.nozzle_throat_area = (pi / 4) * (p.nozzle_throat_diameter^2)
    #again, we need something here that determines propellant_isp based on oxidizer_to_fuel_ratio, chamber_pressure, and exit pressure (which we don't really know yet)
end