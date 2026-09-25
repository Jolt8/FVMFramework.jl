using Revise
using Unitful
using ComponentArrays
using Clapeyron
using NonlinearSolve
using ForwardDiff
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using OrderedCollections
using OrdinaryDiffEq
using Roots
using SparseConnectivityTracer
using Dates
using CSV
using DataFrames
using Sparspak

using FVMFramework

Revise.includet(joinpath(@__DIR__, "overloads/clapeyron_tracer_overloads.jl"))

Revise.includet(joinpath(@__DIR__, "design_helper_functions.jl"))

oxidizer_model = PR(["nitrous oxide"])
chamber_model = ReidIdeal(["nitrous oxide", "ethylene"])

mass_density(oxidizer_model, 100000, 298.15, [1.0])
mass_density(chamber_model, 100000, 298.15, [1.0, 1.0])

function adjustable_valve_flow!(du, u, p, t)
    if p.tank_pressure - p.mid_section_pressure <= 0
        return 0.0
    end
    
    volumetric_flow = p.valve_flow_capacity_factor * p.valve_opening * sqrt((p.tank_pressure - p.mid_section_pressure) / p.tank_density)
    #we might need a different equation that activates when the tank starts spitting out gas

    if true == false
        @show p.valve_flow_capacity_factor
        @show p.valve_opening
        @show p.tank_pressure
        @show p.mid_section_pressure
        @show p.tank_density
        @show volumetric_flow
    end

    mass_flow = volumetric_flow * p.tank_density

    du.tank_oxidizer_mass -= mass_flow
    du.tank_oxidizer_internal_energy -= mass_flow * p.tank_specific_enthalpy
    du.mid_section_mass += mass_flow
    du.mid_section_internal_energy += mass_flow * p.tank_specific_enthalpy

    return (p.tank_pressure - p.mid_section_pressure)
end

function injector_valve_flow!(du, u, p, t)
    if p.mid_section_pressure - p.chamber_pressure <= 0
        return 0.0, 0.0
    end

    volumetric_flow = p.injector_discharge_coefficient * p.injector_orifice_area * sqrt((2 * (p.mid_section_pressure - p.chamber_pressure)) / p.mid_section_density)

    if true == false
        println("")
        @show p.injector_discharge_coefficient
        @show p.injector_orifice_area
        @show p.mid_section_pressure
        @show p.chamber_pressure
        @show p.mid_section_density
        @show volumetric_flow
    end

    oxidizer_mass_flow = volumetric_flow * p.mid_section_density

    du.mid_section_mass -= oxidizer_mass_flow
    du.mid_section_internal_energy -= oxidizer_mass_flow * p.mid_section_specific_enthalpy
    du.chamber_gas_mass += oxidizer_mass_flow
    du.chamber_gas_internal_energy += oxidizer_mass_flow * p.mid_section_specific_enthalpy

    return (p.mid_section_pressure - p.chamber_pressure), oxidizer_mass_flow
end

function regression_rate!(du, u, p, t, oxidizer_mass_flow)
    oxidizer_mass_flux = oxidizer_mass_flow / p.fuel_grain_average_cross_sectional_area 
    #oh no, how do we get the oxidizer_mass_flow here
    #since we don't have access to the previous mass flow rate or previous timestamp, we can't get this
    #hmmm....

    oxidizer_mass_flux_g_per_cm2_s = 0.1 * oxidizer_mass_flux #convert to kg/(m^2*s)

    regression_mm_per_s = p.fuel_regression_coeff_a * oxidizer_mass_flux_g_per_cm2_s^p.fuel_regression_coeff_n

    regression_m_per_s = regression_mm_per_s * 0.001

    du.port_diameter += 2 * regression_m_per_s

    fuel_mass_flow = p.fuel_density * p.fuel_grain_burning_surface_area * regression_m_per_s

    du.chamber_gas_mass += fuel_mass_flow

    return fuel_mass_flow
end

function update_mass_fractions!(du, u, p, t, oxidizer_mass_flow, fuel_mass_flow)
    total_mass_flow_in = oxidizer_mass_flow + fuel_mass_flow

    du.chamber_nitrous_oxide_mass_fraction +=
        (oxidizer_mass_flow - u.chamber_nitrous_oxide_mass_fraction * total_mass_flow_in) /
        u.chamber_gas_mass

    du.chamber_hdpe_mass_fraction +=
        (fuel_mass_flow - u.chamber_hdpe_mass_fraction * total_mass_flow_in) /
        u.chamber_gas_mass
end


function nozzle_outlet!(du, u, p, t)
    chamber_gas_mass_flow_out = p.nozzle_discharge_coefficient * ((p.chamber_pressure * p.nozzle_throat_area) / p.propellant_characteristic_velocity)

    du.chamber_gas_mass -= chamber_gas_mass_flow_out

    return chamber_gas_mass_flow_out
end

function combustion_zone_energy_conservation!(du, u, p, t, model, oxidizer_mass_flow, fuel_mass_flow, chamber_gas_mass_flow_out)
    burning_fuel_mass_flow = min(
        fuel_mass_flow,
        oxidizer_mass_flow / p.stoichiometric_oxidizer_fuel_ratio
    )

    combustion_heat_release = p.combustion_efficiency * burning_fuel_mass_flow * p.fuel_lowering_heating_value

    pyrolysis_heat_absorption = burning_fuel_mass_flow * p.fuel_heat_of_pyrolysis

    fuel_specific_enthalpy = mass_enthalpy(model, p.chamber_pressure, p.fuel_surface_temperature, [0.0, 1.0], phase = :vapor)

    ejected_internal_energy = chamber_gas_mass_flow_out * p.chamber_specific_enthalpy

    du.chamber_gas_internal_energy += 
        oxidizer_mass_flow * p.mid_section_specific_enthalpy + 
        fuel_mass_flow * fuel_specific_enthalpy +
        combustion_heat_release - 
        pyrolysis_heat_absorption - 
        ejected_internal_energy -
        p.wall_heat_loss
end

u0 = ComponentVector(
    tank_oxidizer_mass = 0.0u"kg",
    tank_oxidizer_internal_energy = 0.0u"kJ",

    mid_section_mass = 1e-6u"kg",
    mid_section_internal_energy = 0.0u"kJ",

    chamber_gas_mass = 1e-6u"kg",
    chamber_gas_internal_energy = 0.0u"kJ",
    chamber_nitrous_oxide_mass_fraction = 1.0,
    chamber_hdpe_mass_fraction = 0.0,

    port_diameter = 2.0u"cm"
)

function update_u0!(u, p, t, oxidizer_model, chamber_model)
    #Tank
    u.tank_oxidizer_mass = p.u0_tank_oxidizer_mass

    tank_oxidizer_moles = p.u0_tank_oxidizer_mass / p.nitrous_oxide_molecular_weight

    p.tank_volume = volume(oxidizer_model, p.tank_pressure, p.tank_temperature, [tank_oxidizer_moles])

    u.tank_oxidizer_internal_energy = Clapeyron.VT0.internal_energy(oxidizer_model, p.tank_volume, p.tank_temperature, [tank_oxidizer_moles])

    #Mid Section
    u0_mid_section_density = mass_density(oxidizer_model, p.mid_section_pressure, p.mid_section_temperature, [tank_oxidizer_moles])
    u.mid_section_mass = u0_mid_section_density * p.mid_section_volume
    
    mid_section_oxidizer_moles = u.mid_section_mass / p.nitrous_oxide_molecular_weight
    u.mid_section_internal_energy = Clapeyron.VT0.internal_energy(oxidizer_model, p.mid_section_volume, p.mid_section_temperature, [mid_section_oxidizer_moles])

    #Chamber
    u0_chamber_gas_density = mass_density(chamber_model, p.chamber_pressure, p.chamber_temperature, [tank_oxidizer_moles, 0.0])
    u.chamber_gas_mass = u0_chamber_gas_density * p.chamber_volume

    chamber_gas_oxidizer_moles = (u.chamber_gas_mass * p.u0_chamber_nitrous_oxide_mass_fraction) / p.nitrous_oxide_molecular_weight
    chamber_gas_fuel_moles = (u.chamber_gas_mass * p.u0_chamber_hdpe_mass_fraction) / p.ethylene_molecular_weight
    u.chamber_gas_internal_energy = Clapeyron.VT0.internal_energy(chamber_model, p.chamber_volume, p.chamber_temperature, [chamber_gas_oxidizer_moles, chamber_gas_fuel_moles])

    #Fuel Grain
    u.port_diameter = p.u0_fuel_grain_void_diameter

    return nothing
end

function valve_opening_at_t(t)
    return 1.0
end

function valve_flow_capacity_factor(valve_opening, p)
    return p.valve_flow_capacity_factor 
    #since we don't have a physical valve that we've collected experimental data on yet, we're going to optimize the valve_flow_capacity_factor
    #to help us choose which electronic valve we should purchase
end

Revise.includet(joinpath(@__DIR__, "feed_system_update_state_for_ode.jl"))

#=
function update_state!(du, u, p, t, oxidizer_model, chamber_model)
    du .= 0.0 
    
    tank_n_moles = u.tank_oxidizer_mass / p.nitrous_oxide_molecular_weight

    show_debug = true

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

    if p.tank_vapor_fraction <= 0.999
        p.tank_density = mass_density(oxidizer_model, result, 1) #get liquid density because we drawing from the bottom of the tank
        p.tank_specific_enthalpy = mass_enthalpy(oxidizer_model, result, 1)
    else
        p.tank_density = mass_density(oxidizer_model, result) #otherwise, we will be drawing from the remaining vapor in the tank
        p.tank_specific_enthalpy = mass_enthalpy(oxidizer_model, result)
    end

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
    
    #Adjustable Valve
    p.valve_opening = valve_opening_at_t(t)
    p.valve_flow_capacity_factor = valve_flow_capacity_factor(p.valve_opening, p)

    #Injector
    #p.injector_orifice_area = p.injector_number_of_orifices * (pi / 4) * (p.injector_orifice_diameter^2)

    #Fuel Grain
    p.fuel_grain_average_cross_sectional_area = pi * (u.port_diameter / 2)^2
    
    p.fuel_grain_burning_surface_area = pi * u.port_diameter * p.fuel_grain_length
    
    p.fuel_mass = p.fuel_density * (pi * (p.final_fuel_grain_void_diameter / 2)^2 - pi * (u.port_diameter / 2)^2) * p.fuel_grain_length

    #Nozzle
    p.nozzle_throat_area = (pi / 4) * (p.nozzle_throat_diameter^2)
    #again, we need something here that determines propellant_isp based on oxidizer_to_fuel_ratio, chamber_pressure, and exit pressure (which we don't really know yet)
end
=#

function system_ode!(du_vec, u_vec, p_vec, t, oxidizer_model, chamber_model, p_axes, u_axes)
    du = ComponentVector(du_vec, u_axes)
    u = ComponentVector(u_vec, u_axes)
    p = ComponentVector(eltype(u).(p_vec), p_axes)

    update_state!(du, u, p, t, oxidizer_model, chamber_model)

    if true == false
        println("before")
        @show du
        @show u
        @show p
        println("")
    end

    adjustable_valve_pressure_drop = adjustable_valve_flow!(du, u, p, t)

    injector_valve_pressure_drop, oxidizer_mass_flow = injector_valve_flow!(du, u, p, t)
    
    fuel_mass_flow = regression_rate!(du, u, p, t, oxidizer_mass_flow)

    chamber_gas_mass_flow_out = nozzle_outlet!(du, u, p, t)

    update_mass_fractions!(du, u, p, t, oxidizer_mass_flow, fuel_mass_flow)

    combustion_zone_energy_conservation!(du, u, p, t, chamber_model, oxidizer_mass_flow, fuel_mass_flow, chamber_gas_mass_flow_out)

    if true == false
        println("after")
        @show du
        @show u
        @show p
        println("")
    end
end

properties = ComponentVector(
    #meta parameters
    simulation_time = 10.0u"s",

    #Overall rocket properties
    oxidizer_to_fuel_ratio = 7.6,
    desired_oxidizer_to_fuel_ratio = 7.6,
    propellant_isp = 220.0u"s",
    burn_time = 0.0u"s", #for viewing the optimized burn time 
    desired_burn_time = 5.0u"s",
    average_thrust = 0.0u"N", #for viewing the optimized averge thrust
    desired_average_thrust = 1000.0u"N",
    cummulative_impulse = 0.0u"N*s", #for viewing the optimized cumulative impulse
    desired_impulse = 0.0u"N*s", #this will be calculated based on desired_average thrust and desired_burn_time later
    #we would likely want a better correlation in the future that will return the specific impulse for a given oxidizer to fuel ratio and exit pressure
    gravity = 9.81u"m/s^2",
    target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio = 2.0, #not an optimized parameter, but the ideal choice is hard to know
    fuel_grain_max_diameter = (12.0u"inch" |> u"cm"),
    
    #Tank
    u0_tank_oxidizer_mass = 1.80u"kg", #optimized, should actually be 2.05kg to get the impulse we need, but I want to see if the optimizer will produce 2.05kg for debugging reasons
    tank_pressure = 71.0u"bar", #no longer optimized, commercial tanks determine this
    tank_temperature = 21.0u"°C",
    tank_vapor_fraction = 0.0u"m^3",
    tank_density = 0.0u"kg/m^3",
    tank_volume = 50.0u"cm^3",
    tank_specific_enthalpy = 0.0u"J/kg",
    
    #Nitrous oxide parameters
    nitrous_oxide_molecular_weight = 44.013u"g/mol",
    ethylene_molecular_weight = 28.05u"g/mol",

    #adjustable valve
    valve_flow_capacity_factor = 1e-5, #a result of the final valve we choose #we should optimize this to get a good ideal of what we want
    valve_opening = 1.0,

    #Mid section
    mid_section_pressure = 20.0u"bar", #this determines the initial kg of nitrous oxide in the mid_section
    mid_section_temperature = 21.0u"°C",
    mid_section_volume = 100.0u"cm^3", #Changed from 10u"cm" to 100u"cm" because I wanted to decrease solver stiffness for debugging purposes
    mid_section_density = 0.0u"kg/m^3",
    mid_section_vapor_fraction = 0.0,
    mid_section_specific_enthalpy = 0.0u"J/kg",

    #Injector properties
    injector_discharge_coefficient = 0.7,
    #injector_number_of_orifices = 20, #optimized
    #injector_orifice_diameter = 1.5u"mm", #optimized
    injector_orifice_area = 0.1u"cm^2", #optimized
    target_injector_velocity = 50.0u"m/s",

    #Chamber
    u0_chamber_nitrous_oxide_mass_fraction = 1.0,
    u0_chamber_hdpe_mass_fraction = 0.0,
    chamber_pressure = 10.0u"bar", #this determines the initial kg of nitrous oxide in the chamber
    chamber_temperature = 21.0u"°C",
    chamber_volume = 212.0u"cm^3",
    chamber_density = 0.0u"kg/m^3",
    wall_heat_loss = 0.0u"W",
    chamber_vapor_fraction = 0.0,
    chamber_specific_enthalpy = 0.0u"J/kg",


    #Fuel grain
    u0_fuel_grain_void_diameter = 3.0u"cm", #optimized
    fuel_mass = 0.0u"kg", #this will be derived by substracting the volume of the cylinder formed by the u0_fuel_grain_void_diameter by the final_fuel_grain_void_diameter and then multiplying by the fuel density
    fuel_density = 950.0u"kg/m^3",
    #fuel_regression_rate = 0.5u"mm/s",
    additional_fuel_grain_void_diameter = 1.0u"cm", #optimized
    final_fuel_grain_void_diameter = 0.0u"cm",
    fuel_grain_length = 30.0u"cm", #optimized
    fuel_grain_average_cross_sectional_area = 0.0u"m^2",
    fuel_grain_burning_surface_area = 0.0u"m^2",

    #Fuel Grain Empirical Parameters
    #These are for when fuel_regression is measured in mm/s and oxidizer_mass_flux is measured in g/(cm^2*s)
    fuel_regression_coeff_a = 0.248,
    fuel_regression_coeff_n = 0.331,

    #Fuel properties
    fuel_lowering_heating_value = 47.42u"MJ/kg",
    stoichiometric_oxidizer_fuel_ratio = 9.41,
    fuel_heat_of_pyrolysis = 2.3u"MJ/kg",
    fuel_surface_temperature = 1300.0u"K",
    combustion_efficiency = 0.9,

    #Propellant properties
    propellant_characteristic_velocity = 1500.0u"m/s", 

    #Nozzle
    nozzle_throat_diameter = 1.6u"cm", #optimized
    nozzle_throat_area = 0.0u"m^2",
    nozzle_discharge_coefficient = 0.9,
)

struct OptimizedParameter
    name::Symbol
    lb::Number
    ub::Number
end

optimized_properties = [
    #Overall rocket properties

    #Tank
    OptimizedParameter(:u0_tank_oxidizer_mass, 0.5u"kg", 5.0u"kg",),
    #OptimizedParameter(:tank_pressure, 50.0u"bar", 71.0u"bar",),

    #adjustable valve
    OptimizedParameter(:valve_flow_capacity_factor, 1e-7u"m^2", 1e-2u"m^2"),

    #Mid section

    #Injector properties
    OptimizedParameter(:injector_orifice_area, 0.01u"cm^2", 10.0u"cm^2"),

    #Chamber

    #Fuel grain
    OptimizedParameter(:u0_fuel_grain_void_diameter, 1.0u"cm", 10.0u"cm"),
    #OptimizedParameter(:final_fuel_grain_void_diameter, 1.0u"cm", 20.0u"cm"),
    OptimizedParameter(:additional_fuel_grain_void_diameter, 0.1u"cm", 3.0u"cm"),
    OptimizedParameter(:fuel_grain_length, 5.0u"cm", 50.0u"cm"),

    #Fuel Grain Empirical Parameters

    #Propellant properties

    #Nozzle
    OptimizedParameter(:nozzle_throat_diameter, 0.2u"cm", 2.5u"cm")
]

theta_guess, theta_lb, theta_ub, theta_axes, u_axes, p_axes, theta_to_u_map, theta_to_p_map, p_to_u_map = create_theta_guess(optimized_properties, u0, properties)

function trainsient_system_design_loss(theta, u0, p, theta_axes, u_axes, p_axes, oxidizer_model, chamber_model, theta_to_u_map, theta_to_p_map, p_to_u_map, append_optimized_parameters!, update_u0!, isoutofdomain_feedsystem, port_diameter_termination_cb, problem_template)
    theta = ComponentVector(theta, theta_axes)
    
    # Promote p to the type of theta (which will be Dual during ForwardDiff) so it can accept Duals
    u0 = eltype(theta).(u0)
    p = ComponentVector(eltype(theta).(p), p_axes)

    append_optimized_parameters!(theta, u0, p, theta_to_u_map, theta_to_p_map, p_to_u_map)
    update_u0!(u0, p, 0.0, oxidizer_model, chamber_model)

    #=
    prob = get!(task_local_storage(), :implicit_prob) do    
        # Build a new closure bound to the thread-isolated copies
        f_closure = (du, u, p, t) -> system_ode!(du, u, p, t, oxidizer_model, chamber_model, p_axes)

        ode_func = ODEFunction(f_closure)#, jac_prototype = float.(jac_sparsity))

        tspan = (0.0, p.simulation_time)

        prob = ODEProblem(ode_func, u0, tspan, p)
    end
    =#
    prob = remake(problem_template; u0 = Vector(u0), p = Vector(p), tspan = (0.0, p.simulation_time))

    sol = 0

    @show "before ode solve"
    #try
    sol = solve(prob,
        saveat = (p.simulation_time / 200),
        callback = CallbackSet(port_diameter_termination_cb, approximate_time_to_finish_cb),
        isoutofdomain = isoutofdomain_feedsystem,
        maxiters = 1000,
        dtmin = 1e-5
    )
    #catch e
        #println(e)
        #return theta, p, 1e10
    #end
    @show "after ode solve"

    @show sol.retcode

    if !(sol.retcode == SciMLBase.ReturnCode.Success || sol.retcode == SciMLBase.ReturnCode.Terminated)
        return theta, p, (1e10 + p_named.final_fuel_grain_void_diameter - p_named.u0_fuel_grain_void_diameter) * 1e8 #we want the solver to prioritize runs that got closer to expending all the fuel even if they failed
    end


    #Losses updated every iteration
    injector_velocity_loss = 0.0
    pressure_drop_ratio_loss = 0.0
    oxidizer_to_fuel_ratio_loss = 0.0
    thrust_loss = 0.0
    
    #Loss updated once
    unburned_fuel_loss = 0.0
    unutilized_oxidizer_loss = 0.0

    #Performance Metrics
    cummulative_impulse = 0.0
    depletion_time = 0.0

    du_temporary = ComponentVector(similar(sol.u[1]), u_axes)

    found_depletion_time = false

    last_fuel_mass = p.fuel_mass

    for i in eachindex(sol.u)
        curr_t = sol.t[i]
        u_named = ComponentVector(sol.u[i], u_axes)

        update_state!(du_temporary, u_named, p, curr_t, oxidizer_model, chamber_model)

        if (u_named.port_diameter >= p.final_fuel_grain_void_diameter) && found_depletion_time == false
            depletion_time = curr_t
            found_depletion_time = true
            unutilized_oxidizer_loss = 0.001 * abs2(u_named.tank_oxidizer_mass)
            break #stop evaluating loss if the fuel has burned out
        end

        if (u_named.tank_oxidizer_mass <= 1e-6) && found_depletion_time == false
            depletion_time = curr_t
            found_depletion_time = true
            unburned_fuel_loss = 0.001 * abs2(p.fuel_mass)
            break #stop evaluating loss if the oxidizer has burned out
        end

        adjustable_valve_pressure_drop = adjustable_valve_flow!(du_temporary, u_named, p, curr_t)

        injector_valve_pressure_drop, oxidizer_mass_flow = injector_valve_flow!(du_temporary, u_named, p, curr_t)

        #@show p.target_injector_velocity
        #@show oxidizer_mass_flow / (p.mid_section_density * p.injector_orifice_area)
        injector_velocity_loss += 0.00001 * (1 / length(sol.u)) * abs2(p.target_injector_velocity - oxidizer_mass_flow / (p.mid_section_density * p.injector_orifice_area))

        #@show p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio
        #@show injector_valve_pressure_drop / adjustable_valve_pressure_drop
        pressure_drop_ratio_loss += 0.00001 * (1 / length(sol.u)) * abs2(p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio - (injector_valve_pressure_drop / adjustable_valve_pressure_drop))

        #OBSERVATION: it seems like we're going to have to weigh the importance of injector velocity against adjustable valve authority

        fuel_mass_flow = regression_rate!(du_temporary, u_named, p, curr_t, oxidizer_mass_flow)

        #=
        chamber_gas_mass_flow_out = nozzle_outlet!(du_temporary, u_named, p, curr_t)

        update_mass_fractions!(du_temporary, u_named, p, curr_t, oxidizer_mass_flow, fuel_mass_flow)

        combustion_zone_energy_conservation!(du_temporary, u_named, p, curr_t, chamber_model, oxidizer_mass_flow, fuel_mass_flow, chamber_gas_mass_flow_out)
        =#

        #if we want, we can get the derivative of the cummulative_impulse over time to plot the thrust profile of the engine!
        if i > 1
            u_named_prev = ComponentVector(sol.u[i-1], u_axes)
            dt = sol.t[i] - sol.t[i-1]

            #Cummulative impulse loss calcs
            oxidizer_used = -(u_named.tank_oxidizer_mass - u_named_prev.tank_oxidizer_mass)
            fuel_used = -(p.fuel_mass - last_fuel_mass)
            last_fuel_mass = p.fuel_mass

            propellant_used = oxidizer_used + fuel_used

            cummulative_impulse += p.propellant_isp * p.gravity * propellant_used

            #Oxidizer to fuel ratio loss
            oxidizer_to_fuel_ratio_loss += (1 / length(sol.u)) * 0.001 * abs2(p.desired_oxidizer_to_fuel_ratio - (oxidizer_used / fuel_used))

            #Thrust loss calcs
            oxidizer_mass_flow = oxidizer_used / dt
            fuel_mass_flow = fuel_used / dt

            propallant_mass_flow = oxidizer_mass_flow + fuel_mass_flow

            thrust_produced = p.propellant_isp * p.gravity * propallant_mass_flow
            
            p.average_thrust += (1 / (length(sol.t) - 1)) * thrust_produced
            thrust_loss += (1 / (length(sol.t) - 1)) * 0.001 * abs2(p.desired_average_thrust - thrust_produced)
        end

        if i == length(sol.u)
            if depletion_time == 0.0 && found_depletion_time == false
                @warn "the simulation did no run long enough to completely burn out the fuel grain"
                @show remaining_fuel_grain = (u_named.port_diameter - p.final_fuel_grain_void_diameter)

                @show u_named.port_diameter
                @show p.final_fuel_grain_void_diameter
            end
        end
    end

    p.cummulative_impulse = cummulative_impulse
    impulse_loss = 0.0001 * abs2(p.desired_average_thrust - cummulative_impulse)

    p.burn_time = depletion_time
    burn_time_loss = 0.01 * abs2(p.desired_burn_time - depletion_time)

    #above_max_fuel_grain_diameter_loss = 0.01 * abs2(p.u0_fuel_grain_void_diameter + p.additional_fuel_grain_void_diameter - p.fuel_grain_max_diameter)
    #we'll just enforce a bound on final_fuel_grain_void_diameter

    all_losses = ComponentVector(
        #Updated every solver iteration
        injector_velocity_loss = injector_velocity_loss,
        pressure_drop_ratio_loss = pressure_drop_ratio_loss,
        oxidizer_to_fuel_ratio_loss = oxidizer_to_fuel_ratio_loss,
        thrust_loss = thrust_loss,

        #Updated once at the end of the simulation
        impulse_loss = impulse_loss,
        burn_time_loss = burn_time_loss,
        unburned_fuel_loss = unburned_fuel_loss,
        unutilized_oxidizer_loss = unutilized_oxidizer_loss
    )

    @show sum(all_losses)

    return theta, p, all_losses
end

theta_guess, theta_lb, theta_ub, theta_axes, u_axes, p_axes, theta_to_u_map, theta_to_p_map, p_to_u_map = create_theta_guess(optimized_properties, u0, properties)

theta_to_u_map
theta_to_p_map
p_to_u_map

#convert    all to base SI and then strip away units
theta_guess_unitless = ustrip.(upreferred.(theta_guess))
theta_lb_unitless = ustrip.(upreferred.(theta_lb))
theta_ub_unitless = ustrip.(upreferred.(theta_ub))
u0_unitless = ustrip.(upreferred.(u0))
du_test = deepcopy(u0_unitless)
properties_unitless = ustrip.(upreferred.(properties))

p_axes = getaxes(properties)
u_axes = getaxes(u0)

f_closure = let
    om = oxidizer_model
    cm = chamber_model
    p_axes_local = p_axes
    u_axes_local = u_axes

    (du, u, p, t) -> system_ode!(du, u, p, t, om, cm, p_axes_local, u_axes_local)
end


detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, properties_unitless, 0.0), Vector(du_test), Vector(u0_unitless), detector
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 10.0
tspan = (t0, tMax)

append_optimized_parameters!(Vector(theta_guess_unitless), u0_unitless, properties_unitless, theta_to_u_map, theta_to_p_map, p_to_u_map)
update_u0!(u0_unitless, properties_unitless, 0.0, oxidizer_model, chamber_model)

implicit_prob = ODEProblem(ode_func, Vector(u0_unitless), tspan, Vector(properties_unitless))

system_ode!(du_test, Vector(u0_unitless), Vector(properties_unitless), 0.0, oxidizer_model, chamber_model, p_axes, u_axes)

f_closure(du_test, u0_unitless, properties_unitless, 0.0)

function isoutofdomain_feedsystem_expanded(u, p, t, u_axes)
    u_named = ComponentVector(u, u_axes)

    #=
    if u_named.tank_oxidizer_mass < 0.0 || u_named.mid_section_mass < 0.0 || u_named.chamber_gas_mass < 0.0
        @show "caught out of domain"
    end
    =#

    return (
        u_named.tank_oxidizer_mass < 0.0 ||
        u_named.mid_section_mass < 0.0 ||
        u_named.chamber_gas_mass < 0.0
    )
end

isoutofdomain_feedsystem = let
    u_axes_local = u_axes
    (u, p, t) -> isoutofdomain_feedsystem_expanded(u, p, t, u_axes_local)
end

# Stop the burn when the fuel port reaches the outside diameter of the grain.
# Only the positive crossing is active because regression increases port_diameter.
function port_diameter_limit(u, t, integrator, u_axes, p_axes)
    u_named = ComponentVector(u, u_axes)
    p_named = ComponentVector(integrator.p, p_axes)

    update_state!([0.0], u_named, p_named, t, oxidizer_model, chamber_model)

    #@show u_named.port_diameter - p_named.final_fuel_grain_void_diameter

    #if it gets within 1e-5 m of the final grain diameter, terminate to allow stiff problems to still solve
    u_named.port_diameter - p_named.final_fuel_grain_void_diameter + 1e-5 #0.01 mm 
end

port_diameter_termination_cb = ContinuousCallback(
    (u, t, integrator) -> port_diameter_limit(u, t, integrator, u_axes, p_axes),
    terminate!,
    nothing;
    save_positions = (true, false),
)

sol = solve(
    implicit_prob,
    #FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, ),
    #FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true),
    #AutoTsit5(FBDF(linsolve = SparspakFactorization(), autodiff = AutoFiniteDiff())),
    callback = CallbackSet(
        port_diameter_termination_cb,
        approximate_time_to_finish_cb,
    ),
    isoutofdomain = isoutofdomain_feedsystem
)

viewable_system_design_loss_closure = let 
    u_local = u0_unitless
    theta_axes_local = theta_axes
    u_axes_local = u_axes
    p_axes_local = p_axes
    oxidizer_model_local = oxidizer_model
    chamber_model_local = chamber_model
    theta_to_u_map_local = theta_to_u_map
    theta_to_p_map_local = theta_to_p_map
    p_to_u_map_local = p_to_u_map
    append_optimized_parameters_local! = append_optimized_parameters!
    update_u0_local! = update_u0!
    port_diameter_termination_cb_local = port_diameter_termination_cb
    isoutofdomain_feedsystem_local = isoutofdomain_feedsystem
    problem_template_local = implicit_prob
    
    (theta, p) -> trainsient_system_design_loss(
        theta, u_local, p, 
        theta_axes_local, u_axes_local, p_axes_local, 
        oxidizer_model_local, chamber_model_local, 
        theta_to_u_map_local, theta_to_p_map_local, p_to_u_map_local, 
        append_optimized_parameters_local!, update_u0_local!, 
        isoutofdomain_feedsystem_local, port_diameter_termination_cb_local,
        problem_template_local
    )
end

function viewable_system_design_loss(theta, p)
    theta, p, all_losses = viewable_system_design_loss_closure(theta, p)
    @show theta
    println("")
    @show p
    println("")
    @show all_losses
    println("")
    return sum(all_losses)
end

system_design_loss_closure = let 
    u_local = u0_unitless
    theta_axes_local = theta_axes
    u_axes_local = u_axes
    p_axes_local = p_axes
    oxidizer_model_local = oxidizer_model
    chamber_model_local = chamber_model
    theta_to_u_map_local = theta_to_u_map
    theta_to_p_map_local = theta_to_p_map
    p_to_u_map_local = p_to_u_map
    append_optimized_parameters_local! = append_optimized_parameters!
    update_u0_local! = update_u0!
    isoutofdomain_feedsystem_local = isoutofdomain_feedsystem
    port_diameter_termination_cb_local = port_diameter_termination_cb
    problem_template_local = implicit_prob

    (theta, p) -> trainsient_system_design_loss(
        theta, u_local, p, 
        theta_axes_local, u_axes_local, p_axes_local, 
        oxidizer_model_local, chamber_model_local, 
        theta_to_u_map_local, theta_to_p_map_local, p_to_u_map_local, 
        append_optimized_parameters_local!, update_u0_local!, 
        isoutofdomain_feedsystem_local, port_diameter_termination_cb_local,
        problem_template_local
    )
end

pure_system_design_loss_closure = let
    u_initial = u0_unitless
    θ_axes = theta_axes
    state_axes = u_axes
    parameter_axes = p_axes
    om = oxidizer_model
    cm = chamber_model
    θ_to_u = theta_to_u_map
    θ_to_p = theta_to_p_map
    p_to_u = p_to_u_map
    append_parameters! = append_optimized_parameters!
    initialize_state! = update_u0!
    isoutofdomain_feedsystem_local = isoutofdomain_feedsystem
    port_diameter_termination_cb_local = port_diameter_termination_cb
    problem_template_local = implicit_prob

    function (theta, p)
        _, _, losses = trainsient_system_design_loss(
            theta, u_initial, p,
            θ_axes, state_axes, parameter_axes,
            om, cm,
            θ_to_u, θ_to_p, p_to_u,
            append_parameters!, initialize_state!,
            isoutofdomain_feedsystem_local, port_diameter_termination_cb_local,
            problem_template_local
        )

        return sum(losses)
    end
end

opt_f = OptimizationFunction(pure_system_design_loss_closure, Optimization.AutoFiniteDiff())
opt_prob = OptimizationProblem(opt_f, Vector(theta_guess_unitless), Vector(properties_unitless), lb = Vector(theta_lb_unitless), ub = Vector(theta_ub_unitless))

LOSS = Float64[]
PARS = []

# Ensure the directory exists and use a filename-safe date format (colons are invalid on Windows)
mkpath(joinpath(@__DIR__, "optimization_results"))
results_path = joinpath(@__DIR__, "optimization_results", "optimization_results_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).csv")

# Create the file and manually write the header string using propertynames
open(results_path, "w") do io
    header_str = "loss," * join(string.(propertynames(theta_guess_unitless)), ",")
    println(io, header_str)
end

const cb_lock = ReentrantLock()

cb = function (state, l)
    display(l)
    display(state.u)
    
    lock(cb_lock) do
        push!(LOSS, l)
        push!(PARS, state.u)
        
        # Convert state.u to a named tuple using your p_axes so the CSV has nice column headers
        theta_named = NamedTuple(ComponentVector(state.u, theta_axes))
        row = merge((loss = l, ), theta_named)
        
        # CSV.write with append=true automatically opens, appends, and closes (flushes) the file
        CSV.write(results_path, DataFrame([row]), append=true)
    end
    
    false
end

pure_system_design_loss_closure(theta_guess_unitless, properties_unitless)

sol = solve(opt_prob, LBFGS(), callback = cb, reltol = 1e-4)

sol = solve(opt_prob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(),
    callback = cb,
    PoulationSize = 100,
    #maxiters = 1,
    #maxtime = 60.0,
    Method = :RandomSearcher,
    verbose = true
    #Method = :SepReal
)

losses = Float64[]
successful_parameters = []

parameter_sweep_steps = 5
sweep_parameters = Vector(theta_guess_unitless)

for parameter_index in eachindex(sweep_parameters)
    starting_value = sweep_parameters[parameter_index]

    for bound in (theta_lb_unitless[parameter_index], theta_ub_unitless[parameter_index])
        bound == starting_value && continue

        parameter_values = range(
            starting_value,
            stop = bound,
            length = parameter_sweep_steps + 1,
        )

        for parameter_value in Iterators.drop(parameter_values, 1)
            sweep_parameters[parameter_index] = parameter_value
            l = pure_system_design_loss_closure(
                sweep_parameters,
                properties_unitless,
            )

            if l <= 1e9
                push!(losses, l)
                push!(successful_parameters, copy(sweep_parameters))
            end
        end
    end

    sweep_parameters[parameter_index] = starting_value
end

@show losses
@show successful_parameters

for i in 1:1000
    θ = theta_lb_unitless .+
        rand(length(theta_lb_unitless)) .*
        (theta_ub_unitless .- theta_lb_unitless)

    @show i θ

    l = pure_system_design_loss_closure(
        θ,
        properties_unitless,
    )

    if l <= 1e9
        push!(losses, l)
        push!(successful_parameters, θ)
    end
end

@show losses
@show successful_parameters

#viewable_system_design_loss(sol.u, properties_unitless)
viewable_system_design_loss(new_theta_guess, properties_unitless)

final_properties = ComponentVector(sol.u, u_axes)
