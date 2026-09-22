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

using FVMFramework

Revise.includet(joinpath(@__DIR__, "design_helper_functions.jl"))

model = PR(["nitrous oxide"])

function adjustable_valve_flow!(du, u, p, t, model)
    volumetric_flow = p.valve_flow_capacity_factor * p.valve_opening * sqrt((p.tank_pressure - p.mid_section_pressure) / nitrous_oxide_density)
    #we might need a different equation that activates when the tank starts spitting out gas

    mass_flow = volumetric_flow * u.mid_section_density

    du.tank_oxidizer_mass -= mass_flow
    du.tank_oxidizer_internal_energy -= mass_flow * p.tank_specific_enthalpy
    du.mid_section_mass += mass_flow
    du.mid_section_internal_energy += mass_flow * p.tank_specific_enthalpy

    return (p.tank_pressure - p.mid_section_pressure)
end

function injector_valve_flow!(du, u, p, t, model)
    pre_injector_density = mass_density(model, p.mid_section_pressure, p.mid_section_temperature, [1.0])

    volumetric_flow = p.injector_discharge_coefficient * p.injector_orifice_area * sqrt((2 * (p.mid_section_pressure - p.chamber_pressure)) / pre_injector_density)

    mass_flow = volumetric_flow * p.chamber_density

    du.mid_esction_mass -= mass_flow
    du.mid_section_internal_energy -= mass_flow * p.tank_specific_enthalpy
    du.chamber_gas_mass += mass_flow
    du.chamber_gas_internal += mass_flow * p.tank_specific_enthalpy

    return (p.mid_section_pressure - p.chamber_pressure)
end

function regression_rate!(du, u, p, t, model)
    oxidizer_mass_flux = p.oxidizer_mass_flow / p.fuel_grain_average_cross_sectional_area 
    #oh no, how do we get the oxidizer_mass_flow here
    #since we don't have access to the previous mass flow rate or previous timestamp, we can't get this
    #hmmm....

    oxidizer_mass_flux_g_per_cm2_s = 0.1 * oxidizer_mass_flux #convert to kg/(m^2*s)

    regression_mm_per_s = p.fuel_regression_coeff_a * oxidizer_mass_flux_g_per_cm2_s^p.fuel_regression_coeff_n

    regression_m_per_s = regression_mm_per_s * 0.001

    du.port_diameter -= regression_m_per_s

    fuel_flow = p.fuel_density * p.fuel_grain_burning_surface_area * p.fuel_regression_rate

    return fuel_flow
end

function nozzle_outlet!(du, u, p, t, model)
    chamber_gas_mass_flow_out = p.nozzle_discharge_coefficient * ((p.chamber_pressure * p.nozzle_throat_area) / p.propellant_characteristic_velocity)

    du.chamber_gas_mass -= chamber_gas_mass_flow_out
    du.chamber_gas_internal_energy -= chamber_gas_mass_flow_out * p.propellant_enthalpy
end

u0 = ComponentVector(
    tank_oxidizer_mass = 0.0u"kg",
    tank_oxidizer_internal_energy = 0.0u"kJ",

    mid_section_mass = 1e-6u"kg",
    mid_section_internal_energy = 0.0u"kJ",

    chamber_gas_mass = 1e-6u"kg",
    chamber_gas_internal_energy = 0.0u"kJ",

    port_radius = 2.0u"cm"
)

function update_u0!(du, u, p, t, model)
    u.tank_oxidizer_mass = p.u0_tank_oxidizer_mass

    u.tank_oxidizer_internal_energy = Clapeyron.VT0.internal_energy(model, p.tank_volume, p.tank_temperature, [1.0])

    u.mid_section_mass = p.u0_mid_section_mass

    u.mid_section_internal_energy = Clapeyron.VT0.internal_energy(model, p.mid_section_volume, p.mid_section_temperature, [1.0])

    u.chamber_gas_mass = p.u0_chamber_gas_mass

    u.chamber_gas_internal_energy = Clapeyron.VT0.internal_energy(model, p.chamber_volume, p.chamber_temperature, [1.0])

    u.port_radius = p.u0_fuel_grain_void_diameter
end

function valve_opening_at_t(t)
    return 1.0
end

function valve_flow_capacity_factor(valve_opening, p)
    return p.valve_flow_capacity_factor 
    #since we don't have a physical valve that we've collected experimental data on yet, we're going to optimize the valve_flow_capacity_factor
    #to help us choose which electronic valve we should purchase
end

function update_state!(du, u, p, t, model)
    tank_n_moles = u.tank_oxidizer_mass / u.nitrous_oxide_molecular_weight

    result = uv_flash(model, u.tank_oxidizer_internal_energy, u.tank_volume, [tank_n_moles])

    p.tank_temperature = result.data.T
    p.tank_pressure = pressure(model, result)
    p.tank_vapor_fraction = result.amounts[1]

    if p.tank_vapor_fraction <= 0.999
        p.tank_density = mass_density(model, result, 1) #get liquid density because we drawing from the bottom of the tank
        p.tank_specific_enthalpy = mass_enthalpy(model, result, 1)
    else
        p.tank_density = mass_density(model, result) #otherwise, we will be drawing from the remaining vapor in the tank
        p.tank_specific_enthalpy = mass_enthalpy(model, result)
    end

    mid_section_n_moles = u.mid_section_mass / u.nitrous_oxide_molecular_weight

    result = uv_flash(model, u.mid_section_internal_energy, u.mid_section_volume, [mid_section_n_moles])

    p.mid_section_temperature = result.data.T
    p.mid_section_pressure = pressure(model, result)
    p.mid_section_vapor_fraction = result.amounts[1]
    p.mid_section_density = mass_density(model, result)
    p.mid_section_specific_enthalpy = mass_enthalpy(model, result)

    chamber_n_moles = u.chamber_gas_mass / u.nitrous_oxide_molecular_weight

    result = uv_flash(model, u.chamber_gas_internal_energy, u.chamber_volume, [chamber_n_moles])
    #Should we try to calculate chamber volume based on other parameters or should we just optimize it?

    p.chamber_temperature = result.data.T
    p.chamber_pressure = pressure(model, result)
    p.chamber_vapor_fraction = result.amounts[1]
    p.chamber_density = mass_density(model, result)
    p.chamber_specific_enthalpy = mass_enthalpy(model, result)

    #Other state updates:
    #Adjustable Valve
    p.valve_opening = valve_opening_at_t(t)
    p.valve_flow_capacity_factor = valve_flow_capacity_factor(p.valve_opening, p)

    #Injector
    #p.injector_orifice_area = p.injector_number_of_orifices * (pi / 4) * (p.injector_orifice_diameter^2)

    #Fuel Grain
    p.fuel_grain_average_cross_sectional_area = pi * (u.port_diameter / 2)^2
    
    p.fuel_grain_burning_surface_area = pi * u.port_diameter * u.fuel_grain_length
    
    p.fuel_mass = p.fuel_density * (pi * (u.final_fuel_grain_void_diameter / 2)^2 - pi * (u.port_diameter / 2)^2)

    #Nozzle
    p.nozzle_throat_area = (pi / 4) * (p.nozzle_throat_diameter^2)
    #again, we need something here that determines propellant_isp based on oxidizer_to_fuel_ratio, chamber_pressure, and exit pressure (which we don't really know yet)
end

properties = ComponentVector(
    #Overall rocket properties
    oxidizer_mass_flow = 0.0u"kg/s", #this will be a cache for the previous oxidizer mass flow
    fuel_mass_flow = 0.0u"kg/s", #this will be a cache for the previous fuel mass flow
    oxidizer_to_fuel_ratio = 7.6,
    desired_oxidizer_to_fuel_ratio = 7.6,
    propellant_isp = 220.0u"s",
    burn_time = 0.0u"s", 
    desired_burn_time = 5.0u"s",
    average_thrust = 0.0u"N", #determined by specific_impulse * (oxidizer_mass_flow * fuel_mass_flow) * gravity
    desired_average_thrust = 1000.0u"N",
    #we would likely want a better correlation in the future that will return the specific impulse for a given oxidizer to fuel ratio and exit pressure
    gravity = 9.81u"m/s^2",
    target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio = 2.0, #not an optimized parameter, but the ideal choice is hard to know
    fuel_grain_max_diameter = (12.0u"inch" |> u"cm"),

    #Tank
    tank_pressure = 70u"bar", #optimized
    tank_temperature = 21.0u"°C",
    u0_tank_oxidizer_mass = 5.0u"kg", #optimized
    tank_vapor_fraction = 0.0u"m^3",
    tank_density = 0.0u"kg/m^3",
    tank_volume = 50.0u"cm^3",
    
    #Nitrous oxide parameters
    nitrous_oxide_molecular_weight = 44.013u"g/mol",

    #adjustable valve
    valve_flow_capacity_factor = 1e-5, #a result of the final valve we choose #we should optimize this to get a good ideal of what we want
    valve_opening = 1.0,

    #Mid section
    mid_section_pressure = 40.0u"bar",
    mid_section_temperature = 21.0u"°C",
    mid_section_volume = 10.0u"cm^3",
    u0_mid_section_mass = 1e-6u"kg",
    mid_section_density = 0.0u"kg/m^3",

    #Injector properties
    injector_discharge_coefficient = 0.7,
    #injector_number_of_orifices = 20, #optimized
    #injector_orifice_diameter = 1.5u"mm", #optimized
    injector_orifice_area = 0.1u"cm^2", #optimized
    target_injector_velocity = 50.0u"m/s",

    #Chamber
    chamber_pressure = 30.0u"bar",
    chamber_temperature = 0.0u"°C",
    chamber_volume = 50.0u"cm^3",
    u0_chamber_gas_mass = 1e-6u"kg",
    chamber_density = 0.0u"kg/m^3",

    #Fuel grain
    fuel_mass = 0.0u"kg", #this will be derived by substracting the volume of the cylinder formed by the u0_fuel_grain_void_diameter by the final_fuel_grain_void_diameter and then multiplying by the fuel density
    fuel_density = 950.0u"kg/m^3",
    fuel_regression_rate = 0.5u"mm/s",
    u0_fuel_grain_void_diameter = 3.0u"cm", #optimized
    final_fuel_grain_void_diameter = 0.0u"cm", #optimized
    fuel_grain_length = 30.0u"cm", #optimized
    fuel_grain_average_cross_sectional_area = 0.0u"m^2",
    fuel_grain_burning_surface_area = 0.0u"m^2",

    #Fuel Grain Empirical Parameters
    #These are for when fuel_regression is measured in mm/s and oxidizer_mass_flux is measured in g/(cm^2*s)
    fuel_regression_coeff_a = 0.248,
    fuel_regression_coeff_n = 0.331,

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
    OptimizedParameter(:oxidizer_mass_flow, 1e-6u"kg/s", properties.oxidizer_mass_flow * 10.0),

    #Tank
    OptimizedParameter(:u0_tank_oxidizer_mass, 1e-6u"kg", 50.0u"kg",),
    OptimizedParameter(:tank_volume, 10.0u"cm^3", 30.0u"L",),

    #adjustable valve
    OptimizedParameter(:valve_flow_capacity_factor, 1e-7u"m^2", 1e-2u"m^2"),

    #Mid section

    #Injector properties
    OptimizedParameter(:injector_orifice_area, 0.01u"cm^2", 10.0u"cm^2"),

    #Chamber

    #Fuel grain
    OptimizedParameter(:u0_fuel_grain_void_diameter, 1.0u"cm", 10.0u"cm"),
    OptimizedParameter(:final_fuel_grain_void_diameter, 1.0u"cm", 20.0u"cm"),
    OptimizedParameter(:fuel_grain_length, 5.0u"cm", 200.0u"cm"),

    #Fuel Grain Empirical Parameters

    #Propellant properties

    #Nozzle
    OptimizedParameter(:nozzle_throat_diameter, 0.2u"cm", 2.5u"cm")
]

function system_ode!(du, u, p, t, model, p_axes)
    p = ComponentVector(eltype(u).(p), p_axes)

    update_state!(du, u, p, t, model)

    adjustable_valve_pressure_drop = adjustable_valve_flow(du, u, p, t, model)

    injector_valve_pressure_drop = injector_valve_flow(du, u, p, t, model)

    fuel_flow = regression_rate(du, u, p, t, model)
end

#=
prob = get!(task_local_storage(), :dry_run_implicit_prob) do    
    # Build a new closure bound to the thread-isolated copies
    f_closure = (du, u, p, t) -> system_ode!(du, u, p, t, model, p_axes)
    
    #=
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
    )
    =#

    ode_func = ODEFunction(f_closure)#, jac_prototype = float.(jac_sparsity))

    t0 = 0.0
    tMax = ustrip(upreferred(thermocouple_data.timestamps[end]))
    tspan = (t0, tMax)

    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)
end
=#

function system_design_loss(theta, p, theta_axes, u_axes, p_axes, model, theta_to_p_map, theta_to_u_map, p_to_u_map, append_optimized_parameters!, update_properties!)
    theta = ComponentVector(theta, theta_axes)
    
    # Promote p to the type of theta (which will be Dual during ForwardDiff) so it can accept Duals
    p = ComponentVector(eltype(theta).(p), p_axes)

    append_optimized_parameters!(theta, p, theta_to_p_map, theta_to_u_map, p_to_u_map)
    update_u0!(du, u, p, t, model)

    #we run the ODEProblem here

    #Losses
    injector_velocity_loss = 0.0
    pressure_drop_ratio_loss = 0.0

    #Performance Metrics
    cummulative_impulse = 0.0
    fuel_depletion_time = 0.0
    oxidizer_depletion_time = 0.0

    du_temporary = similar(sol.u[1])

    for i in eachindex(sol.u)
        curr_t = sol.t[i]
        u_named = ComponentVector(sol.u[i], u_axes)

        update_state!(du_temporary, u_named, p, curr_t, model)

        adjustable_valve_pressure_drop = adjustable_valve_flow(du_temporary, u_named, p, curr_t, model)

        injector_valve_pressure_drop = injector_valve_flow(du_temporary, u_named, p, curr_t, model)
        injector_velocity_loss = 0.00001 * abs2(p.target_injector_velocity - injector_valve_oxidizer_mass_flow / (pre_injector_density * p.injector_orifice_area))

        pressure_drop_ratio_loss = 0.00001 * abs2(p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio - (injector_valve_pressure_drop / adjustable_valve_pressure_drop))

        #OBSERVATION: it seems like we're going to have to weigh the importance of injector velocity against adjustable valve authority

        if (abs(u.port_diameter - u.final_fuel_grain_void_diameter) <= 1e-6)
            fuel_depletion_time = curr_t
        end

        if (u.tank_oxidizer_mass <= 1e-6)
            oxidizer_depletion_time = curr_t
        end

        #if we want, we can get the derivative of the cummulative_impulse over time to plot the thrust profile of the engine!
        if i > 1
            u_named_prev = ComponentVector(sol.u[i-1], u_axes)

            change_in_oxidizer_mass = u_named.tank_oxidizer_mass - u_named_prev.tank_oxidizer_mass
            change_in_fuel_mass = u_named.fuel_mass - u_named_prev.fuel_mass

            change_in_propellant_mass = change_in_oxidizer_mass + change_in_fuel_mass

            cummulative_impulse += p.propallant_isp * p.gravity * change_in_propellant_mass
        end
    end

    burn_time_loss = 0.0001 * abs2(p.desired_burn_time - fuel_depletion_time)
    average_thrust_loss = 0.0001 * abs2(p.desired_average_thrust - cummulative_impulse / fuel_depletion_time)

    #above_max_fuel_grain_diameter_loss = 0.01 * abs2(p.u0_fuel_grain_void_diameter + p.additional_fuel_grain_void_diameter - p.fuel_grain_max_diameter)
    #we'll just enforce a bound on final_fuel_grain_void_diameter

    all_losses = ComponentVector(
        injector_velocity_loss = injector_velocity_loss,
        pressure_drop_ratio_loss = pressure_drop_ratio_loss,
        burn_time_loss = burn_time_loss,
        average_thrust_loss = average_thrust_loss
    )

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
properties_unitless = ustrip.(upreferred.(properties))

p_axes = getaxes(properties)
u_axes = getaxes(u0)

function viewable_system_design_loss(theta, p)
    theta, p, all_losses = system_design_loss(theta, p, theta_axes, u_axes, p_axes, model, theta_to_p_map, theta_to_u_map, p_to_u_map, append_optimized_parameters!, update_properties!)
    @show theta
    println("")
    @show p
    println("")
    @show all_losses
    println("")
    return sum(all_losses)
end

function system_design_loss_closure(theta, p)
    theta, p, all_losses = system_design_loss(theta, p, theta_axes, u_axes, p_axes, model, theta_to_p_map, theta_to_u_map, p_to_u_map, append_optimized_parameters!, update_properties!)

    return sum(all_losses)
end

f_closure = (du, u, p, t) -> system_ode!(du, u, p, t, model, p_axes)

ode_func = ODEFunction(f_closure)

t0 = 0.0
tMax = 1000.0
tspan = (t0, tMax)

du_test = deepcopy(u0_unitless)

implicit_prob = ODEProblem(ode_func, u0_unitless, tspan, properties_unitless)

append_optimized_parameters!(Vector(theta_guess_unitless), u0_unitless, properties_unitless, theta_to_u_map, theta_to_p_map, p_to_u_map)
update_u0!(du_test, u0_unitless, properties_unitless, 0.0, model)

sol = solve(implicit_prob, Tsit5(), callback = approximate_time_to_finish_cb)

opt_f = OptimizationFunction(system_design_loss_closure, Optimization.AutoForwardDiff())
opt_prob = OptimizationProblem(opt_f, Vector(theta_guess_unitless), Vector(properties_unitless), lb = Vector(theta_lb_unitless), ub = Vector(theta_ub_unitless))

cb = function (state, l)
    display(l)
    display(state.u)
    false
end

@time sol = solve(opt_prob, callback = cb, LBFGS(), reltol = 1e-4, maxiters = 1000)

viewable_system_design_loss(sol.u, properties_unitless)

final_properties = ComponentVector(sol.u, u_axes)