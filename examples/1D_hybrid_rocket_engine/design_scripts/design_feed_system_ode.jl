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
using Plots

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
    
    volumetric_flow = p.adjusted_valve_flow_capacity_factor * sqrt((p.tank_pressure - p.mid_section_pressure) / p.tank_density)
    #we might need a different equation that activates when the tank starts spitting out gas

    if true == false
        @show p.adjusted_valve_flow_capacity_factor
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

    #@show "injector contribution"
    #@show du.chamber_gas_mass

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
    #@show "regression contribution"
    #@show du.chamber_gas_mass

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
    #@show "nozzle contribution"
    #@show du.chamber_gas_mass

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
    u.port_diameter = p.final_fuel_grain_void_diameter - p.additional_fuel_grain_void_diameter
    #p.u0_fuel_grain_void_diameter

    return nothing
end

Revise.includet(joinpath(@__DIR__, "feed_system_update_state_for_ode.jl"))

function system_ode!(du_vec, u_vec, p_vec, t, oxidizer_model, chamber_model, p_axes, u_axes)
    du = ComponentVector(du_vec, u_axes)
    u = ComponentVector(u_vec, u_axes)
    p = ComponentVector(eltype(u).(p_vec), p_axes)

    if t <= 0.1 #reset the p.burned_out tracker at the start of every simulation or else it stays at 1.0 for the next solve
        p.burned_out = 0.0
    end

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
    
    fuel_mass_flow = 0.0

    #if the fuel hasn't completely burned up yet, still allow it to combust
    #otherwise the port_diameter could go lower than the final_fuel_grain_void_diameter
    #if p.burned_out != 1.0
    if u.port_diameter < p.final_fuel_grain_void_diameter
        fuel_mass_flow = regression_rate!(du, u, p, t, oxidizer_mass_flow)
    end
    
    p.oxidizer_to_fuel_ratio = oxidizer_mass_flow / max(fuel_mass_flow, 1e-9)
    p.propellant_characteristic_velocity = cstar_interpolator_Pa(p.chamber_pressure, p.oxidizer_to_fuel_ratio)

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
    simulation_time = 30.0u"s",

    #Overall rocket properties
    oxidizer_to_fuel_ratio = 0.0,
    desired_oxidizer_to_fuel_ratio = 7.6,
    propellant_isp = 0.0u"s", #220.0u"s", This is now based on a table in CEA_lookup_table.jl
    propellant_characteristic_velocity = 0.0u"m/s", #1500.0u"m/s", This is now based on a table in CEA_lookup_table.jl
    burn_time = 0.0u"s", #for viewing the optimized burn time 
    desired_burn_time = 5.0u"s",
    average_thrust = 0.0u"N", #for viewing the optimized averge thrust
    desired_average_thrust = 600.0u"N",
    cummulative_impulse = 0.0u"N*s", #for viewing the optimized cumulative impulse
    desired_impulse = 0.0u"N*s", #this will be calculated based on desired_average thrust and desired_burn_time later
    #we would likely want a better correlation in the future that will return the specific impulse for a given oxidizer to fuel ratio and exit pressure
    gravity = 9.81u"m/s^2",
    target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio = 2.0, #not an optimized parameter, but the ideal choice is hard to know
    fuel_grain_max_diameter = (12.0u"inch" |> u"cm"),
    burned_out = 0.0, #this switches to 1 with a callback whenever the fuel has burned out
    burn_out_time = 0.0, #this gets set to the t in which 
    
    #Tank
    u0_tank_oxidizer_mass = 1.134u"kg", #optimized, should actually be 2.05kg to get the impulse we need, but I want to see if the optimizer will produce 2.05kg for debugging reasons
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
    valve_flow_capacity_factor = 10.0u"mm^2", #a result of the final valve we choose #we should optimize this to get a good ideal of what we want
    adjusted_valve_flow_capacity_factor = 0.0u"mm^2",
    valve_opening = 1.0,

    #Mid section
    mid_section_pressure = 1.0u"atm", #this determines the initial kg of nitrous oxide in the mid_section
    mid_section_temperature = 21.0u"°C",
    mid_section_volume = 100.0u"cm^3", #Changed from 10u"cm" to 100u"cm" because I wanted to decrease solver stiffness for debugging purposes
    mid_section_density = 0.0u"kg/m^3",
    mid_section_vapor_fraction = 0.0,
    mid_section_specific_enthalpy = 0.0u"J/kg",

    #Injector properties
    injector_discharge_coefficient = 0.7,
    #injector_number_of_orifices = 20, #optimized
    #injector_orifice_diameter = 1.5u"mm", #optimized
    injector_orifice_area = 12.0u"mm^2", #optimized
    target_injector_velocity = 50.0u"m/s",

    #Chamber
    u0_chamber_nitrous_oxide_mass_fraction = 1.0,
    u0_chamber_hdpe_mass_fraction = 0.0,
    chamber_pressure = 1.0u"atm", #this determines the initial kg of nitrous oxide in the chamber
    chamber_temperature = 21.0u"°C",
    chamber_volume = 212.0u"cm^3",
    chamber_density = 0.0u"kg/m^3",
    wall_heat_loss = 0.0u"W",
    chamber_vapor_fraction = 0.0,
    chamber_specific_enthalpy = 0.0u"J/kg",

    #Fuel grain
    #u0_fuel_grain_void_diameter = 3.0u"cm",
    fuel_mass = 0.0u"kg", #this will be derived by substracting the volume of the cylinder formed by the u0_fuel_grain_void_diameter by the final_fuel_grain_void_diameter and then multiplying by the fuel density
    fuel_density = 950.0u"kg/m^3",
    #fuel_regression_rate = 0.5u"mm/s",
    additional_fuel_grain_void_diameter = 1.0u"cm", #optimized
    final_fuel_grain_void_diameter = 23.8u"mm", 
    #not optimized, we're going to be using a COTS phenolic liner, so we're just going to use the ID that the manufactuerer specifies
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
    #OptimizedParameter(:u0_tank_oxidizer_mass, 1.0u"kg", 1.5u"kg",),
    #since we're going to be using a 2.5lbs/1.134kg tank COTS tank, we're no longer going to optimize it and treat it as fixed
    #OptimizedParameter(:tank_pressure, 50.0u"bar", 71.0u"bar",),

    #adjustable valve
    OptimizedParameter(:valve_flow_capacity_factor, 5.0u"mm^2", 20.0u"mm^2"),

    #Mid section

    #Injector properties
    OptimizedParameter(:injector_orifice_area, 7.0u"mm^2", 25.0u"mm^2"),

    #Chamber

    #Fuel grain
    #OptimizedParameter(:u0_fuel_grain_void_diameter, 1.0u"cm", 5.0u"cm"),
    #OptimizedParameter(:final_fuel_grain_void_diameter, 1.0u"cm", 20.0u"cm"),
    OptimizedParameter(:additional_fuel_grain_void_diameter, 0.5u"cm", 2.0u"cm"),
    OptimizedParameter(:fuel_grain_length, 20.0u"cm", 100.0u"cm"),

    #Fuel Grain Empirical Parameters

    #Propellant properties

    #Nozzle
    OptimizedParameter(:nozzle_throat_diameter, 10.0u"mm", 20.0u"mm")
]

Revise.includet(joinpath(@__DIR__, "CEA_lookup_table.jl"))
Revise.includet(joinpath(@__DIR__, "design_feed_sytem_ode_loss.jl"))

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

ode_func = ODEFunction(f_closure)#, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 8.0
tspan = (t0, tMax)

append_optimized_parameters!(Vector(theta_guess_unitless), u0_unitless, properties_unitless, theta_to_u_map, theta_to_p_map, p_to_u_map)
update_u0!(u0_unitless, properties_unitless, 0.0, oxidizer_model, chamber_model)

implicit_prob = ODEProblem(ode_func, Vector(u0_unitless), tspan, Vector(properties_unitless))

system_ode!(du_test, Vector(u0_unitless), Vector(properties_unitless), 0.0, oxidizer_model, chamber_model, p_axes, u_axes)

f_closure(du_test, u0_unitless, properties_unitless, 0.0)

Revise.includet(joinpath(@__DIR__, "internals/out_of_domain_funcs.jl"))
Revise.includet(joinpath(@__DIR__, "internals/callbacks.jl"))

#we're no longer stopping the simulation right after the fuel is burnt, this misses some of the residual thrust 
cb_set = CallbackSet(
    approximate_time_to_finish_cb,
    #port_diameter_termination_cb,
    chamber_pressure_termination_cb,
    #fuel_burnout_cb
);

sol = solve(
    implicit_prob,
    #FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, ),
    #FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true),
    #AutoTsit5(FBDF(linsolve = SparspakFactorization(), autodiff = AutoFiniteDiff())),
    callback = cb_set,
    isoutofdomain = isoutofdomain_feedsystem
)

function plot_sol_states(sol, u_axes, p_axes, oxidizer_model, chamber_model)
    u_named_vec = []
    p_named_vec = []
    thrust_vec = []

    for i in eachindex(sol.t)
        du_named = ComponentVector(similar(sol.u[i]), u_axes)
        u_named = ComponentVector(deepcopy(sol.u[i]), u_axes)
        p_named = ComponentVector(deepcopy(sol.prob.p), p_axes)
        
        update_state!(du_named, u_named, p_named, 0.0, oxidizer_model, chamber_model)

        @show p_named.chamber_pressure

        push!(u_named_vec, u_named)
        push!(p_named_vec, p_named)

        if i > 1
            u_named_prev = ComponentVector(sol.u[i-1], u_axes)
            dt = sol.t[i] - sol.t[i-1]

            #Cummulative impulse loss calcs
            oxidizer_used = u_named_prev.tank_oxidizer_mass - u_named.tank_oxidizer_mass
            fuel_used = p_named.fuel_density * (π / 4) * (u_named_prev.port_diameter^2 - u_named.port_diameter^2) * p_named.fuel_grain_length
            
            oxidizer_mass_flow = oxidizer_used / dt
            fuel_mass_flow = fuel_used / dt

            propellant_used = oxidizer_used + fuel_used
            propellant_mass_flow = oxidizer_mass_flow + fuel_mass_flow

            oxidizer_to_fuel_ratio = oxidizer_mass_flow / max(fuel_mass_flow, 1e-9)

            p_named.propellant_isp = isp_interpolator_Pa(p_named.chamber_pressure, oxidizer_to_fuel_ratio)

            p_named.propellant_characteristic_velocity = cstar_interpolator_Pa(p_named.chamber_pressure, oxidizer_to_fuel_ratio)

            chamber_gas_mass_flow_out = p_named.nozzle_discharge_coefficient * ((p_named.chamber_pressure * p_named.nozzle_throat_area) / p_named.propellant_characteristic_velocity)

            thrust_produced = p_named.propellant_isp * p_named.gravity * chamber_gas_mass_flow_out
        else
            chamber_gas_mass_flow_out = 0.0
            thrust_produced = 0.0
        end

        push!(thrust_vec, thrust_produced)
    end

    port_diameter_plot = plot(sol.t, [u_named_vec[i].port_diameter for i in eachindex(sol.t)], title = "port diameter")
    plot!(port_diameter_plot, sol.t, fill(p_named_vec[1].final_fuel_grain_void_diameter, length(sol.t)))
    display(port_diameter_plot)

    tank_oxidizer_mass_plot = plot(sol.t, [u_named_vec[i].tank_oxidizer_mass for i in eachindex(sol.t)], title = "tank oxidizer mass")
    display(tank_oxidizer_mass_plot)

    tank_pressure_plot = plot(sol.t, [p_named_vec[i].tank_pressure for i in eachindex(sol.t)], title = "tank pressure")
    display(tank_pressure_plot)

    mid_section_mass_plot = plot(sol.t, [u_named_vec[i].mid_section_mass for i in eachindex(sol.t)], title = "mid section mass")
    display(mid_section_mass_plot)

    mid_section_pressure_plot = plot(sol.t, [p_named_vec[i].mid_section_pressure for i in eachindex(sol.t)], title = "mid section pressure")
    display(mid_section_pressure_plot)

    chamber_mass_plot = plot(sol.t, [u_named_vec[i].chamber_gas_mass for i in eachindex(sol.t)], title = "chamber mass")
    display(chamber_mass_plot)

    chamber_pressure_plot = plot(sol.t, [p_named_vec[i].chamber_pressure for i in eachindex(sol.t)], title = "chamber pressure")
    display(chamber_pressure_plot)

    chamber_temperature_plot = plot(sol.t, [p_named_vec[i].chamber_temperature for i in eachindex(sol.t)], title = "chamber temperature")
    display(chamber_temperature_plot)

    thrust_plot = plot(sol.t, [thrust_vec[i] for i in eachindex(sol.t)], title = "thrust")
    display(thrust_plot)

    @show u_named_vec[end-2].chamber_gas_mass
    @show u_named_vec[end-1].chamber_gas_mass
    @show u_named_vec[end].chamber_gas_mass

    @show p_named_vec[end].chamber_temperature
end

plot_sol_states(sol, u_axes, p_axes, oxidizer_model, chamber_model)

#plot(sol.t, [ComponentVector(sol.u[i], u_axes).port_diameter for i in eachindex(sol.t)])

optimized_cb_set = CallbackSet(
    approximate_time_to_finish_cb,
    #port_diameter_termination_cb,
    chamber_pressure_termination_cb,
    #fuel_burnout_cb
);

Revise.includet(joinpath(@__DIR__, "internals/loss_closures.jl"))

viewable_system_design_loss(
    ComponentVector(
        #u0_tank_oxidizer_mass = 1.1306761278408962,
        valve_flow_capacity_factor = 8.088661048535489e-6,
        injector_orifice_area = 0.9e-5,
        additional_fuel_grain_void_diameter = 0.01,
        fuel_grain_length = 0.48,
        nozzle_throat_diameter = 0.021
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        #u0_tank_oxidizer_mass = 1.1306761278408962,
        valve_flow_capacity_factor = 8.088661048535489e-6,
        injector_orifice_area = 8.135451620577921e-6,
        additional_fuel_grain_void_diameter = 0.01,
        fuel_grain_length = 0.48,
        nozzle_throat_diameter = 0.021
    ), properties_unitless
)
#=
viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 1.0106761278408962,
        valve_flow_capacity_factor = 8.088661048535489e-6,
        injector_orifice_area = 8.135451620577921e-6,
        additional_fuel_grain_void_diameter = 0.009824607761510886,
        fuel_grain_length = 0.5669914705724927, 
        nozzle_throat_diameter = 0.015820023622496553
    ), properties_unitless
)
#=
viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 1.23,
        valve_flow_capacity_factor = 1.5e-5,
        injector_orifice_area = 2.3e-5,
        #u0_fuel_grain_void_diameter = 0.020,
        additional_fuel_grain_void_diameter = 0.01,
        fuel_grain_length = 0.40, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.0113
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 1.33,
        valve_flow_capacity_factor = 1.5e-5,
        injector_orifice_area = 2.3e-5,
        #u0_fuel_grain_void_diameter = 0.016,
        additional_fuel_grain_void_diameter = 0.0093,
        fuel_grain_length = 0.45, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.0113
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 1.3338466690919596,
        valve_flow_capacity_factor = 1.4864625820186499e-5,
        injector_orifice_area = 2.307584272651544e-5,
        #u0_fuel_grain_void_diameter = 0.015622001356552163,
        additional_fuel_grain_void_diameter = 0.009261793472995626,
        fuel_grain_length = 0.4516855289918622, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.011327930260216029
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 1.3338466690919596,
        valve_flow_capacity_factor = 1.4864625820186499e-5,
        injector_orifice_area = 2.307584272651544e-5,
        #u0_fuel_grain_void_diameter = 0.015622001356552163,
        additional_fuel_grain_void_diameter = 0.009061793472995626,
        fuel_grain_length = 0.6616855289918622, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.011127930260216029
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 2.034594773122675,
        valve_flow_capacity_factor = 6.6906403746981306e-6,
        injector_orifice_area = 2.0728250081324405e-5,
        #u0_fuel_grain_void_diameter = 0.03403845868495702,
        additional_fuel_grain_void_diameter = 0.009822880727716907,
        fuel_grain_length = 0.6997446134380847, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.014563969006102736
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 2.034594773122675,
        valve_flow_capacity_factor = 6.6906403746981306e-6,
        injector_orifice_area = 2.0728250081324405e-5,
        #u0_fuel_grain_void_diameter = 0.03403845868495702,
        additional_fuel_grain_void_diameter = 0.009822880727716907,
        fuel_grain_length = 0.6997446134380847, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.014563969006102736
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 2.06,
        valve_flow_capacity_factor = 1.0e-4,
        injector_orifice_area = 6.0e-4,
        #u0_fuel_grain_void_diameter = 0.01,
        additional_fuel_grain_void_diameter = 0.0132,
        fuel_grain_length = 0.6, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.0233
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 2.06,
        valve_flow_capacity_factor = 3.0e-6,
        injector_orifice_area = 7.0e-6,
        #u0_fuel_grain_void_diameter = 0.01,
        additional_fuel_grain_void_diameter = 0.013,
        fuel_grain_length = 0.6, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.015
    ), properties_unitless
)

viewable_system_design_loss(
    ComponentVector(
        u0_tank_oxidizer_mass = 5.0,
        valve_flow_capacity_factor = 1.0e-4,
        injector_orifice_area = 1.0e-4,
        #u0_fuel_grain_void_diameter = 0.03,
        additional_fuel_grain_void_diameter = 0.009,
        fuel_grain_length = 0.30, 
        #hmm, increasing the fuel grain length doesn't seem to change the burn time or impulse at all which shouldn't happen
        nozzle_throat_diameter = 0.0156
    ), properties_unitless
)

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

#opt_sol = solve(opt_prob, LBFGS(), callback = cb, reltol = 1e-4)

opt_sol = solve(opt_prob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(),
    callback = cb,
    PoulationSize = 1000,
    #maxiters = 1,
    #maxtime = 60.0,
    Method = :RandomSearcher,
    verbose = true
    #Method = :SepReal
)
#=
new_opt_f = OptimizationFunction(pure_system_design_loss_closure, Optimization.AutoFiniteDiff())
new_opt_prob = OptimizationProblem(new_opt_f, Vector(opt_sol.u), Vector(properties_unitless), lb = Vector(theta_lb_unitless), ub = Vector(theta_ub_unitless))

opt_sol = solve(new_opt_prob, LBFGS(), callback = cb)

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

min(successful_parameters)

best_parameters_idx = argmin(losses)
best_parameters = successful_parameters[best_parameters_idx]

viewable_system_design_loss(best_parameters, properties_unitless)
viewable_system_design_loss(theta_guess_unitless, properties_unitless)

test = (theta_guess_unitless, theta_axes)

#Note: for some reason increasing the additional_fuel-grain_void_diameter by just 0.005 causes the optimizer to fail due to the solver getting way too stiff
#No idea why. 
#It also seems that as long as the difference between the additional fuel graind void diameter and the initial fuel grain void diameter is kept smaller than 0.01, the solver works
#determining whether or not this is a stiffness problem that's actually physical and some error in the solver is going to be hell

#Also, I noticed that increasing the initial port diameter too much leads to the solver failin again due to too much stiffness

#Interestingly, decreasing the valve flow capacity factor or the injector orifice area 
#to 1.0e-6 while also increasing the additional_fuel_grain_void_diameter allows the solver to solve instantly

#Also, increasing the initial oxidizer mass to 5.0 makes the solver significantly less stiff

#So I think this is just an issue of choosing values that are physically meaningful and do not result in too much stiffness

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
