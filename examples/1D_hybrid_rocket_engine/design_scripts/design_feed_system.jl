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

Revise.includet(joinpath(@__DIR__, "../common_data/overall_rocket_data.jl"))

overall_rocket_data = get_overall_rocket_data()

tank_pressure = overall_rocket_data.tank_initial_pressure
desired_chamber_pressure = 30u"bar"

overshoot_factor = 1.2

post_adjustable_valve_minimum_pressure = desired_chamber_pressure * overshoot_factor

#NOTE: chamber pressure is determined by propellant mass flow, nozzle throat area, and characteristic velocity of the propellants

model = PR(["nitrous oxide"])

function adjustable_valve_flow(p, model)
    nitrous_oxide_density = mass_density(model, p.tank_pressure, p.tank_temperature, [1.0])

    if p.mid_section_pressure > p.tank_pressure
        return 1e10, 1e10
    end
    
    volumetric_flow = p.valve_flow_capacity_factor * p.valve_opening * sqrt((p.tank_pressure - p.mid_section_pressure) / nitrous_oxide_density)

    h_initial = enthalpy(model, p.tank_pressure, p.tank_temperature, [1.0])

    flash_result = ph_flash(model, p.mid_section_pressure, h_initial, [1.0])

    p.mid_section_temperature = flash_result.data.T

    mid_section_density = mass_density(model, p.mid_section_pressure, p.mid_section_temperature, [1.0])

    mass_flow = volumetric_flow * mid_section_density

    return mass_flow, (p.tank_pressure - p.mid_section_pressure)
end

function injector_valve_flow(p, model)
    pre_injector_density = mass_density(model, p.mid_section_pressure, p.mid_section_temperature, [1.0])

    if p.chamber_pressure > p.mid_section_pressure
        return 1e10, 1e10, 1e10
    end
    
    volumetric_flow = p.injector_discharge_coefficient * p.injector_orifice_area * sqrt((2 * (p.mid_section_pressure - p.chamber_pressure)) / pre_injector_density)

    h_initial = enthalpy(model, p.mid_section_pressure, p.mid_section_temperature, [1.0])

    flash_result = ph_flash(model, p.chamber_pressure, h_initial, [1.0]).data

    chamber_density = mass_density(model, p.chamber_pressure, flash_result.T, [1.0])

    mass_flow = volumetric_flow * chamber_density

    return mass_flow, (p.mid_section_pressure - p.chamber_pressure), pre_injector_density
end

function regression_rate(p, model)
    oxidizer_mass_flux = p.oxidizer_mass_flow / p.fuel_grain_average_cross_sectional_area

    oxidizer_mass_flux_g_per_cm2_s = 0.1 * oxidizer_mass_flux #convert to kg/(m^2*s)

    regression_mm_per_s = p.fuel_regression_coeff_a * oxidizer_mass_flux_g_per_cm2_s^p.fuel_regression_coeff_n

    return regression_mm_per_s * 0.001 #convert to m/s
end

function fuel_flow_with_regression_model(p, model)
    return p.fuel_density * p.fuel_grain_burning_surface_area * p.fuel_regression_rate
end

function oxidizer_to_fuel_ratio(p, model)
    return p.oxidizer_mass_flow / p.fuel_mass_flow
end

function get_chamber_pressure(p, model)
    return (p.oxidizer_mass_flow + p.fuel_mass_flow) * p.propellant_characteristic_velocity / (p.nozzle_throat_area * p.nozzle_discharge_coefficient)
end

properties = ComponentVector(
    #Overall rocket properties
    oxidizer_mass_flow = 0.39u"kg/s",
    fuel_mass_flow = 0.051u"kg/s",
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
    tank_kg_of_nitrous_oxide = 5.0u"kg", #optimized
    required_tank_volume = 0.0u"m^3", #we will try to minimize this
    
    #Nitrous oxide parameters
    molecular_weight = 44.013u"g/mol",

    #adjustable valve
    valve_flow_capacity_factor = 1e-5, #a result of the final valve we choose #we should optimize this to get a good ideal of what we want
    valve_opening = 1.0,

    #Mid section
    mid_section_pressure = 40.0u"bar",
    mid_section_temperature = 21.0u"°C",

    #Injector properties
    injector_discharge_coefficient = 0.7,
    #injector_number_of_orifices = 20, #optimized
    #injector_orifice_diameter = 1.5u"mm", #optimized
    injector_orifice_area = 0.1u"cm^2", #optimized
    target_injector_velocity = 50.0u"m/s",

    #Chamber
    chamber_pressure = 30.0u"bar",

    #Fuel grain
    fuel_mass = 0.0u"kg", #this will be derived by substracting the volume of the cylinder formed by the initial_fuel_grain_void_diameter by the final_fuel_grain_void_diameter and then multiplying by the fuel density
    fuel_density = 950.0u"kg/m^3",
    fuel_regression_rate = 0.5u"mm/s",
    initial_fuel_grain_void_diameter = 3.0u"cm", #optimized
    additional_fuel_grain_void_diameter = 5.0u"cm", #optimied
    final_fuel_grain_void_diameter = 0.0u"cm", #determined by initial + additional
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
    OptimizedParameter(
        :oxidizer_mass_flow, 
        1e-6u"kg/s", 
        properties.oxidizer_mass_flow * 10.0
    ),
    OptimizedParameter(
        :fuel_mass_flow, 
        1e-6u"kg/s", 
        properties.fuel_mass_flow * 10.0
    ),
    #=
    OptimizedParameter(
        :oxidizer_to_fuel_ratio, 
        4.0, 
        8.5
    ),
    =#

    #Tank
    OptimizedParameter(
        :tank_pressure, 
        51.0u"bar", 
        72.0u"bar"
    ),
    OptimizedParameter(
        :tank_kg_of_nitrous_oxide, 
        1e-6u"kg", 
        50.0u"kg",
    ),

    #Nitrous oxide parameters

    #adjustable valve
    OptimizedParameter(
        :valve_flow_capacity_factor, 
        1e-7u"m^2", 
        1e-2u"m^2"
    ),

    #Mid section
    OptimizedParameter(
        :mid_section_pressure, 
        11.0u"bar", 
        69.0u"bar"
    ),

    #Injector properties
    #=
    OptimizedParameter(
        :injector_orifice_diameter, 
        1.0u"mm", 
        10.0u"mm"
    ),
    OptimizedParameter(
        :injector_number_of_orifices, 
        1, 
        50
    ),
    =#
    OptimizedParameter(
        :injector_orifice_area, 
        0.01u"cm^2", 
        10.0u"cm^2"
    ),


    #Chamber
    OptimizedParameter(
        :chamber_pressure, 
        10.0u"bar", 
        69.0u"bar"
    ),

    #Fuel grain
    OptimizedParameter(
        :fuel_regression_rate, 
        0.0u"mm/s", 
        5.0u"mm/s"
    ),
    OptimizedParameter(
        :initial_fuel_grain_void_diameter, 
        1.0u"cm", 
        5.0u"cm"
    ),
    OptimizedParameter(
        :additional_fuel_grain_void_diameter, 
        0.1u"cm", 
        50.0u"cm"
    ),
    OptimizedParameter(
        :fuel_grain_length, 
        5.0u"cm",
        300.0u"cm"
    ),

    #Fuel Grain Empirical Parameters

    #Propellant properties

    #Nozzle
    OptimizedParameter(
        :nozzle_throat_diameter,
        0.2u"cm",
        2.5u"cm",
    )
]

function update_properties!(p, model)    
    #Tank
    nitrous_oxide_density = mass_density(model, p.tank_pressure, p.tank_temperature, [1.0])
    p.required_tank_volume = p.tank_kg_of_nitrous_oxide / nitrous_oxide_density
    
    #Fuel grain
    p.final_fuel_grain_void_diameter = p.initial_fuel_grain_void_diameter + p.additional_fuel_grain_void_diameter

    initial_void_radius = p.initial_fuel_grain_void_diameter / 2
    final_void_radius = p.final_fuel_grain_void_diameter / 2

    p.fuel_grain_average_cross_sectional_area = (pi / 3) * (initial_void_radius^2 + initial_void_radius*final_void_radius + final_void_radius^2)

    p.fuel_grain_burning_surface_area = pi * (initial_void_radius + final_void_radius) * sqrt(p.fuel_grain_length^2 + (final_void_radius - initial_void_radius)^2)

    p.fuel_mass = p.fuel_density * (pi * (final_void_radius^2 - initial_void_radius^2) * p.fuel_grain_length)

    #Injector
    #p.injector_orifice_area = p.injector_number_of_orifices * (pi / 4) * (p.injector_orifice_diameter^2)

    #Nozzle
    p.nozzle_throat_area = (pi / 4) * (p.nozzle_throat_diameter^2)

    #Overall rocket properties
    
    #again, we need something here that determines propellant_isp based on oxidizer_to_fuel_ratio, chamber_pressure, and exit pressure (which we don't really know yet)
end

function system_design_loss(u, p, u_axes, p_axes, model, u_to_p_map, append_optimized_parameters!, update_properties!)
    u = ComponentVector(u, u_axes)
    
    # Promote p to the type of u (which will be Dual during ForwardDiff) so it can accept Duals
    p = ComponentVector(eltype(u).(p), p_axes)

    append_optimized_parameters!(u, p, u_to_p_map)
    update_properties!(p, model)

    if any(iszero, p)
        @error "It seems like p is missing a value"
    end

    adjustable_valve_oxidizer_mass_flow, adjustable_valve_pressure_drop = adjustable_valve_flow(p, model)
    adjustable_valve_oxidizer_mass_flow_loss = abs2(p.oxidizer_mass_flow - adjustable_valve_oxidizer_mass_flow)

    injector_valve_oxidizer_mass_flow, injector_valve_pressure_drop, pre_injector_density = injector_valve_flow(p, model)
    injector_valve_oxidizer_mass_flow_loss = abs2(p.oxidizer_mass_flow - injector_valve_oxidizer_mass_flow)

    regression_rate_val = regression_rate(p, model)
    fuel_regression_rate_loss = abs2(p.fuel_regression_rate - regression_rate_val)

    fuel_mass_flow_val = fuel_flow_with_regression_model(p, model)
    fuel_mass_flow_loss = abs2(p.fuel_mass_flow - fuel_mass_flow_val)

    oxidizer_to_fuel_ratio_val = oxidizer_to_fuel_ratio(p, model)
    p.oxidizer_to_fuel_ratio = oxidizer_to_fuel_ratio_val
    oxidizer_to_fuel_ratio_loss = abs2(p.oxidizer_to_fuel_ratio - oxidizer_to_fuel_ratio_val)
    #TODO: find a way to actually estimate Isp based on the oxidizer_to_fuel_ratio

    chamber_pressure_val = get_chamber_pressure(p, model)
    chamber_pressure_loss = abs2(p.chamber_pressure - chamber_pressure_val)

    #Overall rocket performance
    total_mass_flow = p.fuel_mass_flow + p.oxidizer_mass_flow
    
    p.average_thrust = p.propellant_isp * p.gravity * total_mass_flow
    
    p.burn_time = min(p.tank_kg_of_nitrous_oxide / p.oxidizer_mass_flow, p.fuel_mass / p.fuel_mass_flow)
    #end of overall rocket_performance_calculations

    required_tank_volume_loss = 0.01 * p.required_tank_volume

    average_thrust_loss = 0.01 * abs2(p.desired_average_thrust - p.average_thrust)

    burn_time_loss = 0.01 * abs2(p.desired_burn_time - p.burn_time)

    if true == false
        println("pressure drop ratio loss")
        @show p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio
        @show adjustable_valve_pressure_drop
        @show injector_valve_pressure_drop
        @show injector_valve_pressure_drop / adjustable_valve_pressure_drop
        println("")
    end

    pressure_drop_ratio_loss = 0.00001 * abs2(p.target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio - (injector_valve_pressure_drop / adjustable_valve_pressure_drop))

    if true == false
        println("-----")
        println("Injector Velocity Loss:")
        @show p.target_injector_velocity
        @show injector_valve_oxidizer_mass_flow
        @show pre_injector_density
        @show p.injector_orifice_area
        @show injector_valve_oxidizer_mass_flow / (pre_injector_density * p.injector_orifice_area)
        println("")
    end
    

    injector_velocity_loss = 0.00001 * abs2(p.target_injector_velocity - injector_valve_oxidizer_mass_flow / (pre_injector_density * p.injector_orifice_area))

    #=
    if p.initial_fuel_grain_void_diameter + p.additional_fuel_grain_void_diameter >= p.fuel_grain_max_diameter
        above_max_fuel_grain_diameter_loss = 1000.0
    else
        above_max_fuel_grain_diameter_loss = 0.0
    end
    =#

    above_max_fuel_grain_diameter_loss = 0.01 * abs2(p.initial_fuel_grain_void_diameter + p.additional_fuel_grain_void_diameter - p.fuel_grain_max_diameter)



    #above_max_fuel_grain_diameter_loss = 0.0

    #OBSERVATION: it seems like we're going to have to weigh the importance of injector velocity against adjustable valve authority

    #=
    @show p.target_injector_velocity
    @show injector_valve_oxidizer_mass_flow
    @show pre_injector_density
    @show p.injector_orifice_area
    @show injector_valve_oxidizer_mass_flow / (pre_injector_density * p.injector_orifice_area)
    =#

    all_losses = ComponentVector(
        adjustable_valve_oxidizer_mass_flow_loss = adjustable_valve_oxidizer_mass_flow_loss,
        injector_valve_oxidizer_mass_flow_loss = injector_valve_oxidizer_mass_flow_loss,
        fuel_regression_rate_loss = fuel_regression_rate_loss,
        fuel_mass_flow_loss = fuel_mass_flow_loss,
        oxidizer_to_fuel_ratio_loss = oxidizer_to_fuel_ratio_loss,
        chamber_pressure_loss = chamber_pressure_loss,
        required_tank_volume_loss = required_tank_volume_loss,
        average_thrust_loss = average_thrust_loss,
        burn_time_loss = burn_time_loss,
        pressure_drop_ratio_loss = pressure_drop_ratio_loss,
        injector_velocity_loss = injector_velocity_loss,
        above_max_fuel_grain_diameter_loss = above_max_fuel_grain_diameter_loss
    )

    #@show all_losses

    return u, p, all_losses

    #return sum(all_losses)
end

function create_u_guess(optimized_properties, properties)
    u_guess_dict = OrderedDict{Symbol, Number}()
    u_lb_dict = OrderedDict{Symbol, Number}()
    u_ub_dict = OrderedDict{Symbol, Number}()

    for optimized_property in optimized_properties
        u_guess_dict[optimized_property.name] = ustrip(upreferred(getproperty(properties, optimized_property.name)))
        u_lb_dict[optimized_property.name] = ustrip(upreferred(optimized_property.lb))
        u_ub_dict[optimized_property.name] = ustrip(upreferred(optimized_property.ub))
    end

    u_guess = ComponentVector(u_guess_dict)
    u_lb = ComponentVector(u_lb_dict)
    u_ub = ComponentVector(u_ub_dict)

    u_axes = getaxes(u_guess)

    u_to_p_map = Tuple{Int, Int}[]

    for name in propertynames(u_guess)
        push!(u_to_p_map, (label2index(properties, string(name))[1], label2index(u_guess, string(name))[1]))
    end

    return u_guess, u_lb, u_ub, u_axes, u_to_p_map
end

function append_optimized_parameters!(u, p, u_to_p_map)
    for (p_idx, u_idx) in u_to_p_map
        p[p_idx] = u[u_idx]
    end
end

u_guess, u_lb, u_ub, u_axes, u_to_p_map = create_u_guess(optimized_properties, properties)

#convert all to base SI and then strip away units
u_guess_unitless = ustrip.(upreferred.(u_guess))
u_lb_unitless = ustrip.(upreferred.(u_lb))
u_ub_unitless = ustrip.(upreferred.(u_ub))
properties_unitless = ustrip.(upreferred.(properties))

p_axes = getaxes(properties)

function viewable_system_design_loss(u, p)
    u, p, all_losses = system_design_loss(u, p, u_axes, p_axes, model, u_to_p_map, append_optimized_parameters!, update_properties!)
    @show u
    println("")
    @show p
    println("")
    @show all_losses
    println("")
    return sum(all_losses)
end

function system_design_loss_closure(u, p)
    u, p, all_losses = system_design_loss(u, p, u_axes, p_axes, model, u_to_p_map, append_optimized_parameters!, update_properties!)

    return sum(all_losses)
end

opt_f = OptimizationFunction(system_design_loss_closure, Optimization.AutoForwardDiff())
opt_prob = OptimizationProblem(opt_f, Vector(u_guess_unitless), Vector(properties_unitless), lb = Vector(u_lb_unitless), ub = Vector(u_ub_unitless))

cb = function (state, l)
    display(l)
    display(state.u)
    false
end

@time sol = solve(opt_prob, callback = cb, LBFGS(), reltol = 1e-4, maxiters = 1000)

viewable_system_design_loss(sol.u, properties_unitless)

final_properties = ComponentVector(sol.u, u_axes)