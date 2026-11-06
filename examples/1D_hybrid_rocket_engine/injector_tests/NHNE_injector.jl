using Clapeyron
using Unitful
using BenchmarkTools
using ComponentArrays
using Optim
using Plots

model = PR(["nitrous oxide"])

function get_spi_mass_flow(
    injector_discharge_coefficient, injector_area,
    inlet_pressure, chamber_pressure,
    nitrous_oxide_liquid_density, 
)
    spi_mass_flow = injector_discharge_coefficient * injector_area * sqrt(2 * nitrous_oxide_liquid_density * (inlet_pressure - chamber_pressure)) |> u"kg/s"

    return spi_mass_flow |> u"kg/s"
end

function get_hem_mass_flow(
    injector_discharge_coefficient, injector_area,
    inlet_pressure, chamber_pressure, 
    starting_entropy, pre_injector_enthalpy
)
    chamber_pressures = collect(range(inlet_pressure, chamber_pressure, 30))
    mass_fluxes = []

    for (i, chamber_pressure_test) in enumerate(chamber_pressures)
        post_injector_result = ps_flash(model, chamber_pressure_test, starting_entropy)
        post_injector_density = mass_density(model, post_injector_result) * u"kg/m^3"
        post_injector_enthalpy = mass_enthalpy(model, post_injector_result) * u"J/kg"

        mass_flux = post_injector_density * sqrt(2 * (pre_injector_enthalpy - post_injector_enthalpy)) |> u"kg/(m^2*s)"

        if i > 1
            if mass_flux < mass_fluxes[i-1]
                break
            end
        end

        push!(mass_fluxes, mass_flux) 
    end

    return injector_discharge_coefficient * injector_area * mass_fluxes[end] |> u"kg/s"

    #plot(ustrip.(upreferred.(range(chamber_pressure, inlet_pressure, 30)))[1:length(mass_fluxes)], ustrip.(upreferred.(mass_fluxes)))
end

function get_nhne_mass_flow(p)
    inlet_temperature = p.inlet_temperature
    inlet_pressure = p.inlet_pressure
    chamber_pressure = p.chamber_pressure
    injector_discharge_coefficient = p.injector_discharge_coefficient
    injector_area = p.injector_area

    inlet_pressure = saturation_pressure(model, inlet_temperature)[1]u"Pa"

    starting_entropy = mass_entropy(model, inlet_pressure, inlet_temperature)

    pre_injector_result = ps_flash(model, inlet_pressure, starting_entropy)
    nitrous_oxide_liquid_density = mass_density(model, pre_injector_result.data.p*u"Pa", pre_injector_result.data.T*u"K", phase=:liquid)
    pre_injector_enthalpy = mass_enthalpy(model, pre_injector_result) * u"J/kg"

    spi_mass_flow = get_spi_mass_flow(
        injector_discharge_coefficient, injector_area,
        inlet_pressure, chamber_pressure,
        nitrous_oxide_liquid_density,
    )

    hem_mass_flow = get_hem_mass_flow(
        injector_discharge_coefficient, injector_area, 
        inlet_pressure, chamber_pressure, 
        starting_entropy, pre_injector_enthalpy
    )

    non_eq_parameter = sqrt((inlet_pressure - chamber_pressure) / (saturation_pressure(model, inlet_temperature)[1]u"Pa" - chamber_pressure))

    nhne_mass_flow = (non_eq_parameter / (1 + non_eq_parameter)) * spi_mass_flow + (1/(1+non_eq_parameter)) * hem_mass_flow

    return spi_mass_flow, hem_mass_flow, nhne_mass_flow
end

p = ComponentVector(
    inlet_temperature = 21.0u"°C",
    inlet_pressure = 50.0u"bar",
    chamber_pressure = 20.0u"bar",
    injector_discharge_coefficient = 0.7,
    injector_area = 4.00u"mm^2"
)

function graph_varied_parameter(parameter_name::Symbol, parameter_vec, p)
    spi_mass_flows = typeof(1.0u"kg/s")[]
    hem_mass_flows = typeof(1.0u"kg/s")[]
    nhne_mass_flows = typeof(1.0u"kg/s")[]

    new_p = deepcopy(p)

    for parameter_val in parameter_vec
        setproperty!(new_p, parameter_name, parameter_val)

        spi_mass_flow, hem_mass_flow, nhne_mass_flow = get_nhne_mass_flow(new_p)

        push!(spi_mass_flows, spi_mass_flow)
        push!(hem_mass_flows, hem_mass_flow)
        push!(nhne_mass_flows, nhne_mass_flow)
    end

    return spi_mass_flows, hem_mass_flows, nhne_mass_flows
end


#parameter_sym = :chamber_pressure
#parameter_range = collect(range(1.0u"atm", 20.0u"bar", 50))

#parameter_sym = :inlet_pressure
#parameter_range = collect(range(20.0u"bar", 100.0u"bar", 50))

#parameter_sym = :injector_area
#parameter_range = collect(range(1.0u"mm^2", 10.0u"mm^2", 50))

parameter_sym = :inlet_temperature
parameter_range = collect(range((p.inlet_temperature-35.0u"K"), (p.inlet_temperature), 50))

saturation_pressure(model, p.inlet_temperature)[1]u"Pa" |> u"bar"

spi_mass_flows, hem_mass_flows, nhne_mass_flows = graph_varied_parameter(parameter_sym, parameter_range, p)

plot(parameter_range, spi_mass_flows, label = "SPI")
plot!(parameter_range, hem_mass_flows, label = "HEM")
plot!(parameter_range, nhne_mass_flows, label = "NHNE")

#=
post_injector_result = ps_flash(model, chamber_pressure, starting_entropy)
post_injector_density = mass_density(model, post_injector_result) * u"kg/m^3"
post_injector_enthalpy = mass_enthalpy(model, post_injector_result) * u"J/kg"

mass_flux = post_injector_density * sqrt(2 * (pre_injector_enthalpy - post_injector_enthalpy))


mass_flow_HEM = injector_discharge_coefficient * injector_area * post_injector_density * sqrt(2 * (pre_injector_enthalpy - post_injector_enthalpy)) |> u"kg/s"
=#