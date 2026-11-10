using Clapeyron
using Unitful
using BenchmarkTools
using ComponentArrays
using Optim
using Plots

model = PR(["nitrous oxide"])

inlet_temperature = 21.0u"°C"
inlet_pressure = 50.0u"bar"
chamber_pressure = 20.0u"bar"

starting_entropy = mass_entropy(model, inlet_pressure, inlet_temperature)

@time ps_flash(model, 6.0u"bar", starting_entropy)

saturation_pressure(model, inlet_temperature)

injector_discharge_coefficient = 0.7

injector_area = 4.00u"mm^2"

pre_injector_result = ps_flash(model, inlet_pressure, starting_entropy)
nitrous_oxide_liquid_density = mass_density(model, pre_injector_result.data.p*u"Pa", pre_injector_result.data.T*u"K", phase=:liquid)
pre_injector_enthalpy = mass_enthalpy(model, pre_injector_result) * u"J/kg"

spi_mass_flow = injector_discharge_coefficient * injector_area * sqrt(2 * nitrous_oxide_liquid_density * (inlet_pressure - chamber_pressure)) |> u"kg/s"

function get_hem_mass_flux()
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

    return mass_fluxes[end]

    #plot(ustrip.(upreferred.(range(chamber_pressure, inlet_pressure, 30)))[1:length(mass_fluxes)], ustrip.(upreferred.(mass_fluxes)))
end

hem_mass_flow = injector_discharge_coefficient * injector_area * get_hem_mass_flux() |> u"kg/s"

non_eq_parameter = sqrt((inlet_pressure - chamber_pressure) / (saturation_pressure(model, inlet_temperature)[1]u"Pa" - chamber_pressure))

nhne_mass_flow = (non_eq_parameter / (1 + non_eq_parameter)) * spi_mass_flow + (1/(1+non_eq_parameter)) * hem_mass_flow


#=
post_injector_result = ps_flash(model, chamber_pressure, starting_entropy)
post_injector_density = mass_density(model, post_injector_result) * u"kg/m^3"
post_injector_enthalpy = mass_enthalpy(model, post_injector_result) * u"J/kg"

mass_flux = post_injector_density * sqrt(2 * (pre_injector_enthalpy - post_injector_enthalpy))


mass_flow_HEM = injector_discharge_coefficient * injector_area * post_injector_density * sqrt(2 * (pre_injector_enthalpy - post_injector_enthalpy)) |> u"kg/s"
=#