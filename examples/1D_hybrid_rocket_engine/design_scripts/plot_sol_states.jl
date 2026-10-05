
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

            p_named.propellant_isp = isp_interpolator_Pa(p_named.chamber_pressure, oxidizer_to_fuel_ratio, p_named.expansion_ratio)

            p_named.propellant_characteristic_velocity = cstar_interpolator_Pa(p_named.chamber_pressure, oxidizer_to_fuel_ratio, p_named.expansion_ratio)

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

    mid_section_mass_plot = plot(sol.t, [u_named_vec[i].mid_section_mass for i in eachindex(sol.t)], title = "mid section mass")
    display(mid_section_mass_plot)

    mid_section_pressure_plot = plot(sol.t, [p_named_vec[i].mid_section_pressure for i in eachindex(sol.t)], title = "mid section pressure")
    display(mid_section_pressure_plot)

    chamber_mass_plot = plot(sol.t, [u_named_vec[i].chamber_gas_mass for i in eachindex(sol.t)], title = "chamber mass")
    display(chamber_mass_plot)

    pressure_plot = plot(sol.t, [p_named_vec[i].tank_pressure for i in eachindex(sol.t)], title = "System Pressures", label = "tank pressure", xlabel = "time (s)", ylabel = "pressure (Pa)")
    plot!(pressure_plot, sol.t, [p_named_vec[i].mid_section_pressure for i in eachindex(sol.t)], label = "mid section pressure")
    plot!(pressure_plot, sol.t, [p_named_vec[i].chamber_pressure for i in eachindex(sol.t)], label = "chamber pressure")
    display(pressure_plot)

    chamber_temperature_plot = plot(sol.t, [p_named_vec[i].chamber_temperature for i in eachindex(sol.t)], title = "chamber temperature", xlabel = "time (s)", ylabel = "temperature (K)")
    display(chamber_temperature_plot)

    thrust_plot = plot(sol.t, [thrust_vec[i] for i in eachindex(sol.t)], title = "thrust")
    display(thrust_plot)

    @show u_named_vec[end-2].chamber_gas_mass
    @show u_named_vec[end-1].chamber_gas_mass
    @show u_named_vec[end].chamber_gas_mass

    @show p_named_vec[end].chamber_temperature
end
