
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
    optimized_cb_set_local = optimized_cb_set
    isoutofdomain_feedsystem_local = isoutofdomain_feedsystem
    problem_template_local = implicit_prob
    
    (theta, p) -> trainsient_system_design_loss(
        theta, u_local, p, 
        theta_axes_local, u_axes_local, p_axes_local, 
        oxidizer_model_local, chamber_model_local, 
        theta_to_u_map_local, theta_to_p_map_local, p_to_u_map_local, 
        append_optimized_parameters_local!, update_u0_local!, 
        isoutofdomain_feedsystem_local, optimized_cb_set_local,
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
    optimized_cb_set_local = optimized_cb_set
    problem_template_local = implicit_prob

    (theta, p) -> trainsient_system_design_loss(
        theta, u_local, p, 
        theta_axes_local, u_axes_local, p_axes_local, 
        oxidizer_model_local, chamber_model_local, 
        theta_to_u_map_local, theta_to_p_map_local, p_to_u_map_local, 
        append_optimized_parameters_local!, update_u0_local!, 
        isoutofdomain_feedsystem_local, optimized_cb_set_local,
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
    optimized_cb_set_local = optimized_cb_set
    problem_template_local = implicit_prob

    function (theta, p)
        _, _, losses = trainsient_system_design_loss(
            theta, u_initial, p,
            θ_axes, state_axes, parameter_axes,
            om, cm,
            θ_to_u, θ_to_p, p_to_u,
            append_parameters!, initialize_state!,
            isoutofdomain_feedsystem_local, optimized_cb_set_local,
            problem_template_local
        )

        return sum(losses)
    end
end