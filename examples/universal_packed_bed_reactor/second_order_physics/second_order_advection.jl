function first_order_upwind(du, u, idx_a, idx_b, face_idx, var_a, var_b)
    if ustrip(du.mass_face[idx_a, face_idx]) < 0.0
        return var_a
    else
        return var_b
    end
end

function second_order_upwind(
    du, u, 
    idx_a, idx_b, face_idx, 
    u_a, u_b,
    variable_gradients,
    u_max, u_min,
    cell_face_normal_a, 
    cell_volume_a, cell_volume_b,
    limiter_function
)
    if ustrip(du.mass_face[idx_a, face_idx]) < 0.0
        phi = limiter_function(
            u_a, u_b,
            variable_gradients,
            u_max, u_min,
            cell_face_normal_a,
            cell_volume_a,
            idx_a
        )
        return u_a + phi * dot(view(variable_gradients, idx_a, :), cell_face_normal_a)
    else
        phi = limiter_function(
            u_a, u_b,
            variable_gradients,
            u_max, u_min,
            cell_face_normal_a,
            cell_volume_b,
            idx_b
        )
        return u_b + phi * dot(view(variable_gradients, idx_b, :), cell_face_normal_a)
    end
end

function second_order_all_species_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist, 
    cell_face_normal_a, 
    cell_volume_a, cell_volume_b,
    limiter_function
)
    upwinded_methanol_mass_fraction = second_order_upwind(
        du, u, 
        idx_a, idx_b, face_idx, 
        u.mass_fractions.methanol[idx_a], u.mass_fractions.methanol[idx_b],
        u.mass_fractions_grad.methanol,
        u.mass_fractions_max.methanol[idx_a], u.mass_fractions_min.methanol[idx_a],
        cell_face_normal_a,
        cell_volume_a, cell_volume_b,
        limiter_function
    )

    du.species_masses.methanol[idx_a] += (du.mass_face[idx_a, face_idx] * upwinded_methanol_mass_fraction)

    upwinded_water_mass_fraction = second_order_upwind(
        du, u, 
        idx_a, idx_b, face_idx, 
        u.mass_fractions.water[idx_a], u.mass_fractions.water[idx_b],
        u.mass_fractions_grad.water,
        u.mass_fractions_max.water[idx_a], u.mass_fractions_min.water[idx_a],
        cell_face_normal_a,
        cell_volume_a, cell_volume_b,
        limiter_function
    )

    du.species_masses.water[idx_a] += (du.mass_face[idx_a, face_idx] * upwinded_water_mass_fraction)

    upwinded_carbon_monoxide_mass_fraction = second_order_upwind(
        du, u, 
        idx_a, idx_b, face_idx, 
        u.mass_fractions.carbon_monoxide[idx_a], u.mass_fractions.carbon_monoxide[idx_b],
        u.mass_fractions_grad.carbon_monoxide,
        u.mass_fractions_max.carbon_monoxide[idx_a], u.mass_fractions_min.carbon_monoxide[idx_a],
        cell_face_normal_a,
        cell_volume_a, cell_volume_b,
        limiter_function
    )

    du.species_masses.carbon_monoxide[idx_a] += (du.mass_face[idx_a, face_idx] * upwinded_carbon_monoxide_mass_fraction)

    upwinded_hydrogen_mass_fraction = second_order_upwind(
        du, u, 
        idx_a, idx_b, face_idx, 
        u.mass_fractions.hydrogen[idx_a], u.mass_fractions.hydrogen[idx_b],
        u.mass_fractions_grad.hydrogen,
        u.mass_fractions_max.hydrogen[idx_a], u.mass_fractions_min.hydrogen[idx_a],
        cell_face_normal_a,
        cell_volume_a, cell_volume_b,
        limiter_function
    )

    du.species_masses.hydrogen[idx_a] += (du.mass_face[idx_a, face_idx] * upwinded_hydrogen_mass_fraction)

    upwinded_carbon_dioxide_mass_fraction = second_order_upwind(
        du, u, 
        idx_a, idx_b, face_idx, 
        u.mass_fractions.carbon_dioxide[idx_a], u.mass_fractions.carbon_dioxide[idx_b],
        u.mass_fractions_grad.carbon_dioxide,
        u.mass_fractions_max.carbon_dioxide[idx_a], u.mass_fractions_min.carbon_dioxide[idx_a],
        cell_face_normal_a,
        cell_volume_a, cell_volume_b,
        limiter_function
    )

    du.species_masses.carbon_dioxide[idx_a] += (du.mass_face[idx_a, face_idx] * upwinded_carbon_dioxide_mass_fraction)

    upwinded_air_mass_fraction = second_order_upwind(
        du, u, 
        idx_a, idx_b, face_idx, 
        u.mass_fractions.air[idx_a], u.mass_fractions.air[idx_b],
        u.mass_fractions_grad.air,
        u.mass_fractions_max.air[idx_a], u.mass_fractions_min.air[idx_a],
        cell_face_normal_a,
        cell_volume_a, cell_volume_b,
        limiter_function
    )

    du.species_masses.air[idx_a] += (du.mass_face[idx_a, face_idx] * upwinded_air_mass_fraction)
end

function second_order_enthalpy_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist, 
    cell_face_normal_a, 
    cell_volume_a, cell_volume_b,
    limiter_function
)
    #=cp_upwinded = second_order_upwind(
        du, u, idx_a, idx_b, face_idx, 
        u.cp[idx_a], u.cp[idx_b],
        u.grad_cp[idx_a], u.grad_cp[idx_b],
        cell_face_normal_a, cell_face_normal_b
    )=#

    cp_upwinded = first_order_upwind(
        du, u, idx_a, idx_b, face_idx, 
        u.cp[idx_a], u.cp[idx_b]
    )
    
    temp_upwinded = second_order_upwind(
        du, u, 
        idx_a, idx_b, face_idx, 
        u.temp[idx_a], u.temp[idx_b],
        u.temp_grad,
        u.temp_max[idx_a], u.temp_min[idx_a],
        cell_face_normal_a, 
        cell_volume_a, cell_volume_b,
        limiter_function
    )

    energy_flux = du.mass_face[idx_a, face_idx] * cp_upwinded * temp_upwinded
    
    du.heat[idx_a] += energy_flux
end