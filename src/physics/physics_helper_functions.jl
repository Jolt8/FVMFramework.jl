
"""
    upwind(du, u, idx_a, )

Performs upwinding for a given variable at a face based on the direction of du.mass_face.

Parameters:
    - du: State derivatives vector
    - u: State vector
    - idx_a: Index of the upstream cell
    - idx_b: Index of the downstream cell
    - face_idx: Index of the face
    - var_a: Value of the variable in the upstream cell
    - var_b: Value of the variable in the downstream cell
    
Returns:
    - Upwind value of the variable
"""
function upwind(
    du, u,
    idx_a, idx_b, face_idx,
    var_a, var_b
)
    if ustrip(du.mass_face[idx_a, face_idx]) < 0.0
        return var_a
    else
        return var_b
    end
end

#NOTE to future developers, don't try to use a tanh based upwind, it's much slower and offers no real benefits
#=
function upwind(
    du, u,
    idx_a, idx_b, face_idx,
    var_a, var_b
)
    #wow, this is a lot slower than the basic upwind function, I thought it would allow the solver to take larger steps, guess not
    k = 10
    return ((var_a - var_b) / 2) * tanh(k * du.mass_face[idx_a, face_idx]) + ((var_a + var_b) / 2)
end
=#

"""
    harmonic_mean(val_a, val_b)
    
    Calculates the harmonic mean of two values.
    
    # Arguments
    - `val_a::Number`: The first value
    - `val_b::Number`: The second value
    
    # Returns
    - `Number`: The harmonic mean of `val_a` and `val_b`
"""
function harmonic_mean(val_a, val_b)
    return 2 * val_a * val_b / (val_a + val_b)
end


"""
    mw_avg!(u, cell_id)
    
    Calculates the average molecular weight of the gas in a cell and writes to `u.mw_avg[cell_id]`.
    
    # Arguments
    - `u::ComponentArray`: State vector
    - `cell_id::Int`: Cell ID
    
    # Returns
    - `nothing`: Modifies `u` in-place
"""
function mw_avg!(u, cell_id)
    inv_mw = 0.0

    for_fields!(u.mass_fractions, u.molecular_weights) do species, mass_fractions, molecular_weights
        inv_mw += mass_fractions[species[cell_id]] / molecular_weights[species[cell_id]]
    end

    u.mw_avg[cell_id] = 1.0 / inv_mw
end

"""
    rho_ideal!(u, cell_id)
    
    Calculates the density of the gas in a cell using the ideal gas law and writes to `u.rho[cell_id]`.
    
    # Arguments
    - `u::ComponentArray`: State vector
    - `cell_id::Int`: Cell ID
    
    # Returns
    - `nothing`: Modifies `u` in-place
"""
function rho_ideal!(u, cell_id)
    u.rho[cell_id] = (u.pressure[cell_id] * u.mw_avg[cell_id]) / (u.R_gas[cell_id] * u.temp[cell_id])
end

"""
    molar_concentrations!(u, cell_id)
    
    Calculates the molar concentrations of each species in a cell and writes to `u.molar_concentrations[species[cell_id]]`.
    
    # Arguments
    - `u::ComponentArray`: State vector
    - `cell_id::Int`: Cell ID
    
    # Returns
    - `nothing`: Modifies `u` in-place
"""
function molar_concentrations!(u, cell_id)
    for_fields!(u.mass_fractions, u.molecular_weights, u.molar_concentrations) do species, mass_fractions, molecular_weights, molar_concentrations
        molar_concentrations[species[cell_id]] = (u.rho[cell_id] * mass_fractions[species[cell_id]]) / molecular_weights[species[cell_id]]
    end
end

"""
    molar_fractions!(u, cell_id)
    
    Calculates the molar fractions of each species in a cell and writes to `u.molar_fractions[species[cell_id]]`.
    
    # Arguments
    - `u::ComponentArray`: State vector
    - `cell_id::Int`: Cell ID
    
    # Returns
    - `nothing`: Modifies `u` in-place
"""
function molar_fractions!(u, cell_id)
    for_fields!(u.mass_fractions, u.molecular_weights, u.molar_fractions) do species, mass_fractions, molecular_weights, molar_fractions
        molar_fractions[species[cell_id]] = (mass_fractions[species[cell_id]] / molecular_weights[species[cell_id]]) * u.mw_avg[cell_id]
    end
end

"""
    cp_avg!(u, cell_id)
    
    Calculates the average specific heat capacity of the gas in a cell and writes to `u.cp_avg[cell_id]`.
    
    # Arguments
    - `u::ComponentArray`: State vector
    - `cell_id::Int`: Cell ID
    
    # Returns
    - `nothing`: Modifies `u` in-place
"""
function cp_avg!(u, cell_id)
    u.cp_avg[cell_id] *= 0.0
    for_fields!(u.mass_fractions, u.species_cps) do species, mass_fractions, species_cps
        u.cp_avg[species[cell_id]] += mass_fractions[species[cell_id]] * species_cps[species[cell_id]]
    end
end

"""
    van_t_hoff(A, dH, T, R)
    
    Calculates the equilibrium constant using the Van't Hoff equation.
    
    # Arguments
    - `A::Number`: Pre-exponential factor
    - `dH::Number`: Heat of reaction
    - `T::Number`: Temperature
    - `R::Number`: Gas constant
    
    # Returns
    - `Number`: Equilibrium constant
    
    # Notes
    - This function is used for calculating the equilibrium constant
    - It is a helper function for the reaction flux calculation
    
    # Example
    ```julia
    van_t_hoff(1, 2, 3, 4)
    ```
"""


function K_gibbs_free(u, cell_id, reaction)
    K_ref = exp(-reaction.ref_delta_G[cell_id] / (u.R_gas[cell_id] * reaction.ref_temp[cell_id]))

    ln_K_ratio = (-reaction.heat_of_reaction[cell_id] / u.R_gas[cell_id]) * (1 / u.temp[cell_id] - 1 / reaction.ref_temp[cell_id])

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

#probably unecessary, but whatever
function arrenhius_k(A, Ea, T, R)
    return A * exp(-Ea / (R * T))
end
