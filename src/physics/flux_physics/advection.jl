
"""
    all_specise_advection!(du, u, idx_a, idx_b, face_idx, area, norm, dist, vol_a, vol_b)

This function computes the advection of all species for a single face given an upwinded mass flux across the face.

Parameters:
    - du: State derivatives vector 
    - u: State vector 
    - idx_a: Index of 'a' cell 
    - idx_b: Index of 'b' cell 
    - face_idx: Index of face 
    - area: Area of face 
    - norm: Normal vector of face 
    - dist: Distance from face to center of 'a' cell 
    - vol_a: Volume of 'a' cell 
    - vol_b: Volume of 'b' cell 

Returns:
    None
"""
function all_species_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    for_fields!(u.mass_fractions, du.species_masses) do species, mass_fractions, species_masses
        upwinded_mass_fraction = upwind(du, u, idx_a, idx_b, face_idx, mass_fractions[species[idx_a]], mass_fractions[species[idx_b]])
        species_masses[species[idx_a]] += (du.mass_face[idx_a, face_idx] * upwinded_mass_fraction)
    end
end

"""
    enthalpy_advection!(du, u, idx_a, idx_b, face_idx, area, norm, dist, vol_a, vol_b)

This function computes the advection of enthalpy for a single face given an upwinded mass flux, specific heat capacity, and temperature across the face.
Note that this function assumes a constant specific heat capacity, but if you want to change this just copy this and get the actual heat capacities with Clapeyron.jl

Parameters:
    - du: State derivatives vector 
    - u: State vector 
    - idx_a: Index of 'a' cell 
    - idx_b: Index of 'b' cell 
    - face_idx: Index of face 
    - area: Area of face 
    - norm: Normal vector of face 
    - dist: Distance from face to center of 'a' cell 
    - vol_a: Volume of 'a' cell 
    - vol_b: Volume of 'b' cell 

Returns:
    None
"""
function enthalpy_advection!(
    du, u,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    cp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.fluid_cp[idx_a], u.fluid_cp[idx_b]) #TODO: change back to u.cp later

    temp_upwinded = upwind(du, u, idx_a, idx_b, face_idx, u.temp[idx_a], u.temp[idx_b])

    energy_flux = du.mass_face[idx_a, face_idx] * cp_upwinded * temp_upwinded
    
    du.heat[idx_a] += energy_flux
end


