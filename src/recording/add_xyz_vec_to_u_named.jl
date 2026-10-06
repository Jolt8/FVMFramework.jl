

"""
    add_xyz_vec_to_u_named!(u_named, field_name::Symbol, field_x::Symbol, field_y::Symbol, field_z::Symbol)

Adds a new field to the u_named object that is a vector of x, y, and z components.

# Arguments
- `u_named`: Named vector of state vectors (created by using regenerate_fvm_state)
- `field_name`: Name of the new field where all three components will be stored
- `field_x`: Name of the x component
- `field_y`: Name of the y component
- `field_z`: Name of the z component

# Example
```julia
add_xyz_vec_to_u_named!(u_named, :velocity, :vel_u, :vel_v, :vel_w)
```
"""
function add_xyz_vec_to_u_named!(u_named, field_name::Symbol, field_x::Symbol, field_y::Symbol, field_z::Symbol)
    for i in eachindex(u_named)
        field_vec = SVector{3, Float64}[
            SVector{3, Float64}(
                #why must we do this cursedness just to get it to look right?
                view(u_named[i], field_z)[c],
                -view(u_named[i], field_y)[c], 
                view(u_named[i], field_x)[c], 
            ) for c in 1:n_cells
        ]

        field_matrix = reduce(hcat, field_vec)
        field_matrix = transpose(field_matrix)
        
        u_named[i] = merge_properties(u_named[i], ComponentVector(field_name = field_matrix,))
    end
end