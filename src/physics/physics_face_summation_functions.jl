"""
    sum_mass_flux_face_to_cell!(du, u, cell_id, vol)
    
    Sum mass flux from faces to cell.
    
    # Arguments
    - `du::ComponentArray`: Residual for mass
    - `u::ComponentArray`: Solution
    - `cell_id::Int`: Cell ID
    - `vol::ComponentArray`: Cell volumes
    
    # Returns
    - nothing
"""
function sum_mass_flux_face_to_cell!(du, u, cell_id, vol)
    du.mass[cell_id] += sum(view(du.mass_face, cell_id, :))
end