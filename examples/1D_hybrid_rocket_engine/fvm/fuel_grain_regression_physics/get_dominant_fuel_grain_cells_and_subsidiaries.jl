
function get_dominant_fuel_grain_cells_and_subsidiaries(grid, geo)
    fuel_wall_patch_cells = unique(face.idx[1] for face in getfacetset(grid, "fuel_wall_patch")) #contains fluid domain and fuel grain cells

    dominant_fuel_grain_cells = []
    dominant_fuel_grain_closest_fluid_domain_cell_map = Dict{Int, Int}()
    dominant_fuel_grain_cell_centroids_z = []
    for cell_id in fuel_wall_patch_cells
        if cell_id in getcellset(grid, "fuel_grain_volume")
            push!(dominant_fuel_grain_cells, cell_id)
            push!(dominant_fuel_grain_cell_centroids_z, geo.cell_centroids[cell_id][3])
        end
    end

    for cell_id in fuel_wall_patch_cells
        if cell_id in getcellset(grid, "fluid_domain")
            cell_z = geo.cell_centroids[cell_id][3]
            idx = findfirst(x -> abs(x - cell_z) < 1e-10, dominant_fuel_grain_cell_centroids_z)
            if !isnothing(idx)
                dominant_fuel_grain_closest_fluid_domain_cell_map[dominant_fuel_grain_cells[idx]] = cell_id
            end
        end
    end

    fuel_grain_subsidiary_cell_ids_to_dominant_cell_id_map = Dict{Int, Int}() #fuel grain cells that are behind the cells that form a face with the fluid domain
    fuel_grain_dominant_cell_id_to_subsidiaries_map = Dict(
        cell_id => Int[] for cell_id in dominant_fuel_grain_cells
    ) #fuel grain cells that are behind the cells that form a face with the fluid domain

    for cell_id in getcellset(grid, "fuel_grain_volume")
        if cell_id in dominant_fuel_grain_cells
            continue
        end

        cell_z = geo.cell_centroids[cell_id][3]

        found_dominant_fuel_grain_index = findfirst((x) -> abs(x - cell_z) < 1e-10 , dominant_fuel_grain_cell_centroids_z)

        if !isnothing(found_dominant_fuel_grain_index)
            fuel_grain_subsidiary_cell_ids_to_dominant_cell_id_map[cell_id] = dominant_fuel_grain_cells[found_dominant_fuel_grain_index]
            dominant_cell_id = dominant_fuel_grain_cells[found_dominant_fuel_grain_index]
            push!(fuel_grain_dominant_cell_id_to_subsidiaries_map[dominant_cell_id], cell_id)
        else
            #this shouldn't happen
        end
    end
    
    fluid_domain_subsidiary_cell_ids_to_dominant_cell_id_map = Dict{Int, Int}() #fluid domain cells that are next to a dominant fuel grain cell
    fluid_domain_dominant_cell_id_to_subsidiaries_map = Dict(
        cell_id => Int[] for cell_id in dominant_fuel_grain_cells
    ) #fluid domain cells that are behind the cells that form a face with the fluid domain
    

    for cell_id in getcellset(grid, "fluid_domain")
        if cell_id in dominant_fuel_grain_cells
            continue
        end

        cell_z = geo.cell_centroids[cell_id][3]

        found_dominant_fuel_grain_index = findfirst((x) -> abs(x - cell_z) < 1e-10 , dominant_fuel_grain_cell_centroids_z)

        if !isnothing(found_dominant_fuel_grain_index)
            fluid_domain_subsidiary_cell_ids_to_dominant_cell_id_map[cell_id] = dominant_fuel_grain_cells[found_dominant_fuel_grain_index]
            dominant_cell_id = dominant_fuel_grain_cells[found_dominant_fuel_grain_index]
            push!(fluid_domain_dominant_cell_id_to_subsidiaries_map[dominant_cell_id], cell_id)
        else
            #this shouldn't happen because some fluid_domain cells are not next to a fuel grain
        end
    end

    return (
        dominant_fuel_grain_cells, 
        fuel_grain_subsidiary_cell_ids_to_dominant_cell_id_map,
        fuel_grain_dominant_cell_id_to_subsidiaries_map,
        fluid_domain_subsidiary_cell_ids_to_dominant_cell_id_map,
        fluid_domain_dominant_cell_id_to_subsidiaries_map
    )
end
