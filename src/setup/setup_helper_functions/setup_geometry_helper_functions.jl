function get_cellset_volume(config, set_name::String)
    total_volume = 0.0
    for cell_id in get_cellset(config.geo, set_name)
        total_volume += config.geo.cell_volumes[cell_id]
    end

    return total_volume
end

function get_facetset_area(config, set_name::String)
    #a cell can have several faces in the set (e.g. a stepped contact), so every listed face is added
    total_area = 0.0

    for (cell_id, neighbor_list) in get_facetset(config.geo, set_name)
        for (neighbor_id, facet_idx, neighbor_facet_idx) in neighbor_list
            face_area = config.geo.cell_face_areas[cell_id, facet_idx]
            total_area += face_area
        end
    end

    return total_area
end

#this style of function is not needed for cell volumes because it's very easy to access cell_volumes with just [cell_id]
function get_facetset_area_per_cell(config, set_name::String)
    n_cells = config.geo.n_cells
    cell_ids_respective_areas = zeros(Float64, n_cells)

    for (cell_id, neighbor_list) in get_facetset(config.geo, set_name)
        for (neighbor_id, facet_idx, neighbor_facet_idx) in neighbor_list
            face_area = config.geo.cell_face_areas[cell_id, facet_idx]
            cell_ids_respective_areas[cell_id] += face_area #we += because a cell might have 2 faces in this facet_set
        end
    end

    return cell_ids_respective_areas
end
