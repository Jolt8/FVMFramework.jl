#this could probably also be handled by dynamic dispatch for facets, but it helps the user know a different routine is happening
#this needs to be updated, it currently doesn't work
function get_neighboring_cell_and_face_idx_from_face_idx(cell_id, face_idx, top) #returning the neighbor_face_idx is not really necessary
    neighbor_info = top.face_face_neighbor[cell_id, face_idx]

    if !isempty(neighbor_info)
        face_idx_neighbor = first(neighbor_info)
        neighbor_id = face_idx_neighbor.idx[1]
        neighbor_face_idx = face_idx_neighbor.idx[2]
    else
        return nothing, nothing
    end

    return neighbor_id, neighbor_face_idx
end

function add_patch!(
    config, name;
    properties,
    patch_function,
)
    cell_ids_and_face_idxs = [cell_id_facet_idx.idx for cell_id_facet_idx in keys(config.grid.facetsets[name].dict)]
    #[(cell_id, face_idx), (cell_id, face_idx), ...]

    n_cells = length(config.grid.cells)
    
    cell_neighbors = [(cell_id, Vector{Tuple{Int, Int, Int}}()) for cell_id in 1:n_cells]

    for (cell_id, face_idx_a) in cell_ids_and_face_idxs
        neighbor_id, neighbor_face_idx = get_neighboring_cell_and_face_idx_from_face_idx(cell_id, face_idx_a, config.top)
        if !isnothing(neighbor_id)
            push!(cell_neighbors[cell_id][2], (neighbor_id, face_idx_a, face_idx_b))
            #push!(cell_neighbors[neighbor_id][2], (cell_id, neighbor_face_idx)) #This is not necessary
        else
            push!(cell_neighbors[cell_id][2], (0, face_idx_a, 0))
        end
    end

    
    filter!(conn -> !(isempty(conn[2])), cell_neighbors)
    #get rid of empty connections

    patch = PatchSetupInfo(name, properties, patch_function, cell_neighbors)

    all_patch_names = [patch.name for patch in config.patches]

    if name in all_patch_names
        existing_patch_idx = findfirst(x -> x == name, all_patch_names)

        config.patches[existing_patch_idx] = patch
    else   
        push!(config.patches, patch)
    end
    return 
end

#this is for updating the state vector (u0) when doing different experimental trials
function update_patch!(config, name) 
    all_patch_names = [patch.name for patch in config.patches]
    existing_patch_idx = findfirst(x -> x == name, all_patch_names)
    
    patch = config.patches[existing_patch_idx]

    cell_ids_and_face_idxs = [cell_id_facet_idx.idx for cell_id_facet_idx in keys(config.grid.facetsets[name].dict)]
    #[(cell_id, face_idx), (cell_id, face_idx), ...]

    n_cells = length(config.grid.cells)
    
    cell_neighbors = [(cell_id, Vector{Tuple{Int, Int, Int}}()) for cell_id in 1:n_cells]

    for (cell_id, face_idx_a) in cell_ids_and_face_idxs
        neighbor_id, neighbor_face_idx = get_neighboring_cell_and_face_idx_from_face_idx(cell_id, face_idx_a, config.top)
        if !isnothing(neighbor_id)
            push!(cell_neighbors[cell_id][2], (neighbor_id, face_idx_a, face_idx_b))
            #push!(cell_neighbors[neighbor_id][2], (cell_id, neighbor_face_idx)) #This is not necessary
        else
            push!(cell_neighbors[cell_id][2], (0, face_idx_a, 0))
        end
    end

    config.regions[existing_region_idx] = region
end