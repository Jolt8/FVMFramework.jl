"""
    patch_faces_from_geo(geo, name, owner_region)

The faces of face set `name` (a boundary patch, or "<region_a>__<region_b>"),
in the `cell_neighbors` layout, keeping only cells of `owner_region` unless
it is "" (all cells).
"""
function patch_faces_from_geo(geo, name, owner_region)
    facetset = get_facetset(geo, name)

    if owner_region == ""
        owner_region_cell_ids = nothing
    else
        owner_region_cell_ids = Set(get_cellset(geo, owner_region))
    end

    cell_neighbors = Tuple{Int, Vector{Tuple{Int, Int, Int}}}[]
    for (idx_a, neighbor_list) in facetset
        if !isnothing(owner_region_cell_ids) && !(idx_a in owner_region_cell_ids)
            continue
        end
        push!(cell_neighbors, (idx_a, copy(neighbor_list)))
    end
    return cell_neighbors
end

function add_patch!(
    config, name;
    owner_region = nothing,
    properties = ComponentVector(),
    patch_function,
)
    owner_region_name = isnothing(owner_region) ? "" : owner_region
    cell_neighbors = patch_faces_from_geo(config.geo, name, owner_region_name)

    patch = PatchSetupInfo(name, properties, patch_function, cell_neighbors, owner_region_name)

    all_patch_names = [patch.name for patch in config.patches]

    if name in all_patch_names
        existing_patch_idx = findfirst(x -> x == name, all_patch_names)

        config.patches[existing_patch_idx] = patch
    else
        push!(config.patches, patch)
    end
    return
end

function add_patch!(
    config;
    owner = nothing,
    connection = nothing,
    properties = ComponentVector(),
    patch_function,
)
    cell_neighbors = get_facetset(config.geo, owner, connection)
    patch_name = "$(owner)__$(connection)"

    patch = PatchSetupInfo(patch_name, properties, patch_function, cell_neighbors, owner)

    all_patch_names = [patch.name for patch in config.patches]

    if patch_name in all_patch_names
        existing_patch_idx = findfirst(x -> x == patch_name, all_patch_names)

        config.patches[existing_patch_idx] = patch
    else
        push!(config.patches, patch)
    end
    return
end

#rebuild a patch's faces from geo (keeps its properties, function and owner region)
function update_patch!(config, name)
    all_patch_names = [patch.name for patch in config.patches]
    existing_patch_idx = findfirst(x -> x == name, all_patch_names)
    if isnothing(existing_patch_idx)
        error("no patch named \"$name\"")
    end

    patch = config.patches[existing_patch_idx]
    patch.cell_neighbors = patch_faces_from_geo(config.geo, name, patch.owner_region)
    return
end
