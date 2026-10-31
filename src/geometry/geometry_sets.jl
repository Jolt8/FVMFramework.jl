"""
    build_named_facetset(cell_face_pairs, cell_neighbors)

Turn a list of `(cell_id, face_idx)` into the `cell_neighbors` layout:
`[(idx_a, [(idx_b, face_idx_a, face_idx_b), ...]), ...]`, sorted by cell id.
The neighbour of each face is looked up in `cell_neighbors`; a face with no
neighbour (a boundary face) becomes `(0, face_idx_a, 0)`.
"""
function build_named_facetset(cell_face_pairs, cell_neighbors)
    neighbor_of_face = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()
    for (idx_a, neighbor_list) in cell_neighbors
        for (idx_b, face_idx_a, face_idx_b) in neighbor_list
            neighbor_of_face[(idx_a, face_idx_a)] = (idx_b, face_idx_b)
        end
    end

    faces_of_cell = OrderedDict{Int, Vector{NTuple{3, Int}}}()
    for (idx_a, face_idx_a) in sort(collect(cell_face_pairs))
        idx_b, face_idx_b = get(neighbor_of_face, (idx_a, face_idx_a), (0, 0))
        push!(get!(faces_of_cell, idx_a, NTuple{3, Int}[]), (idx_b, face_idx_a, face_idx_b))
    end

    return [(idx_a, face_list) for (idx_a, face_list) in faces_of_cell]
end

"""
    get_cellset(geo, name)

Cell ids of a named cell set (a region).
"""
function get_cellset(geo, name::AbstractString)
    if !haskey(geo.cellsets, name)
        error("no cell set named \"$name\". Available: $(join(keys(geo.cellsets), ", "))")
    end
    return geo.cellsets[name]
end

"""
    get_facetset(geo, name)

A named face set (e.g. a boundary patch), in the `cell_neighbors` layout.
A name of the form "<region_a>__<region_b>" that is not stored (meshes do
not store region pairs) is derived with `get_facetset(geo, region_a, region_b)`.
"""
function get_facetset(geo, name::AbstractString)
    if haskey(geo.facetsets, name)
        return geo.facetsets[name]
    end

    if occursin("__", name)
        region_a, region_b = split(name, "__"; limit = 2)
        if haskey(geo.cellsets, region_a) && haskey(geo.cellsets, region_b)
            return get_facetset(geo, String(region_a), String(region_b))
        end
    end

    error("no face set named \"$name\". Available: $(join(keys(geo.facetsets), ", "))")
end

"""
    get_facetset(geo, region_a, region_b)

The faces between two regions, listed from `region_a`'s side: every entry is
`(cell in region_a, [(cell in region_b, face_idx_a, face_idx_b), ...])`.
Each shared face appears once, so looping over the result applies a flux
once per face. Swapping the arguments gives the same faces seen from the
other side.

Where the faces come from:
  - a stored set "<region_a>__<region_b>" (ROM exports store these),
  - a stored "<region_b>__<region_a>", flipped to region_a's side,
  - otherwise derived from `geo.cellsets` and `geo.cell_neighbors` (meshes).
For region_a == region_b this is the faces inside the region, each listed
once from the lower cell id (or the stored "<region>__<region>" set).
"""
function get_facetset(geo, region_a::AbstractString, region_b::AbstractString)
    forward_name = "$(region_a)__$(region_b)"
    reverse_name = "$(region_b)__$(region_a)"

    if haskey(geo.facetsets, forward_name)
        return geo.facetsets[forward_name]
    end

    if region_a != region_b && haskey(geo.facetsets, reverse_name)
        # Same faces, seen from the other cell.
        flipped_faces_of_cell = OrderedDict{Int, Vector{NTuple{3, Int}}}()
        for (idx_b, neighbor_list) in geo.facetsets[reverse_name]
            for (idx_a, face_idx_b, face_idx_a) in neighbor_list
                push!(get!(flipped_faces_of_cell, idx_a, NTuple{3, Int}[]), (idx_b, face_idx_a, face_idx_b))
            end
        end
        sorted_cell_ids = sort(collect(keys(flipped_faces_of_cell)))
        return [(idx_a, flipped_faces_of_cell[idx_a]) for idx_a in sorted_cell_ids]
    end

    # Derive from the cell sets: faces of region_a cells whose neighbour is in region_b.
    cells_in_region_b = Set(get_cellset(geo, region_b))
    region_pair_faces = Tuple{Int, Vector{NTuple{3, Int}}}[]
    for idx_a in get_cellset(geo, region_a)
        faces_to_region_b = NTuple{3, Int}[]
        for (idx_b, face_idx_a, face_idx_b) in geo.cell_neighbors[idx_a][2]
            if !(idx_b in cells_in_region_b)
                continue
            end
            # Inside one region every face is seen from both cells; keep one.
            if region_a == region_b && idx_b < idx_a
                continue
            end
            push!(faces_to_region_b, (idx_b, face_idx_a, face_idx_b))
        end
        if !isempty(faces_to_region_b)
            push!(region_pair_faces, (idx_a, faces_to_region_b))
        end
    end
    return region_pair_faces
end
