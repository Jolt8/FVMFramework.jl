#=
Loading a reduced-order-model geometry exported by the ROMGeometry FreeCAD
workbench (format_version 1, see freecad-rom-geometry/docs/export_format.md).

A ROM is a very coarse unstructured mesh whose cells have unusual shapes. The
export already holds every volume, area, normal and distance, so the
FVMGeometry is filled straight from the file: no Ferrite grid, no nodes.

    rom_data = read_rom_geometry_json(path)
    add_virtual_cell!(rom_data; name = "recirculation_zone", volume = ..., faces = [...])
    config = create_fvm_config_from_rom(rom_data, u_proto)

u_proto needs rom_cell_count(rom_data) entries per field (count the virtual
cells too, so build u_proto after adding them).
=#

const ROM_GEOMETRY_FORMAT_VERSION = 1

const ROMCoordinate = SVector{3, Float64}

rom_coordinate(value::Nothing) = zero(ROMCoordinate)
rom_coordinate(value) = ROMCoordinate(Float64(value[1]), Float64(value[2]), Float64(value[3]))

"""
    read_rom_geometry_json(path)

Read a ROMGeometry export. Returns the parsed file (an ordered dictionary, so
region and face-set order are kept) for `add_virtual_cell!` and
`create_fvm_config_from_rom`.
"""
function read_rom_geometry_json(path)
    rom_data = JSON.parsefile(path; dicttype = OrderedDict{String, Any})

    format_version = get(rom_data, "format_version", nothing)
    if format_version != ROM_GEOMETRY_FORMAT_VERSION
        error("$path has format_version $(format_version); this FVMFramework reads format_version $(ROM_GEOMETRY_FORMAT_VERSION)")
    end

    for (cell_position, cell) in enumerate(rom_data["cells"])
        if cell["id"] != cell_position
            error("$path: cell number $cell_position has id $(cell["id"]); ids must be 1, 2, 3, ... in order")
        end
    end

    return rom_data
end

"""
    rom_cell_count(rom_data)

Number of cells, including virtual cells added so far.
"""
rom_cell_count(rom_data) = length(rom_data["cells"])

"""
    rom_region_cell_id(rom_data, region_name, axial_index)

Cell id of one cell of a region: `axial_index` is an integer (1 = lowest z),
`:first` or `:last`.
"""
function rom_region_cell_id(rom_data, region_name, axial_index)
    if !haskey(rom_data["cellsets"], region_name)
        error("no region named \"$region_name\". Available: $(join(keys(rom_data["cellsets"]), ", "))")
    end
    region_cell_ids = rom_data["cellsets"][region_name]

    if axial_index == :first
        return Int(first(region_cell_ids))
    elseif axial_index == :last
        return Int(last(region_cell_ids))
    end

    if !(1 <= axial_index <= length(region_cell_ids))
        error("region \"$region_name\" has $(length(region_cell_ids)) cell(s); axial_index $axial_index is out of range")
    end
    return Int(region_cell_ids[axial_index])
end

"""
    add_virtual_cell!(rom_data; name, volume, faces, centroid = (0, 0, 0))

Add a cell with no CAD volume (e.g. `recirculation_zone`) and connect it to
real cells. Must be called BEFORE the geometry is built: the half-face
layout (CellFaceLayout offsets) is fixed after that.

Each entry of `faces` is a NamedTuple:
  - `neighbor_region`: the real region it touches, e.g. "fluid_domain"
  - `neighbor_axial_index`: which cell of that region: an integer, `:first` or `:last`
  - `area`: shared face area (m^2)
  - `half_cell_distance`: from the virtual cell's node to the face (m)
  - `neighbor_half_cell_distance`: from the real cell's node to the face (m)
  - optional `normal`: unit normal pointing out of the virtual cell (default zero vector)

For every face this adds the matching half-face to the real cell, so the
interface is reciprocal like any exported one. It also adds the cell set
`name` and, per touched region, the face set "<real_region>__<name>"
(holding the real cells' half-faces, like any "<a>__<b>" set).
"""
function add_virtual_cell!(rom_data; name, volume, faces, centroid = (0.0, 0.0, 0.0))
    if haskey(rom_data["cellsets"], name)
        error("a cell set named \"$name\" already exists")
    end

    virtual_cell_id = rom_cell_count(rom_data) + 1
    virtual_cell_faces = Any[]

    for (virtual_face_idx, face_spec) in enumerate(faces)
        neighbor_region = face_spec.neighbor_region
        neighbor_cell_id = rom_region_cell_id(rom_data, neighbor_region, face_spec.neighbor_axial_index)
        neighbor_cell = rom_data["cells"][neighbor_cell_id]
        neighbor_face_idx = length(neighbor_cell["faces"]) + 1

        virtual_normal = hasproperty(face_spec, :normal) ? collect(Float64, face_spec.normal) : nothing
        neighbor_normal = isnothing(virtual_normal) ? nothing : -virtual_normal

        # Series distance between the two nodes, as for exported interfaces.
        neighbor_distance = face_spec.half_cell_distance + face_spec.neighbor_half_cell_distance

        push!(virtual_cell_faces, OrderedDict{String, Any}(
            "local_idx" => virtual_face_idx, "kind" => "interface",
            "neighbor_cell" => neighbor_cell_id, "neighbor_local_idx" => neighbor_face_idx,
            "area" => face_spec.area, "surface" => "virtual", "normal" => virtual_normal, "centroid" => nothing,
            "half_cell_distance" => face_spec.half_cell_distance, "neighbor_distance" => neighbor_distance,
            "distance_method" => "virtual", "notes" => "added by add_virtual_cell!",
        ))

        push!(neighbor_cell["faces"], OrderedDict{String, Any}(
            "local_idx" => neighbor_face_idx, "kind" => "interface",
            "neighbor_cell" => virtual_cell_id, "neighbor_local_idx" => virtual_face_idx,
            "area" => face_spec.area, "surface" => "virtual", "normal" => neighbor_normal, "centroid" => nothing,
            "half_cell_distance" => face_spec.neighbor_half_cell_distance, "neighbor_distance" => neighbor_distance,
            "distance_method" => "virtual", "notes" => "face to virtual cell \"$name\"",
        ))

        face_set_name = "$(neighbor_region)__$(name)"
        push!(get!(rom_data["facetsets"], face_set_name, Any[]), Any[neighbor_cell_id, neighbor_face_idx])
    end

    push!(rom_data["cells"], OrderedDict{String, Any}(
        "id" => virtual_cell_id, "region" => name, "axial_index" => 1,
        "volume" => volume, "centroid" => collect(Float64, centroid), "faces" => virtual_cell_faces,
    ))
    rom_data["cellsets"][name] = Any[virtual_cell_id]
    push!(rom_data["regions"], OrderedDict{String, Any}("name" => name, "virtual" => true, "cell_ids" => Any[virtual_cell_id]))

    return virtual_cell_id
end

"""
    add_virtual_face!(rom_data; region_a, axial_index_a, region_b, axial_index_b,
                      area, half_cell_distance_a, half_cell_distance_b, normal = nothing)

Connect two existing cells that do not touch in the CAD (e.g. pre_chamber and
mid_section) with a face of your own. Must be called BEFORE the geometry is
built: the half-face layout (CellFaceLayout offsets) is fixed after that.

  - `region_a`, `axial_index_a`: the first cell; `axial_index_a` is an integer
    (1 = lowest z), `:first` or `:last`. Same for `region_b`, `axial_index_b`.
  - `area`: face area (m^2), e.g. a flow area from a closure
  - `half_cell_distance_a`: from cell a's node to the face (m); same for b.
    Their sum is the node-to-node distance the flux functions use.
  - optional `normal`: unit normal pointing from cell a into cell b (default
    zero vector, as for a full 360 degree cylinder)

Each cell gets one new half-face pointing at the other, so the result is a
normal reciprocal interface: it shows up in `cell_neighbors`, in the
connection group of the two regions (finish_fvm_config), and in
`get_facetset(geo, region_a, region_b)` from either side. It is stored in
the region pair's existing face set if there is one ("<a>__<b>" or
"<b>__<a>"), otherwise in a new "<region_a>__<region_b>". Returns the two
new local face indices `(face_idx_a, face_idx_b)`.
"""
function add_virtual_face!(
    rom_data;
    region_a, axial_index_a,
    region_b, axial_index_b,
    area, half_cell_distance_a, half_cell_distance_b,
    normal = nothing,
)
    cell_id_a = rom_region_cell_id(rom_data, region_a, axial_index_a)
    cell_id_b = rom_region_cell_id(rom_data, region_b, axial_index_b)
    if cell_id_a == cell_id_b
        error("add_virtual_face!: both sides are cell $cell_id_a; a face needs two different cells")
    end

    cell_a = rom_data["cells"][cell_id_a]
    cell_b = rom_data["cells"][cell_id_b]
    face_idx_a = length(cell_a["faces"]) + 1
    face_idx_b = length(cell_b["faces"]) + 1

    normal_a = isnothing(normal) ? nothing : collect(Float64, normal)
    normal_b = isnothing(normal_a) ? nothing : -normal_a

    # Series distance between the two nodes, as for exported interfaces.
    neighbor_distance = half_cell_distance_a + half_cell_distance_b

    push!(cell_a["faces"], OrderedDict{String, Any}(
        "local_idx" => face_idx_a, "kind" => "interface",
        "neighbor_cell" => cell_id_b, "neighbor_local_idx" => face_idx_b,
        "area" => area, "surface" => "virtual", "normal" => normal_a, "centroid" => nothing,
        "half_cell_distance" => half_cell_distance_a, "neighbor_distance" => neighbor_distance,
        "distance_method" => "virtual", "notes" => "virtual face to $region_b, added by add_virtual_face!",
    ))

    push!(cell_b["faces"], OrderedDict{String, Any}(
        "local_idx" => face_idx_b, "kind" => "interface",
        "neighbor_cell" => cell_id_a, "neighbor_local_idx" => face_idx_a,
        "area" => area, "surface" => "virtual", "normal" => normal_b, "centroid" => nothing,
        "half_cell_distance" => half_cell_distance_b, "neighbor_distance" => neighbor_distance,
        "distance_method" => "virtual", "notes" => "virtual face to $region_a, added by add_virtual_face!",
    ))

    # One face set per region pair: get_facetset(geo, a, b) looks up
    # "<a>__<b>" first, so a second "<b>__<a>" set would hide these faces.
    # Inside one region the set holds the lower cell id's half-face.
    forward_set_name = "$(region_a)__$(region_b)"
    reverse_set_name = "$(region_b)__$(region_a)"
    if region_a == region_b
        if cell_id_a < cell_id_b
            face_set_entry = Any[cell_id_a, face_idx_a]
        else
            face_set_entry = Any[cell_id_b, face_idx_b]
        end
        push!(get!(rom_data["facetsets"], forward_set_name, Any[]), face_set_entry)
    elseif haskey(rom_data["facetsets"], forward_set_name) || !haskey(rom_data["facetsets"], reverse_set_name)
        push!(get!(rom_data["facetsets"], forward_set_name, Any[]), Any[cell_id_a, face_idx_a])
    else
        push!(rom_data["facetsets"][reverse_set_name], Any[cell_id_b, face_idx_b])
    end

    return face_idx_a, face_idx_b
end

"""
    build_rom_fvm_geometry(rom_data)

Fill an `FVMGeometry` from a ROM export. Every face of every cell becomes one
half-face, in the export's local order. Node-based fields are left empty
(nothing in the setup or solver path reads them), `null` normals (full 360
degree cylinders and cones, and virtual faces without a normal) become zero
vectors, and `update_fvm_geometry!` must NOT be called on the result.
"""
function build_rom_fvm_geometry(rom_data)
    rom_cells = rom_data["cells"]
    n_cells = length(rom_cells)

    cell_n_faces = [length(rom_cell["faces"]) for rom_cell in rom_cells]
    cell_face_layout = CellFaceLayout(cell_n_faces)
    cell_face_map = [collect(cell_face_layout[cell_id]) for cell_id in 1:n_cells]

    cell_volumes = zeros(Float64, n_cells)
    cell_centroids = zeros(ROMCoordinate, n_cells)
    cell_neighbor_normals = cell_face_layout_zeros(ROMCoordinate, cell_face_layout)
    cell_neighbor_distances = cell_face_layout_zeros(Float64, cell_face_layout)
    cell_face_areas = cell_face_layout_zeros(Float64, cell_face_layout)
    cell_face_normals = cell_face_layout_zeros(ROMCoordinate, cell_face_layout)
    cell_face_distances = cell_face_layout_zeros(Float64, cell_face_layout)
    cell_face_centroids = cell_face_layout_zeros(ROMCoordinate, cell_face_layout)

    cell_neighbors = [(cell_id, NTuple{3, Int}[]) for cell_id in 1:n_cells]
    unconnected_cell_face_map = Tuple{Int, Int}[]

    for cell_id in 1:n_cells
        rom_cell = rom_cells[cell_id]
        cell_volumes[cell_id] = rom_cell["volume"]
        cell_centroids[cell_id] = rom_coordinate(rom_cell["centroid"])

        for (face_idx, rom_face) in enumerate(rom_cell["faces"])
            if rom_face["local_idx"] != face_idx
                error("cell $cell_id: face number $face_idx has local_idx $(rom_face["local_idx"])")
            end

            face_normal = rom_coordinate(rom_face["normal"])
            cell_face_areas[cell_id, face_idx] = rom_face["area"]
            cell_face_normals[cell_id, face_idx] = face_normal
            cell_face_distances[cell_id, face_idx] = rom_face["half_cell_distance"]
            cell_face_centroids[cell_id, face_idx] = rom_coordinate(rom_face["centroid"])

            # Ignored interfaces and boundary faces have no neighbour.
            if rom_face["kind"] != "interface"
                push!(unconnected_cell_face_map, (cell_id, face_idx))
                continue
            end

            neighbor_cell_id = rom_face["neighbor_cell"]
            neighbor_face_idx = rom_face["neighbor_local_idx"]
            push!(cell_neighbors[cell_id][2], (neighbor_cell_id, face_idx, neighbor_face_idx))

            # The outward normal of an interface face points into the
            # neighbour, so it is also the neighbour direction.
            cell_neighbor_normals[cell_id, face_idx] = face_normal
            cell_neighbor_distances[cell_id, face_idx] = rom_face["neighbor_distance"]
        end
    end

    cellsets = OrderedDict{String, Vector{Int}}()
    for (set_name, cell_ids) in rom_data["cellsets"]
        cellsets[set_name] = Int.(cell_ids)
    end

    facetsets = OrderedDict{String, Vector{Tuple{Int, Vector{NTuple{3, Int}}}}}()
    for (set_name, cell_face_pairs) in rom_data["facetsets"]
        facetsets[set_name] = build_named_facetset([(Int(pair[1]), Int(pair[2])) for pair in cell_face_pairs], cell_neighbors)
    end

    n_total_faces = total_faces(cell_face_layout)
    return FVMGeometry(
        #Topology information (no nodes in a ROM)
        Vector{Int}[],
        cell_neighbors,
        [[Int[] for face_idx in 1:cell_n_faces[cell_id]] for cell_id in 1:n_cells],
        cell_face_map,
        [Int[] for packed_face_id in 1:n_total_faces],
        cell_n_faces,
        unconnected_cell_face_map,
        n_cells,

        #Geometry information
        cell_face_layout,
        cell_volumes,
        cell_centroids,
        cell_neighbor_normals,
        cell_neighbor_distances,
        cell_face_areas,
        cell_face_normals,
        cell_face_distances,
        cell_face_centroids,
        ROMCoordinate[],

        #Named sets
        cellsets,
        facetsets,
    )
end

"""
    create_fvm_config_from_rom(rom_data, u_proto)

`SimulationConfigInfo` for a ROM. After this, `add_region!` (with the region
names of the export), `add_patch!` (with a patch name or "<region_a>__<region_b>")
and `finish_fvm_config` work exactly as for a mesh.
"""
function create_fvm_config_from_rom(rom_data, u_proto)
    geo = build_rom_fvm_geometry(rom_data)
    return create_fvm_config_from_geo(geo, u_proto)
end
