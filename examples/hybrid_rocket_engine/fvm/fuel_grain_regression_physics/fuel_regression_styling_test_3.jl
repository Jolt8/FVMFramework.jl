"""
Geometric information for one cell face lying on an axial cross-section of the
fuel port.

The two radial face indices are retained because the surface temperature model
needs the cell-to-wall distance on the face nearest the burning fuel surface.
"""
struct CellFaceOnAxialCrossSection
    cell_id::Int
    axial_face_idx::Int
    face_idx_toward_fuel_wall::Int
    face_idx_away_from_fuel_wall::Int
end

"""
Description of the fluid and solid geometry at one axial fuel-port station.

One fuel thickness state is stored at every station. During regression, the
inner fuel surface moves outward while the outside of the fuel grain remains
fixed.
"""
struct FuelPortStation
    axial_coordinate::Float64
    upstream_fluid_cross_section::Vector{CellFaceOnAxialCrossSection}
    downstream_fluid_cross_section::Vector{CellFaceOnAxialCrossSection}
    upstream_fuel_cross_section::Vector{CellFaceOnAxialCrossSection}
    downstream_fuel_cross_section::Vector{CellFaceOnAxialCrossSection}
    initial_axis_coordinate::Float64
    initial_port_radius::Float64
    initial_fuel_thickness::Float64
    fixed_fuel_outer_radius::Float64
    fuel_thickness_state_idx::Int
end

"""
Affine radial motion assigned to one mesh node.

`fraction_of_wall_displacement` is one at the fuel surface and zero at either
the symmetry axis or the fixed outside surface of the grain. Intermediate
nodes move linearly between those boundaries.
"""
struct FuelRegressionMeshNode
    node_id::Int
    fuel_port_station_id::Int
    initial_radial_coordinate::Float64
    fraction_of_wall_displacement::Float64
end

"""A fluid/solid face pair on the regressing fuel wall."""
struct FuelWallInterface
    port_fluid_cell_id::Int
    fluid_face_idx::Int
    solid_fuel_cell_id::Int
    fuel_face_idx::Int
end

struct FuelRegressionWorkspace{NV, RR}
    mesh_node_velocities::NV
    fuel_surface_regression_rates::RR
end

"""
All geometry and reusable numerical bookkeeping required by the fuel
regression model.

The model assumes an axisymmetric meridional mesh whose radial coordinate is
`x` and whose motor-axis coordinate is `z`.
"""
struct FuelRegressionModel{W}
    fuel_port_stations::Vector{FuelPortStation}
    moving_mesh_nodes::Vector{FuelRegressionMeshNode}
    fuel_wall_interfaces::Vector{FuelWallInterface}
    fluid_cells_affected_by_wall_motion::Vector{Int}
    is_fluid_cell::BitVector
    neighbor_cell_id_by_face::Vector{Vector{Int}}
    initial_node_coordinates::Vector{SVector{3, Float64}}
    workspaces::W
end


# -----------------------------------------------------------------------------
# Describe the fuel grain, port, and fluid/solid interface
# -----------------------------------------------------------------------------

function sorted_unique_coordinates(coordinates; absolute_tolerance)
    sorted_coordinates = sort!(collect(coordinates))
    if isempty(sorted_coordinates)
        return Float64[]
    end

    unique_coordinates = Float64[first(sorted_coordinates)]
    for coordinate in @view sorted_coordinates[2:end]
        if abs(coordinate - last(unique_coordinates)) > absolute_tolerance
            push!(unique_coordinates, coordinate)
        end
    end

    return unique_coordinates
end

function node_ids_on_cell_face(geo, cell_id, face_idx)
    packed_face_id = geo.cell_face_layout[cell_id, face_idx]
    return geo.cell_faces_node_ids[packed_face_id]
end

function find_fuel_port_station_id(station_axial_coordinates, node_axial_coordinate; absolute_tolerance)
    for (station_id, station_axial_coordinate) in enumerate(station_axial_coordinates)
        if abs(station_axial_coordinate - node_axial_coordinate) <= absolute_tolerance
            return station_id
        end
    end

    error("No fuel-port station was found at z = $node_axial_coordinate")
end

function find_cell_face_on_axial_plane(geo, cell_id, axial_coordinate; absolute_tolerance)
    for face_idx in 1:geo.cell_n_faces[cell_id]
        face_node_ids = node_ids_on_cell_face(geo, cell_id, face_idx)
        all_face_nodes_are_on_plane = true

        for node_id in face_node_ids
            node_axial_coordinate = geo.node_coordinates[node_id][3]
            if abs(node_axial_coordinate - axial_coordinate) > absolute_tolerance
                all_face_nodes_are_on_plane = false
                break
            end
        end

        if all_face_nodes_are_on_plane
            return face_idx
        end
    end

    return nothing
end


function identify_radial_faces(geo, cell_id, material_region::Symbol)
    radial_face_indices = Int[]

    for face_idx in 1:geo.cell_n_faces[cell_id]
        outward_face_normal = geo.cell_face_normals[cell_id, face_idx]
        if abs(outward_face_normal[1]) > 0.5
            push!(radial_face_indices, face_idx)
        end
    end

    if length(radial_face_indices) != 2
        error("Expected two radial faces for cell $cell_id, found $(length(radial_face_indices))")
    end

    sort!(radial_face_indices; by = face_idx -> geo.cell_face_centroids[cell_id, face_idx][1])
    inner_radial_face_idx = first(radial_face_indices)
    outer_radial_face_idx = last(radial_face_indices)

    if material_region === :port_fluid
        face_idx_toward_fuel_wall = outer_radial_face_idx
        face_idx_away_from_fuel_wall = inner_radial_face_idx
    elseif material_region === :solid_fuel
        face_idx_toward_fuel_wall = inner_radial_face_idx
        face_idx_away_from_fuel_wall = outer_radial_face_idx
    else
        error("Unknown material region: $material_region")
    end

    return face_idx_toward_fuel_wall, face_idx_away_from_fuel_wall
end

function describe_cell_face_on_cross_section(geo, cell_id, axial_face_idx, material_region)
    face_idx_toward_fuel_wall, face_idx_away_from_fuel_wall = identify_radial_faces(geo, cell_id, material_region)

    return CellFaceOnAxialCrossSection(
        cell_id,
        axial_face_idx,
        face_idx_toward_fuel_wall,
        face_idx_away_from_fuel_wall,
    )
end

function collect_cells_touching_cross_section(
        geo,
        candidate_cell_ids,
        axial_coordinate,
        side_of_cross_section::Symbol,
        material_region::Symbol;
        absolute_tolerance,
    )
    cells_on_cross_section = CellFaceOnAxialCrossSection[]

    for cell_id in candidate_cell_ids
        cell_axial_coordinate = geo.cell_centroids[cell_id][3]

        if side_of_cross_section === :upstream
            cell_is_on_requested_side = cell_axial_coordinate < axial_coordinate - absolute_tolerance
        elseif side_of_cross_section === :downstream
            cell_is_on_requested_side = cell_axial_coordinate > axial_coordinate + absolute_tolerance
        else
            error("Unknown side of cross-section: $side_of_cross_section")
        end

        if !cell_is_on_requested_side
            continue
        end

        axial_face_idx = find_cell_face_on_axial_plane(geo, cell_id, axial_coordinate; absolute_tolerance)
        if isnothing(axial_face_idx)
            continue
        end

        cross_section_face = describe_cell_face_on_cross_section(geo, cell_id, axial_face_idx, material_region)
        push!(cells_on_cross_section, cross_section_face)
    end

    # Sorting from the motor axis toward the wall makes the first and last
    # entries have an immediate geometric meaning later in the model.
    sort!(cells_on_cross_section; by = entry -> geo.cell_centroids[entry.cell_id][1])

    return cells_on_cross_section
end

function identify_fuel_wall_interfaces(geo, port_fluid_cell_ids, solid_fuel_cell_ids)
    solid_fuel_cell_set = Set(solid_fuel_cell_ids)
    fuel_wall_interfaces = FuelWallInterface[]

    for port_fluid_cell_id in port_fluid_cell_ids
        for (neighbor_cell_id, fluid_face_idx, neighbor_face_idx) in geo.cell_neighbors[port_fluid_cell_id][2]
            if neighbor_cell_id in solid_fuel_cell_set
                push!(fuel_wall_interfaces, FuelWallInterface(
                    port_fluid_cell_id,
                    fluid_face_idx,
                    neighbor_cell_id,
                    neighbor_face_idx,
                ))
            end
        end
    end

    return fuel_wall_interfaces
end

function map_neighbor_cells_to_faces(geo)
    neighbor_cell_id_by_face = [zeros(Int, geo.cell_n_faces[cell_id]) for cell_id in 1:geo.n_cells]

    for (cell_id, cell_neighbors) in geo.cell_neighbors
        for (neighbor_cell_id, cell_face_idx, _) in cell_neighbors
            neighbor_cell_id_by_face[cell_id][cell_face_idx] = neighbor_cell_id
        end
    end

    return neighbor_cell_id_by_face
end

"""
    create_fuel_regression_model(grid, geo; absolute_tolerance = 1e-10)

Inspect the initial mesh and construct a physical description of the fuel port.

Modeling assumptions:

- `x` is the radial coordinate and `z` is the motor-axis coordinate.
- One fuel-thickness degree of freedom is used at every axial node plane.
- The port wall is initially the largest fluid radius at each station.
- The outside surface of the fuel grain is fixed.
- Interior mesh nodes move linearly between a fixed boundary and the port wall.
"""
function create_fuel_regression_model(grid, geo; absolute_tolerance = 1e-10)
    solid_fuel_cell_ids = sort!(collect(getcellset(grid, "fuel_grain_volume")))
    all_fluid_cell_ids = sort!(collect(getcellset(grid, "fluid_domain")))

    if isempty(solid_fuel_cell_ids)
        error("The fuel_grain_volume cell set is empty")
    end

    fuel_cell_axial_coordinates = sorted_unique_coordinates(
        (geo.cell_centroids[cell_id][3] for cell_id in solid_fuel_cell_ids);
        absolute_tolerance,
    )

    port_fluid_cell_ids = Int[]
    for fluid_cell_id in all_fluid_cell_ids
        fluid_cell_axial_coordinate = geo.cell_centroids[fluid_cell_id][3]

        for fuel_cell_axial_coordinate in fuel_cell_axial_coordinates
            if abs(fluid_cell_axial_coordinate - fuel_cell_axial_coordinate) <= absolute_tolerance
                push!(port_fluid_cell_ids, fluid_cell_id)
                break
            end
        end
    end

    solid_fuel_node_ids = Int[]
    for solid_fuel_cell_id in solid_fuel_cell_ids
        append!(solid_fuel_node_ids, geo.cell_node_ids[solid_fuel_cell_id])
    end
    sort!(solid_fuel_node_ids)
    unique!(solid_fuel_node_ids)

    station_axial_coordinates = sorted_unique_coordinates(
        (geo.node_coordinates[node_id][3] for node_id in solid_fuel_node_ids);
        absolute_tolerance,
    )

    all_node_ids_at_station = [Set{Int}() for _ in station_axial_coordinates]
    fluid_node_ids_at_station = [Set{Int}() for _ in station_axial_coordinates]
    fuel_node_ids_at_station = [Set{Int}() for _ in station_axial_coordinates]

    for port_fluid_cell_id in port_fluid_cell_ids
        for node_id in geo.cell_node_ids[port_fluid_cell_id]
            node_axial_coordinate = geo.node_coordinates[node_id][3]
            station_id = find_fuel_port_station_id(station_axial_coordinates, node_axial_coordinate; absolute_tolerance)

            push!(all_node_ids_at_station[station_id], node_id)
            push!(fluid_node_ids_at_station[station_id], node_id)
        end
    end

    for solid_fuel_cell_id in solid_fuel_cell_ids
        for node_id in geo.cell_node_ids[solid_fuel_cell_id]
            node_axial_coordinate = geo.node_coordinates[node_id][3]
            station_id = find_fuel_port_station_id(station_axial_coordinates, node_axial_coordinate; absolute_tolerance)

            push!(all_node_ids_at_station[station_id], node_id)
            push!(fuel_node_ids_at_station[station_id], node_id)
        end
    end

    fuel_port_stations = FuelPortStation[]
    moving_mesh_nodes = FuelRegressionMeshNode[]

    for (station_id, axial_coordinate) in enumerate(station_axial_coordinates)
        fluid_node_ids = collect(fluid_node_ids_at_station[station_id])
        solid_fuel_node_ids_at_this_station = collect(fuel_node_ids_at_station[station_id])

        if isempty(fluid_node_ids)
            error("No port-fluid nodes were found at z = $axial_coordinate")
        end
        if isempty(solid_fuel_node_ids_at_this_station)
            error("No solid-fuel nodes were found at z = $axial_coordinate")
        end

        initial_axis_coordinate = minimum(node_id -> geo.node_coordinates[node_id][1], fluid_node_ids)
        initial_port_radius = maximum(node_id -> geo.node_coordinates[node_id][1], fluid_node_ids)
        fixed_fuel_outer_radius = maximum(node_id -> geo.node_coordinates[node_id][1], solid_fuel_node_ids_at_this_station)
        initial_fuel_thickness = fixed_fuel_outer_radius - initial_port_radius

        upstream_fluid_cross_section = collect_cells_touching_cross_section(
            geo, port_fluid_cell_ids, axial_coordinate, :upstream, :port_fluid; absolute_tolerance,
        )
        downstream_fluid_cross_section = collect_cells_touching_cross_section(
            geo, port_fluid_cell_ids, axial_coordinate, :downstream, :port_fluid; absolute_tolerance,
        )
        upstream_fuel_cross_section = collect_cells_touching_cross_section(
            geo, solid_fuel_cell_ids, axial_coordinate, :upstream, :solid_fuel; absolute_tolerance,
        )
        downstream_fuel_cross_section = collect_cells_touching_cross_section(
            geo, solid_fuel_cell_ids, axial_coordinate, :downstream, :solid_fuel; absolute_tolerance,
        )

        push!(fuel_port_stations, FuelPortStation(
            axial_coordinate,
            upstream_fluid_cross_section,
            downstream_fluid_cross_section,
            upstream_fuel_cross_section,
            downstream_fuel_cross_section,
            initial_axis_coordinate,
            initial_port_radius,
            initial_fuel_thickness,
            fixed_fuel_outer_radius,
            station_id,
        ))

        station_node_ids = sort!(collect(all_node_ids_at_station[station_id]))
        for node_id in station_node_ids
            initial_radial_coordinate = geo.node_coordinates[node_id][1]

            if initial_radial_coordinate <= initial_port_radius + absolute_tolerance
                distance_from_axis = initial_radial_coordinate - initial_axis_coordinate
                initial_port_span = initial_port_radius - initial_axis_coordinate
                fraction_of_wall_displacement = distance_from_axis / initial_port_span
            else
                distance_from_initial_wall = initial_radial_coordinate - initial_port_radius
                fraction_through_solid_fuel = distance_from_initial_wall / initial_fuel_thickness
                fraction_of_wall_displacement = 1 - fraction_through_solid_fuel
            end

            push!(moving_mesh_nodes, FuelRegressionMeshNode(
                node_id,
                station_id,
                initial_radial_coordinate,
                fraction_of_wall_displacement,
            ))
        end
    end

    moving_node_ids = Set{Int}()
    for mesh_node in moving_mesh_nodes
        if abs(mesh_node.fraction_of_wall_displacement) > absolute_tolerance
            push!(moving_node_ids, mesh_node.node_id)
        end
    end

    fluid_cells_affected_by_wall_motion = Int[]
    for fluid_cell_id in all_fluid_cell_ids
        cell_moves = false
        for node_id in geo.cell_node_ids[fluid_cell_id]
            if node_id in moving_node_ids
                cell_moves = true
                break
            end
        end

        if cell_moves
            push!(fluid_cells_affected_by_wall_motion, fluid_cell_id)
        end
    end

    is_fluid_cell = falses(geo.n_cells)
    is_fluid_cell[all_fluid_cell_ids] .= true

    initial_node_coordinates = SVector{3, Float64}.(geo.node_coordinates)
    number_of_mesh_nodes = length(geo.node_coordinates)
    number_of_fuel_port_stations = length(fuel_port_stations)

    workspaces = GeneralLazyBufferCache() do fuel_thickness_seed
        zero_rate = zero(fuel_thickness_seed)
        zero_mesh_velocity = SVector(zero_rate, zero_rate, zero_rate)

        FuelRegressionWorkspace(
            fill(zero_mesh_velocity, number_of_mesh_nodes),
            fill(zero_rate, number_of_fuel_port_stations),
        )
    end

    fuel_wall_interfaces = identify_fuel_wall_interfaces(geo, port_fluid_cell_ids, solid_fuel_cell_ids)
    neighbor_cell_id_by_face = map_neighbor_cells_to_faces(geo)

    return FuelRegressionModel(
        fuel_port_stations,
        moving_mesh_nodes,
        fuel_wall_interfaces,
        fluid_cells_affected_by_wall_motion,
        is_fluid_cell,
        neighbor_cell_id_by_face,
        initial_node_coordinates,
        workspaces,
    )
end


# -----------------------------------------------------------------------------
# Initialize and move the regressing fuel geometry
# -----------------------------------------------------------------------------

function initialize_fuel_grain_thickness!(fuel_grain_thickness, fuel_regression_model::FuelRegressionModel)
    thickness_unit = oneunit(first(fuel_grain_thickness))

    for fuel_port_station in fuel_regression_model.fuel_port_stations
        state_idx = fuel_port_station.fuel_thickness_state_idx
        initial_fuel_thickness = fuel_port_station.initial_fuel_thickness * thickness_unit
        fuel_grain_thickness[state_idx] = initial_fuel_thickness
    end

    return fuel_grain_thickness
end

function move_fuel_wall_and_interior_mesh!(u, fuel_regression_model::FuelRegressionModel, geo)
    # Always reconstruct from the initial mesh. This avoids accumulating mesh
    # motion error when an implicit integrator evaluates the same time repeatedly.
    copyto!(geo.node_coordinates, fuel_regression_model.initial_node_coordinates)

    for mesh_node in fuel_regression_model.moving_mesh_nodes
        fuel_port_station = fuel_regression_model.fuel_port_stations[mesh_node.fuel_port_station_id]
        fuel_thickness = u.fuel_grain_thickness[fuel_port_station.fuel_thickness_state_idx]
        fuel_wall_displacement = fuel_port_station.initial_fuel_thickness - fuel_thickness

        radial_node_displacement = mesh_node.fraction_of_wall_displacement * fuel_wall_displacement
        new_radial_coordinate = mesh_node.initial_radial_coordinate + radial_node_displacement

        old_node_coordinate = geo.node_coordinates[mesh_node.node_id]
        geo.node_coordinates[mesh_node.node_id] = setindex(old_node_coordinate, new_radial_coordinate, 1)
    end

    return nothing
end

function update_fuel_regression_geometry!(u, fuel_regression_model::FuelRegressionModel, geo, grid)
    move_fuel_wall_and_interior_mesh!(u, fuel_regression_model, geo)

    update_fvm_geometry!(geo, geo.node_coordinates)

    # The hybrid-motor mesh needs its axisymmetric geometric corrections
    # reapplied after every wall movement.
    correct_hybrid_geo!(geo, grid)

    return nothing
end


# -----------------------------------------------------------------------------
# Evaluate the local regression law at every axial fuel-port station
# -----------------------------------------------------------------------------

function choose_cross_section_for_regression_measurements(fuel_port_station::FuelPortStation)
    if !isempty(fuel_port_station.upstream_fluid_cross_section)
        fluid_cross_section = fuel_port_station.upstream_fluid_cross_section
        fuel_cross_section = fuel_port_station.upstream_fuel_cross_section
    else
        # The inlet station has no cells upstream, so use its downstream face.
        fluid_cross_section = fuel_port_station.downstream_fluid_cross_section
        fuel_cross_section = fuel_port_station.downstream_fuel_cross_section
    end

    return fluid_cross_section, fuel_cross_section
end

function calculate_cross_sectional_oxidizer_mass_flux(u, geo, fluid_cross_section)
    reference_fluid_face = first(fluid_cross_section)
    reference_fluid_cell_id = reference_fluid_face.cell_id
    reference_face_area = geo.cell_face_areas[reference_fluid_cell_id, reference_fluid_face.axial_face_idx]
    reference_axial_velocity = u.momentum_density_w[reference_fluid_cell_id] / oneunit(u.density[reference_fluid_cell_id])

    nitrous_oxide_mass_flow = zero(reference_axial_velocity * u.species_densities.nitrous_oxide[reference_fluid_cell_id] * reference_face_area)
    port_cross_sectional_area = zero(reference_face_area)

    for fluid_face in fluid_cross_section
        fluid_cell_id = fluid_face.cell_id
        axial_face_idx = fluid_face.axial_face_idx

        fluid_density = u.density[fluid_cell_id]
        if iszero(fluid_density)
            axial_fluid_velocity = zero(u.momentum_density_w[fluid_cell_id] / oneunit(fluid_density))
        else
            axial_fluid_velocity = u.momentum_density_w[fluid_cell_id] / fluid_density
        end

        nitrous_oxide_density = u.species_densities.nitrous_oxide[fluid_cell_id]
        axial_face_area = geo.cell_face_areas[fluid_cell_id, axial_face_idx]
        local_nitrous_oxide_mass_flow = abs(axial_fluid_velocity) * nitrous_oxide_density * axial_face_area

        nitrous_oxide_mass_flow += local_nitrous_oxide_mass_flow
        port_cross_sectional_area += axial_face_area
    end

    if iszero(port_cross_sectional_area)
        return zero(nitrous_oxide_mass_flow / oneunit(port_cross_sectional_area))
    end

    nitrous_oxide_mass_flux = nitrous_oxide_mass_flow / port_cross_sectional_area
    return nitrous_oxide_mass_flux
end

function estimate_fuel_surface_temperature(u, geo, fluid_cross_section, fuel_cross_section)
    fluid_cell_next_to_wall = last(fluid_cross_section)
    fuel_cell_next_to_wall = first(fuel_cross_section)

    fluid_cell_id = fluid_cell_next_to_wall.cell_id
    solid_fuel_cell_id = fuel_cell_next_to_wall.cell_id

    fluid_distance_to_wall = geo.cell_face_distances[fluid_cell_id, fluid_cell_next_to_wall.face_idx_toward_fuel_wall]
    fuel_distance_to_wall = geo.cell_face_distances[solid_fuel_cell_id, fuel_cell_next_to_wall.face_idx_toward_fuel_wall]

    fluid_thermal_conductance_per_area = u.k[fluid_cell_id] / fluid_distance_to_wall
    fuel_thermal_conductance_per_area = u.k[solid_fuel_cell_id] / fuel_distance_to_wall
    total_thermal_conductance_per_area = fluid_thermal_conductance_per_area + fuel_thermal_conductance_per_area

    # This is a two-resistance conductive interface temperature. It does not
    # yet include a flame-side convection model, radiation, or a surface energy
    # balance for pyrolysis. Those omissions are important model limitations.
    if iszero(total_thermal_conductance_per_area)
        fuel_surface_temperature = (u.temperature[fluid_cell_id] + u.temperature[solid_fuel_cell_id]) / 2
    else
        conductance_weighted_temperature = fluid_thermal_conductance_per_area * u.temperature[fluid_cell_id] + fuel_thermal_conductance_per_area * u.temperature[solid_fuel_cell_id]
        fuel_surface_temperature = conductance_weighted_temperature / total_thermal_conductance_per_area
    end

    return fuel_surface_temperature
end

function calculate_fuel_regression_rates!(du, u, p, fuel_regression_model::FuelRegressionModel, workspace, geo)
    fill!(workspace.fuel_surface_regression_rates, zero(eltype(workspace.fuel_surface_regression_rates)))

    for (station_id, fuel_port_station) in enumerate(fuel_regression_model.fuel_port_stations)
        fluid_cross_section, fuel_cross_section = choose_cross_section_for_regression_measurements(fuel_port_station)

        if isempty(fluid_cross_section) || isempty(fuel_cross_section)
            continue
        end

        nitrous_oxide_mass_flux = calculate_cross_sectional_oxidizer_mass_flux(u, geo, fluid_cross_section)
        fuel_surface_temperature = estimate_fuel_surface_temperature(u, geo, fluid_cross_section, fuel_cross_section)

        # The empirical law is enabled only after the solid surface reaches the
        # specified pyrolysis threshold. The correlation is r_dot = a G_ox^n.
        if fuel_surface_temperature > p.min_fuel_regression_temperature
            nonnegative_nitrous_oxide_mass_flux = max(nitrous_oxide_mass_flux, zero(nitrous_oxide_mass_flux))
            fuel_surface_regression_rate = p.regression_coefficient_a * nonnegative_nitrous_oxide_mass_flux ^ p.regression_exponent_n

            fuel_thickness_state_idx = fuel_port_station.fuel_thickness_state_idx
            du.fuel_grain_thickness[fuel_thickness_state_idx] = -fuel_surface_regression_rate
            workspace.fuel_surface_regression_rates[station_id] = fuel_surface_regression_rate
        end
    end

    return nothing
end


# -----------------------------------------------------------------------------
# Convert fuel regression into mesh motion
# -----------------------------------------------------------------------------

function calculate_mesh_node_velocities!(fuel_regression_model::FuelRegressionModel, workspace)
    fill!(workspace.mesh_node_velocities, zero(eltype(workspace.mesh_node_velocities)))

    for mesh_node in fuel_regression_model.moving_mesh_nodes
        station_regression_rate = workspace.fuel_surface_regression_rates[mesh_node.fuel_port_station_id]
        radial_mesh_velocity = mesh_node.fraction_of_wall_displacement * station_regression_rate

        workspace.mesh_node_velocities[mesh_node.node_id] = SVector(
            radial_mesh_velocity,
            zero(radial_mesh_velocity),
            zero(radial_mesh_velocity),
        )
    end

    return nothing
end

function calculate_face_mesh_velocity(workspace, geo, cell_id, face_idx)
    face_node_ids = node_ids_on_cell_face(geo, cell_id, face_idx)
    face_mesh_velocity = zero(eltype(workspace.mesh_node_velocities))

    for node_id in face_node_ids
        face_mesh_velocity += workspace.mesh_node_velocities[node_id]
    end

    face_mesh_velocity /= length(face_node_ids)
    return face_mesh_velocity
end


# -----------------------------------------------------------------------------
# Conserve fluid quantities while the computational cells change volume
# -----------------------------------------------------------------------------

function calculate_moving_mesh_flux_for_conserved_quantity(
        conserved_quantity,
        fuel_regression_model,
        workspace,
        geo,
        fluid_cell_id,
        cell_volume_rate,
    )
    moving_mesh_flux = zero(conserved_quantity[fluid_cell_id] * cell_volume_rate)

    for face_idx in 1:geo.cell_n_faces[fluid_cell_id]
        face_mesh_velocity = calculate_face_mesh_velocity(workspace, geo, fluid_cell_id, face_idx)
        outward_face_normal = geo.cell_face_normals[fluid_cell_id, face_idx]
        face_area = geo.cell_face_areas[fluid_cell_id, face_idx]
        face_swept_volume_rate = dot(face_mesh_velocity, outward_face_normal) * face_area

        neighbor_cell_id = fuel_regression_model.neighbor_cell_id_by_face[fluid_cell_id][face_idx]
        face_has_fluid_neighbor = neighbor_cell_id != 0 && fuel_regression_model.is_fluid_cell[neighbor_cell_id]

        if face_has_fluid_neighbor
            face_conserved_quantity = (conserved_quantity[fluid_cell_id] + conserved_quantity[neighbor_cell_id]) / 2
            moving_mesh_flux += face_conserved_quantity * face_swept_volume_rate
        end
    end

    # The second term accounts for dilution/concentration caused by the local
    # control-volume change and enforces the geometric conservation law.
    control_volume_change = conserved_quantity[fluid_cell_id] * cell_volume_rate
    return moving_mesh_flux - control_volume_change
end

function apply_moving_mesh_conservation!(du, u, fuel_regression_model::FuelRegressionModel, workspace, geo)
    for fluid_cell_id in fuel_regression_model.fluid_cells_affected_by_wall_motion
        representative_mesh_speed = first(workspace.fuel_surface_regression_rates)
        representative_face_area = geo.cell_face_areas[fluid_cell_id, 1]
        cell_volume_rate = zero(representative_mesh_speed * representative_face_area)

        for face_idx in 1:geo.cell_n_faces[fluid_cell_id]
            face_mesh_velocity = calculate_face_mesh_velocity(workspace, geo, fluid_cell_id, face_idx)
            outward_face_normal = geo.cell_face_normals[fluid_cell_id, face_idx]
            face_area = geo.cell_face_areas[fluid_cell_id, face_idx]
            face_swept_volume_rate = dot(face_mesh_velocity, outward_face_normal) * face_area

            cell_volume_rate += face_swept_volume_rate
        end

        inverse_cell_volume = inv(geo.cell_volumes[fluid_cell_id])

        density_mesh_flux = calculate_moving_mesh_flux_for_conserved_quantity(u.density, fuel_regression_model, workspace, geo, fluid_cell_id, cell_volume_rate)
        radial_momentum_mesh_flux = calculate_moving_mesh_flux_for_conserved_quantity(u.momentum_density_u, fuel_regression_model, workspace, geo, fluid_cell_id, cell_volume_rate)
        transverse_momentum_mesh_flux = calculate_moving_mesh_flux_for_conserved_quantity(u.momentum_density_v, fuel_regression_model, workspace, geo, fluid_cell_id, cell_volume_rate)
        axial_momentum_mesh_flux = calculate_moving_mesh_flux_for_conserved_quantity(u.momentum_density_w, fuel_regression_model, workspace, geo, fluid_cell_id, cell_volume_rate)
        total_energy_mesh_flux = calculate_moving_mesh_flux_for_conserved_quantity(u.volumetric_energy, fuel_regression_model, workspace, geo, fluid_cell_id, cell_volume_rate)

        du.density[fluid_cell_id] += density_mesh_flux * inverse_cell_volume
        du.momentum_density_u[fluid_cell_id] += radial_momentum_mesh_flux * inverse_cell_volume
        du.momentum_density_v[fluid_cell_id] += transverse_momentum_mesh_flux * inverse_cell_volume
        du.momentum_density_w[fluid_cell_id] += axial_momentum_mesh_flux * inverse_cell_volume
        du.volumetric_energy[fluid_cell_id] += total_energy_mesh_flux * inverse_cell_volume

        for species_name in propertynames(u.species_densities)
            species_density = getproperty(u.species_densities, species_name)
            species_density_rate = getproperty(du.species_densities, species_name)
            species_density_mesh_flux = calculate_moving_mesh_flux_for_conserved_quantity(species_density, fuel_regression_model, workspace, geo, fluid_cell_id, cell_volume_rate)

            species_density_rate[fluid_cell_id] += species_density_mesh_flux * inverse_cell_volume
        end

        if hasproperty(u, :turbulent_kinetic_energy_density)
            turbulent_kinetic_energy_mesh_flux = calculate_moving_mesh_flux_for_conserved_quantity(u.turbulent_kinetic_energy_density, fuel_regression_model, workspace, geo, fluid_cell_id, cell_volume_rate)
            specific_dissipation_rate_mesh_flux = calculate_moving_mesh_flux_for_conserved_quantity(u.specific_dissipation_rate_density, fuel_regression_model, workspace, geo, fluid_cell_id, cell_volume_rate)

            du.turbulent_kinetic_energy_density[fluid_cell_id] += turbulent_kinetic_energy_mesh_flux * inverse_cell_volume
            du.specific_dissipation_rate_density[fluid_cell_id] += specific_dissipation_rate_mesh_flux * inverse_cell_volume
        end
    end

    return nothing
end


# -----------------------------------------------------------------------------
# Transfer pyrolyzed solid fuel into the port flow
# -----------------------------------------------------------------------------

function inject_pyrolyzed_fuel_into_port!(du, u, p, fuel_regression_model::FuelRegressionModel, workspace, geo)
    for fuel_wall_interface in fuel_regression_model.fuel_wall_interfaces
        port_fluid_cell_id = fuel_wall_interface.port_fluid_cell_id
        solid_fuel_cell_id = fuel_wall_interface.solid_fuel_cell_id
        fluid_face_idx = fuel_wall_interface.fluid_face_idx

        fuel_surface_velocity = calculate_face_mesh_velocity(workspace, geo, port_fluid_cell_id, fluid_face_idx)
        fluid_outward_face_normal = geo.cell_face_normals[port_fluid_cell_id, fluid_face_idx]
        fuel_wall_face_area = geo.cell_face_areas[port_fluid_cell_id, fluid_face_idx]

        volume_of_fuel_consumed_per_time = dot(fuel_surface_velocity, fluid_outward_face_normal) * fuel_wall_face_area
        nonnegative_fuel_volume_flow = max(volume_of_fuel_consumed_per_time, zero(volume_of_fuel_consumed_per_time))
        solid_fuel_density = u.density[solid_fuel_cell_id]
        pyrolyzed_fuel_mass_flow = solid_fuel_density * nonnegative_fuel_volume_flow

        # The present chemistry treats all regressed fuel as gaseous ethylene.
        # The injected momentum follows the receding surface velocity; a more
        # complete model could prescribe a separate wall-normal blowing speed.
        du.density_flow[port_fluid_cell_id] += pyrolyzed_fuel_mass_flow
        du.species_density_flow.ethylene[port_fluid_cell_id] += pyrolyzed_fuel_mass_flow
        du.momentum_density_u_flow[port_fluid_cell_id] += pyrolyzed_fuel_mass_flow * fuel_surface_velocity[1]
        du.momentum_density_v_flow[port_fluid_cell_id] += pyrolyzed_fuel_mass_flow * fuel_surface_velocity[2]
        du.momentum_density_w_flow[port_fluid_cell_id] += pyrolyzed_fuel_mass_flow * fuel_surface_velocity[3]

        solid_fuel_specific_energy = u.volumetric_energy[solid_fuel_cell_id] / solid_fuel_density
        pyrolysis_gas_specific_energy = solid_fuel_specific_energy - p.fuel_heat_of_pyrolysis
        pyrolysis_energy_flow = pyrolyzed_fuel_mass_flow * pyrolysis_gas_specific_energy

        du.volumetric_energy_flow[port_fluid_cell_id] += pyrolysis_energy_flow
    end

    return nothing
end


# -----------------------------------------------------------------------------
# Fuel-regression physics evaluated during one right-hand-side call
# -----------------------------------------------------------------------------

function regress_fuel_layers!(du, u, p, t, system, geo)
    fuel_regression_model = system.fuel_regression_model
    first_fuel_thickness_state_idx = first(fuel_regression_model.fuel_port_stations).fuel_thickness_state_idx
    fuel_thickness_seed = u.fuel_grain_thickness[first_fuel_thickness_state_idx]
    workspace = get_tmp(fuel_regression_model.workspaces, fuel_thickness_seed)

    # 1. Measure oxidizer flow and surface temperature, then evaluate r_dot = a G_ox^n.
    calculate_fuel_regression_rates!(du, u, p, fuel_regression_model, workspace, geo)

    # 2. Convert the local regression rate into a velocity for every moving mesh node.
    calculate_mesh_node_velocities!(fuel_regression_model, workspace)

    # 3. Conserve the existing port-fluid state as the ALE control volumes expand.
    apply_moving_mesh_conservation!(du, u, fuel_regression_model, workspace, geo)

    # 4. Add the solid mass consumed at the wall as gaseous fuel, momentum, and energy.
    inject_pyrolyzed_fuel_into_port!(du, u, p, fuel_regression_model, workspace, geo)

    return nothing
end
