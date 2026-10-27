"""
Fuel-grain regression model for the hybrid rocket port.

Coordinate convention used by this model:

  * `z` is the motor axis and the positive `z` direction is downstream.
  * `x` is the resolved radial direction from the port toward the fuel grain.
  * Each axial mesh plane owns one fuel-thickness state.

The calculation performed during each residual evaluation is:

  1. Evaluate nitrous-oxide mass flux and fuel-surface temperature at every
     axial station.
  2. Apply the empirical regression law wherever the surface is hot enough.
  3. Convert the regression rates into radial mesh-node velocities.
  4. Apply ALE conservation terms to fluid cells deformed by the moving wall.
  5. Inject regressed fuel mass, momentum, and energy into the port flow.

The mesh-description types below are constructed once. `FuelRegressionWorkspace`
holds scalar-type-dependent temporary arrays used during residual evaluations.
"""


# -----------------------------------------------------------------------------
# Physical mesh description
# -----------------------------------------------------------------------------

"""A cell face lying on an axial (`z = constant`) plane."""
struct AxialFace
    node_ids_on_z_plane::Vector{Int}
    cell_id::Int
    face_idx_on_z_plane::Int
    face_idx_facing_toward_fuel_wall::Int
    face_idx_facing_away_from_fuel_wall::Int
    closest_node_ids_from_fuel_wall::Vector{Int}
    furthest_node_ids_from_fuel_wall::Vector{Int}
end

"""
Mesh and fuel-thickness information at one axial station.

The upstream face set is normally used to measure the flow at the station. The
downstream set exists so that the inlet end of the grain, which has no upstream
cell at its own `z` plane, can use the adjacent downstream cells instead.
"""
struct AxialPosition
    z_pos::Float64
    sorted_upstream_fluid_cell_entries::Vector{AxialFace}
    sorted_downstream_fluid_cell_entries::Vector{AxialFace}
    sorted_upstream_fuel_cell_entries::Vector{AxialFace}
    sorted_downstream_fuel_cell_entries::Vector{AxialFace}
    min_fluid_x_bound::Float64
    initial_fuel_grain_x_bound::Float64
    initial_fuel_grain_thickness::Float64
    max_fuel_grain_x_bound::Float64
    node_ids::Vector{Int}
    thickness_state_idx::Int
end

"""
Affine radial-motion data for one mesh node.

`port_radius_multiplier` is one at the initial fuel surface. It decreases
linearly to zero at the fixed port-side and outer-fuel boundaries. This spreads
the fuel-wall displacement through the neighboring mesh without introducing a
separate mesh-deformation solve.
"""
struct FuelRegressionNode
    node_id::Int
    axial_position_id::Int
    initial_x::Float64
    port_radius_multiplier::Float64
end

"""A fluid/fuel face pair on the regressing fuel surface."""
struct FuelWallFace
    fluid_cell_id::Int
    fluid_face_idx::Int
    fuel_cell_id::Int
    fuel_face_idx::Int
end

struct FuelRegressionWorkspace{NV, PR}
    node_velocities::NV
    port_radius_rates::PR
end

struct FuelRegressionModel{W}
    axial_positions::Vector{AxialPosition}
    regression_nodes::Vector{FuelRegressionNode}
    wall_faces::Vector{FuelWallFace}
    affected_fluid_cell_ids::Vector{Int}
    fluid_cell_mask::BitVector
    neighbor_cell_by_face::Vector{Vector{Int}}
    initial_node_coordinates::Vector{SVector{3, Float64}}
    workspaces::W
end


# -----------------------------------------------------------------------------
# One-time interpretation of the fuel and port mesh
# -----------------------------------------------------------------------------

function unique_sorted_coordinates(values; absolute_tolerance)
    sorted_coordinates = sort!(collect(values))
    isempty(sorted_coordinates) && return Float64[]

    unique_coordinates = Float64[first(sorted_coordinates)]
    for coordinate in @view sorted_coordinates[2:end]
        previous_coordinate = last(unique_coordinates)
        if abs(coordinate - previous_coordinate) > absolute_tolerance
            push!(unique_coordinates, coordinate)
        end
    end

    return unique_coordinates
end

@inline function face_node_ids(geo, cell_id, face_idx)
    packed_face_id = geo.cell_face_layout[cell_id, face_idx]
    return geo.cell_faces_node_ids[packed_face_id]
end

function find_axial_position_id(axial_z_positions, node_z; absolute_tolerance)
    position_id = findfirst(axial_z_positions) do axial_z
        abs(axial_z - node_z) <= absolute_tolerance
    end

    isnothing(position_id) && error(
        "No fuel-grain axial position was found at z = $node_z",
    )
    return position_id
end

function find_faces_toward_and_away_from_fuel_wall(geo, cell_id, material::Symbol)
    radial_face_indices = Int[]

    for face_idx in 1:geo.cell_n_faces[cell_id]
        outward_face_normal = geo.cell_face_normals[cell_id, face_idx]
        face_is_radial = abs(outward_face_normal[1]) > 0.5
        if face_is_radial
            push!(radial_face_indices, face_idx)
        end
    end

    length(radial_face_indices) == 2 || error(
        "Expected two radial faces for cell $cell_id, found $(length(radial_face_indices))",
    )

    sort!(
        radial_face_indices;
        by = face_idx -> geo.cell_face_centroids[cell_id, face_idx][1],
    )
    lower_x_face_idx = first(radial_face_indices)
    higher_x_face_idx = last(radial_face_indices)

    # The fuel surface is the high-x face of a port-fluid cell and the low-x
    # face of a fuel cell.
    if material === :fluid
        return higher_x_face_idx, lower_x_face_idx
    elseif material === :fuel
        return lower_x_face_idx, higher_x_face_idx
    end

    error("Unknown cell material: $material")
end

function describe_axial_face(geo, cell_id, axial_face_idx, material::Symbol)
    face_idx_facing_toward_fuel_wall,
    face_idx_facing_away_from_fuel_wall =
        find_faces_toward_and_away_from_fuel_wall(geo, cell_id, material)

    return AxialFace(
        copy(face_node_ids(geo, cell_id, axial_face_idx)),
        cell_id,
        axial_face_idx,
        face_idx_facing_toward_fuel_wall,
        face_idx_facing_away_from_fuel_wall,
        copy(face_node_ids(geo, cell_id, face_idx_facing_toward_fuel_wall)),
        copy(face_node_ids(geo, cell_id, face_idx_facing_away_from_fuel_wall)),
    )
end

function collect_cells_touching_axial_plane(
    geo,
    candidate_cell_ids,
    axial_z,
    side_of_plane::Symbol,
    material::Symbol;
    absolute_tolerance,
)
    axial_faces = AxialFace[]

    for cell_id in candidate_cell_ids
        cell_centroid_z = geo.cell_centroids[cell_id][3]
        cell_is_on_requested_side = if side_of_plane === :upstream
            cell_centroid_z < axial_z - absolute_tolerance
        elseif side_of_plane === :downstream
            cell_centroid_z > axial_z + absolute_tolerance
        else
            error("Unknown side of axial plane: $side_of_plane")
        end

        if !cell_is_on_requested_side
            continue
        end

        face_idx_on_axial_plane = nothing
        for candidate_face_idx in 1:geo.cell_n_faces[cell_id]
            candidate_face_node_ids = face_node_ids(geo, cell_id, candidate_face_idx)
            all_nodes_are_on_axial_plane = all(candidate_face_node_ids) do node_id
                node_z = geo.node_coordinates[node_id][3]
                abs(node_z - axial_z) <= absolute_tolerance
            end

            if all_nodes_are_on_axial_plane
                face_idx_on_axial_plane = candidate_face_idx
                break
            end
        end

        if !isnothing(face_idx_on_axial_plane)
            push!(
                axial_faces,
                describe_axial_face(geo, cell_id, face_idx_on_axial_plane, material),
            )
        end
    end

    # Radial ordering lets the wall-temperature calculation select the fluid
    # and fuel cells immediately adjacent to the fuel surface.
    sort!(axial_faces; by = face -> geo.cell_centroids[face.cell_id][1])
    return axial_faces
end

function build_neighbor_cell_by_face(geo)
    neighbor_cell_by_face = [
        zeros(Int, geo.cell_n_faces[cell_id]) for cell_id in 1:geo.n_cells
    ]

    for (cell_id, neighbors) in geo.cell_neighbors
        for (neighbor_cell_id, face_idx, _) in neighbors
            neighbor_cell_by_face[cell_id][face_idx] = neighbor_cell_id
        end
    end

    return neighbor_cell_by_face
end

function find_fuel_wall_faces(geo, port_fluid_cell_ids, fuel_cell_ids)
    fuel_cell_set = Set(fuel_cell_ids)
    fuel_wall_faces = FuelWallFace[]

    for fluid_cell_id in port_fluid_cell_ids
        for (
            neighbor_cell_id,
            fluid_face_idx,
            neighbor_face_idx,
        ) in geo.cell_neighbors[fluid_cell_id][2]
            if neighbor_cell_id in fuel_cell_set
                push!(
                    fuel_wall_faces,
                    FuelWallFace(
                        fluid_cell_id,
                        fluid_face_idx,
                        neighbor_cell_id,
                        neighbor_face_idx,
                    ),
                )
            end
        end
    end

    return fuel_wall_faces
end

"""
    create_fuel_regression_model(grid, geo; atol = 1e-10)

Interpret the port and fuel-grain cell sets as a sequence of axial fuel
stations. The first `length(model.axial_positions)` entries of
`u.fuel_grain_thickness` are used as the remaining fuel thicknesses.

Mesh motion is prescribed rather than solved. At every axial station, the fuel
surface moves by the consumed thickness. Motion decreases linearly to zero at
the inner fluid boundary and at the outer fuel boundary.
"""
function create_fuel_regression_model(grid, geo; atol = 1e-10)
    fuel_cell_ids = sort!(collect(getcellset(grid, "fuel_grain_volume")))
    all_fluid_cell_ids = sort!(collect(getcellset(grid, "fluid_domain")))
    isempty(fuel_cell_ids) && error("The fuel_grain_volume cell set is empty")

    # Restrict the moving fluid region to axial segments occupied by fuel. This
    # prevents injector and nozzle fluid cells from being included merely
    # because they belong to the global fluid-domain cell set.
    fuel_segment_centroid_z_positions = unique_sorted_coordinates(
        (geo.cell_centroids[cell_id][3] for cell_id in fuel_cell_ids);
        absolute_tolerance = atol,
    )
    port_fluid_cell_ids = [
        cell_id for cell_id in all_fluid_cell_ids if any(
            fuel_segment_centroid_z_positions,
        ) do fuel_segment_z
            fluid_cell_z = geo.cell_centroids[cell_id][3]
            abs(fluid_cell_z - fuel_segment_z) <= atol
        end
    ]

    fuel_node_ids = unique!(sort!(reduce(
        vcat,
        (geo.cell_node_ids[cell_id] for cell_id in fuel_cell_ids),
    )))
    axial_z_positions = unique_sorted_coordinates(
        (geo.node_coordinates[node_id][3] for node_id in fuel_node_ids);
        absolute_tolerance = atol,
    )

    # Gather the fluid and fuel nodes belonging to each axial plane. The union
    # is the portion of the mesh deformed by that plane's thickness state.
    all_node_ids_at_axial_position = [Set{Int}() for _ in axial_z_positions]
    fluid_node_ids_at_axial_position = [Set{Int}() for _ in axial_z_positions]
    fuel_node_ids_at_axial_position = [Set{Int}() for _ in axial_z_positions]

    for fluid_cell_id in port_fluid_cell_ids
        for node_id in geo.cell_node_ids[fluid_cell_id]
            node_z = geo.node_coordinates[node_id][3]
            axial_position_id = find_axial_position_id(
                axial_z_positions,
                node_z;
                absolute_tolerance = atol,
            )
            push!(all_node_ids_at_axial_position[axial_position_id], node_id)
            push!(fluid_node_ids_at_axial_position[axial_position_id], node_id)
        end
    end

    for fuel_cell_id in fuel_cell_ids
        for node_id in geo.cell_node_ids[fuel_cell_id]
            node_z = geo.node_coordinates[node_id][3]
            axial_position_id = find_axial_position_id(
                axial_z_positions,
                node_z;
                absolute_tolerance = atol,
            )
            push!(all_node_ids_at_axial_position[axial_position_id], node_id)
            push!(fuel_node_ids_at_axial_position[axial_position_id], node_id)
        end
    end

    axial_positions = AxialPosition[]
    regression_nodes = FuelRegressionNode[]

    for (axial_position_id, axial_z) in enumerate(axial_z_positions)
        fluid_node_ids = collect(
            fluid_node_ids_at_axial_position[axial_position_id],
        )
        fuel_node_ids_on_plane = collect(
            fuel_node_ids_at_axial_position[axial_position_id],
        )
        isempty(fluid_node_ids) && error("No fluid nodes were found at z = $axial_z")
        isempty(fuel_node_ids_on_plane) && error(
            "No fuel nodes were found at z = $axial_z",
        )

        minimum_port_fluid_x = minimum(
            node_id -> geo.node_coordinates[node_id][1],
            fluid_node_ids,
        )
        initial_fuel_surface_x = maximum(
            node_id -> geo.node_coordinates[node_id][1],
            fluid_node_ids,
        )
        fixed_outer_fuel_x = maximum(
            node_id -> geo.node_coordinates[node_id][1],
            fuel_node_ids_on_plane,
        )
        initial_fuel_grain_thickness =
            fixed_outer_fuel_x - initial_fuel_surface_x

        upstream_fluid_faces = collect_cells_touching_axial_plane(
            geo,
            port_fluid_cell_ids,
            axial_z,
            :upstream,
            :fluid;
            absolute_tolerance = atol,
        )
        downstream_fluid_faces = collect_cells_touching_axial_plane(
            geo,
            port_fluid_cell_ids,
            axial_z,
            :downstream,
            :fluid;
            absolute_tolerance = atol,
        )
        upstream_fuel_faces = collect_cells_touching_axial_plane(
            geo,
            fuel_cell_ids,
            axial_z,
            :upstream,
            :fuel;
            absolute_tolerance = atol,
        )
        downstream_fuel_faces = collect_cells_touching_axial_plane(
            geo,
            fuel_cell_ids,
            axial_z,
            :downstream,
            :fuel;
            absolute_tolerance = atol,
        )

        node_ids_on_axial_plane = sort!(collect(
            all_node_ids_at_axial_position[axial_position_id],
        ))
        push!(
            axial_positions,
            AxialPosition(
                axial_z,
                upstream_fluid_faces,
                downstream_fluid_faces,
                upstream_fuel_faces,
                downstream_fuel_faces,
                minimum_port_fluid_x,
                initial_fuel_surface_x,
                initial_fuel_grain_thickness,
                fixed_outer_fuel_x,
                node_ids_on_axial_plane,
                axial_position_id,
            ),
        )

        # Distribute the fuel-surface displacement linearly into both adjacent
        # materials. Nodes at the fuel surface move with the full regression
        # displacement; nodes at either fixed radial boundary do not move.
        for node_id in node_ids_on_axial_plane
            initial_node_x = geo.node_coordinates[node_id][1]

            if initial_node_x <= initial_fuel_surface_x + atol
                distance_from_fixed_fluid_boundary =
                    initial_node_x - minimum_port_fluid_x
                initial_fluid_region_width =
                    initial_fuel_surface_x - minimum_port_fluid_x
                port_radius_multiplier =
                    distance_from_fixed_fluid_boundary / initial_fluid_region_width
            else
                distance_into_fuel = initial_node_x - initial_fuel_surface_x
                fraction_through_fuel =
                    distance_into_fuel / initial_fuel_grain_thickness
                port_radius_multiplier = 1 - fraction_through_fuel
            end

            push!(
                regression_nodes,
                FuelRegressionNode(
                    node_id,
                    axial_position_id,
                    initial_node_x,
                    port_radius_multiplier,
                ),
            )
        end
    end

    moving_node_ids = Set(
        node.node_id for node in regression_nodes
        if abs(node.port_radius_multiplier) > atol
    )
    affected_fluid_cell_ids = [
        fluid_cell_id for fluid_cell_id in all_fluid_cell_ids if any(
            node_id -> node_id in moving_node_ids,
            geo.cell_node_ids[fluid_cell_id],
        )
    ]

    fluid_cell_mask = falses(geo.n_cells)
    fluid_cell_mask[all_fluid_cell_ids] .= true

    initial_node_coordinates = SVector{3, Float64}.(geo.node_coordinates)
    number_of_mesh_nodes = length(geo.node_coordinates)
    number_of_axial_positions = length(axial_positions)

    workspaces = GeneralLazyBufferCache() do scalar_seed
        scalar_zero = zero(scalar_seed)
        stationary_node_velocity = SVector(
            scalar_zero,
            scalar_zero,
            scalar_zero,
        )
        FuelRegressionWorkspace(
            fill(stationary_node_velocity, number_of_mesh_nodes),
            fill(scalar_zero, number_of_axial_positions),
        )
    end

    fuel_wall_faces = find_fuel_wall_faces(
        geo,
        port_fluid_cell_ids,
        fuel_cell_ids,
    )
    neighbor_cell_by_face = build_neighbor_cell_by_face(geo)

    return FuelRegressionModel(
        axial_positions,
        regression_nodes,
        fuel_wall_faces,
        affected_fluid_cell_ids,
        fluid_cell_mask,
        neighbor_cell_by_face,
        initial_node_coordinates,
        workspaces,
    )
end


# -----------------------------------------------------------------------------
# Fuel-thickness state and moving geometry
# -----------------------------------------------------------------------------

function initialize_fuel_grain_thickness!(
    fuel_grain_thickness,
    model::FuelRegressionModel,
)
    thickness_unit = oneunit(first(fuel_grain_thickness))

    for axial_position in model.axial_positions
        thickness_state_idx = axial_position.thickness_state_idx
        initial_thickness = axial_position.initial_fuel_grain_thickness
        fuel_grain_thickness[thickness_state_idx] =
            initial_thickness * thickness_unit
    end

    return fuel_grain_thickness
end

function move_fuel_wall!(u, model::FuelRegressionModel, geo)
    # Always reconstruct the current mesh from its initial coordinates. This
    # makes displacement a direct function of the thickness state and avoids
    # accumulating mesh-motion error across residual evaluations.
    copyto!(geo.node_coordinates, model.initial_node_coordinates)

    for regression_node in model.regression_nodes
        axial_position = model.axial_positions[regression_node.axial_position_id]

        remaining_fuel_thickness = u.fuel_grain_thickness[axial_position.thickness_state_idx]

        consumed_fuel_thickness = axial_position.initial_fuel_grain_thickness - remaining_fuel_thickness

        radial_node_displacement = regression_node.port_radius_multiplier * consumed_fuel_thickness

        current_node_x = regression_node.initial_x + radial_node_displacement

        initial_coordinate = geo.node_coordinates[regression_node.node_id]

        geo.node_coordinates[regression_node.node_id] = setindex(initial_coordinate, current_node_x, 1)
    end

    return nothing
end

function update_fuel_regression_geometry!(u, model::FuelRegressionModel, geo, grid)
    move_fuel_wall!(u, model, geo)
    update_fvm_geometry!(geo, geo.node_coordinates)

    # Reapply the hybrid-motor-specific geometry corrections after the generic
    # finite-volume geometry has been rebuilt from the moved nodes.
    correct_hybrid_geo!(geo, grid)
    return nothing
end


# -----------------------------------------------------------------------------
# Regression-driving oxidizer flux and fuel-surface temperature
# -----------------------------------------------------------------------------

function calculate_axial_nitrous_oxide_mass_flux(u, geo, fluid_axial_faces)
    nitrous_oxide_mass_flow_rate = 0.0
    port_cross_sectional_area = 0.0

    for axial_face in fluid_axial_faces
        fluid_cell_id = axial_face.cell_id
        axial_face_idx = axial_face.face_idx_on_z_plane

        mixture_density = u.density[fluid_cell_id]
        axial_momentum_density = u.momentum_density_w[fluid_cell_id]
        axial_velocity = if iszero(mixture_density)
            zero(mixture_density)
        else
            axial_momentum_density / mixture_density
        end

        nitrous_oxide_density =
            u.species_densities.nitrous_oxide[fluid_cell_id]
        axial_face_area = geo.cell_face_areas[fluid_cell_id, axial_face_idx]
        nitrous_oxide_mass_flow_through_face =
            abs(axial_velocity) * nitrous_oxide_density * axial_face_area

        nitrous_oxide_mass_flow_rate += nitrous_oxide_mass_flow_through_face
        port_cross_sectional_area += axial_face_area
    end

    if iszero(port_cross_sectional_area)
        return zero(nitrous_oxide_mass_flow_rate)
    end

    return nitrous_oxide_mass_flow_rate / port_cross_sectional_area
end

function calculate_fuel_grain_face_temperature(
    u,
    geo,
    fluid_cell_nearest_fuel_wall::AxialFace,
    fuel_cell_nearest_port::AxialFace,
)
    fluid_cell_id = fluid_cell_nearest_fuel_wall.cell_id
    fuel_cell_id = fuel_cell_nearest_port.cell_id

    fluid_face_idx =
        fluid_cell_nearest_fuel_wall.face_idx_facing_toward_fuel_wall
    fuel_face_idx = fuel_cell_nearest_port.face_idx_facing_toward_fuel_wall
    fluid_centroid_to_wall_distance =
        geo.cell_face_distances[fluid_cell_id, fluid_face_idx]
    fuel_centroid_to_wall_distance =
        geo.cell_face_distances[fuel_cell_id, fuel_face_idx]

    # Treat the fluid-side and fuel-side half cells as two one-dimensional
    # thermal resistances meeting at a common interface temperature. Radiation
    # and a separate wall heat-transfer correlation are not included here.
    fluid_thermal_conductance_per_area =
        u.k[fluid_cell_id] / fluid_centroid_to_wall_distance
    fuel_thermal_conductance_per_area =
        u.k[fuel_cell_id] / fuel_centroid_to_wall_distance
    total_thermal_conductance_per_area =
        fluid_thermal_conductance_per_area + fuel_thermal_conductance_per_area

    if iszero(total_thermal_conductance_per_area)
        return (u.temperature[fluid_cell_id] + u.temperature[fuel_cell_id]) / 2
    end

    fluid_temperature_contribution =
        fluid_thermal_conductance_per_area * u.temperature[fluid_cell_id]
    fuel_temperature_contribution =
        fuel_thermal_conductance_per_area * u.temperature[fuel_cell_id]
    return (
        fluid_temperature_contribution + fuel_temperature_contribution
    ) / total_thermal_conductance_per_area
end

function get_axial_oxidizer_mass_flux_and_fuel_grain_face_temperature(
    u,
    geo,
    axial_position::AxialPosition,
)
    # Prefer the cells immediately upstream of the axial plane so mass flux is
    # sampled consistently. At the grain inlet, use the downstream cells
    # because there is no upstream fuel/port segment on that plane.
    use_downstream_cells =
        isempty(axial_position.sorted_upstream_fluid_cell_entries)

    fluid_axial_faces = if use_downstream_cells
        axial_position.sorted_downstream_fluid_cell_entries
    else
        axial_position.sorted_upstream_fluid_cell_entries
    end
    fuel_axial_faces = if use_downstream_cells
        axial_position.sorted_downstream_fuel_cell_entries
    else
        axial_position.sorted_upstream_fuel_cell_entries
    end

    if isempty(fluid_axial_faces) || isempty(fuel_axial_faces)
        return 0.0, 0.0
    end

    nitrous_oxide_mass_flux = calculate_axial_nitrous_oxide_mass_flux(
        u,
        geo,
        fluid_axial_faces,
    )

    # Faces are sorted from low x to high x. Therefore the last fluid cell and
    # first fuel cell are the pair touching the fuel surface.
    fluid_cell_nearest_fuel_wall = last(fluid_axial_faces)
    fuel_cell_nearest_port = first(fuel_axial_faces)
    fuel_grain_face_temperature = calculate_fuel_grain_face_temperature(
        u,
        geo,
        fluid_cell_nearest_fuel_wall,
        fuel_cell_nearest_port,
    )

    return nitrous_oxide_mass_flux, fuel_grain_face_temperature
end

function calculate_fuel_regression_rates!(
    du,
    u,
    p,
    t,
    system,
    model::FuelRegressionModel,
    workspace,
    geo,
)
    for (axial_position_id, axial_position) in enumerate(model.axial_positions)
        nitrous_oxide_mass_flux,
        fuel_grain_face_temperature =
            get_axial_oxidizer_mass_flux_and_fuel_grain_face_temperature(
                u,
                geo,
                axial_position,
            )

        nonnegative_nitrous_oxide_mass_flux = max(
            nitrous_oxide_mass_flux,
            zero(nitrous_oxide_mass_flux),
        )

        # Doran et al. correlation: r_dot = a * G_ox^n. Regression is disabled
        # below the specified surface-temperature threshold so cold fuel does
        # not release ethylene before ignition.
        fuel_is_regressing =
            fuel_grain_face_temperature > p.min_fuel_regression_temperature
        fuel_surface_regression_rate = if fuel_is_regressing
            p.regression_coefficient_a *
            nonnegative_nitrous_oxide_mass_flux^p.regression_exponent_n
        else
            zero(nonnegative_nitrous_oxide_mass_flux)
        end

        thickness_state_idx = axial_position.thickness_state_idx
        du.fuel_grain_thickness[thickness_state_idx] =
            -fuel_surface_regression_rate
        workspace.port_radius_rates[axial_position_id] =
            fuel_surface_regression_rate
    end

    return nothing
end


# -----------------------------------------------------------------------------
# Mesh velocity induced by fuel regression
# -----------------------------------------------------------------------------

function calculate_regression_mesh_node_velocities!(
    model::FuelRegressionModel,
    workspace,
)
    stationary_node_velocity = zero(eltype(workspace.node_velocities))
    fill!(workspace.node_velocities, stationary_node_velocity)

    for regression_node in model.regression_nodes
        local_fuel_surface_regression_rate = workspace.port_radius_rates[
            regression_node.axial_position_id
        ]
        radial_mesh_velocity =
            regression_node.port_radius_multiplier *
            local_fuel_surface_regression_rate

        workspace.node_velocities[regression_node.node_id] = SVector(
            radial_mesh_velocity,
            zero(radial_mesh_velocity),
            zero(radial_mesh_velocity),
        )
    end

    return nothing
end

@inline function calculate_face_mesh_velocity(
    model::FuelRegressionModel,
    workspace,
    geo,
    cell_id,
    face_idx,
)
    node_ids_on_face = face_node_ids(geo, cell_id, face_idx)
    face_mesh_velocity = zero(eltype(workspace.node_velocities))

    for node_id in node_ids_on_face
        face_mesh_velocity += workspace.node_velocities[node_id]
    end

    return face_mesh_velocity / length(node_ids_on_face)
end

function calculate_cell_volume_rate(
    model::FuelRegressionModel,
    workspace,
    geo,
    cell_id,
)
    cell_volume_rate = 0.0

    for face_idx in 1:geo.cell_n_faces[cell_id]
        face_mesh_velocity = calculate_face_mesh_velocity(
            model,
            workspace,
            geo,
            cell_id,
            face_idx,
        )
        outward_face_normal = geo.cell_face_normals[cell_id, face_idx]
        face_area = geo.cell_face_areas[cell_id, face_idx]
        outward_swept_volume_rate =
            dot(face_mesh_velocity, outward_face_normal) * face_area
        cell_volume_rate += outward_swept_volume_rate
    end

    return cell_volume_rate
end


# -----------------------------------------------------------------------------
# Arbitrary-Lagrangian-Eulerian conservation for the deforming port mesh
# -----------------------------------------------------------------------------

function calculate_ale_transport_for_conserved_density(
    conserved_density,
    model::FuelRegressionModel,
    workspace,
    geo,
    fluid_cell_id,
    cell_volume_rate,
)
    conserved_quantity_flow_due_to_mesh_motion = zero(
        conserved_density[fluid_cell_id] * cell_volume_rate,
    )

    for face_idx in 1:geo.cell_n_faces[fluid_cell_id]
        face_mesh_velocity = calculate_face_mesh_velocity(
            model,
            workspace,
            geo,
            fluid_cell_id,
            face_idx,
        )
        outward_face_normal = geo.cell_face_normals[fluid_cell_id, face_idx]
        face_area = geo.cell_face_areas[fluid_cell_id, face_idx]
        outward_swept_volume_rate =
            dot(face_mesh_velocity, outward_face_normal) * face_area

        neighbor_cell_id =
            model.neighbor_cell_by_face[fluid_cell_id][face_idx]
        face_has_fluid_neighbor =
            neighbor_cell_id != 0 && model.fluid_cell_mask[neighbor_cell_id]

        if face_has_fluid_neighbor
            face_conserved_density = (
                conserved_density[fluid_cell_id] +
                conserved_density[neighbor_cell_id]
            ) / 2
            conserved_quantity_flow_due_to_mesh_motion +=
                face_conserved_density * outward_swept_volume_rate
        end
    end

    # This geometric-conservation term accounts for the changing denominator
    # when the solver stores a conserved quantity per unit cell volume.
    changing_cell_volume_correction =
        conserved_density[fluid_cell_id] * cell_volume_rate
    return (
        conserved_quantity_flow_due_to_mesh_motion -
        changing_cell_volume_correction
    )
end

function apply_fuel_regression_ale!(
    du,
    u,
    model::FuelRegressionModel,
    workspace,
    geo,
)
    for fluid_cell_id in model.affected_fluid_cell_ids
        cell_volume_rate = calculate_cell_volume_rate(
            model,
            workspace,
            geo,
            fluid_cell_id,
        )
        inverse_cell_volume = inv(geo.cell_volumes[fluid_cell_id])

        density_ale_transport = calculate_ale_transport_for_conserved_density(
            u.density,
            model,
            workspace,
            geo,
            fluid_cell_id,
            cell_volume_rate,
        )
        axial_momentum_ale_transport =
            calculate_ale_transport_for_conserved_density(
                u.momentum_density_w,
                model,
                workspace,
                geo,
                fluid_cell_id,
                cell_volume_rate,
            )
        radial_momentum_ale_transport =
            calculate_ale_transport_for_conserved_density(
                u.momentum_density_u,
                model,
                workspace,
                geo,
                fluid_cell_id,
                cell_volume_rate,
            )
        transverse_momentum_ale_transport =
            calculate_ale_transport_for_conserved_density(
                u.momentum_density_v,
                model,
                workspace,
                geo,
                fluid_cell_id,
                cell_volume_rate,
            )
        energy_ale_transport = calculate_ale_transport_for_conserved_density(
            u.volumetric_energy,
            model,
            workspace,
            geo,
            fluid_cell_id,
            cell_volume_rate,
        )

        du.density[fluid_cell_id] +=
            density_ale_transport * inverse_cell_volume
        du.momentum_density_u[fluid_cell_id] +=
            radial_momentum_ale_transport * inverse_cell_volume
        du.momentum_density_v[fluid_cell_id] +=
            transverse_momentum_ale_transport * inverse_cell_volume
        du.momentum_density_w[fluid_cell_id] +=
            axial_momentum_ale_transport * inverse_cell_volume
        du.volumetric_energy[fluid_cell_id] +=
            energy_ale_transport * inverse_cell_volume

        # Every species density experiences the same mesh motion. The generic
        # field loop removes bookkeeping repetition without hiding which
        # physical family of conserved quantities is being transported.
        for_fields!(
            du.species_densities,
            u.species_densities,
        ) do species, du_species_densities, species_densities
            species_ale_transport =
                calculate_ale_transport_for_conserved_density(
                    species_densities[species],
                    model,
                    workspace,
                    geo,
                    fluid_cell_id,
                    cell_volume_rate,
                )
            du_species_densities[species[fluid_cell_id]] +=
                species_ale_transport * inverse_cell_volume
        end

        if hasproperty(u, :turbulent_kinetic_energy_density)
            turbulent_kinetic_energy_ale_transport =
                calculate_ale_transport_for_conserved_density(
                    u.turbulent_kinetic_energy_density,
                    model,
                    workspace,
                    geo,
                    fluid_cell_id,
                    cell_volume_rate,
                )
            specific_dissipation_rate_ale_transport =
                calculate_ale_transport_for_conserved_density(
                    u.specific_dissipation_rate_density,
                    model,
                    workspace,
                    geo,
                    fluid_cell_id,
                    cell_volume_rate,
                )

            du.turbulent_kinetic_energy_density[fluid_cell_id] +=
                turbulent_kinetic_energy_ale_transport * inverse_cell_volume
            du.specific_dissipation_rate_density[fluid_cell_id] +=
                specific_dissipation_rate_ale_transport * inverse_cell_volume
        end
    end

    return nothing
end


# -----------------------------------------------------------------------------
# Fuel mass, momentum, and energy released at the regressing wall
# -----------------------------------------------------------------------------

function apply_fuel_regression_sources!(
    du,
    u,
    p,
    t,
    system,
    model::FuelRegressionModel,
    workspace,
    geo,
)
    for fuel_wall_face in model.wall_faces
        fluid_cell_id = fuel_wall_face.fluid_cell_id
        fuel_cell_id = fuel_wall_face.fuel_cell_id
        fluid_face_idx = fuel_wall_face.fluid_face_idx

        fuel_surface_velocity = calculate_face_mesh_velocity(
            model,
            workspace,
            geo,
            fluid_cell_id,
            fluid_face_idx,
        )
        outward_fluid_face_normal =
            geo.cell_face_normals[fluid_cell_id, fluid_face_idx]
        fuel_wall_face_area =
            geo.cell_face_areas[fluid_cell_id, fluid_face_idx]
        volume_rate_opened_by_regressing_fuel =
            dot(fuel_surface_velocity, outward_fluid_face_normal) *
            fuel_wall_face_area

        # Only wall motion that enlarges the fluid cell releases fuel into it.
        # The fuel is represented as ethylene in the gas-phase species model.
        positive_released_fuel_volume_rate = max(
            volume_rate_opened_by_regressing_fuel,
            zero(volume_rate_opened_by_regressing_fuel),
        )
        solid_fuel_density = u.density[fuel_cell_id]
        pyrolyzed_fuel_mass_flow_rate =
            solid_fuel_density * positive_released_fuel_volume_rate

        du.density_flow[fluid_cell_id] += pyrolyzed_fuel_mass_flow_rate
        du.species_density_flow.ethylene[fluid_cell_id] +=
            pyrolyzed_fuel_mass_flow_rate

        radial_fuel_momentum_flow =
            pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[1]
        transverse_fuel_momentum_flow =
            pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[2]
        axial_fuel_momentum_flow =
            pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[3]
        du.momentum_density_u_flow[fluid_cell_id] +=
            radial_fuel_momentum_flow
        du.momentum_density_v_flow[fluid_cell_id] +=
            transverse_fuel_momentum_flow
        du.momentum_density_w_flow[fluid_cell_id] +=
            axial_fuel_momentum_flow

        solid_fuel_specific_energy =
            u.volumetric_energy[fuel_cell_id] / solid_fuel_density
        injected_fuel_specific_energy =
            solid_fuel_specific_energy - p.fuel_heat_of_pyrolysis
        pyrolyzed_fuel_energy_flow =
            pyrolyzed_fuel_mass_flow_rate * injected_fuel_specific_energy
        du.volumetric_energy_flow[fluid_cell_id] +=
            pyrolyzed_fuel_energy_flow
    end

    return nothing
end

solid_fuel_specific_energy = u.volumetric_energy[fuel_cell_id] / solid_fuel_density

injected_fuel_specific_energy = solid_fuel_specific_energy - p.fuel_heat_of_pyrolysis

pyrolyzed_fuel_energy_flow = pyrolyzed_fuel_mass_flow_rate * injected_fuel_specific_energy

du.volumetric_energy_flow[fluid_cell_id] += pyrolyzed_fuel_energy_flow


# -----------------------------------------------------------------------------
# Fuel-regression sequence called by the coupled rocket residual
# -----------------------------------------------------------------------------

function regress_fuel_layers!(du, u, p, t, system, geo)
    fuel_regression_model = system.fuel_regression_model
    first_axial_position = first(fuel_regression_model.axial_positions)
    first_thickness_state_idx = first_axial_position.thickness_state_idx
    scalar_type_seed = u.fuel_grain_thickness[first_thickness_state_idx]
    fuel_regression_workspace = get_tmp(
        fuel_regression_model.workspaces,
        scalar_type_seed,
    )

    calculate_fuel_regression_rates!(
        du,
        u,
        p,
        t,
        system,
        fuel_regression_model,
        fuel_regression_workspace,
        geo,
    )
    calculate_regression_mesh_node_velocities!(
        fuel_regression_model,
        fuel_regression_workspace,
    )
    apply_fuel_regression_ale!(
        du,
        u,
        fuel_regression_model,
        fuel_regression_workspace,
        geo,
    )
    apply_fuel_regression_sources!(
        du,
        u,
        p,
        t,
        system,
        fuel_regression_model,
        fuel_regression_workspace,
        geo,
    )

    return nothing
end
