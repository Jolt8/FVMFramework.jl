"""
Fuel-grain regression model for the hybrid rocket port.

The model uses `z` as the motor axis and `x` as the resolved radial direction.
Each axial mesh plane has one fuel-thickness state. During a residual evaluation
the model:

  1. calculates oxidizer mass flux and fuel-surface temperature;
  2. evaluates the empirical fuel-regression law;
  3. converts fuel regression into mesh-node velocities;
  4. applies ALE conservation terms in the deforming fluid cells; and
  5. injects pyrolyzed fuel mass, momentum, and energy into the port flow.

The mesh geometry is described once by `FuelRegressionModel`. Temporary arrays
whose scalar type must follow the ODE state are stored in
`FuelRegressionWorkspace`.
"""


# -----------------------------------------------------------------------------
# Fuel-port geometry
# -----------------------------------------------------------------------------

"""A fluid or fuel cell face that lies on an axial (`z = constant`) plane."""
struct AxialFace
    node_ids_on_z_plane::Vector{Int}
    cell_id::Int
    face_idx_on_z_plane::Int
    face_idx_facing_toward_fuel_wall::Int
    face_idx_facing_away_from_fuel_wall::Int
    closest_node_ids_from_fuel_wall::Vector{Int}
    furthest_node_ids_from_fuel_wall::Vector{Int}
end

"""The port, fuel, and thickness-state description at one axial station."""
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
Prescribed radial motion for a node affected by fuel regression.

`port_radius_multiplier` is one at the fuel surface and falls linearly to zero
at the fixed inner-fluid and outer-fuel boundaries.
"""
struct FuelRegressionNode
    node_id::Int
    axial_position_id::Int
    initial_x::Float64
    port_radius_multiplier::Float64
end

"""A coincident fluid/fuel face pair on the regressing fuel surface."""
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
# One-time interpretation of the CFD mesh
# -----------------------------------------------------------------------------

function _unique_sorted_with_tolerance(values; atol)
    sorted_values = sort!(collect(values))
    isempty(sorted_values) && return Float64[]

    unique_values = Float64[first(sorted_values)]

    for value in @view sorted_values[2:end]
        value_is_distinct = abs(value - last(unique_values)) > atol

        if value_is_distinct
            push!(unique_values, value)
        end
    end

    return unique_values
end

@inline function _face_node_ids(geo, cell_id, face_idx)
    packed_face_id = geo.cell_face_layout[cell_id, face_idx]
    return geo.cell_faces_node_ids[packed_face_id]
end

function _find_axial_position_id(axial_z_positions, node_z; atol)
    axial_position_id = findfirst(z_position -> abs(z_position - node_z) <= atol, axial_z_positions)

    isnothing(axial_position_id) && error("No fuel-grain axial position was found at z = $node_z")

    return axial_position_id
end

function _find_faces_toward_and_away_from_fuel_wall(geo, cell_id, material::Symbol)
    radial_face_indices = Int[]

    for face_idx in 1:geo.cell_n_faces[cell_id]
        outward_face_normal = geo.cell_face_normals[cell_id, face_idx]
        face_is_radial = abs(outward_face_normal[1]) > 0.5

        if face_is_radial
            push!(radial_face_indices, face_idx)
        end
    end

    length(radial_face_indices) == 2 || error( #UGH, I don't like this syntax, explicit if statements are always better
        "Expected two radial faces for cell $cell_id, found $(length(radial_face_indices))",
    )

    sort!(radial_face_indices; by = face_idx -> geo.cell_face_centroids[cell_id, face_idx][1])

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

function _describe_axial_face(geo, cell_id, axial_face_idx, material::Symbol)
    toward_fuel_wall_face_idx, away_from_fuel_wall_face_idx =
        _find_faces_toward_and_away_from_fuel_wall(geo, cell_id, material)

    return AxialFace(
        copy(_face_node_ids(geo, cell_id, axial_face_idx)),
        cell_id,
        axial_face_idx,
        toward_fuel_wall_face_idx,
        away_from_fuel_wall_face_idx,
        copy(_face_node_ids(geo, cell_id, toward_fuel_wall_face_idx)),
        copy(_face_node_ids(geo, cell_id, away_from_fuel_wall_face_idx)),
    )
end

function _collect_axial_faces(geo, cell_ids, axial_z, side_of_plane::Symbol, material::Symbol; atol)
    axial_faces = AxialFace[]

    for cell_id in cell_ids
        cell_centroid_z = geo.cell_centroids[cell_id][3]

        if side_of_plane === :upstream
            cell_is_on_requested_side = cell_centroid_z < axial_z - atol
        elseif side_of_plane === :downstream
            cell_is_on_requested_side = cell_centroid_z > axial_z + atol
        else
            error("Unknown side of axial plane: $side_of_plane")
        end

        cell_is_on_requested_side || continue

        face_idx_on_axial_plane = nothing

        for candidate_face_idx in 1:geo.cell_n_faces[cell_id]
            candidate_face_node_ids = _face_node_ids(geo, cell_id, candidate_face_idx)
            all_nodes_are_on_axial_plane = all(
                node_id -> abs(geo.node_coordinates[node_id][3] - axial_z) <= atol,
                candidate_face_node_ids,
            )

            if all_nodes_are_on_axial_plane
                face_idx_on_axial_plane = candidate_face_idx
                break
            end
        end

        isnothing(face_idx_on_axial_plane) && continue

        push!(axial_faces, _describe_axial_face(geo, cell_id, face_idx_on_axial_plane, material))
    end

    # Radial ordering identifies the fluid and fuel cells touching the wall.
    sort!(axial_faces; by = axial_face -> geo.cell_centroids[axial_face.cell_id][1])

    return axial_faces
end

function _build_neighbor_cell_by_face(geo)
    neighbor_cell_by_face = [zeros(Int, geo.cell_n_faces[cell_id]) for cell_id in 1:geo.n_cells]

    for (cell_id, neighbors) in geo.cell_neighbors
        for (neighbor_cell_id, face_idx, _) in neighbors
            neighbor_cell_by_face[cell_id][face_idx] = neighbor_cell_id
        end
    end

    return neighbor_cell_by_face
end

function _find_fuel_wall_faces(geo, port_fluid_cell_ids, fuel_cell_ids)
    fuel_cell_set = Set(fuel_cell_ids)
    fuel_wall_faces = FuelWallFace[]

    for fluid_cell_id in port_fluid_cell_ids
        for (neighbor_cell_id, fluid_face_idx, fuel_face_idx) in geo.cell_neighbors[fluid_cell_id][2]
            neighbor_is_fuel = neighbor_cell_id in fuel_cell_set

            if neighbor_is_fuel
                push!(
                    fuel_wall_faces,
                    FuelWallFace(fluid_cell_id, fluid_face_idx, neighbor_cell_id, fuel_face_idx),
                )
            end
        end
    end

    return fuel_wall_faces
end

"""
    create_fuel_regression_model(grid, geo; atol = 1e-10)

Describe the fuel port as a sequence of axial stations and construct the affine
node-motion map used during regression. The first
`length(model.axial_positions)` entries of `u.fuel_grain_thickness` hold the
remaining fuel thickness at those stations.

Mesh deformation is prescribed rather than solved. The fuel surface moves by
the consumed fuel thickness. The displacement decreases linearly to zero at the
inner fluid boundary and the outer fuel boundary.
"""
function create_fuel_regression_model(grid, geo; atol = 1e-10)
    fuel_cell_ids = sort!(collect(getcellset(grid, "fuel_grain_volume")))
    all_fluid_cell_ids = sort!(collect(getcellset(grid, "fluid_domain")))

    isempty(fuel_cell_ids) && error("The fuel_grain_volume cell set is empty")

    # Limit the moving fluid region to axial segments occupied by the grain.
    # This excludes injector and nozzle cells from the port deformation model.
    fuel_segment_z_positions = _unique_sorted_with_tolerance(
        (geo.cell_centroids[cell_id][3] for cell_id in fuel_cell_ids);
        atol,
    )
    port_fluid_cell_ids = [
        cell_id for cell_id in all_fluid_cell_ids
        if any(z -> abs(geo.cell_centroids[cell_id][3] - z) <= atol, fuel_segment_z_positions)
    ]

    #SHOULD BE:
    port_fluid_cell_ids = Int[]

    for cell_id in all_fluid_cell_ids
        if any(z -> abs(geo.cell_centroids[cell_id][3] - z) <= atol, fuel_segment_z_positions) #write a comment that says what this does because it isn't exactly clear
            push!(port_fluid_cell_ids, cell_id)
        end
    end
    #END

    fuel_node_ids = unique!(sort!(reduce(vcat, (geo.cell_node_ids[cell_id] for cell_id in fuel_cell_ids))))
    axial_z_positions = _unique_sorted_with_tolerance(
        (geo.node_coordinates[node_id][3] for node_id in fuel_node_ids);
        atol,
    )

    # Separate node sets preserve which material occupies each axial plane.
    all_node_ids_at_axial_position = [Set{Int}() for _ in axial_z_positions]
    fluid_node_ids_at_axial_position = [Set{Int}() for _ in axial_z_positions]
    fuel_node_ids_at_axial_position = [Set{Int}() for _ in axial_z_positions]

    for fluid_cell_id in port_fluid_cell_ids
        for node_id in geo.cell_node_ids[fluid_cell_id]
            node_z = geo.node_coordinates[node_id][3]
            axial_position_id = _find_axial_position_id(axial_z_positions, node_z; atol)

            push!(all_node_ids_at_axial_position[axial_position_id], node_id)
            push!(fluid_node_ids_at_axial_position[axial_position_id], node_id)
        end
    end

    for fuel_cell_id in fuel_cell_ids
        for node_id in geo.cell_node_ids[fuel_cell_id]
            node_z = geo.node_coordinates[node_id][3]
            axial_position_id = _find_axial_position_id(axial_z_positions, node_z; atol)

            push!(all_node_ids_at_axial_position[axial_position_id], node_id)
            push!(fuel_node_ids_at_axial_position[axial_position_id], node_id)
        end
    end

    axial_positions = AxialPosition[]
    regression_nodes = FuelRegressionNode[]

    for (axial_position_id, axial_z) in enumerate(axial_z_positions)
        fluid_node_ids = collect(fluid_node_ids_at_axial_position[axial_position_id])
        fuel_node_ids = collect(fuel_node_ids_at_axial_position[axial_position_id])

        isempty(fluid_node_ids) && error("No fluid nodes were found at z = $axial_z")
        isempty(fuel_node_ids) && error("No fuel nodes were found at z = $axial_z")

        fixed_inner_fluid_x = minimum(node_id -> geo.node_coordinates[node_id][1], fluid_node_ids)
        initial_fuel_surface_x = maximum(node_id -> geo.node_coordinates[node_id][1], fluid_node_ids)
        fixed_outer_fuel_x = maximum(node_id -> geo.node_coordinates[node_id][1], fuel_node_ids)
        initial_fuel_grain_thickness = fixed_outer_fuel_x - initial_fuel_surface_x

        upstream_fluid_faces = _collect_axial_faces(geo, port_fluid_cell_ids, axial_z, :upstream, :fluid; atol)
        downstream_fluid_faces = _collect_axial_faces(geo, port_fluid_cell_ids, axial_z, :downstream, :fluid; atol)
        upstream_fuel_faces = _collect_axial_faces(geo, fuel_cell_ids, axial_z, :upstream, :fuel; atol)
        downstream_fuel_faces = _collect_axial_faces(geo, fuel_cell_ids, axial_z, :downstream, :fuel; atol)

        node_ids_on_axial_plane = sort!(collect(all_node_ids_at_axial_position[axial_position_id]))

        axial_position = AxialPosition(
            axial_z,
            upstream_fluid_faces,
            downstream_fluid_faces,
            upstream_fuel_faces,
            downstream_fuel_faces,
            fixed_inner_fluid_x,
            initial_fuel_surface_x,
            initial_fuel_grain_thickness,
            fixed_outer_fuel_x,
            node_ids_on_axial_plane,
            axial_position_id,
        )
        push!(axial_positions, axial_position)

        # Spread the fuel-surface displacement linearly through the neighboring
        # fluid and fuel mesh. The two outer bounds remain fixed.
        for node_id in node_ids_on_axial_plane
            initial_node_x = geo.node_coordinates[node_id][1]

            if initial_node_x <= initial_fuel_surface_x + atol
                distance_from_fixed_fluid_boundary = initial_node_x - fixed_inner_fluid_x
                initial_fluid_region_width = initial_fuel_surface_x - fixed_inner_fluid_x
                port_radius_multiplier = distance_from_fixed_fluid_boundary / initial_fluid_region_width
            else
                distance_into_fuel = initial_node_x - initial_fuel_surface_x
                fraction_through_fuel = distance_into_fuel / initial_fuel_grain_thickness
                port_radius_multiplier = 1 - fraction_through_fuel
            end

            regression_node = FuelRegressionNode(
                node_id,
                axial_position_id,
                initial_node_x,
                port_radius_multiplier,
            )
            push!(regression_nodes, regression_node)
        end
    end

    moving_node_ids = Set(
        node.node_id for node in regression_nodes if abs(node.port_radius_multiplier) > atol
    )
    affected_fluid_cell_ids = [cell_id for cell_id in all_fluid_cell_ids if any(node_id -> node_id in moving_node_ids, geo.cell_node_ids[cell_id])]

    #SHOULD BE:
    moving_node_ids = Int[]

    for node in regression_nodes
        if abs(node.port_radius_multiplier) > atol
            push!(moving_node_ids, node.node_id)
        end
    end

    affected_fluid_cell_ids = Int[]

    for cell_id in all_fluid_cell_ids
        for node_id in geo.cell_node_ids[cell_id]
            if node_id in moving_node_ids
                push!(affected_fluid_cell_ids, cell_id)
                break
            end
        end
    end
    #END

    fluid_cell_mask = falses(geo.n_cells)
    fluid_cell_mask[all_fluid_cell_ids] .= true

    initial_node_coordinates = SVector{3, Float64}.(geo.node_coordinates)
    number_of_mesh_nodes = length(geo.node_coordinates)
    number_of_axial_positions = length(axial_positions)

    workspaces = GeneralLazyBufferCache() do scalar_seed
        scalar_zero = zero(scalar_seed)
        stationary_node_velocity = SVector(scalar_zero, scalar_zero, scalar_zero)

        FuelRegressionWorkspace(
            fill(stationary_node_velocity, number_of_mesh_nodes),
            fill(scalar_zero, number_of_axial_positions),
        )
    end

    fuel_wall_faces = _find_fuel_wall_faces(geo, port_fluid_cell_ids, fuel_cell_ids)
    neighbor_cell_by_face = _build_neighbor_cell_by_face(geo)

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

function initialize_fuel_grain_thickness!(fuel_grain_thickness, model::FuelRegressionModel)
    thickness_unit = oneunit(first(fuel_grain_thickness))

    for axial_position in model.axial_positions
        thickness_state_idx = axial_position.thickness_state_idx
        initial_fuel_grain_thickness = axial_position.initial_fuel_grain_thickness * thickness_unit

        fuel_grain_thickness[thickness_state_idx] = initial_fuel_grain_thickness
    end

    return fuel_grain_thickness
end

function move_fuel_wall!(u, model::FuelRegressionModel, geo)
    # Reconstruct the current mesh from its original coordinates every time.
    # This avoids accumulating mesh-displacement error between residual calls.
    copyto!(geo.node_coordinates, model.initial_node_coordinates)

    for regression_node in model.regression_nodes
        axial_position = model.axial_positions[regression_node.axial_position_id]
        remaining_fuel_thickness = u.fuel_grain_thickness[axial_position.thickness_state_idx]
        consumed_fuel_thickness = axial_position.initial_fuel_grain_thickness - remaining_fuel_thickness

        radial_node_displacement = regression_node.port_radius_multiplier * consumed_fuel_thickness
        current_node_x = regression_node.initial_x + radial_node_displacement

        current_node_coordinate = geo.node_coordinates[regression_node.node_id]
        geo.node_coordinates[regression_node.node_id] = setindex(current_node_coordinate, current_node_x, 1)
    end

    return nothing
end

function update_fuel_regression_geometry!(u, model::FuelRegressionModel, geo, grid)
    move_fuel_wall!(u, model, geo)

    update_fvm_geometry!(geo, geo.node_coordinates)

    # Reapply the hybrid-motor-specific corrections after the generic
    # finite-volume geometry has been rebuilt from the moved nodes.
    correct_hybrid_geo!(geo, grid)

    return nothing
end

#SHOULD BE:
function update_fuel_regression_geometry!(u, model::FuelRegressionModel, geo, grid)
    move_fuel_wall!(u, model, geo)
    update_fvm_geometry!(geo, geo.node_coordinates)

    # Reapply the hybrid-motor-specific corrections after the generic
    # finite-volume geometry has been rebuilt from the moved nodes.
    correct_hybrid_geo!(geo, grid)
end
#END


# -----------------------------------------------------------------------------
# Oxidizer mass flux and fuel-surface temperature
# -----------------------------------------------------------------------------

function calculate_axial_nitrous_oxide_mass_flux(u, geo, fluid_axial_faces)
    nitrous_oxide_mass_flow_rate = 0.0
    port_cross_sectional_area = 0.0

    for axial_face in fluid_axial_faces
        fluid_cell_id = axial_face.cell_id
        axial_face_idx = axial_face.face_idx_on_z_plane

        mixture_density = u.density[fluid_cell_id]
        axial_momentum_density = u.momentum_density_w[fluid_cell_id]
        axial_velocity = iszero(mixture_density) ? zero(mixture_density) : axial_momentum_density / mixture_density #I hate this syntax, if statements always look nicer to me

        nitrous_oxide_density = u.species_densities.nitrous_oxide[fluid_cell_id]
        axial_face_area = geo.cell_face_areas[fluid_cell_id, axial_face_idx]
        nitrous_oxide_mass_flow_through_face = abs(axial_velocity) * nitrous_oxide_density * axial_face_area

        nitrous_oxide_mass_flow_rate += nitrous_oxide_mass_flow_through_face
        port_cross_sectional_area += axial_face_area
    end

    if iszero(port_cross_sectional_area)
        return zero(nitrous_oxide_mass_flow_rate)
    end

    nitrous_oxide_mass_flux = nitrous_oxide_mass_flow_rate / port_cross_sectional_area
    return nitrous_oxide_mass_flux
end

function calculate_fuel_grain_face_temperature(u, geo, fluid_cell_nearest_wall::AxialFace, fuel_cell_nearest_wall::AxialFace)
    fluid_cell_id = fluid_cell_nearest_wall.cell_id
    fuel_cell_id = fuel_cell_nearest_wall.cell_id

    fluid_wall_face_idx = fluid_cell_nearest_wall.face_idx_facing_toward_fuel_wall
    fuel_wall_face_idx = fuel_cell_nearest_wall.face_idx_facing_toward_fuel_wall
    fluid_centroid_to_wall_distance = geo.cell_face_distances[fluid_cell_id, fluid_wall_face_idx]
    fuel_centroid_to_wall_distance = geo.cell_face_distances[fuel_cell_id, fuel_wall_face_idx]

    # Treat the fluid-side and fuel-side half cells as one-dimensional thermal
    # resistances meeting at a common interface temperature. Radiation and a
    # separate wall heat-transfer correlation are not included.
    fluid_thermal_conductance_per_area = u.k[fluid_cell_id] / fluid_centroid_to_wall_distance
    fuel_thermal_conductance_per_area = u.k[fuel_cell_id] / fuel_centroid_to_wall_distance
    total_thermal_conductance_per_area = fluid_thermal_conductance_per_area + fuel_thermal_conductance_per_area

    if iszero(total_thermal_conductance_per_area)
        average_adjacent_cell_temperature = (u.temperature[fluid_cell_id] + u.temperature[fuel_cell_id]) / 2
        return average_adjacent_cell_temperature
    end

    fluid_temperature_contribution = fluid_thermal_conductance_per_area * u.temperature[fluid_cell_id]
    fuel_temperature_contribution = fuel_thermal_conductance_per_area * u.temperature[fuel_cell_id]
    fuel_grain_face_temperature = (fluid_temperature_contribution + fuel_temperature_contribution) / total_thermal_conductance_per_area

    return fuel_grain_face_temperature
end

function get_axial_oxidizer_mass_flux_and_fuel_grain_face_temperature(u, geo, axial_position::AxialPosition)
    # Use cells immediately upstream of each plane. At the grain inlet there is
    # no upstream grain segment, so use the cells immediately downstream.
    use_downstream_cells = isempty(axial_position.sorted_upstream_fluid_cell_entries)

    if use_downstream_cells
        fluid_axial_faces = axial_position.sorted_downstream_fluid_cell_entries
        fuel_axial_faces = axial_position.sorted_downstream_fuel_cell_entries
    else
        fluid_axial_faces = axial_position.sorted_upstream_fluid_cell_entries
        fuel_axial_faces = axial_position.sorted_upstream_fuel_cell_entries
    end

    if isempty(fluid_axial_faces) || isempty(fuel_axial_faces)
        return 0.0, 0.0
    end

    nitrous_oxide_mass_flux = calculate_axial_nitrous_oxide_mass_flux(u, geo, fluid_axial_faces)

    # Axial faces are sorted from low x to high x. The last fluid cell and first
    # fuel cell are therefore the pair adjacent to the fuel surface.
    fluid_cell_nearest_wall = last(fluid_axial_faces)
    fuel_cell_nearest_wall = first(fuel_axial_faces)
    fuel_grain_face_temperature = calculate_fuel_grain_face_temperature(
        u,
        geo,
        fluid_cell_nearest_wall,
        fuel_cell_nearest_wall,
    )

    #SHOULD BE:
    fuel_grain_face_temperature = calculate_fuel_grain_face_temperature(u, geo, fluid_cell_nearest_wall, fuel_cell_nearest_wall)
    #END

    return nitrous_oxide_mass_flux, fuel_grain_face_temperature
end


# -----------------------------------------------------------------------------
# Empirical fuel-regression law
# -----------------------------------------------------------------------------

function calculate_fuel_regression_rates!(du, u, p, t, system, model::FuelRegressionModel, workspace, geo)
    for (axial_position_id, axial_position) in enumerate(model.axial_positions)
        nitrous_oxide_mass_flux, fuel_grain_face_temperature =
            get_axial_oxidizer_mass_flux_and_fuel_grain_face_temperature(u, geo, axial_position)

        #SHOULD BE:
        nitrous_oxide_mass_flux, fuel_grain_face_temperature = get_axial_oxidizer_mass_flux_and_fuel_grain_face_temperature(u, geo, axial_position)
        #END

        nonnegative_nitrous_oxide_mass_flux = max(nitrous_oxide_mass_flux, zero(nitrous_oxide_mass_flux))
        fuel_surface_is_hot_enough = fuel_grain_face_temperature > p.min_fuel_regression_temperature

        # Doran et al. correlation: r_dot = a * G_ox^n. The temperature cutoff
        # prevents cold fuel from releasing ethylene before ignition.
        if fuel_surface_is_hot_enough
            fuel_surface_regression_rate = p.regression_coefficient_a * nonnegative_nitrous_oxide_mass_flux^p.regression_exponent_n
        else
            fuel_surface_regression_rate = zero(nonnegative_nitrous_oxide_mass_flux)
        end

        thickness_state_idx = axial_position.thickness_state_idx
        du.fuel_grain_thickness[thickness_state_idx] = -fuel_surface_regression_rate
        workspace.port_radius_rates[axial_position_id] = fuel_surface_regression_rate
    end

    return nothing
end


# -----------------------------------------------------------------------------
# Mesh velocity induced by fuel regression
# -----------------------------------------------------------------------------

function calculate_regression_mesh_node_velocities!(model::FuelRegressionModel, workspace)
    stationary_node_velocity = zero(eltype(workspace.node_velocities))
    fill!(workspace.node_velocities, stationary_node_velocity)

    for regression_node in model.regression_nodes
        fuel_surface_regression_rate = workspace.port_radius_rates[regression_node.axial_position_id]
        radial_mesh_velocity = regression_node.port_radius_multiplier * fuel_surface_regression_rate

        workspace.node_velocities[regression_node.node_id] = SVector(radial_mesh_velocity,
            zero(radial_mesh_velocity),
            zero(radial_mesh_velocity),
        )

        #SHOULD BE:
        workspace.node_velocities[regression_node.node_id] = SVector(radial_mesh_velocity, zero(radial_mesh_velocity), zero(radial_mesh_velocity))
        #END

        #OR AT LEAST:
        workspace.node_velocities[regression_node.node_id] = SVector(
            radial_mesh_velocity,
            zero(radial_mesh_velocity),
            zero(radial_mesh_velocity),
        )
        #END
    end

    return nothing
end

@inline function calculate_face_mesh_velocity(model::FuelRegressionModel, workspace, geo, cell_id, face_idx)
    node_ids_on_face = _face_node_ids(geo, cell_id, face_idx)
    face_mesh_velocity = zero(eltype(workspace.node_velocities))

    for node_id in node_ids_on_face
        face_mesh_velocity += workspace.node_velocities[node_id]
    end

    face_mesh_velocity /= length(node_ids_on_face)
    return face_mesh_velocity

    #SHOULD BE:
    face_mesh_velocity /= length(node_ids_on_face)

    return face_mesh_velocity
    #END
end

function calculate_cell_volume_rate(model::FuelRegressionModel, workspace, geo, cell_id)
    cell_volume_rate = 0.0

    for face_idx in 1:geo.cell_n_faces[cell_id]
        face_mesh_velocity = calculate_face_mesh_velocity(model, workspace, geo, cell_id, face_idx)
        outward_face_normal = geo.cell_face_normals[cell_id, face_idx]
        face_area = geo.cell_face_areas[cell_id, face_idx]
        outward_swept_volume_rate = dot(face_mesh_velocity, outward_face_normal) * face_area

        cell_volume_rate += outward_swept_volume_rate
    end

    return cell_volume_rate
end


# -----------------------------------------------------------------------------
# ALE conservation in the deforming port mesh
# -----------------------------------------------------------------------------

function calculate_ale_transport_for_conserved_density(
    conserved_density,
    model::FuelRegressionModel,
    workspace,
    geo,
    fluid_cell_id,
    cell_volume_rate,
)

#SHOULD BE:
function calculate_ale_transport_for_conserved_density(conserved_density, model::FuelRegressionModel, workspace, geo, fluid_cell_id, cell_volume_rate)
#END
    conserved_quantity_flow_due_to_mesh_motion = zero(conserved_density[fluid_cell_id] * cell_volume_rate)

    for face_idx in 1:geo.cell_n_faces[fluid_cell_id]
        face_mesh_velocity = calculate_face_mesh_velocity(model, workspace, geo, fluid_cell_id, face_idx)
        outward_face_normal = geo.cell_face_normals[fluid_cell_id, face_idx]
        face_area = geo.cell_face_areas[fluid_cell_id, face_idx]
        outward_swept_volume_rate = dot(face_mesh_velocity, outward_face_normal) * face_area

        neighbor_cell_id = model.neighbor_cell_by_face[fluid_cell_id][face_idx]
        face_has_fluid_neighbor = neighbor_cell_id != 0 && model.fluid_cell_mask[neighbor_cell_id]

        if face_has_fluid_neighbor
            conserved_density_at_face = (conserved_density[fluid_cell_id] + conserved_density[neighbor_cell_id]) / 2
            conserved_quantity_flow_through_face = conserved_density_at_face * outward_swept_volume_rate

            conserved_quantity_flow_due_to_mesh_motion += conserved_quantity_flow_through_face
        end
    end

    # The state stores conserved quantity per unit volume. This term accounts
    # for the changing denominator as the cell volume changes.
    changing_cell_volume_correction = conserved_density[fluid_cell_id] * cell_volume_rate
    net_ale_transport = conserved_quantity_flow_due_to_mesh_motion - changing_cell_volume_correction

    return net_ale_transport
end

function apply_fuel_regression_ale!(du, u, model::FuelRegressionModel, workspace, geo)
    for fluid_cell_id in model.affected_fluid_cell_ids
        cell_volume_rate = calculate_cell_volume_rate(model, workspace, geo, fluid_cell_id)
        inverse_cell_volume = inv(geo.cell_volumes[fluid_cell_id])

        density_ale_transport = calculate_ale_transport_for_conserved_density(
            u.density, model, workspace, geo, fluid_cell_id, cell_volume_rate,
        )
        radial_momentum_ale_transport = calculate_ale_transport_for_conserved_density(
            u.momentum_density_u, model, workspace, geo, fluid_cell_id, cell_volume_rate,
        )
        transverse_momentum_ale_transport = calculate_ale_transport_for_conserved_density(
            u.momentum_density_v, model, workspace, geo, fluid_cell_id, cell_volume_rate,
        )
        axial_momentum_ale_transport = calculate_ale_transport_for_conserved_density(
            u.momentum_density_w, model, workspace, geo, fluid_cell_id, cell_volume_rate,
        )
        energy_ale_transport = calculate_ale_transport_for_conserved_density(
            u.volumetric_energy, model, workspace, geo, fluid_cell_id, cell_volume_rate,
        )

        #SHOULD BE:
        density_ale_transport = calculate_ale_transport_for_conserved_density(u.density, model, workspace, geo, fluid_cell_id, cell_volume_rate)

        u_momentum_ale_transport = calculate_ale_transport_for_conserved_density(u.momentum_density_u, model, workspace, geo, fluid_cell_id, cell_volume_rate)
        v_momentum_ale_transport = calculate_ale_transport_for_conserved_density(u.momentum_density_v, model, workspace, geo, fluid_cell_id, cell_volume_rate)
        w_momentum_ale_transport = calculate_ale_transport_for_conserved_density(u.momentum_density_w, model, workspace, geo, fluid_cell_id, cell_volume_rate)

        energy_ale_transport = calculate_ale_transport_for_conserved_density(u.volumetric_energy, model, workspace, geo, fluid_cell_id, cell_volume_rate)
        #END

        du.density[fluid_cell_id] += density_ale_transport * inverse_cell_volume
        du.momentum_density_u[fluid_cell_id] += radial_momentum_ale_transport * inverse_cell_volume
        du.momentum_density_v[fluid_cell_id] += transverse_momentum_ale_transport * inverse_cell_volume
        du.momentum_density_w[fluid_cell_id] += axial_momentum_ale_transport * inverse_cell_volume
        du.volumetric_energy[fluid_cell_id] += energy_ale_transport * inverse_cell_volume

        #SHOULD BE:
        du.density[fluid_cell_id] += density_ale_transport * inverse_cell_volume

        du.momentum_density_u[fluid_cell_id] += radial_momentum_ale_transport * inverse_cell_volume
        du.momentum_density_v[fluid_cell_id] += transverse_momentum_ale_transport * inverse_cell_volume
        du.momentum_density_w[fluid_cell_id] += axial_momentum_ale_transport * inverse_cell_volume

        du.volumetric_energy[fluid_cell_id] += energy_ale_transport * inverse_cell_volume
        #END
        

        # All species undergo the same numerical mesh transport. Keeping this
        # loop generic avoids repeating identical ALE bookkeeping for each gas.
        for_fields!(du.species_densities, u.species_densities) do species, du_species, u_species
            species_ale_transport = calculate_ale_transport_for_conserved_density(
                u_species[species], model, workspace, geo, fluid_cell_id, cell_volume_rate,
            )

            #SHOULD BE:
            species_ale_transport = calculate_ale_transport_for_conserved_density(u_species[species], model, workspace, geo, fluid_cell_id, cell_volume_rate)
            #END

            du_species[species[fluid_cell_id]] += species_ale_transport * inverse_cell_volume
        end

        if hasproperty(u, :turbulent_kinetic_energy_density)
            turbulent_kinetic_energy_ale_transport = calculate_ale_transport_for_conserved_density(
                u.turbulent_kinetic_energy_density, model, workspace, geo, fluid_cell_id, cell_volume_rate,
            )
            specific_dissipation_rate_ale_transport = calculate_ale_transport_for_conserved_density(
                u.specific_dissipation_rate_density, model, workspace, geo, fluid_cell_id, cell_volume_rate,
            )

            #SHOULD BE:
            turbulent_kinetic_energy_ale_transport = calculate_ale_transport_for_conserved_density(u.turbulent_kinetic_energy_density, model, workspace, geo, fluid_cell_id, cell_volume_rate)
            specific_dissipation_rate_ale_transport = calculate_ale_transport_for_conserved_density(u.specific_dissipation_rate_density, model, workspace, geo, fluid_cell_id, cell_volume_rate)
            #END

            du.turbulent_kinetic_energy_density[fluid_cell_id] += turbulent_kinetic_energy_ale_transport * inverse_cell_volume
            du.specific_dissipation_rate_density[fluid_cell_id] += specific_dissipation_rate_ale_transport * inverse_cell_volume
        end
    end

    return nothing
end


# -----------------------------------------------------------------------------
# Fuel mass, momentum, and energy released at the regressing wall
# -----------------------------------------------------------------------------

function apply_fuel_regression_sources!(
    du, u, p, t, system, model::FuelRegressionModel, workspace, geo,
)
    for fuel_wall_face in model.wall_faces
        fluid_cell_id = fuel_wall_face.fluid_cell_id
        fuel_cell_id = fuel_wall_face.fuel_cell_id
        fluid_face_idx = fuel_wall_face.fluid_face_idx

        fuel_surface_velocity = calculate_face_mesh_velocity(model, workspace, geo, fluid_cell_id, fluid_face_idx)
        outward_fluid_face_normal = geo.cell_face_normals[fluid_cell_id, fluid_face_idx]
        fuel_wall_face_area = geo.cell_face_areas[fluid_cell_id, fluid_face_idx]
        volume_rate_opened_by_regressing_fuel = dot(fuel_surface_velocity, outward_fluid_face_normal) * fuel_wall_face_area

        # Only wall motion that enlarges the fluid cell releases fuel. The
        # gas-phase chemistry model represents the pyrolyzed fuel as ethylene.
        released_fuel_volume_rate = max(volume_rate_opened_by_regressing_fuel, zero(volume_rate_opened_by_regressing_fuel))
        solid_fuel_density = u.density[fuel_cell_id]
        pyrolyzed_fuel_mass_flow_rate = solid_fuel_density * released_fuel_volume_rate
        #I think this the equations here don't have enough new lines between them, it's getting a bit hard to read

        du.density_flow[fluid_cell_id] += pyrolyzed_fuel_mass_flow_rate
        du.species_density_flow.ethylene[fluid_cell_id] += pyrolyzed_fuel_mass_flow_rate

        radial_fuel_momentum_flow = pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[1]
        transverse_fuel_momentum_flow = pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[2]
        axial_fuel_momentum_flow = pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[3]

        du.momentum_density_u_flow[fluid_cell_id] += radial_fuel_momentum_flow
        du.momentum_density_v_flow[fluid_cell_id] += transverse_fuel_momentum_flow
        du.momentum_density_w_flow[fluid_cell_id] += axial_fuel_momentum_flow

        #SHOULD BE:
        du.momentum_density_u_flow[fluid_cell_id] += pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[1]
        du.momentum_density_v_flow[fluid_cell_id] += pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[2]
        du.momentum_density_w_flow[fluid_cell_id] += pyrolyzed_fuel_mass_flow_rate * fuel_surface_velocity[3]
        #END

        solid_fuel_specific_energy = u.volumetric_energy[fuel_cell_id] / solid_fuel_density
        injected_fuel_specific_energy = solid_fuel_specific_energy - p.fuel_heat_of_pyrolysis
        pyrolyzed_fuel_energy_flow = pyrolyzed_fuel_mass_flow_rate * injected_fuel_specific_energy

        du.volumetric_energy_flow[fluid_cell_id] += pyrolyzed_fuel_energy_flow
    end

    return nothing
end


# -----------------------------------------------------------------------------
# Fuel-regression sequence called by the coupled rocket residual
# -----------------------------------------------------------------------------

function regress_fuel_layers!(du, u, p, t, system, geo)
    fuel_regression_model = system.fuel_regression_model

    first_axial_position = first(fuel_regression_model.axial_positions)
    first_thickness_state_idx = first_axial_position.thickness_state_idx
    scalar_type_seed = u.fuel_grain_thickness[first_thickness_state_idx]
    fuel_regression_workspace = get_tmp(fuel_regression_model.workspaces, scalar_type_seed)

    calculate_fuel_regression_rates!(
        du, u, p, t, system, fuel_regression_model, fuel_regression_workspace, geo,
    )

    #SHOULD BE:
    calculate_fuel_regression_rates!(du, u, p, t, system, fuel_regression_model, fuel_regression_workspace, geo)
    #END

    calculate_regression_mesh_node_velocities!(fuel_regression_model, fuel_regression_workspace)

    apply_fuel_regression_ale!(du, u, fuel_regression_model, fuel_regression_workspace, geo)

    apply_fuel_regression_sources!(
        du, u, p, t, system, fuel_regression_model, fuel_regression_workspace, geo,
    )

    #SHOULD BE:
    apply_fuel_regression_sources!(du, u, p, t, system, fuel_regression_model, fuel_regression_workspace, geo)
    #END

    return nothing
end
