struct AxialFace
    node_ids_on_z_plane::Vector{Int}
    cell_id::Int
    face_idx_on_z_plane::Int
    face_idx_facing_toward_fuel_wall::Int
    face_idx_facing_away_from_fuel_wall::Int
    closest_node_ids_from_fuel_wall::Vector{Int}
    furthest_node_ids_from_fuel_wall::Vector{Int}
end

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

struct FuelRegressionNode
    node_id::Int
    axial_position_id::Int
    initial_x::Float64
    port_radius_multiplier::Float64
end

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

function _unique_sorted_with_tolerance(values; atol)
    sorted_values = sort!(collect(values))
    isempty(sorted_values) && return Float64[]

    unique_values = Float64[first(sorted_values)]
    for value in @view sorted_values[2:end]
        if abs(value - last(unique_values)) > atol
            push!(unique_values, value)
        end
    end
    return unique_values
end

@inline function _face_node_ids(geo, cell_id, face_idx)
    packed_face_id = geo.cell_face_layout[cell_id, face_idx]
    return geo.cell_faces_node_ids[packed_face_id]
end

function _find_axial_position_id(axial_z_positions, z; atol)
    position_id = findfirst(z_position -> abs(z_position - z) <= atol, axial_z_positions)
    isnothing(position_id) && error("No fuel-grain axial position was found at z = $z")
    return position_id
end

function _find_axial_face_idx(geo, cell_id, z; atol)
    for face_idx in 1:geo.cell_n_faces[cell_id]
        face_node_ids = _face_node_ids(geo, cell_id, face_idx)
        if all(node_id -> abs(geo.node_coordinates[node_id][3] - z) <= atol, face_node_ids)
            return face_idx
        end
    end
    error("Cell $cell_id does not have a face on the axial plane z = $z")
end

function _find_radial_face_indices(geo, cell_id, cell_kind::Symbol)
    radial_face_indices = Int[]
    for face_idx in 1:geo.cell_n_faces[cell_id]
        normal = geo.cell_face_normals[cell_id, face_idx]
        if abs(normal[1]) > 0.5
            push!(radial_face_indices, face_idx)
        end
    end
    length(radial_face_indices) == 2 || error(
        "Expected two radial faces for cell $cell_id, found $(length(radial_face_indices))",
    )

    sort!(radial_face_indices; by = face_idx -> geo.cell_face_centroids[cell_id, face_idx][1])
    if cell_kind === :fluid
        return last(radial_face_indices), first(radial_face_indices)
    elseif cell_kind === :fuel
        return first(radial_face_indices), last(radial_face_indices)
    end
    error("Unknown cell kind: $cell_kind")
end

function _make_axial_face(geo, cell_id, axial_face_idx, cell_kind::Symbol)
    toward_wall_face_idx, away_from_wall_face_idx =
        _find_radial_face_indices(geo, cell_id, cell_kind)
    return AxialFace(
        copy(_face_node_ids(geo, cell_id, axial_face_idx)),
        cell_id,
        axial_face_idx,
        toward_wall_face_idx,
        away_from_wall_face_idx,
        copy(_face_node_ids(geo, cell_id, toward_wall_face_idx)),
        copy(_face_node_ids(geo, cell_id, away_from_wall_face_idx)),
    )
end

function _collect_axial_faces(geo, cell_ids, z, side::Symbol, cell_kind::Symbol; atol)
    entries = AxialFace[]
    for cell_id in cell_ids
        centroid_z = geo.cell_centroids[cell_id][3]
        is_on_side = side === :upstream ? centroid_z < z - atol : centroid_z > z + atol
        is_on_side || continue

        face_idx = nothing
        for candidate_face_idx in 1:geo.cell_n_faces[cell_id]
            face_node_ids = _face_node_ids(geo, cell_id, candidate_face_idx)
            if all(node_id -> abs(geo.node_coordinates[node_id][3] - z) <= atol, face_node_ids)
                face_idx = candidate_face_idx
                break
            end
        end
        isnothing(face_idx) && continue
        push!(entries, _make_axial_face(geo, cell_id, face_idx, cell_kind))
    end
    sort!(entries; by = entry -> geo.cell_centroids[entry.cell_id][1])
    return entries
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

function _build_wall_faces(geo, port_fluid_cell_ids, fuel_cell_ids)
    fuel_cell_set = Set(fuel_cell_ids)
    wall_faces = FuelWallFace[]
    for fluid_cell_id in port_fluid_cell_ids
        for (neighbor_cell_id, fluid_face_idx, neighbor_face_idx) in geo.cell_neighbors[fluid_cell_id][2]
            if neighbor_cell_id in fuel_cell_set
                push!(wall_faces, FuelWallFace(
                    fluid_cell_id,
                    fluid_face_idx,
                    neighbor_cell_id,
                    neighbor_face_idx,
                ))
            end
        end
    end
    return wall_faces
end

"""
    create_fuel_regression_model(grid, geo; kwargs...)

Build the axial fuel-wall description and the affine node maps used during fuel
regression. The thickness states occupy the first `length(axial_positions)`
entries of `u.fuel_grain_thickness`.
"""
function create_fuel_regression_model(grid, geo; atol = 1e-10)
    fuel_cell_ids = sort!(collect(getcellset(grid, "fuel_grain_volume")))
    all_fluid_cell_ids = sort!(collect(getcellset(grid, "fluid_domain")))
    isempty(fuel_cell_ids) && error("The fuel_grain_volume cell set is empty")

    fuel_segment_z_positions = _unique_sorted_with_tolerance(
        (geo.cell_centroids[cell_id][3] for cell_id in fuel_cell_ids);
        atol,
    )
    port_fluid_cell_ids = [
        cell_id for cell_id in all_fluid_cell_ids
        if any(z -> abs(geo.cell_centroids[cell_id][3] - z) <= atol, fuel_segment_z_positions)
    ]

    fuel_node_ids = unique!(sort!(reduce(vcat, (geo.cell_node_ids[cell_id] for cell_id in fuel_cell_ids))))
    axial_z_positions = _unique_sorted_with_tolerance(
        (geo.node_coordinates[node_id][3] for node_id in fuel_node_ids);
        atol,
    )

    plane_node_sets = [Set{Int}() for _ in axial_z_positions]
    fluid_plane_node_sets = [Set{Int}() for _ in axial_z_positions]
    fuel_plane_node_sets = [Set{Int}() for _ in axial_z_positions]

    for cell_id in port_fluid_cell_ids
        for node_id in geo.cell_node_ids[cell_id]
            position_id = _find_axial_position_id(
                axial_z_positions,
                geo.node_coordinates[node_id][3];
                atol,
            )
            push!(plane_node_sets[position_id], node_id)
            push!(fluid_plane_node_sets[position_id], node_id)
        end
    end
    for cell_id in fuel_cell_ids
        for node_id in geo.cell_node_ids[cell_id]
            position_id = _find_axial_position_id(
                axial_z_positions,
                geo.node_coordinates[node_id][3];
                atol,
            )
            push!(plane_node_sets[position_id], node_id)
            push!(fuel_plane_node_sets[position_id], node_id)
        end
    end

    axial_positions = AxialPosition[]
    regression_nodes = FuelRegressionNode[]
    for (position_id, z) in enumerate(axial_z_positions)
        fluid_node_ids = collect(fluid_plane_node_sets[position_id])
        fuel_plane_node_ids = collect(fuel_plane_node_sets[position_id])
        isempty(fluid_node_ids) && error("No fluid nodes were found at z = $z")
        isempty(fuel_plane_node_ids) && error("No fuel nodes were found at z = $z")

        min_fluid_x = minimum(node_id -> geo.node_coordinates[node_id][1], fluid_node_ids)
        initial_fuel_x = maximum(node_id -> geo.node_coordinates[node_id][1], fluid_node_ids)
        max_fuel_x = maximum(node_id -> geo.node_coordinates[node_id][1], fuel_plane_node_ids)
        initial_thickness = max_fuel_x - initial_fuel_x

        upstream_fluid_faces = _collect_axial_faces(
            geo, port_fluid_cell_ids, z, :upstream, :fluid; atol,
        )
        downstream_fluid_faces = _collect_axial_faces(
            geo, port_fluid_cell_ids, z, :downstream, :fluid; atol,
        )
        upstream_fuel_faces = _collect_axial_faces(
            geo, fuel_cell_ids, z, :upstream, :fuel; atol,
        )
        downstream_fuel_faces = _collect_axial_faces(
            geo, fuel_cell_ids, z, :downstream, :fuel; atol,
        )

        node_ids = sort!(collect(plane_node_sets[position_id]))
        push!(axial_positions, AxialPosition(
            z,
            upstream_fluid_faces,
            downstream_fluid_faces,
            upstream_fuel_faces,
            downstream_fuel_faces,
            min_fluid_x,
            initial_fuel_x,
            initial_thickness,
            max_fuel_x,
            node_ids,
            position_id,
        ))

        for node_id in node_ids
            initial_x = geo.node_coordinates[node_id][1]
            if initial_x <= initial_fuel_x + atol
                multiplier = (initial_x - min_fluid_x) / (initial_fuel_x - min_fluid_x)
            else
                fraction_through_fuel =
                    (initial_x - initial_fuel_x) / (max_fuel_x - initial_fuel_x)
                multiplier = 1 - fraction_through_fuel
            end
            push!(regression_nodes, FuelRegressionNode(
                node_id,
                position_id,
                initial_x,
                multiplier,
            ))
        end
    end

    moving_node_ids = Set(
        node.node_id for node in regression_nodes if abs(node.port_radius_multiplier) > atol
    )
    affected_fluid_cell_ids = [
        cell_id for cell_id in all_fluid_cell_ids
        if any(node_id -> node_id in moving_node_ids, geo.cell_node_ids[cell_id])
    ]
    fluid_cell_mask = falses(geo.n_cells)
    fluid_cell_mask[all_fluid_cell_ids] .= true

    initial_node_coordinates = SVector{3, Float64}.(geo.node_coordinates)
    n_nodes = length(geo.node_coordinates)
    n_axial_positions = length(axial_positions)
    workspaces = GeneralLazyBufferCache() do seed
        zero_value = zero(seed)
        zero_velocity = SVector(zero_value, zero_value, zero_value)
        FuelRegressionWorkspace(
            fill(zero_velocity, n_nodes),
            fill(zero_value, n_axial_positions),
        )
    end

    return FuelRegressionModel(
        axial_positions,
        regression_nodes,
        _build_wall_faces(geo, port_fluid_cell_ids, fuel_cell_ids),
        affected_fluid_cell_ids,
        fluid_cell_mask,
        _build_neighbor_cell_by_face(geo),
        initial_node_coordinates,
        workspaces,
    )
end

function initialize_fuel_grain_thickness!(fuel_grain_thickness, model::FuelRegressionModel)
    thickness_unit = oneunit(first(fuel_grain_thickness))
    for axial_position in model.axial_positions
        fuel_grain_thickness[axial_position.thickness_state_idx] =
            axial_position.initial_fuel_grain_thickness * thickness_unit
    end
    return fuel_grain_thickness
end

function move_fuel_wall!(u, model::FuelRegressionModel, geo)
    copyto!(geo.node_coordinates, model.initial_node_coordinates)
    for regression_node in model.regression_nodes
        axial_position = model.axial_positions[regression_node.axial_position_id]
        wall_displacement = axial_position.initial_fuel_grain_thickness -
                            u.fuel_grain_thickness[axial_position.thickness_state_idx]
        new_x = regression_node.initial_x +
                regression_node.port_radius_multiplier * wall_displacement
        coordinate = geo.node_coordinates[regression_node.node_id]
        geo.node_coordinates[regression_node.node_id] = setindex(coordinate, new_x, 1)
    end
    return nothing
end

function update_fuel_regression_geometry!(u, model::FuelRegressionModel, geo, grid)
    move_fuel_wall!(u, model, geo)
    update_fvm_geometry!(geo, geo.node_coordinates)
    correct_hybrid_geo!(geo, grid)
    return nothing
end

function get_axial_oxidizer_mass_flux_and_fuel_grain_face_temperature(u, geo, axial_position::AxialPosition)
    use_downstream_entries = isempty(axial_position.sorted_upstream_fluid_cell_entries)
    face_entries = use_downstream_entries ?
        axial_position.sorted_downstream_fluid_cell_entries :
        axial_position.sorted_upstream_fluid_cell_entries
    fuel_entries = use_downstream_entries ?
        axial_position.sorted_downstream_fuel_cell_entries :
        axial_position.sorted_upstream_fuel_cell_entries
    isempty(face_entries) && return 0.0, 0.0
    isempty(fuel_entries) && return 0.0, 0.0

    oxidizer_mass_flow = 0.0
    cross_sectional_area = 0.0
    for face_entry in face_entries
        cell_id = face_entry.cell_id
        face_idx = face_entry.face_idx_on_z_plane
        density = u.density[cell_id]
        axial_velocity = iszero(density) ? zero(density) : u.momentum_density_w[cell_id] / density
        face_area = geo.cell_face_areas[cell_id, face_idx]
        oxidizer_mass_flow += abs(axial_velocity) * u.species_densities.nitrous_oxide[cell_id] * face_area
        cross_sectional_area += face_area
    end

    fluid_entry = last(face_entries)
    fuel_entry = first(fuel_entries)
    fluid_cell_id = fluid_entry.cell_id
    fuel_cell_id = fuel_entry.cell_id
    fluid_face_distance = geo.cell_face_distances[
        fluid_cell_id, fluid_entry.face_idx_facing_toward_fuel_wall,
    ]
    fuel_face_distance = geo.cell_face_distances[
        fuel_cell_id, fuel_entry.face_idx_facing_toward_fuel_wall,
    ]
    fluid_conductance = u.k[fluid_cell_id] / fluid_face_distance
    fuel_conductance = u.k[fuel_cell_id] / fuel_face_distance
    total_conductance = fluid_conductance + fuel_conductance
    fuel_grain_face_temperature = iszero(total_conductance) ?
        (u.temperature[fluid_cell_id] + u.temperature[fuel_cell_id]) / 2 :
        (fluid_conductance * u.temperature[fluid_cell_id] +
         fuel_conductance * u.temperature[fuel_cell_id]) / total_conductance

    oxidizer_mass_flux = iszero(cross_sectional_area) ?
        zero(oxidizer_mass_flow) : oxidizer_mass_flow / cross_sectional_area
    return oxidizer_mass_flux, fuel_grain_face_temperature
end

function calculate_fuel_regression_rates!(du, u, p, t, system, model::FuelRegressionModel, workspace, geo)
    for (position_id, axial_position) in enumerate(model.axial_positions)
        oxidizer_mass_flux, fuel_grain_face_temperature =
            get_axial_oxidizer_mass_flux_and_fuel_grain_face_temperature(u, geo, axial_position)
        if fuel_grain_face_temperature > p.min_fuel_regression_temperature
            regression_rate = p.regression_coefficient_a * max(oxidizer_mass_flux, zero(oxidizer_mass_flux)) ^ p.regression_exponent_n
            du.fuel_grain_thickness[axial_position.thickness_state_idx] = -regression_rate
            workspace.port_radius_rates[position_id] = regression_rate
        end
    end
    return nothing
end

function _update_node_velocities!(model::FuelRegressionModel, workspace)
    fill!(workspace.node_velocities, zero(eltype(workspace.node_velocities)))
    for regression_node in model.regression_nodes
        radial_velocity = regression_node.port_radius_multiplier *
                          workspace.port_radius_rates[regression_node.axial_position_id]
        workspace.node_velocities[regression_node.node_id] = SVector(
            radial_velocity,
            zero(radial_velocity),
            zero(radial_velocity),
        )
    end
    return nothing
end

@inline function _face_mesh_velocity(model::FuelRegressionModel, workspace, geo, cell_id, face_idx)
    face_node_ids = _face_node_ids(geo, cell_id, face_idx)
    velocity = zero(eltype(workspace.node_velocities))
    for node_id in face_node_ids
        velocity += workspace.node_velocities[node_id]
    end
    return velocity / length(face_node_ids)
end

function _mesh_flux_for_scalar(u_values, model, workspace, geo, cell_id, volume_rate)
    mesh_flux = zero(u_values[cell_id] * volume_rate)
    for face_idx in 1:geo.cell_n_faces[cell_id]
        face_velocity = _face_mesh_velocity(model, workspace, geo, cell_id, face_idx)
        swept_volume_rate = dot(face_velocity, geo.cell_face_normals[cell_id, face_idx]) *
                            geo.cell_face_areas[cell_id, face_idx]
        neighbor_cell_id = model.neighbor_cell_by_face[cell_id][face_idx]
        if neighbor_cell_id != 0 && model.fluid_cell_mask[neighbor_cell_id]
            face_value = (u_values[cell_id] + u_values[neighbor_cell_id]) / 2
            mesh_flux += face_value * swept_volume_rate
        end
    end
    return mesh_flux - u_values[cell_id] * volume_rate
end

function apply_fuel_regression_ale!(du, u, model::FuelRegressionModel, workspace, geo)
    for cell_id in model.affected_fluid_cell_ids
        volume_rate = 0.0
        for face_idx in 1:geo.cell_n_faces[cell_id]
            face_velocity = _face_mesh_velocity(model, workspace, geo, cell_id, face_idx)
            volume_rate += dot(face_velocity, geo.cell_face_normals[cell_id, face_idx]) *
                           geo.cell_face_areas[cell_id, face_idx]
        end
        inverse_volume = inv(geo.cell_volumes[cell_id])

        du.density[cell_id] += _mesh_flux_for_scalar(
            u.density, model, workspace, geo, cell_id, volume_rate,
        ) * inverse_volume
        du.momentum_density_u[cell_id] += _mesh_flux_for_scalar(
            u.momentum_density_u, model, workspace, geo, cell_id, volume_rate,
        ) * inverse_volume
        du.momentum_density_v[cell_id] += _mesh_flux_for_scalar(
            u.momentum_density_v, model, workspace, geo, cell_id, volume_rate,
        ) * inverse_volume
        du.momentum_density_w[cell_id] += _mesh_flux_for_scalar(
            u.momentum_density_w, model, workspace, geo, cell_id, volume_rate,
        ) * inverse_volume
        du.volumetric_energy[cell_id] += _mesh_flux_for_scalar(
            u.volumetric_energy, model, workspace, geo, cell_id, volume_rate,
        ) * inverse_volume

        for_fields!(du.species_densities, u.species_densities) do species, du_species, u_species
            du_species[species[cell_id]] += _mesh_flux_for_scalar(
                u_species[species], model, workspace, geo, cell_id, volume_rate,
            ) * inverse_volume
        end

        if hasproperty(u, :turbulent_kinetic_energy_density)
            du.turbulent_kinetic_energy_density[cell_id] += _mesh_flux_for_scalar(
                u.turbulent_kinetic_energy_density, model, workspace, geo, cell_id, volume_rate,
            ) * inverse_volume
            du.specific_dissipation_rate_density[cell_id] += _mesh_flux_for_scalar(
                u.specific_dissipation_rate_density, model, workspace, geo, cell_id, volume_rate,
            ) * inverse_volume
        end
    end
    return nothing
end

function apply_fuel_regression_sources!(
        du, u, p, t, system, model::FuelRegressionModel, workspace, geo,
    )
    for wall_face in model.wall_faces
        fluid_cell_id = wall_face.fluid_cell_id
        fuel_cell_id = wall_face.fuel_cell_id
        face_idx = wall_face.fluid_face_idx
        face_velocity = _face_mesh_velocity(model, workspace, geo, fluid_cell_id, face_idx)
        swept_volume_rate = dot(face_velocity, geo.cell_face_normals[fluid_cell_id, face_idx]) *
                            geo.cell_face_areas[fluid_cell_id, face_idx]
        fuel_mass_flow = u.density[fuel_cell_id] * max(swept_volume_rate, 0.0)

        du.density_flow[fluid_cell_id] += fuel_mass_flow
        du.species_density_flow.ethylene[fluid_cell_id] += fuel_mass_flow
        du.momentum_density_u_flow[fluid_cell_id] += fuel_mass_flow * face_velocity[1]
        du.momentum_density_v_flow[fluid_cell_id] += fuel_mass_flow * face_velocity[2]
        du.momentum_density_w_flow[fluid_cell_id] += fuel_mass_flow * face_velocity[3]

        fuel_specific_energy = u.volumetric_energy[fuel_cell_id] / u.density[fuel_cell_id]
        injected_specific_energy = fuel_specific_energy - p.fuel_heat_of_pyrolysis
        du.volumetric_energy_flow[fluid_cell_id] += fuel_mass_flow * injected_specific_energy
    end
    return nothing
end

function regress_fuel_layers!(du, u, p, t, system, geo)
    model = system.fuel_regression_model
    seed = u.fuel_grain_thickness[first(model.axial_positions).thickness_state_idx]
    workspace = get_tmp(model.workspaces, seed)
    calculate_fuel_regression_rates!(du, u, p, t, system, model, workspace, geo)
    _update_node_velocities!(model, workspace)
    apply_fuel_regression_ale!(du, u, model, workspace, geo)
    apply_fuel_regression_sources!(du, u, p, t, system, model, workspace, geo)
    return nothing
end
