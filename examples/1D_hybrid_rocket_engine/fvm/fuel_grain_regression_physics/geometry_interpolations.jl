function get_min_max_x_fuel_grain_bounds(grid, cells, node_coordinates)
    cells_node_coordinates = []
    for cell_id in cells 
        append!(cells_node_coordinates, [node_coordinates[node_id][1] for node_id in grid.cells[cell_id].nodes])
    end

    sort!(cells_node_coordinates)

    min_x_fuel_grain_bound = cells_node_coordinates[1]
    max_x_fuel_grain_bound = cells_node_coordinates[end]

    return min_x_fuel_grain_bound, max_x_fuel_grain_bound
end

struct FuelGrainRegressionGeometryInterpolations
    cell_volumes_interps::Vector
    cell_centroids_interps::Vector
    cell_neighbor_normals_interps::Vector
    cell_neighbor_distances_interps::Vector
    cell_face_areas_interps::Vector
    cell_face_normals_interps::Vector
    cell_face_distances_interps::Vector
    cell_face_centroids_interps::Vector
end

function make_face_scalar_interpolations(values_over_steps, cell_id, n_faces, interpolation_points)
    return SVector{n_faces}(ntuple(
        face_idx -> LinearInterpolation(
            [values[cell_id, face_idx] for values in values_over_steps],
            interpolation_points,
        ),
        n_faces,
    ))
end

function make_face_vector_interpolations(values_over_steps, cell_id, n_faces, interpolation_points)
    return SVector{n_faces}(ntuple(
        face_idx -> SVector{3}(ntuple(
            component -> LinearInterpolation(
                [values[cell_id, face_idx][component] for values in values_over_steps],
                interpolation_points,
            ),
            3,
        )),
        n_faces,
    ))
end

evaluate_geometry_interpolation(interpolation::DataInterpolations.LinearInterpolation, x) = interpolation(x)
evaluate_geometry_interpolation(interpolations::SVector, x) = map(
    interpolation -> evaluate_geometry_interpolation(interpolation, x),
    interpolations,
)

function generate_fuel_grain_regression_node_coordinates(
    starting_node_coordinates,
    fluid_node_ids,
    fuel_node_ids,
    min_x_fluid_bound,
    initial_interface_x,
    max_x_fuel_grain_bound,
    n_steps,
)
    interface_x_over_steps = collect(range(
        initial_interface_x,
        max_x_fuel_grain_bound,
        length = n_steps + 1,
    ))[1:end-1]
    node_coordinates_over_steps = []

    for interface_x in interface_x_over_steps
        node_coordinates = deepcopy(starting_node_coordinates)

        for node_id in fluid_node_ids
            starting_x = starting_node_coordinates[node_id][1]
            fluid_fraction = (starting_x - min_x_fluid_bound) /
                (initial_interface_x - min_x_fluid_bound)
            interpolated_x = min_x_fluid_bound +
                fluid_fraction * (interface_x - min_x_fluid_bound)

            if isapprox(starting_x, min_x_fluid_bound; atol = 1e-12, rtol = 0)
                interpolated_x = starting_x
            end

            node_coordinates[node_id] = setindex(node_coordinates[node_id], interpolated_x, 1)
        end

        for node_id in fuel_node_ids
            starting_x = starting_node_coordinates[node_id][1]
            fuel_fraction = (starting_x - initial_interface_x) /
                (max_x_fuel_grain_bound - initial_interface_x)
            interpolated_x = interface_x +
                fuel_fraction * (max_x_fuel_grain_bound - interface_x)

            if isapprox(starting_x, max_x_fuel_grain_bound; atol = 1e-12, rtol = 0)
                interpolated_x = starting_x
            end

            node_coordinates[node_id] = setindex(node_coordinates[node_id], interpolated_x, 1)
        end

        push!(node_coordinates_over_steps, node_coordinates)
    end

    return node_coordinates_over_steps, interface_x_over_steps
end

function precompute_fuel_grain_regression_geometry_interpolations(
    grid, geo, 
    dominant_fuel_grain_cells, 
    fuel_grain_subsidiary_cell_ids_to_dominant_cell_id_map,
    fuel_grain_dominant_cell_id_to_subsidiaries_map,
    fluid_domain_subsidiary_cell_ids_to_dominant_cell_id_map,
    fluid_domain_dominant_cell_id_to_subsidiaries_map,
)
    node_coordinates = get_node_coordinates(grid)

    fuel_grain_cells = vcat(
        dominant_fuel_grain_cells,
        collect(keys(fuel_grain_subsidiary_cell_ids_to_dominant_cell_id_map)),
    )
    fluid_domain_cells = collect(keys(fluid_domain_subsidiary_cell_ids_to_dominant_cell_id_map))

    min_x_fuel_grain_bound, max_x_fuel_grain_bound = get_min_max_x_fuel_grain_bounds(
        grid,
        fuel_grain_cells,
        node_coordinates,
    )
    min_x_fluid_bound, max_x_fluid_bound = get_min_max_x_fuel_grain_bounds(
        grid,
        fluid_domain_cells,
        node_coordinates,
    )
    @assert isapprox(min_x_fuel_grain_bound, max_x_fluid_bound)

    #since every cell is a rectangular prism with 90* corners, we can use some simple logic here to get the interpolations we need
    #However, the geometry of these cells have been modified previously with geometry_editing so this might be more challenging that I think
    #A lazy way would just be to rerun update_fvm_geometry!(geo, node_coordinates)
    #then rerun the geometry_editing_script
    #this might actually be the most efficient way to do this since it's really easy and probably quite fast
    #yep, we're going with this

    #interpolations we need (for all cells)
    n_cells = geo.n_cells
    cell_volumes_interps = Vector{DataInterpolations.LinearInterpolation}(undef, n_cells)
    cell_centroids_interps = Vector{DataInterpolations.LinearInterpolation}(undef, n_cells)

    #   Neighbor information, stored in the same half-face layout
    cell_neighbor_normals_interps = Vector(undef, n_cells)
    cell_neighbor_distances_interps = Vector(undef, n_cells)
    
    #   Face information
    cell_face_areas_interps = Vector(undef, n_cells)
    cell_face_normals_interps = Vector(undef, n_cells)
    cell_face_distances_interps = Vector(undef, n_cells)
    cell_face_centroids_interps = Vector(undef, n_cells)

    cell_volumes_over_steps = []
    cell_centroids_over_steps = []

    cell_neighbor_normals_over_steps = []
    cell_neighbor_distances_over_steps = []

    cell_face_areas_over_steps = []
    cell_face_normals_over_steps = []
    cell_face_distances_over_steps = []
    cell_face_centroids_over_steps = []

    n_steps = 5

    starting_node_coordinates = deepcopy(node_coordinates)
    fluid_node_ids = unique([
        node_id for cell_id in fluid_domain_cells for node_id in grid.cells[cell_id].nodes
    ])
    fuel_node_ids = unique([
        node_id for cell_id in fuel_grain_cells for node_id in grid.cells[cell_id].nodes
    ])
    node_coordinates_over_steps, interface_x_over_steps =
        generate_fuel_grain_regression_node_coordinates(
            starting_node_coordinates,
            fluid_node_ids,
            fuel_node_ids,
            min_x_fluid_bound,
            max_x_fluid_bound,
            max_x_fuel_grain_bound,
            n_steps,
        )
    
    for node_coordinates_at_step in node_coordinates_over_steps
        step_geo = deepcopy(geo)
        update_fvm_geometry!(step_geo, node_coordinates_at_step)
        correct_hybrid_geo!(step_geo, grid)

        push!(cell_volumes_over_steps, deepcopy(step_geo.cell_volumes))
        push!(cell_centroids_over_steps, deepcopy(step_geo.cell_centroids))

        push!(cell_neighbor_normals_over_steps, deepcopy(step_geo.cell_neighbor_normals))
        push!(cell_neighbor_distances_over_steps, deepcopy(step_geo.cell_neighbor_distances))

        push!(cell_face_areas_over_steps, deepcopy(step_geo.cell_face_areas))
        push!(cell_face_normals_over_steps, deepcopy(step_geo.cell_face_normals))
        push!(cell_face_distances_over_steps, deepcopy(step_geo.cell_face_distances))
        push!(cell_face_centroids_over_steps, deepcopy(step_geo.cell_face_centroids))
    end

    for cell_id in 1:geo.n_cells
        cell_volumes = [values[cell_id] for values in cell_volumes_over_steps]
        cell_centroids = [values[cell_id] for values in cell_centroids_over_steps]
        n_faces = geo.cell_n_faces[cell_id]

        cell_volumes_interps[cell_id] = LinearInterpolation(cell_volumes, interface_x_over_steps)
        cell_centroids_interps[cell_id] = LinearInterpolation(cell_centroids, interface_x_over_steps)

        cell_neighbor_normals_interps[cell_id] = make_face_vector_interpolations(
            cell_neighbor_normals_over_steps, cell_id, n_faces, interface_x_over_steps,
        )
        cell_neighbor_distances_interps[cell_id] = make_face_scalar_interpolations(
            cell_neighbor_distances_over_steps, cell_id, n_faces, interface_x_over_steps,
        )

        cell_face_areas_interps[cell_id] = make_face_scalar_interpolations(
            cell_face_areas_over_steps, cell_id, n_faces, interface_x_over_steps,
        )
        cell_face_normals_interps[cell_id] = make_face_vector_interpolations(
            cell_face_normals_over_steps, cell_id, n_faces, interface_x_over_steps,
        )
        cell_face_distances_interps[cell_id] = make_face_scalar_interpolations(
            cell_face_distances_over_steps, cell_id, n_faces, interface_x_over_steps,
        )
        cell_face_centroids_interps[cell_id] = make_face_vector_interpolations(
            cell_face_centroids_over_steps, cell_id, n_faces, interface_x_over_steps,
        )
    end

    return FuelGrainRegressionGeometryInterpolations(
        cell_volumes_interps,
        cell_centroids_interps,
        cell_neighbor_normals_interps,
        cell_neighbor_distances_interps,
        cell_face_areas_interps,
        cell_face_normals_interps,
        cell_face_distances_interps,
        cell_face_centroids_interps,
    )
end

function apply_fuel_grain_regression_geometry_interpolations!(geo, fuel_grain_regression_geometry_interpolations, fuel_port_radius)
    for cell_id in 1:geo.n_cells
        geo.cell_volumes[cell_id] = fuel_grain_regression_geometry_interpolations.cell_volumes_interps[cell_id](fuel_port_radius)
        geo.cell_centroids[cell_id] = fuel_grain_regression_geometry_interpolations.cell_centroids_interps[cell_id](fuel_port_radius)

        geo.cell_neighbor_normals[cell_id, :] .= evaluate_geometry_interpolation(fuel_grain_regression_geometry_interpolations.cell_neighbor_normals_interps[cell_id], fuel_port_radius)
        geo.cell_neighbor_distances[cell_id, :] .= evaluate_geometry_interpolation(fuel_grain_regression_geometry_interpolations.cell_neighbor_distances_interps[cell_id], fuel_port_radius)

        geo.cell_face_areas[cell_id, :] .= evaluate_geometry_interpolation(fuel_grain_regression_geometry_interpolations.cell_face_areas_interps[cell_id], fuel_port_radius)
        geo.cell_face_normals[cell_id, :] .= evaluate_geometry_interpolation(fuel_grain_regression_geometry_interpolations.cell_face_normals_interps[cell_id], fuel_port_radius)
        geo.cell_face_distances[cell_id, :] .= evaluate_geometry_interpolation(fuel_grain_regression_geometry_interpolations.cell_face_distances_interps[cell_id], fuel_port_radius)
        geo.cell_face_centroids[cell_id, :] .= evaluate_geometry_interpolation(fuel_grain_regression_geometry_interpolations.cell_face_centroids_interps[cell_id], fuel_port_radius)
    end
end
