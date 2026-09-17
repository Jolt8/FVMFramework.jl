import LinearAlgebra

"""
    WeightedLeastSquaresStencil

Precomputed coefficients for cell-centred weighted least-squares gradients.
`coefficients[cell_id][:, j]` multiplies the difference between the value in
`neighbor_ids[cell_id][j]` and the value in `cell_id`.
"""
struct WeightedLeastSquaresStencil{T}
    neighbor_ids::Vector{Vector{Int}}
    coefficients::Vector{Matrix{T}}
end

"""
    build_weighted_least_squares_stencil(geo; weight_power=2)

Build a weighted least-squares gradient stencil from cell centroids and
face-neighbour connectivity. A neighbour at distance `r` receives weight
`1 / r^weight_power`.

The weighted displacement matrix is inverted with a pseudoinverse. This makes
the reconstruction usable for lower-dimensional meshes embedded in 3D (for
example, the 100 x 1 x 1 mesh in this example): unresolved gradient components
are set to zero instead of making the normal equations singular.
"""
function build_weighted_least_squares_stencil(geo; weight_power=2)
    weight_power >= 0 || throw(ArgumentError("weight_power must be non-negative"))

    n_cells = length(geo.cell_centroids)
    coordinate_type = float(eltype(eltype(geo.cell_centroids)))
    neighbor_ids = [Int[] for _ in 1:n_cells]
    coefficients = [zeros(coordinate_type, 3, 0) for _ in 1:n_cells]

    for (cell_id, neighbors_with_faces) in geo.cell_neighbors
        ids = Int[neighbor_id for (neighbor_id, _) in neighbors_with_faces if neighbor_id > 0]
        neighbor_ids[cell_id] = ids

        isempty(ids) && continue

        weighted_displacements = zeros(coordinate_type, length(ids), 3)
        square_root_weights = zeros(coordinate_type, length(ids))
        center = geo.cell_centroids[cell_id]

        for (row, neighbor_id) in enumerate(ids)
            displacement = geo.cell_centroids[neighbor_id] - center
            distance = LinearAlgebra.norm(displacement)
            iszero(distance) && throw(ArgumentError(
                "cells $cell_id and $neighbor_id have coincident centroids",
            ))

            square_root_weight = distance^(-weight_power / 2)
            square_root_weights[row] = square_root_weight

            for dimension in 1:3
                weighted_displacements[row, dimension] =
                    square_root_weight * displacement[dimension]
            end
        end

        # If A contains centroid displacements and W contains the weights, the
        # gradient is (sqrt(W) * A)^+ * sqrt(W) * delta_phi.
        coefficients[cell_id] =
            LinearAlgebra.pinv(weighted_displacements) *
            LinearAlgebra.Diagonal(square_root_weights)
    end

    return WeightedLeastSquaresStencil(neighbor_ids, coefficients)
end

"""
    populate_weighted_least_squares_gradient!(gradient, values, stencil)

Populate an `n_cells x 3` gradient cache for one cell-centred scalar field.
The calculation only uses additions and scalar multiplications at runtime so
it remains compatible with automatic differentiation and sparsity tracing.
"""
function populate_weighted_least_squares_gradient!(gradient, values, stencil)
    length(values) == length(stencil.neighbor_ids) || throw(DimensionMismatch(
        "the value field and weighted least-squares stencil have different cell counts",
    ))
    size(gradient) == (length(values), 3) || throw(DimensionMismatch(
        "the gradient cache must have size (n_cells, 3)",
    ))

    for cell_id in eachindex(stencil.neighbor_ids)
        center_value = values[cell_id]
        coefficient_matrix = stencil.coefficients[cell_id]

        gradient_x = zero(center_value)
        gradient_y = zero(center_value)
        gradient_z = zero(center_value)

        for (column, neighbor_id) in enumerate(stencil.neighbor_ids[cell_id])
            value_difference = values[neighbor_id] - center_value
            gradient_x += coefficient_matrix[1, column] * value_difference
            gradient_y += coefficient_matrix[2, column] * value_difference
            gradient_z += coefficient_matrix[3, column] * value_difference
        end

        gradient[cell_id, 1] = gradient_x
        gradient[cell_id, 2] = gradient_y
        gradient[cell_id, 3] = gradient_z
    end

    return gradient
end

#this will be treated as a flux function
function corrected_face_gradient(
    phi_a,
    phi_b,
    grad_a,
    grad_b,
    normal,
    distance
)
    grad_avg = 0.5 * (grad_a + grad_b)

    normal_derivative =
        (phi_b - phi_a) / distance

    return grad_avg +
        (
            normal_derivative -
            dot(grad_avg, normal)
        ) * normal
end

function populate_weighted_least_squares_face_values!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals, cell_face_distances,
    cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    u.grad_vel_u_face[idx_a, face_idx] =
    corrected_face_gradient(
        u.vel_u[idx_a],
        u.vel_u[idx_b],
        u.grad_vel_u[idx_a],
        u.grad_vel_u[idx_b],
        cell_face_normals[idx_a][face_idx],
        cell_face_distances[idx_a][face_idx]
    )

    u.grad_vel_v_face[idx_a, face_idx] =
    corrected_face_gradient(
        u.vel_v[idx_a],
        u.vel_v[idx_b],
        u.grad_vel_v[idx_a],
        u.grad_vel_v[idx_b],
        cell_face_normals[idx_a][face_idx],
        cell_face_distances[idx_a][face_idx]
    )

    u.grad_vel_w_face[idx_a, face_idx] =
    corrected_face_gradient(
        u.vel_w[idx_a],
        u.vel_w[idx_b],
        u.grad_vel_w[idx_a],
        u.grad_vel_w[idx_b],
        cell_face_normals[idx_a][face_idx],
        cell_face_distances[idx_a][face_idx]
    )

    u.grad_temperature_face[idx_a, face_idx] =
    corrected_face_gradient(
        u.temperature[idx_a],
        u.temperature[idx_b],
        u.grad_temperature[idx_a],
        u.grad_temperature[idx_b],
        cell_face_normals[idx_a][face_idx],
        cell_face_distances[idx_a][face_idx]
    )
end

"""
    update_weighted_least_squares_gradients!(u, stencil)

Populate all primitive-variable gradients used by the Navier-Stokes example.
Velocity and temperature must already have been updated for every cell.
"""
function update_weighted_least_squares_gradients!(u, stencil)
    populate_weighted_least_squares_gradient!(u.grad_vel_u, u.vel_u, stencil)
    populate_weighted_least_squares_gradient!(u.grad_vel_v, u.vel_v, stencil)
    populate_weighted_least_squares_gradient!(u.grad_vel_w, u.vel_w, stencil)
    populate_weighted_least_squares_gradient!(u.grad_temperature, u.temperature, stencil)
    return nothing
end
