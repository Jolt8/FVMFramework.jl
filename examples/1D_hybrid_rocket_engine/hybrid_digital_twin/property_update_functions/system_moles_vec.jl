#=
Per-cell species vectors (moles and mole fractions) that Clapeyron needs as plain AbstractVectors.

The values of every cell are stored in one flat vector laid out cell by cell:
    [cell 1: species 1..n_species | cell 2: species 1..n_species | ...]
cell_species_ranges[cell_id] gives the indexes of that cell's species in the flat vector.

The flat vector comes from a DiffCache, so it can hold Float64, ForwardDiff.Dual, or SparseConnectivityTracer numbers.
Once per RHS call (at the top of solve_system!), select_cell_species_buffer! pulls the buffer matching the current
number type out of the DiffCache and stores it in active_flat_values. After that, system.moles_vec[cell_id] is just a
view into that buffer, so writes like system.moles_vec[cell_id][species_idx] = value go straight into it.
=#

mutable struct CellSpeciesVectors
    diff_cache::DiffCache
    cell_species_ranges::Vector{UnitRange{Int}}
    active_flat_values::AbstractVector
end

function CellSpeciesVectors(n_cells::Int, n_species::Int, jacobian_chunk_size::Int)
    flat_values = zeros(n_cells * n_species)

    diff_cache = DiffCache(flat_values, jacobian_chunk_size)

    cell_species_ranges = [((cell_id - 1) * n_species + 1):(cell_id * n_species) for cell_id in 1:n_cells]

    return CellSpeciesVectors(diff_cache, cell_species_ranges, flat_values)
end

#number_sample must have the number type of the current RHS call (e.g. u.density[1] + first(p))
function select_cell_species_buffer!(species_vectors::CellSpeciesVectors, number_sample)
    species_vectors.active_flat_values = get_tmp(species_vectors.diff_cache, number_sample)
end

Base.getindex(species_vectors::CellSpeciesVectors, cell_id::Int) = view(species_vectors.active_flat_values, species_vectors.cell_species_ranges[cell_id])


function update_cell_moles_vec!(du, u, p, t, system, geo, cell_id)
    # species_idx is a Ref because a plain Int that is reassigned inside the do-block closure gets wrapped in
    # an untyped Core.Box, which allocates and dispatches dynamically on every species
    species_idx = Ref(1)

    # u.cell_volume, not geo.cell_volumes, so the port gas cells use their current (growing) volume
    cell_mass = u.density[cell_id] * u.cell_volume[cell_id]

    cell_moles = system.moles_vec[cell_id]

    for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
        cell_moles[species_idx[]] = (u_mass_fractions[species[cell_id]] * cell_mass) / u_molecular_weights[species[cell_id]]

        species_idx[] += 1
    end
end

function update_cell_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    species_idx = Ref(1)

    cell_mole_fractions = system.mole_fractions_vec[cell_id]

    for_fields!(u.mass_fractions, u.molecular_weights) do species, u_mass_fractions, u_molecular_weights
        cell_mole_fractions[species_idx[]] = u_mass_fractions[species[cell_id]] * (u.mw_avg[cell_id] / u_molecular_weights[species[cell_id]])

        species_idx[] += 1
    end
end

function update_cell_moles_and_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
    update_cell_moles_vec!(du, u, p, t, system, geo, cell_id)
    update_cell_mole_fractions_vec!(du, u, p, t, system, geo, cell_id)
end
