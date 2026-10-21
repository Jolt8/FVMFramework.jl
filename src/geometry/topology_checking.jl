function check_cellset_connectivity(grid, set_name)
    #returns a dict of the amount of total connections between the set_name and all other cellsets
    #e.g. Dict("set_name" => n_connections)

    matched_cell_set_times = Dict{String, Int}()

    top = ExclusiveTopology(grid)

    for cell_id in grid.cellsets[set_name]
        cellset_names = collect(keys(grid.cellsets))
        test_set = filter(x -> x != set_name, cellset_names)
        cell_neighbors = getneighborhood(top, grid, CellIndex(cell_id))

        for set_name_j in test_set
            if any(in(cell_neighbors), grid.cellsets[set_name_j])
                if !(set_name_j in keys(matched_cell_set_times))
                    #if we haven't the cellset before, create a new entry
                    push!(matched_cell_set_times, set_name_j => 1) 
                else
                    matched_cell_set_times[set_name_j] += 1
                end
            end
        end
    end
    return matched_cell_set_times
end
    

#=NOTE:
    - When using freecad and exporting a booleanfragments, make sure to add:
        Mesh 3; 
        Coherence Mesh; <--- this one!
        Save "output.msh";
    - because for some reason Ferrite would not recognize the cells as connecting because there were duplicate
      nodes at the same coordinate on the same "shared" surface of two separate physical connection_groups
    - This was such a headache, I wonder if anyone else has struggled with this
    - Also, why is this necessary, what did I do to require doing this 
    - Other things I tried that did not seem to help much
        - switching to CompSolid instead of Standard in freecad booleanfragments
        - putting Coherence in my .geo file
        - Using a compound filter in freecad
    - if you want to check for this in any future gmsh file, just turn on node labels once meshed and check for duplicate node labels
    - I'm also angry because I'm realizing that this is why my initial endeavor with OpenFOAM failed because of this exact issue
=#