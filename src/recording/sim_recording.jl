function write_to_vtk_helper!(vtk, named_step, n_cells, geo; prev_field = "", prefix = "", include_zeros_fields = true)
    for field in keys(named_step)
        data = named_step[field]

        #TODO: we may want to have it so that fields where no change occurs are not included
        #this is especially true for the du field

        if data isa ComponentVector || data isa NamedTuple
            if prev_field == ""
                new_field = string(prefix, field)
            else 
                new_field = string(prev_field, "_", field)
            end
            write_to_vtk_helper!(vtk, data, n_cells, geo, prev_field = new_field, prefix = "")
        elseif data isa Vector || data isa SubArray && length(data) == n_cells 
            #the length check is to handle cases where a field is not present on all cells, which write_to_vtk wouldn't accept
            #this is tricky because we can't exactly write stuff to cells that they're not on
            #one way to handle this is to create an array of length n_cells and then fill it with zeros and then write to the respective cell_ids of the field
            #this doens't work for things that are indexed by controller_id for example, which is what the current logic is handling
            if prev_field == ""
                new_field = string(prefix, field)
            else
                new_field = string(prev_field, "_", field)
            end

            if iszero(data) == false || include_zeros_fields == true
                write_cell_data(vtk, data, new_field)
            end
        elseif data isa AbstractMatrix && size(data, 1) == n_cells #this is to handle per-face values
            if prev_field == ""
                new_field = string(prefix, field)
            else
                new_field = string(prev_field, "_", field)
            end
            
            if iszero(data) == false || include_zeros_fields == true
                pos_x = zeros(n_cells)
                neg_x = zeros(n_cells)
                pos_y = zeros(n_cells)
                neg_y = zeros(n_cells)
                pos_z = zeros(n_cells)
                neg_z = zeros(n_cells)
                
                # Construct a vector field using the magnitude of the face flux and the normal
                vector_data = fill(SVector{3, Float64}(0.0, 0.0, 0.0), n_cells)
                
                for cell_id in 1:n_cells
                    for face_idx in 1:size(data, 2)
                        # Use the absolute face normal to find the nearest axis
                        normal = geo.cell_face_normals[cell_id][face_idx]
                        nx, ny, nz = normal[1], normal[2], normal[3]
                        abs_nx, abs_ny, abs_nz = abs(nx), abs(ny), abs(nz)
                        
                        val = data[cell_id, face_idx]
                        
                        vector_data[cell_id] = (
                            vector_data[cell_id][1] + val * nx,
                            vector_data[cell_id][2] + val * ny,
                            vector_data[cell_id][3] + val * nz
                        )
                        
                        if abs_nx >= abs_ny && abs_nx >= abs_nz
                            if nx > 0
                                pos_x[cell_id] += val
                            else
                                neg_x[cell_id] += val
                            end
                        elseif abs_ny >= abs_nx && abs_ny >= abs_nz
                            if ny > 0
                                pos_y[cell_id] += val
                            else
                                neg_y[cell_id] += val
                            end
                        else
                            if nz > 0
                                pos_z[cell_id] += val
                            else
                                neg_z[cell_id] += val
                            end
                        end
                    end
                end
                if iszero(pos_x) == false || iszero(neg_x) == false || include_zeros_fields == true
                    write_cell_data(vtk, pos_x, string(new_field, "_x_pos"))
                    write_cell_data(vtk, neg_x, string(new_field, "_x_neg"))
                end
                if iszero(pos_y) == false || iszero(neg_y) == false || include_zeros_fields == true
                    write_cell_data(vtk, pos_y, string(new_field, "_y_pos"))
                    write_cell_data(vtk, neg_y, string(new_field, "_y_neg"))
                end
                if iszero(pos_z) == false || iszero(neg_z) == false || include_zeros_fields == true
                    write_cell_data(vtk, pos_z, string(new_field, "_z_pos"))
                    write_cell_data(vtk, neg_z, string(new_field, "_z_neg"))
                end
                write_cell_data(vtk, vector_data, string(new_field, "_vec"))
            end
        else
            #hmm, it appears that's there's no way to write cell face values
            println("not processed: ", field)
        end
    end
end

function sol_to_vtk(sol, du_named, u_named, grid, geo, sim_file, root_dir; include_zeros_fields = true, track_progress = false)
    date_and_time = Dates.format(now(), "I.MM.SS p yyyy-mm-dd")
    #date_and_time = Dates.format(now(), "I.MM.SS p")

    #root_dir = "C://Users//wille//Desktop//Julia_cfd_output_files"

    project_name = replace(basename(sim_file), r".jl" => "")

    sim_folder_name = project_name * " " * date_and_time

    output_dir = joinpath(root_dir, sim_folder_name)

    mkpath(output_dir)

    pvd_filename = joinpath(output_dir, "solution_collection")

    pvd = paraview_collection(pvd_filename)

    step_filename = joinpath(output_dir, "timestep")

    #this step may become a significant bottleneck
    #update: This is a very big bottleneck
    #we get a component array of the simulation's state at this step

    n_cells = length(grid.cells)

    for (step, t) in enumerate(sol.t)
        if track_progress == true
            println("saving timestep $t of $(sol.t[end]), remaining steps = $(length(sol.t)-step)")
        end
        VTKGridFile(step_filename * " $step" * " at $t.vtu", grid) do vtk
            write_to_vtk_helper!(vtk, u_named[step], n_cells, geo; prev_field = "", prefix = "u_", include_zeros_fields = include_zeros_fields)
            write_to_vtk_helper!(vtk, du_named[step], n_cells, geo; prev_field = "", prefix = "du_", include_zeros_fields = include_zeros_fields)
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
end