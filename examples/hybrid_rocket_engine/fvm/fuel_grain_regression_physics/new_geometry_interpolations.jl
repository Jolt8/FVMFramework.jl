struct AxialFace
    node_ids_on_z_plane::Vector{Int}
    fluid_cell_id::Int
    fluid_face_idx_on_z_plane::Int
    fluid_face_idx_facing_toward_fuel_wall::Int #we use this variable to determine moving face velocity
    fluid_face_idx_facing_away_from_fuel_wall::Int
    closest_node_ids_from_fuel_wall::Vector{Int}
    furthest_node_ids_from_fuel_wall::Vector{Int}
end

struct AxialPosition
    z_pos::Float64
    sorted_upstream_fluid_cell_entries::Vector{AxialFace} #sorted from closest to furthest from the z axis/center of the motor (furthest to closest from the fuel grain wall)
    sorted_downstream_fluid_cell_entries::Vector{AxialFace} #sorted from closest to furthest from the z axis/center of the motor (furthest to closest from the fuel grain wall)
    sorted_downstream_fuel_cell_entries::Vector{AxialFace} #sorted from closest to furthest from the z axis/center of the motor (furthest to closest from the fuel grain wall)
    
    min_fluid_x_bound::Float64
    initial_fuel_grain_x_bound::Float64
    initial_fuel_grain_thickness::Float64
    max_fuel_grain_x_bound::Float64
    #initial_port_diameter = max_fuel_grain_x_bound - initial_fuel_grain_x_bound
    #current_port_diameter = max_fuel_grain_x_bound - current_fuel_wall_thickness
    #we can find how much each node needs to be translated based on (current_port_diameter / initial_port_diameter) * (max_fuel_grain_x_bound - initial_node_coordinates_z_pos)

    node_ids::Vector{Int} #all node_ids on this specific axial plane that only includes the fluid_domain and the fuel_grain

    #initial_node_positions::Dict{Int, SVector{3, Float64}} #instead we'll add all the initial_node_coordinates to system

    #fluid_fuel_wall_interface_face_idx::Int #the id of the face of the fluid cell that shares an interface with the fuel_wall
    #fuel_wall_fluid_interface_face_idx::Int #the id of the face of the fuel cell that shares an interface with the fluid domain
    #oh these are both redundant because we can get them by getting the last fluid cells in the sorted list's face id on the z plane or vice versa.
end


function create_axial_position_structs(grid, geo)
end

function get_axial_oxidizer_mass_flux!(du, u, p, t, system, geo, axial_position, axial_position_id)
    nitrous_oxide_mass_flow = 0.0

    port_cross_sectional_area = pi * (axial_position.max_fuel_grain_x_bound - u.fuel_grain_thickness[axial_position_id])^2

    for upstream_face_entry in axial_position.sorted_upstream_fluid_cell_entries
        upstream_cell_id = upstream_face_entry.fluid_cell_id
        upstream_face_idx_on_z_plane = upstream_face_entry.fluid_face_idx_on_z_plane
        
        vel_face = abs(u.vel_w_face[upstream_cell_id, upstream_face_idx_on_z_plane])
        area = geo.cell_face_areas[upstream_cell_id, upstream_face_idx_on_z_plane]
        nitrous_oxide_density = u.species_densities.nitrous_oxide[upstream_cell_id]

        nitrous_oxide_mass_flow += (vel_face * area * nitrous_oxide_density)
    end

    nitrous_oxide_mass_flux = nitrous_oxide_mass_flow / port_cross_sectional_area

    return nitrous_oxide_mass_flux
end

function move_fuel_wall!(du, u, p, t, system, geo)
    for (axial_position_id, axial_position) in enumerate(system.axial_positions)
        initial_port_diameter = axial_position.max_fuel_grain_x_bound - axial_position.initial_fuel_grain_thickness
        current_port_diameter = axial_position.max_fuel_grain_x_bound - u.fuel_grain_thickness[axial_position_id]
        for node_id in axial_position.node_ids
            geo.node_coordinates[node_id][1] = (current_port_diameter / initial_port_diameter) * (axial_position.max_fuel_grain_x_bound - system.initial_node_coordinates[node_id][1])
        end
    end
end

#NOTE: I think we're going to have to store a copy of the last timestep's geo struct

function get_cell_row_fluid_and_fuel_volume(du, u, p, t, system, geo, cell_row_id)
    axial_position_1 = system.axial_positions[cell_row_id] #note that this assumes that each cell_row_id corresponds to the axial_position closer to the top of the motor
    axial_position_2 = system.axial_positions[cell_row_id+1]

    total_fluid_cell_row_volume = 0.0

    for fluid_cell_info in axial_position_1.sorted_downstream_fluid_cell_entries
        fluid_cell_id = fluid_cell_info.fluid_cell_id

        total_fluid_cell_row_volume += geo.cell_volumes[fluid_cell_id]
    end

    total_fuel_cell_row_volume = 0.0

    for fuel_cell_info in axial_position_1.sorted_downstream_fuel_cell_entries
        fuel_cell_id = fuel_cell_info.fluid_cell_id

        total_fuel_cell_row_volume += geo.cell_volumes[fuel_cell_id]
    end

    return total_fluid_cell_row_volume, total_fuel_cell_row_volume
end

function update_previous_cell_volumes!(system, geo)
    n_cell_rows = length(system.axial_positions[1:end-1])

    for cell_row_id in 1:n_cell_rows
        axial_position_1 = system.axial_positions[cell_row_id]
        #axial_position_2 = system.axial_positions[cell_row_id+1]

        for fluid_cell_info in axial_position_1.sorted_downstream_fluid_cell_entries
            fluid_cell_id = fluid_cell_info.fluid_cell_id

            system.previous_cell_volumes[fluid_cell_id] = geo.cell_volumes[fluid_cell_id]
        end

        for fuel_cell_info in axial_position_1.sorted_downstream_fuel_cell_entries
            fuel_cell_id = fuel_cell_info.fluid_cell_id

            system.previous_cell_volumes[fuel_cell_id] = geo.cell_volumes[fuel_cell_id]
        end
    end
end

function update_previous_face_centroid_positions!(system, geo)
    n_cell_rows = length(system.axial_positions[1:end-1])

    for cell_row_id in 1:n_cell_rows
        axial_position_1 = system.axial_positions[cell_row_id]
        #axial_position_2 = system.axial_positions[cell_row_id+1]

        for (radial_position_id, fluid_cell_info) in enumerate(axial_position_1.sorted_downstream_fluid_cell_entries)
            fluid_cell_id = fluid_cell_info.fluid_cell_id

            system.previous_cell_centroid_positions[fluid_cell_id] = geo.face_centroids[fluid_cell_id, fluid_cell_info.fluid_face_idx_facing_toward_fuel_wall]
        end

        for (radial_position_id, fuel_cell_info) in enumerate(axial_position_1.sorted_downstream_fuel_cell_entries)
            fuel_cell_id = fuel_cell_info.fluid_cell_id

            system.previous_cell_centroid_positions[fluid_cell_id] = geo.face_centroids[fuel_cell_id, fuel_cell_info.fuel_face_idx_facing_toward_fluid]
        end
    end
end

function get_closest_and_furthest_face_x_velocities!(du, u, p, t, system, geo, axial_position_id_1, axial_position_id_2, radial_position_id)
    fluid_cell_info = axial_position.sorted_downstream_fluid_cell_entries[radial_position_id]

    axial_position_1_node_id = fluid_cell_info.furthest_node_ids_from_fuel_wall[1]
    axial_position_2_node_id = fluid_cell_info.closest_node_ids_from_fuel_wall[1]

    axial_position_1_velocity = du.fuel_grain_thickness[axial_position_id_1]
    axial_position_2_velocity = du.fuel_grain_thickness[axial_position_id_2]

    cell_node_velocity_1 = (axial_position.max_fuel_grain_x_bound - system.initial_node_coordinates[axial_position_1_node_id][1]) * axial_position_1_velocity
    cell_node_velocity_2 = (axial_position.max_fuel_grain_x_bound - system.initial_node_coordinates[axial_position_2_node_id][1]) * axial_position_2_velocity

    moving_face_vel_x = (cell_node_velocity_1 + cell_node_velocity_2) / 2

    return moving_face_vel_x
end

function get_face_closer_to_fuel_velocity_x(du, u, p, t, system, geo, axial_position_id_1, axial_position_id_2, radial_position_id)
    fluid_cell_info = axial_position.sorted_downstream_fluid_cell_entries[radial_position_id]
    fluid_cell_id = fluid_cell_info.fluid_cell_id

    axial_position_1_node_id = fluid_cell_info.furthest_node_ids_from_fuel_wall[1]
    axial_position_2_node_id = fluid_cell_info.closest_node_ids_from_fuel_wall[1]

    axial_position_1_velocity = du.fuel_grain_thickness[axial_position_id_1]
    axial_position_2_velocity = du.fuel_grain_thickness[axial_position_id_2]

    cell_node_velocity_1 = (axial_position.max_fuel_grain_x_bound - system.initial_node_coordinates[axial_position_1_node_id][1]) * axial_position_1_velocity
    cell_node_velocity_2 = (axial_position.max_fuel_grain_x_bound - system.initial_node_coordinates[axial_position_2_node_id][1]) * axial_position_2_velocity

    moving_face_vel_x = (cell_node_velocity_1 + cell_node_velocity_2) / 2

    return moving_face_vel_x
end


#NOTE: this now requires previous time

function regress_fuel_layers!(du, u, p, t, system, geo)
    for (axial_position_id, axial_position) in enumerate(system.axial_positions)
        nitrous_oxide_mass_flux = get_axial_oxidizer_mass_flux!(du, u, p, t, system, geo, axial_position, axial_position_id)

        regression_rate = p.regression_coeff_a * (nitrous_oxide_mass_flux ^ p.regression_coeff_n)

        du.fuel_grain_thickness[axial_position_id] = -regression_rate
    end

    #note that this requires system to store the previous geometry
    for (axial_position_id_1, axial_position_1) in enumerate(system.axial_positions[1:end-1])
        axial_position_2 = system.axial_positions[axial_position_id_1+1]
    end

    #loop through each axial_position pair and use that for updating the downstream fluid cell
    #note that this requires system to store the previous geometry
    n_cell_rows = length(system.axial_positions[1:end-1])

    dt = t - system.prev_t

    for cell_row_id in 1:n_cell_rows
        axial_position_1 = system.axial_positions[cell_row_id]
        axial_position_2 = system.axial_positions[cell_row_id+1]

        #update the density of cells so that mass is conserved
        for (radial_position_id, fluid_cell_info) in enumerate(axial_position_1.sorted_downstream_fluid_cell_entries)
            fluid_cell_id = fluid_cell_info.fluid_cell_id
            face_idx_facing_toward_fuel_wall = fluid_cell_info.fluid_face_idx_facing_toward_fuel_wall #face_idx facing the fuel grain
            face_idx_facing_away_from_fuel_wall = fluid_cell_info.fluid_face_idx_facing_away_from_fuel_wall

            moving_face_vel_x = get_moving_fluid_face_velocity_x(
                du, u, p, t, system, geo, 
                axial_position_id_1, axial_position_id_2, radial_position_id
            )

            face_volume_flux = moving_face_vel_x * geo.face_normals[fluid_cell_id, face_idx][1] * geo.face_areas[fluid_cell_id, face_idx]

            du.density[fluid_cell_id] -= u.density[fluid_cell_id] * face_volume_flux / geo.cell_volumes[fluid_cell_id]

            du.momentum_density_u[fluid_cell_id] -= u.momentum_density_u[fluid_cell_id] * face_volume_flux
            du.momentum_density_v[fluid_cell_id] -= u.momentum_density_v[fluid_cell_id] * face_volume_flux
            du.momentum_density_w[fluid_cell_id] -= u.momentum_density_w[fluid_cell_id] * face_volume_flux

            du.volumetric_energy[fluid_cell_id] -= u.volumetric_energy[fluid_cell_id] * face_volume_flux

            for_fields!(du.species_densities, u.species_densities) do species, du_species_densities, u_species_densities
                du_species_densities[species[fluid_cell_id]] -= u_species_densities[species[fluid_cell_id]] * face_volume_flux
            end
        end

        #We don't do the same as above to the fuel because we want its mass to decrease (which will be added to the fluid later)

        fuel_cell_id = axial_position_2.downstream_fuel_cell_ids[1] #the fuel cell id corresponding to the fuel grain closest to the top of the motor in this cell row
        #we'll use this cell_id to get the relevant fuel properties

        fluid_info = axial_position_1.sorted_downstream_fluid_cell_entries[1]
        fluid_cell_id = fluid_info.fluid_cell_id #the fluid cell id corresponding to the fluid cell closest to the bottom of the motor in this cell row
        #we'll use this cell_id to get the relevant fluid properties and to apply a flux to this cell_id
        
        current_fluid_cell_row_volume, current_fuel_cell_row_volume = get_cell_row_fluid_and_fuel_volume!(du, u, p, t, system, geo, cell_row_id)

        #then we just apply a flux to this current fluid cell based on the change in volume

        du_fuel_mass = -du.fuel_grain_thickness[fuel_cell_id] * geo.face_areas[fluid_cell_id, fluid_info.fluid_face_idx_facing_toward_fuel_wall] * u.density[fuel_cell_id]

        du_heat = du_fuel_mass * u.fuel_heat_of_pyrolysis[fuel_cell_id] + du_fuel_mass * u.fuel_heat_of_vaporization[fuel_cell_id]

        #This assumption that the heat loss is applied to both is slightly arbitrary although I'm sure it will work fine
        du.volumetric_energy[fuel_cell_id] -= (du_heat / 2) / geo.cell_volumes[fuel_cell_id]
        du.volumetric_energy[fluid_cell_id] -= (du_heat / 2) / geo.cell_volumes[fluid_cell_id]
        
        du.density[fluid_cell_id] += du_fuel_mass / geo.cell_volumes[fluid_cell_id]
        du.species_densities[fluid_cell_id].ethylene += du_fuel_mass / geo.cell_volumes[fluid_cell_id]
    end

    system.previous_total_fluid_cell_row_volumes, system.previous_total_fuel_cell_row_volumes = get_cell_row_fluid_and_fuel_volume(du, u, p, t, system, geo, cell_row_id)
    update_previous_cell_volumes!(system, geo)

    move_fuel_wall!(du, u, p, t, system, geo)
end


        

