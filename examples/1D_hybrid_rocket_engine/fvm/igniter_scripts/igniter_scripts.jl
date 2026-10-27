function get_igniter_total_volume(geo, grid)
    igniter_volume = 0.0
    
    for igniter_cell_id in getcellset(grid, "igniter_volume")
        igniter_volume += geo.cell_volumes[igniter_cell_id]
    end

    return igniter_volume
end

function get_igniter_total_surface_area(geo, grid)
    igniter_total_surface_area = 0.0
    igniter_cell_ids = getcellset(grid, "igniter_volume")
    
    for face in getfacetset(grid, "fluid_domain_to_igniter_puck_interface")
        cell_id = face[1]
        face_idx = face[2]

        if cell_id in igniter_cell_ids
            igniter_total_surface_area += geo.cell_face_areas[cell_id, face_idx]
        end
    end

    return igniter_volume
end

function setup_igniter_properties!(fluid_initial_conditions, fluid_properties, config, input_igniter_properties)
    igniter_puck_initial_condition = deepcopy(fluid_initial_conditions)
    igniter_puck_properties = deepcopy(fluid_properties)

    igniter_puck_initial_conditions.mass_fractions .= 0.0
    igniter_puck_initial_condition.mass_fractions.nitrogen = 1.0

    igniter_puck_initial_conditions.density = ustrip(upreferred(input_igniter_properties.density))

    igniter_volume = get_igniter_total_volume(geo, grid)

    igniter_volume = igniter_volume * u"m^3"

    initial_igniter_mass = (igniter_puck_initial_conditions.density * u"kg/m^3") * igniter_volume

    config.u_proto.initial_igniter_mass[1] = initial_igniter_mass

    igniter_total_surface_area = get_igniter_total_surface_area(geo, grid) #we'll use this to calculate the fraction of the mass flux that each face receives 

    igniter_system_properties = ComponentVector(
        volume = igniter_volume,
        initial_mass = initial_igniter_mass,
        burn_time = input_igniter_properties.burn_time,
        total_surface_area = igniter_total_surface_area
    )

    igniter_system_properties = ustrip.(upreferred.(igniter_system_properties))

    return igniter_puck_initial_condition, igniter_puck_properties, igniter_system_properties
end

function igniter_mass_burn_rate_at_t(du, u, p, t, system, geo)
    if t < system.igniter.burn_time 
        return system.igniter.initial_mass * (1 - system.igniter.burn_time / t)
    else
        return 0.0
    end
end

#this function will be responsible for taking modifying du.igniter_mass which will be used in a later patch function to 
#add mass and energy to the surrounding cells
#I think it's best to just have the igniter
function update_igniter_combustion!(du, u, p, t, system, geo, cell_id)
    du.igniter_mass[1] = -igniter_mass_burn_rate_at_t(du, u, p, t, system, geo)
end