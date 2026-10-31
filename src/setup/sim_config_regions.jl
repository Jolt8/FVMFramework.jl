function add_region!(
    config, name;
    type,
    initial_conditions,
    properties,
    property_update_function,
    region_function
)
    region_cells = collect(getcellset(config.grid, name))

    for cell_id in region_cells
        for property_name in propertynames(initial_conditions)
            var = getproperty(config.u_proto, property_name)
            initial_condition = getproperty(initial_conditions, property_name)
            if var isa ComponentVector
                for sub_name in propertynames(var)
                    getproperty(var, sub_name)[cell_id] = getproperty(initial_condition, sub_name)
                end
            else
                var[cell_id] = initial_condition
            end
        end
    end

    region = RegionSetupInfo(name, type, initial_conditions, properties, property_update_function, region_function, region_cells)

    all_region_names = [region.name for region in config.regions]

    if name in all_region_names
        existing_region_idx = findfirst(x -> x == name, all_region_names)

        config.regions[existing_region_idx] = region
    else   
        push!(config.regions, region)
    end
    return 
end

#this is for updating the state vector (u0) when doing different experimental trials
function update_region!(config, name)
    all_region_names = [region.name for region in config.regions]
    existing_region_idx = findfirst(x -> x == name, all_region_names)
    
    region = config.regions[existing_region_idx]

    region_cells = collect(getcellset(config.grid, name))

    for cell_id in region_cells
        for property_name in propertynames(region.initial_conditions)
            var = getproperty(config.u_proto, property_name)
            initial_condition = getproperty(region.initial_conditions, property_name)
            if var isa ComponentVector
                for sub_name in propertynames(var)
                    getproperty(var, sub_name)[cell_id] = getproperty(initial_condition, sub_name)
                end
            else
                var[cell_id] = initial_condition
            end
        end
    end

    config.regions[existing_region_idx] = region
end