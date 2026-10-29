function _build_blank_dict!(current_dict, properties, n_cells)
    for property_name in propertynames(properties)
        value = getproperty(properties, property_name)
        if value isa ComponentArray || value isa NamedTuple
            if !haskey(current_dict, property_name)
                current_dict[property_name] = OrderedDict{Symbol, Any}()
            end
            _build_blank_dict!(current_dict[property_name], value, n_cells)
        else
            if !haskey(current_dict, property_name)
                if value isa AbstractArray
                    current_dict[property_name] = zeros(n_cells) .* upreferred(unit(value[1]))
                else
                    current_dict[property_name] = zeros(n_cells) .* upreferred(unit(value))
                end
            end
        end
    end
end


function _drill_down_and_fill_properties!(merged_properties, properties, cells)
    for property_name in propertynames(properties)
        if merged_properties[property_name] isa ComponentArray
            _drill_down_and_fill_properties!(getproperty(merged_properties, property_name), getproperty(properties, property_name), cells)
        else
            if merged_properties[property_name] isa AbstractArray 
                if properties[property_name] isa AbstractArray
                    for cell_id in cells
                        view(merged_properties, property_name)[cell_id] = properties[property_name][cell_id]
                    end
                else
                    for cell_id in cells
                        view(merged_properties, property_name)[cell_id] = properties[property_name]
                    end
                end
            elseif merged_properties[property_name] isa Number
                merged_properties[property_name] = properties[property_name]
            else
                error("The property $property_name was not handled")
            end
        end
    end
end

function _drill_down_and_fill_patch_properties!(merged_properties, properties, cells)
    for property_name in propertynames(properties)
        if merged_properties[property_name] isa ComponentArray
            _drill_down_and_fill_patch_properties!(getproperty(merged_properties, property_name), getproperty(properties, property_name), cells)
        else
            if merged_properties[property_name] isa AbstractArray
                if properties[property_name] isa AbstractArray
                    for cell_id in cells
                        view(merged_properties, property_name)[cell_id] = properties[property_name][cell_id]
                    end
                else
                    for cell_id in cells
                        view(merged_properties, property_name)[cell_id] = properties[property_name]
                    end
                end
            elseif merged_properties[property_name] isa Number
                merged_properties[property_name] = properties[property_name]
            else
                error("The property $property_name was not handled")
            end
        end
    end
end
#I don't think this is actually necessary, in fact, it's causing a lot of problems

function _dict_to_namedtuple(dict::OrderedDict{Symbol, Any})
    NamedTuple(k => (v isa OrderedDict{Symbol, Any} ? _dict_to_namedtuple(v) : v) for (k, v) in dict)
end


function merge_region_properties(config)
    n_cells = length(config.geo.cell_volumes)

    prop_dict = OrderedDict{Symbol, Any}()

    #convert component arrays to named tuples
    for region in config.regions
        _build_blank_dict!(prop_dict, region.properties, n_cells)
    end

    for patch in config.patches
        _build_blank_dict!(prop_dict, patch.properties, n_cells)
    end
    
    merged_properties = ComponentVector(_dict_to_namedtuple(prop_dict)) #dict to named tuple seems to be required here, but it breaks caches

    for region in config.regions
        _drill_down_and_fill_properties!(merged_properties, region.properties, region.region_cells)
    end

    for patch in config.patches
        _drill_down_and_fill_patch_properties!(merged_properties, patch.properties, [cell[1] for cell in patch.cell_neighbors])
    end

    return merged_properties
end


function _drill_down_and_fill_caches!(merged_caches, region_cache_syms_and_units, special_caches, merged_properties)
    for property_name in propertynames(region_cache_syms_and_units)
        property_unit = getproperty(region_cache_syms_and_units, property_name)
        
        if property_name in keys(special_caches)
            
            #special caches are already handled during initialization in cache_dict
        elseif merged_caches[property_name] isa ComponentArray
            #recursively drill down both caches and properties if possible
            if hasproperty(merged_properties, property_name)
                _drill_down_and_fill_caches!(getproperty(merged_caches, property_name), property_unit, special_caches, getproperty(merged_properties, property_name))
            else
                #if properties are missing for this group, still drill down for units but with empty properties
                _drill_down_and_fill_caches!(getproperty(merged_caches, property_name), property_unit, special_caches, (_ = nothing,))
            end
        elseif merged_caches[property_name] isa AbstractArray || merged_caches[property_name] isa SubArray
            # since we use the properties as the initial value for the cache, we might not even need properties and could just preemptively merge them
            if hasproperty(merged_properties, property_name)
                view(merged_caches, property_name) .= merged_properties[property_name]
            else
                view(merged_caches, property_name) .= 0.0 .* upreferred(property_unit)
            end
        else
            error("The cache $property_name was not handled")
        end
    end
end

function _build_blank_cache_dict!(current_dict, cache_syms_and_units, special_caches, n_cells)
    for cache_name in propertynames(cache_syms_and_units)
        cache_unit = getproperty(cache_syms_and_units, cache_name)
        
        if cache_name in keys(special_caches)
            if !haskey(current_dict, cache_name)
                current_dict[cache_name] = special_caches[cache_name]
            end
        elseif !haskey(current_dict, cache_name)
            if cache_unit isa ComponentArray || cache_unit isa NamedTuple
                current_dict[cache_name] = OrderedDict{Symbol, Any}()
                _build_blank_cache_dict!(current_dict[cache_name], cache_unit, special_caches, n_cells)
            else
                current_dict[cache_name] = zeros(n_cells) .* upreferred(cache_unit)
            end
        end
    end
end

function merge_region_caches(config, merged_properties)
    n_cells = length(config.geo.cell_volumes)

    cache_dict = OrderedDict{Symbol, Any}()

    for region in config.regions
        _build_blank_cache_dict!(cache_dict, config.cache_syms_and_units, config.special_caches, n_cells)
    end

    merged_caches = ComponentVector(_dict_to_namedtuple(cache_dict))

    # Fill all regions from the merged_properties
    for region in config.regions
        _drill_down_and_fill_caches!(merged_caches, config.cache_syms_and_units, config.special_caches, merged_properties)
    end

    return merged_caches
end