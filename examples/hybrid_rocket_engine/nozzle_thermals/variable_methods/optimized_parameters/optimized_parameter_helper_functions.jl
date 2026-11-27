function add_optimized_parameters_to_properties(properties, optimized_parameters)
    optimized_parameters_dict = OrderedDict()
    
    for name in propertynames(optimized_parameters)
        optimized_parameter = getproperty(optimized_parameters, name)

        optimized_parameters_dict[name] = optimized_parameter.guess
    end

    properties = merge_properties(
        properties,
        ComponentVector(optimized_parameters_dict)
    )

    return properties
end