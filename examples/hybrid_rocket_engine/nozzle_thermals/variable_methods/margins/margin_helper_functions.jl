function add_margins_to_properties(properties, margins)
    margins_dict = OrderedDict()
    
    for name in propertynames(margins)
        margin = getproperty(margins, name)

        if margin.type == :sized_to
            margins_dict[name] = OrderedDict(:desired => margin.desired)
        elseif margin.type == :rated_against
            margins_dict[name] = OrderedDict(:desired => margin.desired, :found => zero(margin.desired))
        else
            error("Margin $(margin) has an invalid type")
        end
    end

    properties = merge_properties(
        properties,
        ComponentVector(margins_dict)
    )

    return properties
end
