
function add_uncertainties_to_properties(properties, uncertainties)
    model_uncertainties = OrderedDict()
    
    for name in propertynames(uncertainties)
        uncertainty = getproperty(uncertainties, name)

        model_uncertainties[name] = uncertainty.nominal
    end

    properties = merge_properties(
        properties,
        ComponentVector(model_uncertainties)
    )

    return properties
end

# modified p functions:
function get_nominal_p(p, uncertainties)
    new_p = deepcopy(p)

    for name in propertynames(uncertainties)
        setproperty!(new_p, name, getproperty(uncertainties, name).nominal)
    end

    return new_p
end

function get_worst_case_p(p, uncertainties)
    new_p = deepcopy(p)

    for name in propertynames(uncertainties)
        if getproperty(uncertainties, name).worst_case == :high
            setproperty!(new_p, name, getproperty(uncertainties, name).high)
        elseif getproperty(uncertainties, name).worst_case == :low
            setproperty!(new_p, name, getproperty(uncertainties, name).low)
        else
            setproperty!(new_p, name, getproperty(uncertainties, name).nominal)
        end
    end

    return new_p
end

function get_best_case_p(p, uncertainties)
    new_p = deepcopy(p)

    for name in propertynames(uncertainties)
        if getproperty(uncertainties, name).worst_case == :high
            setproperty!(new_p, name, getproperty(uncertainties, name).low)
        elseif getproperty(uncertainties, name).worst_case == :low
            setproperty!(new_p, name, getproperty(uncertainties, name).high)
        else
            setproperty!(new_p, name, getproperty(uncertainties, name).nominal)
        end
    end

    return new_p
end