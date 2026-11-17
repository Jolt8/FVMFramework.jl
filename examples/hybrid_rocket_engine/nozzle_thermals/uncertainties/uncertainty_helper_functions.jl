
function add_uncertainties_to_properties(properties, uncertainties)
    model_uncertainties = Dict()
    
    for name in propertynames(uncertainties)
        uncertainty = getproperty(uncertainties, name)
        if uncertainty.type == :model
            model_uncertainties[name] = zero(uncertainty.nominal)
        elseif uncertainty.type == :parameter
            #do nothing here for now
        end
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
        getproperty(new_p, name) = getproperty(uncertainties, name).nominal
    end

    return new_p
end

function get_worst_case_p(p, uncertainties)
    new_p = deepcopy(p)

    for name in propertynames(uncertainties)
        if getproperty(uncertainties, name).worst_case == :high
            getproperty(new_p, name) = getproperty(uncertainties, name).high
        elseif getproperty(uncertainties, name).worst_case == :low
            getproperty(new_p, name) = getproperty(uncertainties, name).low
        else
            getproperty(new_p, name) = getproperty(uncertainties, name).nominal
        end
    end

    return new_p
end

function get_best_case_p(p, uncertainties)
    new_p = deepcopy(p)

    for name in propertynames(uncertainties)
        if getproperty(uncertainties, name).worst_case == :high
            getproperty(new_p, name) = getproperty(uncertainties, name).low
        elseif getproperty(uncertainties, name).worst_case == :low
            getproperty(new_p, name) = getproperty(uncertainties, name).high
        else
            getproperty(new_p, name) = getproperty(uncertainties, name).nominal
        end
    end

    return new_p
end