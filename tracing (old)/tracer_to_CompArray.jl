function inflate_path!(tree::Dict{Symbol, Any}, clean_path::Vector{Symbol}, value)
    curr = tree
    for i in 1:(length(clean_path)-1)
        if !haskey(curr, clean_path[i])
            curr[clean_path[i]] = Dict{Symbol, Any}()
        end
        curr = curr[clean_path[i]]
    end
    curr[clean_path[end]] = value
end

function dict_to_component_vector(d::Dict{Symbol, Any}, n_cells::Int)
    pairs = Pair{Symbol, Any}[]
    for (key, template) in d
        if template isa Dict
            push!(pairs, key => dict_to_component_vector(template, n_cells))
        elseif template == :vector
            push!(pairs, key => zeros(n_cells))
        elseif template == :scalar
            push!(pairs, key => 0.0)
        elseif typeof(template) <: Int
            push!(pairs, key => [zeros(template) for _ in 1:n_cells])
        end
    end
    return ComponentVector(; pairs...)
end

function region_setup(classified_vars::Dict{Vector{Symbol}, VarAccessLog}, region_symbols::Set{Symbol}, n_cells::Int, n_faces::Int)
    tree = Dict{Symbol, Any}()

    face_idxs = [Symbol("face_idx_$i") for i in 1:n_faces]

    for (path, _) in classified_vars
        clean_path = copy(path)
        if any(in.(region_symbols, Ref(path)))
            clean_path = setdiff(path, region_symbols)
            region_symbol = filter(x -> x in region_symbols, path)[1]
            pushfirst!(clean_path, region_symbol)
            template = :scalar
        elseif any(in.(face_idxs, Ref(path)))
            clean_path = setdiff(path, face_idxs)
            template = :scalar
        else
            clean_path = path
            template = :scalar
        end
        inflate_path!(tree, clean_path, template)
    end

    return dict_to_component_vector(tree, n_cells)
end

function build_component_array(classified_vars::Dict{Vector{Symbol}, VarAccessLog}, region_symbols::Set{Symbol}, n_cells::Int, n_faces::Int)
    tree = Dict{Symbol, Any}()

    face_idxs = [Symbol("face_idx_$i") for i in 1:n_faces]

    for (path, _) in classified_vars
        if any(in.(region_symbols, Ref(path)))
            template = :vector
        elseif any(in.(face_idxs, Ref(path)))
            template = n_faces
        else
            template = :scalar
        end
        inflate_path!(tree, path, template)
    end

    return dict_to_component_vector(tree, n_cells)
end


function build_component_array_merge_regions(classified_vars::Dict{Vector{Symbol}, VarAccessLog}, region_symbols::Set{Symbol}, n_cells::Int, n_faces::Int)
    tree = Dict{Symbol, Any}()

    face_idxs = [Symbol("face_idx_$i") for i in 1:n_faces]

    for (path, _) in classified_vars
        clean_path = setdiff(path, region_symbols)
        clean_path = setdiff(clean_path, face_idxs)
        if any(in.(region_symbols, Ref(path)))
            if any(in.(face_idxs, Ref(path)))
                template = n_faces
            else
                template = :vector
            end
        else
            template = :scalar
        end
        inflate_path!(tree, clean_path, template)
    end

    return dict_to_component_vector(tree, n_cells)
end