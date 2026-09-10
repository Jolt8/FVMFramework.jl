#Note that this script was written by AI

"""
    foreach_field_at!(f, cell_id::Int, groups::Vararg{Any, N}) where {N}

Iterates over the elements of each field at a specific cell_id

# Example:
```julia
foreach_field_at!(cell_id, du.mass_fractions, du.molar_concentrations) do species, du_mass_fractions, du_molar_concentrations
    du_mass_fractions[species] += 1.0
    du_molar_concentrations[species] += 1.0
end
```

The reason this function exists is that this causes a ton of allocations and dynamic dispatch:
```julia
for name in propertynames(du.mass_fractions[1])
    getproperty(du.mass_fractions, name)[cell_id] += 1.0 
end
```

However, with how this function works, there are sometimes cases where it does not behave as you would expect.
For example, let's say you have these two vectors and wanted to iterate through both of them:
```julia
mass_fractions = ComponentVector(
    methane = 0.5,
    water = 0.5
)

elemental_compositions = ComponentVector(
    methane = (
        C = 1,
        H = 4,
        O = 0
    ),
    water = (
        C = 0,
        H = 2,
        O = 1
    )
)
```
\n
If you just do:
```julia
foreach_field_at!(cell_id, mass_fractions, elemental_compositions) do species, mass_fractions, elemental_compositions
    mass_fractions[species[cell_id]] += 1.0
    elemental_compositions[species[cell_id]] += 1.0
end
```
**species** will just be [1], [2] etc.
\n
This will not index [elemental_compositions] properly because it needs be indexed as [1:3], [4:6]
thus, **elemental_compositions[species]** will just return the value for **elemental_compositions.methane.C** and nothing else

To make this this doesn't happen, do this:
```julia
mass_fractions_species_idx = 1

foreach_field_at!(cell_id, elemental_compositions) do species, elemental_compositions
    view(mass_fractions, cell_id)[mass_fractions_species_idx] += 1.0
    elemental_compositions[species[cell_id]] += 1.0
    mass_fractions_species_idx += 1
end
```
While this is not the most ideal, it's pretty much the best option avaliable without introducing even more bloated looping functions
"""
@generated function foreach_field_at!(f, cell_id::Int, groups::Vararg{Any, N}) where {N}
    G1 = groups[1]
    local properties
    if G1 <: VirtualFVMArray
        properties = fieldnames(G1.parameters[2])
    elseif G1 <: ComponentArray
        Ax = G1.parameters[4].parameters[1]
        properties = keys(Ax.parameters[1])
    elseif G1 <: SubArray || G1 <: Number
        error("only one field detected in ComponentVector, perhaps you forgot to add a comma after creating a field with only one symbol
            \n
            for example:
            incorrect: 
            ComponentVector(
                mass_fractions = (
                    carbon_dioxide = 1.0
                )
            )
            \n
            correct: 
            ComponentVector(
                mass_fractions = (
                    carbon_dioxide = 1.0,
                )
            )
            "
        )
    else
        properties = ()
    end

    exprs = []
    for n in properties
        view_exprs = []
        for i in 1:N
            if groups[i] <: VirtualFVMArray
                push!(view_exprs, :(_resolve_unified_val(groups[$i], Val{$(QuoteNode(n))}())))
            else
                push!(view_exprs, :(getproperty(groups[$i], $(QuoteNode(n)))))
            end
        end
        push!(exprs, :(f(cell_id, $(view_exprs...))))
    end

    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

"""
    for_fields!(f, groups::Vararg{Any, N}) where {N}

    Iterates over the elements of each field

    # Example:
    ```julia
    for_fields!(du.mass_fractions) do species, du_mass_fractions
        du_mass_fractions[species] += 1.0 
    end
    ```

    Similar to for_fields_at! this function exists because looping using propertynames() causes a significant amount of 
    allocations and dynamic dispatch
"""
@generated function for_fields!(f, groups::Vararg{Any, N}) where {N}
    G1 = groups[1]
    exprs = []
    if G1 <: VirtualFVMArray
        # For VirtualFVMArray, axes are in the NamedTuple 'axes'
        AxType = G1.parameters[2] # NamedTuple type of axes
        for name in fieldnames(AxType)
            # Use getproperty at runtime to get the axis
            push!(exprs, :(f(getfield(getfield(groups[1], :axes), $(QuoteNode(name))), groups...)))
        end
    elseif G1 <: ComponentArray
        # For ComponentArray, the axes are in the 4th type parameter
        # AxTupleType: Tuple{Axis{...}}
        AxType = G1.parameters[4].parameters[1]
        # AxType: Axis{layout}
        layout = AxType.parameters[1] # This is the NamedTuple instance
        for i in eachindex(layout)
            axis_obj = layout[i]
            # Try to get the index/range from the type of the axis object
            # ComponentArrays uses the first type parameter for the index/range in most axes
            T = typeof(axis_obj)
            if hasproperty(T, :parameters) && length(T.parameters) >= 1
                idx = T.parameters[1]
                push!(exprs, :(f($idx, groups...)))
            else
                # Fallback
                push!(exprs, :(f($axis_obj, groups...)))
            end
        end
    elseif G1 <: SubArray || G1 <: Number
        error(
        "Only one field detected in ComponentVector, perhaps you forgot to add a comma after creating a field with only one symbol
        For example:
            Incorrect:
            ComponentVector(
                mass_fractions = (
                    carbon_dioxide = 1.0
                )
            )
            Correct:
            ComponentVector(
                mass_fractions = (
                    carbon_dioxide = 1.0,
                )
            )"
        )
    end

    return quote
        Base.@_inline_meta
        $(exprs...)
        nothing
    end
end

#const foreach_field_at! = for_fields!