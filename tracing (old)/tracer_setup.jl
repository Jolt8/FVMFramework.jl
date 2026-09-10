
# ── Script Generation ───────────────────────────────────────────────────────
# Generates a .jl file that the user edits to specify initial conditions 
# and fixed properties per region. The file, when `include`d, returns 
# the filled-in state and fixed ComponentVectors ready to be populated
# into the merged runtime vectors.

"""
    generate_setup_script(filepath, state_setup_cv, fixed_setup_cv, system)

Writes a Julia file at `filepath` that the user fills in to define:
- Initial values for state variables per region
- Values for fixed parameters per region

The generated file defines `get_initial_state()` and `get_fixed_properties()`
which return the filled-in per-region ComponentVectors.
"""
function generate_setup_script(filepath, state_setup_cv, fixed_setup_cv, system)
    io = IOBuffer()

    println(io, "using ComponentArrays")
    println(io, "using Unitful")
    println(io, "")

    # ── State Variables (Initial Conditions) ──
    println(io, "# ── Initial Conditions (State Variables) ─────────────────────────────")
    println(io, "")

    for reg in system.region_groups
        reg_sym = Symbol(reg.name)
        if hasproperty(state_setup_cv, reg_sym)
            region_state = getproperty(state_setup_cv, reg_sym)
            println(io, "# Region: $(reg.name) ($(length(reg.region_cells)) cells)")
            _write_component_vector_template(io, "state_$(reg.name)", region_state, "")
            println(io, "")
        end
    end

    println(io, "")

    # ── Fixed Properties ──
    println(io, "# ── Fixed Properties ─────────────────────────────────────────────────")
    println(io, "")

    for reg in system.region_groups
        reg_sym = Symbol(reg.name)
        if hasproperty(fixed_setup_cv, reg_sym)
            region_fixed = getproperty(fixed_setup_cv, reg_sym)
            println(io, "# Region: $(reg.name) ($(length(reg.region_cells)) cells)")
            _write_component_vector_template(io, "fixed_$(reg.name)", region_fixed, "")
            println(io, "")
        end
    end

    # ── Return functions ──
    println(io, "")
    println(io, "# ── Do not modify below this line ────────────────────────────────────")
    println(io, "")

    # Build the return dict for state
    println(io, "initial_state = Dict{Symbol, Any}(")
    for (i, reg) in enumerate(system.region_groups)
        reg_sym = Symbol(reg.name)
        if hasproperty(state_setup_cv, reg_sym)
            comma = i < length(system.region_groups) ? "," : ""
            println(io, "    :$(reg.name) => state_$(reg.name)$(comma)")
        end
    end
    println(io, ")")
    println(io, "")

    # Build the return dict for fixed
    println(io, "fixed_properties = Dict{Symbol, Any}(")
    for (i, reg) in enumerate(system.region_groups)
        reg_sym = Symbol(reg.name)
        if hasproperty(fixed_setup_cv, reg_sym)
            comma = i < length(system.region_groups) ? "," : ""
            println(io, "    :$(reg.name) => fixed_$(reg.name)$(comma)")
        end
    end
    println(io, ")")

    # Write to file
    open(filepath, "w") do f
        write(f, String(take!(io)))
    end

    println("Setup script written to: $filepath")
    println("Edit the values, then include() it to get `initial_state` and `fixed_properties` dicts.")
end


"""
Helper: recursively writes a ComponentVector as a Julia code template.
Scalars get placeholder 0.0, nested CVs get nested ComponentVector.
"""
function _write_component_vector_template(io, varname, cv, indent)
    pnames = propertynames(cv)
    println(io, "$(indent)$(varname) = ComponentVector(")
    for (i, pname) in enumerate(pnames)
        val = getproperty(cv, pname)
        comma = i < length(pnames) ? "," : ""
        if val isa ComponentVector
            _write_component_vector_template(io, "$(indent)    $(pname)", val, "$(indent)    ")
            # close with comma
            println(io, "$(indent)    )$(comma)")
        else
            println(io, "$(indent)    $(pname) = 0.0$(comma)")
        end
    end
    if indent == ""
        println(io, ")")
    end
end


# ── Population ──────────────────────────────────────────────────────────────
# Takes the user-filled setup dicts and populates the flat merged ComponentVectors

"""
    populate_merged_vector!(merged_cv, setup_dict, system)

Fills `merged_cv` (the flat, per-cell ComponentVector from `build_component_array_merge_regions`)
using the per-region scalar values from `setup_dict` (the dict returned by the setup script).

For each region, iterates over its cells and sets the corresponding values in the merged vector.
"""
function populate_merged_vector!(merged_cv, setup_dict, system)
    for reg in system.region_groups
        reg_sym = Symbol(reg.name)
        if !haskey(setup_dict, reg_sym)
            continue
        end
        region_values = setup_dict[reg_sym]
        _populate_region_cells!(merged_cv, region_values, reg.region_cells)
    end

    return merged_cv
end


"""
Helper: recursively drills into a ComponentVector template and writes 
scalar values into the merged vector at each cell_id.
"""
function _populate_region_cells!(merged_cv, region_values, region_cells)
    for pname in propertynames(region_values)
        val = getproperty(region_values, pname)
        merged_field = getproperty(merged_cv, pname)

        if val isa ComponentVector && merged_field isa ComponentVector
            # Nested — recurse
            _populate_region_cells!(merged_field, val, region_cells)
        elseif val isa Number && merged_field isa AbstractVector
            # Scalar value → fill into all cells of this region
            for cell_id in region_cells
                merged_field[cell_id] = val
            end
        elseif val isa Number && merged_field isa Number
            # Scalar to scalar
            setproperty!(merged_cv, pname, val)
        end
    end
end