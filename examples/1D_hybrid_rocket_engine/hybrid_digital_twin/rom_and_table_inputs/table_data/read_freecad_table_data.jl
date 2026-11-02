#=
Read freecad_table_data.json (written by export_freecad_table_data.py) into a
NamedTuple of Unitful quantities, in the table's own units (mm, mm^2, °, ...):

    table_data = read_freecad_table_data_json(joinpath(@__DIR__, "freecad_table_data.json"))
    table_data.phenolic_liner_inner_diameter   # 33.32 mm
    table_data.converging_half_angle           # 30.0°
    table_data.flange_n_holes                  # 8 (plain number, no unit)
    uconvert(u"m", table_data.phenolic_liner_inner_diameter)  # 0.03332 m
=#

using JSON
using Unitful

function read_freecad_table_data_json(path)
    table_export = JSON.parsefile(path)
    table_parameters = table_export["parameters"]

    parameter_names = Symbol[]
    parameter_values = Any[]

    for (alias, entry) in table_parameters
        entry_value = entry["value"]
        entry_unit = entry["unit"]

        if entry_unit == ""
            # Counts (flange_n_holes) stay integers, ratios stay floats.
            parameter_value = entry_value
        else
            parameter_value = entry_value * Unitful.uparse(entry_unit)
        end

        push!(parameter_names, Symbol(alias))
        push!(parameter_values, parameter_value)
    end

    return NamedTuple{Tuple(parameter_names)}(Tuple(parameter_values))
end
