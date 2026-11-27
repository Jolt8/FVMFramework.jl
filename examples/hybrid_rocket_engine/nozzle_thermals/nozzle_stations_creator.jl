struct NozzleStation
    top_position::typeof(1.0u"m")
    middle_position::typeof(1.0u"m")
    bottom_position::typeof(1.0u"m")
    top_diameter::typeof(1.0u"m")
    middle_diameter::typeof(1.0u"m")
    bottom_diameter::typeof(1.0u"m")
    slant_length::typeof(1.0u"m")
    half_angle::typeof(1.0u"°")
    wall_area::typeof(1.0u"m^2")
    is_converging::Bool
    is_diverging::Bool
end

function get_nozzle_throat_position(cad_data)
    return cad_data.nozzle_length - cad_data.nozzle_diverging_length
end

function get_diameter_at_nozzle_station(cad_data, nozzle_station)
    nozzle_throat_position = get_nozzle_throat_position(cad_data)

    if nozzle_station < nozzle_throat_position
        return cad_data.nozzle_throat_diameter + 2 * (nozzle_throat_position - nozzle_station) * tand(cad_data.converging_half_angle)
    else
        return cad_data.nozzle_throat_diameter + 2 * (nozzle_station - nozzle_throat_position) * tand(cad_data.diverging_half_angle)
    end
end

function frustum_lateral_area(d_top, d_bottom, dz)
    r_top = d_top / 2
    r_bottom = d_bottom / 2
    slant_length = sqrt((r_bottom - r_top)^2 + dz^2)
    return pi * (r_top + r_bottom) * slant_length
end

function create_nozzle_stations(cad_data; n_converging = 100, n_diverging = 100)
    nozzle_throat_position = get_nozzle_throat_position(cad_data)

    # Separate grids for each cone so a segment edge lands exactly on the throat
    converging_edges = collect(range(0.0u"mm", nozzle_throat_position, n_converging + 1))
    diverging_edges  = collect(range(nozzle_throat_position, cad_data.nozzle_length, n_diverging + 1))
    segment_edges    = vcat(converging_edges, diverging_edges[2:end]) # drop the duplicated throat edge

    nozzle_stations = NozzleStation[]

    for i in 1:(n_converging + n_diverging)
        top_position = segment_edges[i]
        bottom_position = segment_edges[i+1]
        middle_position = (top_position + bottom_position) / 2

        top_diameter = get_diameter_at_nozzle_station(cad_data, top_position)
        middle_diameter = get_diameter_at_nozzle_station(cad_data, middle_position)
        bottom_diameter = get_diameter_at_nozzle_station(cad_data, bottom_position)

        if middle_position < nozzle_throat_position
            is_converging = true
            is_diverging = false
            half_angle = cad_data.converging_half_angle
        else
            is_converging = false
            is_diverging = true
            half_angle = cad_data.diverging_half_angle
        end

        wall_area = frustum_lateral_area(top_diameter, bottom_diameter, bottom_position - top_position)

        slant_length = sqrt((top_diameter / 2 - bottom_diameter / 2)^2 + (bottom_position - top_position)^2)

        push!(nozzle_stations, NozzleStation(
            top_position,
            middle_position,
            bottom_position,
            top_diameter,
            middle_diameter,
            bottom_diameter,
            slant_length,
            half_angle,
            wall_area,
            is_converging,
            is_diverging,
        ))
    end

    return nozzle_stations
end
