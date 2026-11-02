
"""
Typical usage:
vol = cell_geometry(geo, cell_id)
"""
function cell_geometry(geo, cell_id)
    return geo.cell_volumes[cell_id]
end