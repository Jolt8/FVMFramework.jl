"""
    geometry_building!(geo, fvm_topology, node_coordinates)

Compatibility name for rebuilding an already allocated geometry. New code can
call `update_fvm_geometry!` directly.
"""
function geometry_building!(geo, fvm_topology, node_coordinates)
    return update_fvm_geometry!(geo, fvm_topology, node_coordinates)
end
