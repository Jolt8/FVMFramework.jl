mutable struct RegionSetupInfo{P <: AbstractPhysics} #this must be defined before SimulationConfigInfo
    name::String
    type::P
    initial_conditions::ComponentVector
    properties::ComponentVector
    property_update_function::Function
    region_function::Function
    region_cells::Vector{Int}
end

mutable struct PatchSetupInfo #this must be defined before SimulationConfigInfo
    name::String
    properties::ComponentVector
    patch_function::Function
    cell_neighbors::Vector{Tuple{Int, Vector{Tuple{Int, Int, Int}}}}
    owner_region::String #"" = every cell of the face set; kept so update_patch! can rebuild the faces
end

mutable struct SimulationConfigInfo
    # Everything the setup needs (cells, faces, cell sets and face sets) is in
    # geo, so the same config works for a Ferrite mesh and a reduced-order
    # model (create_fvm_config_from_rom). Keep your own `grid` for VTK output.
    geo::FVMGeometry
    regions::Vector{RegionSetupInfo}
    patches::Vector{PatchSetupInfo}

    cache_syms_and_units::NamedTuple
    special_caches::ComponentVector
    second_order_syms::Vector{Symbol}
    optimized_parameters::ComponentVector

    u_proto::ComponentVector

    data::NamedTuple
end


"""
    create_fvm_config(grid, u_proto)

Create a `SimulationConfigInfo` struct with the given grid and 
initial conditions.

# Arguments
- `grid::Ferrite.Grid`: The grid to use for the simulation.
- `u_proto::ComponentVector`: The initial conditions for the simulation.

# Returns
- `SimulationConfigInfo`: The created simulation config struct.

"""
function create_fvm_config(grid, u_proto)
    top = ExclusiveTopology(grid)

    geo = build_fvm_geo_into_struct(grid, top)

    return create_fvm_config_from_geo(geo, u_proto)
end

"""
    create_fvm_config_from_geo(geo, u_proto)

Create a `SimulationConfigInfo` from an already built `FVMGeometry` (from a
mesh via `build_fvm_geo_into_struct`, or from a ROM export).
"""
function create_fvm_config_from_geo(geo, u_proto)
    return SimulationConfigInfo(
        geo,
        RegionSetupInfo[],
        PatchSetupInfo[],
        NamedTuple(), 
        ComponentVector(),
        Symbol[],
        ComponentVector(),
        u_proto,
        ()
    )
end