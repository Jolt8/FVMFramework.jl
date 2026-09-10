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
    cell_neighbors::Vector{Tuple{Int, Vector{Tuple{Int, Int}}}}
end

mutable struct SimulationConfigInfo
    grid::Ferrite.Grid
    geo::FVMGeometry
    top::ExclusiveTopology
    regions::Vector{RegionSetupInfo}
    patches::Vector{PatchSetupInfo}

    cache_syms_and_units::NamedTuple
    special_caches::ComponentVector
    second_order_syms::Vector{Symbol}
    optimized_parameters::ComponentVector

    u_proto::ComponentVector
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

    return SimulationConfigInfo(
        grid,
        geo,
        top, 
        RegionSetupInfo[],
        PatchSetupInfo[],
        ControllerSetupInfo[],
        NamedTuple(), 
        ComponentVector(),
        Symbol[],
        ComponentVector(),
        u_proto
    )
end