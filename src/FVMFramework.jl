module FVMFramework

using Ferrite
using OrdinaryDiffEq
using LinearAlgebra
using SparseArrays
#using SciMLSensitivity
#using Optimization
#using OptimizationPolyalgorithms
#using Zygote
#using Enzyme
using ForwardDiff
using RecursiveArrayTools
#using OptimizationOptimJL
using ILUZero
using NonlinearSolve
using ComponentArrays
using StaticArrays
using ProfileView
using Dates
using WriteVTK
using PreallocationTools
using Polyester
using Unitful

import AlgebraicMultigrid
import SparseConnectivityTracer
import ADTypes
import Logging

#for future me, if you use this project on a new computer, in the FVMFramework environment just do:
#] activate .
#] instantiate

# ----- Geometry -----
include("geometry/geometry_helper_functions.jl")
export get_node_coordinates, get_cell_topology, get_nodes_of_cells, get_face_nodes
export get_cell_neighbors, get_unconnected_map, get_cell_face_map, cross_product

include("geometry/geometry_rebuilding_tetra.jl")
export rebuild_fvm_geometry_tetra!

include("geometry/geometry_rebuilding_hexa.jl")
export rebuild_fvm_geometry_hexa!

include("geometry/geometry_building.jl")
export build_fvm_geo_into_struct, FVMGeometry, FVMGeometryTetra, FVMGeometryHexa

include("geometry/topology_checking.jl")
export check_cellset_connectivity, check_grid_connectivity

# ----- ComponentArray Addons -----
include("component_arrays_addons/simple_merge.jl")
export merge_properties

include("component_arrays_addons/virtual_component_array.jl")
export FaceVectorView, VirtualAxis, VirtualFVMArray, virtual_merge_axes

include("component_arrays_addons/for_fields!_componentarray_looping.jl")
export for_fields!, foreach_field_at!

# ----- Physics ----
#   ---- Physics Types ----
include("physics_types/abstract_physics_types.jl")
export AbstractPhysics

#   ---- Flux Methods ----
#       --- Flux Physics ---
include("physics/flux_physics/advection.jl")
export species_advection!, all_species_advection!, enthalpy_advection!

include("physics/flux_physics/darcy_flow.jl")
export get_darcy_mass_flux, pressure_driven_mass_flux!

include("physics/flux_physics/diffusion.jl")
export species_numerical_flux, mass_fraction_diffusion!

include("physics/flux_physics/heat_transfer.jl")
export get_k_effective, numerical_flux, heat_diffusion!

#   ---- Internal Methods ----
#       --- Internal Capacities ---
include("physics/capacity_helper_functions.jl")
export cap_heat_flux_to_temp_change!, cap_mass_flux_to_pressure_change!, cap_species_mass_flux_to_mass_fraction_change!

#   ---- Helper Functions ----
include("physics/physics_helper_functions.jl")
export upwind, harmonic_mean #for fluxes
export van_t_hoff, arrenhius_k, K_gibbs_free #for chemical reactions
export mw_avg!, rho_ideal!, molar_concentrations!, get_cell_cp #other misc props

#   ---- Summation Functions ----
include("physics/physics_face_summation_functions.jl")
export sum_mass_flux_face_to_cell!

# ----- Setup Methods -----
# ---- Helper Functions ----
include("setup/setup_helper_functions/setup_geometry_helper_functions.jl")
export get_cellset_volume, get_facetset_area, get_facetset_area_per_cell


# ---- Main Setup Methods ----
include("setup/sim_config_base.jl")
export create_fvm_config, RegionSetupInfo, PatchSetupInfo, ControllerSetupInfo, SimulationConfigInfo

include("setup/sim_config_ca_merge_handling.jl")
export merge_region_properties, merge_region_caches

include("setup/sim_config_setup_syms.jl")
export add_setup_syms!

include("setup/sim_config_second_order_setup.jl")
export merge_in_second_order_caches!

include("setup/sim_config_regions.jl")
export add_region!, update_region!

include("setup/sim_config_patches.jl")
export add_patch!, update_patch!

include("setup/sim_config_units.jl")
export run_and_check_units

include("setup/sim_config_finishing.jl")
export finish_fvm_config, FVMSystem

#   ---- Sim Recording ----
include("recording/sim_recording.jl")
export sol_to_vtk

include("recording/state_regenerator.jl")
export regenerate_fvm_state

# ----- Solvers -----
#   ---- Preconditioners ----
include("solvers/preconditioners.jl")
export iluzero, algebraicmultigrid

#   ---- Callbacks ----
#       --- Progress Callbacks
include("solvers/callbacks/progress_callbacks.jl")
export show_t_progress, approximate_time_to_finish_cb

#   ---- Common Operator Methods ----
include("solvers/common_operator_methods/unpack_fvm_state.jl")
export unpack_fvm_state

include("solvers/common_operator_methods/shared_group_functions.jl")
export solve_connection_group!, update_region_group!, solve_region_group!, solve_patch_group!

include("solvers/common_operator_methods/further_simplified_group_functions.jl")
export solve_connection_groups!, update_region_groups!, solve_region_groups!, solve_patch_groups!

#   ---- FVM Operators ----
include("solvers/fvm_operators/fvm_operator.jl")
export fvm_operator!
end
