using Revise
using Logging
using Unitful
using OrdinaryDiffEq
using NonlinearSolve
using Sparspak
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using ILUZero
using LinearAlgebra
using StaticArrays
using FVMFramework
using Clapeyron

mesh_path = joinpath(@__DIR__, "meshes/hybrid_rocket_cfd_3d.msh")
grid = togrid(mesh_path)

n_cells = length(grid.cells)
n_faces = nfacets(grid.cells[1])

compressible_solver_options = (
    use_sst_turbulence = true,
)

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3",
    momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_v = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_w = zeros(n_cells)u"kg/(m^2*s)",
    volumetric_energy = zeros(n_cells)u"J/m^3",
    fuel_grain_thickness = zeros(n_cells)u"m"
)

if compressible_solver_options.use_sst_turbulence == true
    u_proto = merge_properties(u_proto, ComponentVector(
        turbulent_kinetic_energy_density = zeros(n_cells)u"J/m^3",
        specific_dissipation_rate_density = zeros(n_cells)u"kg/(m^3*s)",
    ))
end

config = create_fvm_config(grid, u_proto)

Revise.includet(joinpath(@__DIR__, "geometry_editing/geometry_editing.jl"))

config.geo = correct_hybrid_geo(config.geo, grid)

config.geo.cell_volumes
config.geo.cell_centroids
config.geo.cell_face_normals[1]


