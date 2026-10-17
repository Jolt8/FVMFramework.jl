using Revise
using Logging
using Unitful
using OrdinaryDiffEq
using NonlinearSolve
using Sparspak
using Ferrite
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using ILUZero
using LinearAlgebra
using StaticArrays
using FVMFramework
using Clapeyron

mesh_path = joinpath(@__DIR__, "airfoil_mesh/ramjet/ramjet_mesh-Cut.msh")
grid = togrid(mesh_path)

n_cells = length(grid.cells)
n_faces = nfacets(grid.cells[1])