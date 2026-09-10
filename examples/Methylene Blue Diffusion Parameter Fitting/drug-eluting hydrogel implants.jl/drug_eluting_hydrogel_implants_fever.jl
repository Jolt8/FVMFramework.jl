using Revise

using FVMFramework

using Ferrite
using FerriteGmsh
using OrdinaryDiffEq
using SparseArrays
using ComponentArrays   
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero
using StaticArrays
using PreallocationTools
using ForwardDiff
using Polyester
using DataInterpolations
using DataFrames

using XLSX

using Unitful

mesh_path = joinpath(@__DIR__, "DEHI_meshes/dehi_model_output_fine.msh")

grid = togrid(mesh_path)

grid.cellsets
grid.facetsets

n_cells = length(grid.cells)
u_proto = (
    mass_fractions = (methylene_blue = zeros(n_cells), water = zeros(n_cells)),
)

config = create_fvm_config(grid, u_proto);

function normalize_mass_fractions(mass_fractions)
    total_mass_fractions = 0.0

    for (species_name, mass_fraction) in pairs(mass_fractions)
        total_mass_fractions += mass_fractions[species_name][1]
    end

    for species_name in keys(mass_fractions)
        mass_fractions[species_name][1] /= total_mass_fractions
    end

    return NamedTuple{keys(mass_fractions)}(first.(values(mass_fractions)))
end

#these are just for classifying regions to make sure they do the right connection functions
struct ImplantInterior <: AbstractPhysics end
struct SurroundingTissue <: AbstractPhysics end

dialysis_tubing_initial_mass_fractions = (
    methylene_blue = [0.1],
    water = [0.9]
)

dialysis_tubing_initial_mass_fractions = normalize_mass_fractions(dialysis_tubing_initial_mass_fractions)

add_region!(
    config, "implant_interior";
    type = ImplantInterior(),
    initial_conditions = (
        mass_fractions = dialysis_tubing_initial_mass_fractions,
    ),
    properties = (
        temp = ustrip(40.13u"°C" |> u"K"),
        k = 0.6, # k (W/(m*K))
        cp = 4184, # cp (J/(kg*K))
        rho = 1000, # rho (kg/m^3)
        mu = 1e-3, # mu (Pa*s)
        species_ids = (methylene_blue = 1, water = 2), #we could use mass_fractions for species loops, but this is just more consistent
        diffusion_coefficients = (
            methylene_blue = 1.0e-8,
            water = 1.0e-8
        ), #diffusion coefficients (m^2/s)
        molecular_weights = (
            methylene_blue = 0.31985, 
            water = 0.01802
        ), #species_molecular_weights [kg/mol]
    ), 
    optimized_syms = [],
    cache_syms = [:heat, :molar_concentrations, :mass, :species_masses, :mw_avg, :rho], 
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        #cap_heat_flux_to_temp_change!(du, u, cell_id, vol)

        cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    end
)

surrounding_fluid_initial_mass_fractions = (
    methylene_blue = [0.0],
    water = [1.0]
)

surrounding_fluid_initial_mass_fractions = normalize_mass_fractions(surrounding_fluid_initial_mass_fractions)

surrounding_fluid_properties = (
    temp = ustrip(40.13u"°C" |> u"K"),
    k = 0.6, # k (W/(m*K))
    cp = 4184, # cp (J/(kg*K))
    rho = 1000, # rho (kg/m^3)
    mu = 1e-3, # mu (Pa*s)
    species_ids = (methylene_blue = 1, water = 2), #we could use mass_fractions for species loops, but this is just more consistent
    diffusion_coefficients = (
            methylene_blue = 1.0e-8,
            water = 1.0e-8
        ), #diffusion coefficients (m^2/s)
    molecular_weights = (
        methylene_blue = 0.31985, 
        water = 0.01802
    ), #species_molecular_weights [kg/mol]
)

add_region!(
    config, "surrounding_tissue";
    type = SurroundingTissue(),
    initial_conditions = (
        mass_fractions = surrounding_fluid_initial_mass_fractions,
    ),
    properties = surrounding_fluid_properties,
    optimized_syms = [],
    cache_syms = [:heat, :molar_concentrations, :mass, :species_masses, :mw_avg, :rho], 
    region_function =
    function surrounding_fluid!(du, u, cell_id, vol)
        #property updating/retrieval

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations

        #capacities
        cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    end
)

include("arrenhius_mass_fraction_diffusion_meth_blue_and_water.jl") 

add_patch!(
    config, "implant_surface";
    properties = (
        diffusion_pre_exponential_factor = 0.000139038533,
        diffusion_activation_energy = 28718.5536309222
    ), 
    optimized_syms = [],
    patch_function =
    function membrane_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        
        #zeroing out regular diffusion that would happen across this face
        #du.mass_fractions.methylene_blue[idx_a] *= 0.0
        #du.mass_fractions.water[idx_a] *= 0.0
        #du.mass_fractions.methylene_blue[idx_b] *= 0.0
        #du.mass_fractions.water[idx_b] *= 0.0

        arrenhius_mass_fraction_diffusion_meth_blue_and_water!(
            du, u,
            idx_a, idx_b, face_idx,
            cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
            cell_volumes[idx_a], cell_volumes[idx_b]
        )
    end
)

#Connection functions
function implant_implant_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    mass_fraction_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end

function implant_surrounding_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    #nothing happens because the add_patch! is already taking care of this
end

function surrounding_surrounding_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    mass_fraction_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end


#this is the smallest I could make this function
#I like using <: here because it makes it look nice with syntax highlighting
function connection_map_function(type_a, type_b)
    typeof(type_a) <: ImplantInterior && typeof(type_b) <: ImplantInterior && return implant_implant_flux!
    typeof(type_a) <: SurroundingTissue && typeof(type_b) <: ImplantInterior && return implant_surrounding_flux!
    typeof(type_a) <: ImplantInterior && typeof(type_b) <: SurroundingTissue && return implant_surrounding_flux!
    typeof(type_a) <: SurroundingTissue && typeof(type_b) <: SurroundingTissue && return surrounding_surrounding_flux!
end

n_faces = length(config.geo.cell_neighbor_areas[1])
n_cells = length(config.geo.cell_volumes)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
    species_masses = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names)))
)

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, special_caches);

u_test = (; create_views_inline(u0_vec, system.u_proto_axes)..., create_views_inline(get_tmp(system.u_diff_cache_vec, 0.0), system.u_cache_axes)...
)

du_test = (; create_views_inline(du0_vec, system.du_proto_axes)..., create_views_inline(get_tmp(system.du_diff_cache_vec, 0.0), system.du_cache_axes)...
)

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f_test!(
    du, u, p, t, 

    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,

    system.connection_groups, system.region_groups, system.patch_groups,
    system.merged_properties,

    system.du_diff_cache_vec, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system

p_guess = 0.0
#=
prob = ODEProblem(f_closure_implicit, u0_vec, (0, 0), p_guess)
@time sol = solve(prob, Tsit5(), tspan = (0.0, 1e-12))

t_interval = 3600.0 / 200
@time sol = solve(prob, Tsit5(), tspan = (0.0, 3600.0), saveat = t_interval, callback = approximate_time_to_finish_cb)
=#
detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, system.p_vec, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = ustrip.((1.0u"hr" |> u"s"))#ustrip.((168u"hr" |> u"s"))
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, system.p_vec)

save_interval = tMax / 200

VSCodeServer.@profview sol = solve(
    implicit_prob, 
    FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), 
    callback = approximate_time_to_finish_cb,
    saveat = save_interval
)
#holy moly, KrylovJL_GMRES is so much faster than KLUFactorization for these large problems

u_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

drug_mass_fractions = [u_named[i].mass_fractions.methylene_blue[36103] for i in eachindex(u_named)]

df1 = DataFrame(AA = sol.t, AB = drug_mass_fractions)

output_path = joinpath(@__DIR__, "drug_eluting_hydrogel_implant_fever_natural_diffusion_results.xlsx")

rm(output_path)

XLSX.writetable(output_path, 
    "fever" => df1
)

sim_file = @__FILE__

#sol_to_vtk(sol, u_named, grid, sim_file)
