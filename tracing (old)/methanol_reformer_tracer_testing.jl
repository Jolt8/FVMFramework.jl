#=
using Revise
using Logging

ENV["JULIA_PKG_PRESERVE_TIERED_INSTALLED"] == true
ENV["JULIA_PKG_PRESERVE_TIERED_INSTALLED"] = true

using FVMFramework
#I think this take a long time because Enzyme (despite not being installed) throws a warning from DiffEqBaseEnzyme.jl

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
using StatsBase

using Unitful

grid = togrid("C://Users//wille//Desktop//FreeCad Projects//Methanol Reformer//output.msh")

grid.cellsets

n_cells = length(grid.cells)

struct OptimizedVariable
    guess::Float64
end

function optimize!(var)
    return OptimizedVariable(var)
end

#optimize!() will be inserted into any variable in the setup functions below
#for example, if we wanted to fit kf_A and kf_Ea to experimental data we would just do this:
#then, our tracer function would check if typeof(var) == OptimizedVariable, and if it is, we would make a pointer to a p value 

config = create_fvm_config(grid)

#in the future, we may want to add a method to make all cells that are not under a cell set 
#become part of a default set with no internal physics or variables 

#TODO: Another thing I'd like to implement in my eventual optimizaiton pipeline is a method to extract a very basic
#correlation between the average_temperature measured in the reforming_area cellset to another temp_sensor cellset 
#that could be fed into an arduino to extrapolate sensor data to the reactor's actual internal reforming temp

n_faces = length(config.geo.cell_neighbor_areas[1])

add_region!(
    config, "reforming_area";
    type=Fluid(),
    region_function=
    function reforming_area!(du, u, cell_id, cell_volumes)
        #property updating/retrieval
        molar_concentrations!(u, cell_id)
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)

        #variable summation
        sum_mass_flux_face_to_cell!(du, u, cell_id, n_faces)

        #internal physics

        #we have to figure out if we're going to pass in just a singular cell_volume or all cell_volumes and use cell_id
        #power_law_react_cell!(du, u, cell_id, u.reactions.example_reaction, vol)
        #example of how to do a power law reaction

        PAM_reforming_react_cell!(du, u, cell_id, cell_volumes)

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, cell_volumes)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, cell_volumes)
    end
)

inlet_total_volume = get_cell_set_total_volume(grid, "inlet", config.geo)
#this should be 1e-8 [m^3], which it is

desired_m_dot = ustrip(0.2u"g/s" |> u"kg/s")
corrected_m_dot_per_volume = desired_m_dot / inlet_total_volume

add_region!(
    config, "inlet";
    type=Fluid(),
    region_function=
    function inlet_area!(du, u, cell_id, cell_volumes)
        #property retrieval
        molar_concentrations!(u, cell_id)
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)

        #variable summation
        sum_mass_flux_face_to_cell!(du, u, cell_id, n_faces)

        #internal physics
        #=
        power_law_react_cell!(du, u, cell_id, u.reactions.example_reaction, vol)
        PAM_reforming_react_cell!(du, u, cell_id, vol)
        =#

        #sources
        du.mass[cell_id] += corrected_m_dot_per_volume * cell_volumes[cell_id]
        #ok this is a really weird edge case, we keep track of du.mass on a per face_idx basis, but we need to add mass to the cell
        #we should probably differentiate per face_idx variables like du.mass_face[cell_id][face_idx] 
        #then, we'll just sum up all the du.mass_face[cell_id][face_idx] to get du.mass[cell_id]

        #boundary conditions

        # for some reason, this causes path, tracer_ref, and value to be logged 
        for species_name in propertynames(du.mass_fractions)
            du.mass_fractions[species_name][cell_id] *= 0.0
        end
        #du.mass_fractions[:, cell_id] .= 0.0 # we should create a method for this to actually work, for now, the above works fine

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, cell_volumes)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, cell_volumes)
    end
)

add_region!(
    config, "outlet";
    type=Fluid(),
    region_function=
    function outlet_area!(du, u, cell_id, cell_volumes)
        #property retrieval
        molar_concentrations!(u, cell_id)
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)

        #variable summation
        sum_mass_flux_face_to_cell!(du, u, cell_id, n_faces)

        #internal physics

        #sources

        #boundary conditions
        du.mass[cell_id] *= 0.0

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, cell_volumes)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, cell_volumes)
    end
)

add_region!(
    config, "wall";
    type=Solid(),
    region_function=
    function wall_area!(du, u, cell_id, cell_volumes)
        #property retrieval

        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, cell_volumes)
    end
)

#we might want to add something akin to distribute_over_set_volume!(du += input_wattage)
total_heating_volume = get_cell_set_total_volume(grid, "heating_areas", config.geo)

input_wattage = 100.0 # W
corrected_volumetric_heating = input_wattage / total_heating_volume

add_region!(
    config, "heating_areas";
    type=Solid(),
    region_function=
    function heating_area!(du, u, cell_id, cell_volumes)
        #property retrieval

        #internal physics

        #sources
        du.heat[cell_id] += corrected_volumetric_heating * cell_volumes[cell_id]

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, cell_volumes)
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    continuity_and_momentum_darcy!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )

    diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )

    all_species_advection!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )

    enthalpy_advection!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )

    diffusion_mass_fraction_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function solid_solid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    #hmm, perhaps these physics functions need to be more strictly typed
    #Checking profview, I'm getting some runtime dispatch and GC here, I don't know why 
    diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function fluid_solid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function solid_fluid_flux!(
    du, u, idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    diffusion_temp_exchange!(
        du, u, idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end


#this is the smallest I could make this function
#we use types here because we can create subtypes of types to prevent having to do if a == "reforming_area" && b == "inlet" && return fluid_fluid_flux!
function conneciton_map_function(type_a, type_b)
    type_a <: Fluid && type_b <: Fluid && return fluid_fluid_flux!
    type_a <: Solid && type_b <: Solid && return solid_solid_flux!
    type_a <: Fluid && type_b <: Solid && return fluid_solid_flux!
    type_a <: Solid && type_b <: Fluid && return solid_fluid_flux!
end

geo, system = finish_fvm_config(config, conneciton_map_function)

system.connection_groups
#FIXME: inlet, outlet, and heating_areas are not being recorded in connection groups

#this is required for when the tracer encounters something like for mass_fraction in propertynames(mass_fractions)
n_faces = length(geo.cell_neighbor_areas[1])

ctx_loop_context = Dict(
    :mass_fractions => [:methanol, :water, :carbon_monoxide, :hydrogen, :carbon_dioxide],
    :molar_concentrations => [:methanol, :water, :carbon_monoxide, :hydrogen, :carbon_dioxide],
    :reforming_reactions => [:MSR_rxn, :MD_rxn, :WGS_rxn], :face_idxs => [Symbol("face_idx_$i") for i in 1:n_faces] #this is always going to be required 
)

du_tracer, u_tracer, ctx = create_tracer_context(ctx_loop_context)

for conn in system.connection_groups
    idx_a = 1
    neighbor_list = 1
    idx_b = 1
    face_idx = 1

    i = 1

    while neighbor_list == Tuple{Int64,Int64}[]
        idx_a = conn.cell_neighbors[i][1]
        neighbor_list = conn.cell_neighbors[i][2]
        idx_b = neighbor_list[1][1]
        face_idx = neighbor_list[1][2]
        i += 1
    end

    fake_cell_neighbor_areas = Dict(
        conn.name_a => Dict(
            "face_idx_1" => geo.cell_neighbor_areas[idx_a][face_idx],
            "face_idx_2" => geo.cell_neighbor_areas[idx_b][face_idx],
            "face_idx_3" => geo.cell_neighbor_areas[idx_a][face_idx],
            "face_idx_4" => geo.cell_neighbor_areas[idx_b][face_idx],
            "face_idx_5" => geo.cell_neighbor_areas[idx_a][face_idx],
            "face_idx_6" => geo.cell_neighbor_areas[idx_b][face_idx]
        )
    )

    fake_cell_neighbor_normals = Dict(
        conn.name_a => Dict(
            "face_idx_1" => geo.cell_neighbor_normals[idx_a][face_idx],
            "face_idx_2" => geo.cell_neighbor_normals[idx_b][face_idx],
            "face_idx_3" => geo.cell_neighbor_normals[idx_a][face_idx],
            "face_idx_4" => geo.cell_neighbor_normals[idx_b][face_idx],
            "face_idx_5" => geo.cell_neighbor_normals[idx_a][face_idx],
            "face_idx_6" => geo.cell_neighbor_normals[idx_b][face_idx]
        )
    )

    fake_cell_neighbor_distances = Dict(
        conn.name_a => Dict(
            "face_idx_1" => geo.cell_neighbor_distances[idx_a][face_idx],
            "face_idx_2" => geo.cell_neighbor_distances[idx_b][face_idx],
            "face_idx_3" => geo.cell_neighbor_distances[idx_a][face_idx],
            "face_idx_4" => geo.cell_neighbor_distances[idx_b][face_idx],
            "face_idx_5" => geo.cell_neighbor_distances[idx_a][face_idx],
            "face_idx_6" => geo.cell_neighbor_distances[idx_b][face_idx]
        )
    )

    conn.flux_function!(
        du_tracer, u_tracer, conn.name_a, conn.name_b, "face_idx_$face_idx", #we swap out idx_a and idx_b for conn.name_a and conn.name_b for the tracer
        fake_cell_neighbor_areas, fake_cell_neighbor_normals, fake_cell_neighbor_distances
    )
end

#Internal Physics, Sources, Boundary Conditions, and Capacities Loops 
#oh wait, now we don't even need the other fields for the different regions, we only need the region function
for reg in system.region_groups
    fake_cell_volumes = Dict(reg.name => geo.cell_volumes[reg.region_cells[1]])
    reg.region_function!(
        du_tracer, u_tracer, reg.name, #we swap out cell_id for reg.name for the tracer
        fake_cell_volumes
    )
end

var_access_logs, encountered_paths = merge_trace_results(ctx.access_logs)

state_vars, cache_vars, fixed_vars = classify_variables(var_access_logs)

n_cells = length(grid.cells)

region_symbols = Set([Symbol(reg.name) for reg in system.region_groups])

n_faces = length(geo.cell_neighbor_areas[1])

#NOTE: if you ever get an error like no method matching setindex!(::Symbol, ::Symbol, ::Symbol), it's probably because there's conflicting useage of a variable name 

# ── Step 1: Build the per-region setup templates ─────────────────────────
# region_setup gives us a ComponentVector where each top-level key is a region name
# and the values underneath are scalars that the user fills in
state_setup = region_setup(state_vars, region_symbols, n_cells, n_faces)
fixed_setup = region_setup(fixed_vars, region_symbols, n_cells, n_faces)

# ── Step 2: Build the flat merged vectors (per-cell) ─────────────────────
# These are the actual runtime vectors with zeros(n_cells)
complete_state = build_component_array_merge_regions(state_vars, region_symbols, n_cells, n_faces)
complete_fixed = build_component_array_merge_regions(fixed_vars, region_symbols, n_cells, n_faces)

# ── Step 3: Generate the setup script for the user to fill in ────────────
setup_script_path = joinpath(dirname(@__FILE__), "methanol_reformer_setup_4.jl")
generate_setup_script(setup_script_path, state_setup, fixed_setup, system)

#rm(setup_script_path, force=true)

# ── Step 4: After the user fills in the script, include it ───────────────
include(setup_script_path)
# This gives us `initial_state` and `fixed_properties` dicts

# ── Step 5: Populate the merged vectors from the user's values ───────────
populate_merged_vector!(complete_state, initial_state, system)

populate_merged_vector!(complete_fixed, fixed_properties, system)

N::Int = ForwardDiff.pickchunksize(length(propertynames(complete_state)))
cache_cv = build_component_array_merge_regions(cache_vars, region_symbols, n_cells, n_faces)

complete_fixed = ComponentVector(complete_fixed; NamedTuple(cache_cv)...)

DiffCache(complete_fixed, N) #even though this creates unecessary elements on fixed values, it's good enough for now

du_fixed = complete_fixed .= 0.0
u_fixed = complete_fixed

state_axes = getaxes(complete_state)[1]
u0_flat = Vector(complete_state)
du0_flat = u0_flat .= 0.0

du_full = ComponentVector(complete_state; NamedTuple(du_fixed)...)
u_full = ComponentVector(complete_state; NamedTuple(u_fixed)...)

full_axes = getaxes(du_full)[1]

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f_test!(
    du, u, p, t,
    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,
    system.connection_groups, system.region_groups,
    du_fixed, u_fixed,
    full_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system

t0 = 0.0
tMax = 10.0
tspan = (t0, tMax)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

p_guess = 0.0

test_prob = ODEProblem(f_closure_implicit, u0_flat, (0.0, 0.00000001), p_guess)
VSCodeServer.@profview sol = solve(test_prob, Tsit5(), callback=approximate_time_to_finish_cb)

#somehow adding reactions makes this less stiff

#holy moly, this is soooo stiff, even with just 0.000005 s of sim time, it has to take 1700 steps! It also took 356 seconds 

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_flat, u0_flat, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

implicit_prob = ODEProblem(ode_func, u0_flat, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback=approximate_time_to_finish_cb)
VSCodeServer.@profview sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), callback=approximate_time_to_finish_cb)
#algebraicmultigrid is only better for more than 1e6 cells

record_sol = true

sim_file = @__FILE__

u_named = rebuild_u_named_vel(sol.u, u_proto)

grid.cellsets["reforming_area"]

mass_fractions_beginning = u_named[1].mass_fractions[:, 6389]

mass_fractions_end = u_named[end].mass_fractions[:, 6389]

mass_fractions_beginning == mass_fractions_end
mass_fractions_beginning - mass_fractions_end

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end
=#