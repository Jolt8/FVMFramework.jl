using Revise

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

using Unitful

mesh_path = joinpath(@__DIR__, "output.msh")

grid = togrid(mesh_path)

grid.cellsets

n_cells = length(grid.cells)
u_proto = (
    pressure=zeros(n_cells),
    mass_fractions=(methanol=zeros(n_cells), water=zeros(n_cells), carbon_monoxide=zeros(n_cells), hydrogen=zeros(n_cells), carbon_dioxide=zeros(n_cells)),
    temp=zeros(n_cells)
)

config = create_fvm_config(grid, u_proto)

#in the future, we may want to add a method to make all cells that are not under a cell set 
#become part of a default set with no internal physics or variables 

#controller logic
#We're just going to make this monitor one u field for now 
#we're not going to use this struct, it's just for reference
struct PIDController
    proportional_gain::Float64
    integral_time::Float64
    derivative_time::Float64
    desired_value_comp_vector::NamedTuple
    initial_volumetric_input::Float64
    min_volumetric_input::Float64
    max_volumetric_input::Float64
end

temp_controller = (
    proportional_gain=1.0,
    integral_time=1.0,
    derivative_time=1.0,
    desired_value_comp_vector=(temp = ustrip(270.0u"°C" |> u"K")),
    initial_volumetric_input=1e7,
    min_volumetric_input=0.0,
    max_volumetric_input=1e8
)

#this will probably have to be converted into an integrator internally to get previous temperature and to prevent excessive cache usage

#TODO: Another thing I'd like to implement in my eventual optimizaiton pipeline is a method to extract a very basic
#correlation between the average_temperature measured in the reforming_area cellset to another temp_sensor cellset 
#that could be fed into an arduino to extrapolate sensor data to the reactor's actual internal reforming temp
add_controller!(
    config, "temp_controller";
    controller=temp_controller,
    monitored_cellset="reforming_area",
    affected_cellset="heating_areas",
    controller_function=
    function pid_controller(du, u, monitored_cells, affected_cells, controller_id, controller, cell_volumes)
        field = first(pairs(controller.desired_value_comp_vector))[1]
        desired = controller.desired_value_comp_vector[field] #desired will be a single scalar value
        measured_vec = u[field]
        measured_du_vec = du[field]

        measured_avg = 0.0
        measured_du_avg = 0.0

        @batch for monitored_cell_id in monitored_cells
            measured_avg += measured_vec[monitored_cell_id]
            measured_du_avg += measured_du_vec[monitored_cell_id]
        end

        measured_avg /= length(monitored_cells)
        measured_du_avg /= length(monitored_cells)

        error = measured_avg - desired

        du.integral_error[controller_id] = error

        corrected_volumetric_addition = (
            controller.initial_volumetric_input +
            (controller.proportional_gain * error) +
            (controller.integral_time * u.integral_error[controller_id]) +
            (controller.derivative_time * measured_du_avg)
        )

        corrected_volumetric_addition = clamp(corrected_volumetric_addition, controller.min_volumetric_input, controller.max_volumetric_input)

        @batch for affected_cell_id in affected_cells
            du[field][affected_cell_id] += corrected_volumetric_addition * cell_volumes[affected_cell_id]
        end
    end
)


#for CH3OH, HCOO, OH
van_t_hoff_A = (CH3O=1.7e-6, HCOO=4.74e-13, OH=3.32e-14)
van_t_hoff_dH = (CH3O=-46800.0, HCOO=-115000.0, OH=-110000.0)

MSR_rxn = (
    heat_of_reaction=49500.0, # [J/mol]
    ref_delta_G=-3800.0, # [J/mol]
    ref_temp=298.15, # [K]
    kf_A=1.25e7, # [s^-1]
    kf_Ea=103000.0, # [J/mol]
    reactant_stoich_coeffs=(methanol=1, water=1), # reactant_stoich_coeffs
    product_stoich_coeffs=(carbon_dioxide=1, hydrogen=3),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs=(methanol=-1, water=-1, carbon_monoxide=0, hydrogen=3, carbon_dioxide=1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A=van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH=van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

MD_rxn = (
    heat_of_reaction=90200.0, # [J/mol]
    ref_delta_G=24800.0, # [J/mol]
    ref_temp=298.15, # [K]
    kf_A=1.15e11, # [s^-1]
    kf_Ea=170000.0, # [J/mol]
    reactant_stoich_coeffs=(methanol=1,), # reactant_stoich_coeffs
    product_stoich_coeffs=(carbon_monoxide=1, hydrogen=2), # product_stoich_coeffs: 1 CO + 2 H2
    stoich_coeffs=(methanol=-1, water=0, carbon_monoxide=1, hydrogen=2, carbon_dioxide=0), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A=van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH=van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

WGS_rxn = (
    heat_of_reaction=-41100.0, # [J/mol]
    ref_delta_G=-28600.0, # [J/mol]
    ref_temp=298.15, # [K]
    kf_A=3.65e7, # [s^-1]
    kf_Ea=87500.0, # [J/mol]
    reactant_stoich_coeffs=(carbon_monoxide=1, water=1), # reactant_stoich_coeffs
    product_stoich_coeffs=(carbon_dioxide=1, hydrogen=1), # product_stoich_coeffs: 1 CO2 + 1 H2
    stoich_coeffs=(methanol=0, water=-1, carbon_monoxide=-1, hydrogen=1, carbon_dioxide=1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A=van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH=van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

#since each mass fraction is modified, they have to be vectors
initial_mass_fractions = (
    methanol=[1.0],
    water=[1.3],
    carbon_monoxide=[0.0001],
    hydrogen=[0.02],
    carbon_dioxide=[0.0001]
)

total_mass_fractions = 0.0

for (species_name, mass_fraction) in pairs(initial_mass_fractions)
    total_mass_fractions += initial_mass_fractions[species_name][1]
end

total_mass_fractions

#each mass fraction must be a vector to be mutable
initial_mass_fractions = (
    methanol = [1.0],
    water = [1.3],
    carbon_monoxide = [0.0001],
    hydrogen = [0.02],
    carbon_dioxide = [0.0001]
)

for (species_name, mass_fraction) in pairs(initial_mass_fractions)
    initial_mass_fractions[species_name][1] /= total_mass_fractions
end

initial_mass_fractions = NamedTuple{keys(initial_mass_fractions)}(first.(values(initial_mass_fractions)))
#this is not ideal, but we can't use scalars inside the initial mass_fractions because the values of NamedTuples can't be modified 

reforming_area_properties = (
    k = 237.0, # k (W/(m*K))
    cp = 4.184, # cp (J/(kg*K))
    mu = 1e-5, # mu (Pa*s)
    permeability = 0.6e-11, # permeability (m^2)
    species_ids = (methanol = 1, water = 2, carbon_monoxide = 3, hydrogen = 4, carbon_dioxide = 5), #we could use mass_fractions for species loops, but this is just more consistent
    diffusion_coefficients = (
        methanol = 1e-5,
        water = 1e-5,
        carbon_monoxide = 1e-5,
        hydrogen = 1e-5,
        carbon_dioxide = 1e-5
    ), #diffusion coefficients (m^2/s)
    molecular_weights = (
        methanol = 0.03204,
        water = 0.01802,
        carbon_monoxide = 0.02801,
        hydrogen = 0.00202,
        carbon_dioxide = 0.04401
    ), #species_molecular_weights [kg/mol]
    reactions = (reforming_reactions = (MSR_rxn = MSR_rxn, MD_rxn = MD_rxn, WGS_rxn = WGS_rxn),),
    reactions_kg_cat = (reforming_reactions = (MSR_rxn = 1250.0, MD_rxn = 1250.0, WGS_rxn = 1250.0),), # cell_kg_cat_per_m3_for_each_reaction
)
#remember, for any named tuple with a single field inside it, remember to add a comma to the end 
#ex. (reforming_reactions = (MSR_rxn = MSR_rxn, ),)
#not: (reforming_reaction = (MSR_rxn = MSR_rxn))

#these are just for classifying regions to make sure they do the right connection functions
struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

add_region!(
    config, "reforming_area";
    type = Fluid(),
    initial_conditions = (
        pressure = 100000.0,
        mass_fractions = initial_mass_fractions,
        temp = ustrip(270.0u"°C" |> u"K")
    ),
    properties = reforming_area_properties,
    state_syms = [:pressure, :mass_fractions, :temp],
    cache_syms = [:mass_face, :mass, :heat, :mw_avg, :rho, :molar_concentrations, :net_rates],
    #= 
        Oh no, we're probably going to have to do this eventually at the start of the loop
        However, we have to decide if we should put the cache updates in the region_function or not
        I feel like we should because it doesn't feel consistent for non-flux physics functions that still require cache updates
        However, we'll just not do that for now
    region_cache_update = 
    function reforming_area_cache_update!(du, u, cell_id, vol)
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
    end,
    =#
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        #sources

        #boundary conditions

        #variable summations
        sum_mass_flux_face_to_cell!(du, u, cell_id)

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    end
)

inlet_total_volume = get_cell_set_total_volume(grid, "inlet", config.geo)
#this should be 1e-8 [m^3], which it is

desired_m_dot = ustrip(0.2u"g/s" |> u"kg/s")
corrected_m_dot_per_volume = desired_m_dot / inlet_total_volume

#we should probably create methods to automatically copy the parameters from other already defined regions

add_region!(
    config, "inlet",
    type = Fluid(),
    initial_conditions = (
        pressure = 110000.0,
        mass_fractions = initial_mass_fractions,
        temp = ustrip(270.0u"°C" |> u"K")
    ),
    properties = reforming_area_properties,
    state_syms = [:pressure, :mass_fractions, :temp],
    cache_syms = [:mass_face, :mass, :heat, :mw_avg, :rho, :molar_concentrations, :net_rates],
    #anything else is assumed to be a fixed variable
    region_function =
    function inlet_area!(du, u, cell_id, vol)
        #property retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #sources
        du.mass[cell_id] += corrected_m_dot_per_volume * vol

        #boundary conditions
        for (species_name, value) in pairs(u.species_ids)
            #not sure which one of these is better
            #u.mass_fractions[species_name][cell_id] = initial_mass_fractions[species_name]
            du.mass_fractions[species_name][cell_id] = 0.0
        end

        #variable summations
        sum_mass_flux_face_to_cell!(du, u, cell_id)

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "outlet",
    type = Fluid(),
    initial_conditions = (
        pressure = 100000.0,
        mass_fractions = initial_mass_fractions,
        temp = ustrip(270.0u"°C" |> u"K")
    ),
    properties = reforming_area_properties,
    state_syms = [:pressure, :mass_fractions, :temp],
    cache_syms = [:mass_face, :mass, :heat, :mw_avg, :rho, :molar_concentrations, :net_rates],
    region_function = 
    function outlet_area!(du, u, cell_id, vol)
        #property retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics

        #sources

        #boundary conditions
        du.mass[cell_id] -= corrected_m_dot_per_volume * vol

        #variable summations
        sum_mass_flux_face_to_cell!(du, u, cell_id)

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "wall";
    type = Solid(),
    initial_conditions = (
        temp = ustrip(270.0u"°C" |> u"K"),
    ),
    properties = (
        k = 237.0, # k (W/(m*K))
        rho = 2700.0, # rho (kg/m^3)
        cp = 921.0, # cp (J/(kg*K))
    ),
    state_syms = [:temp],
    cache_syms = [:heat],
    region_function = 
    function wall_area!(du, u, cell_id, vol)
        #property retrieval

        #internal physics

        #sources

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

#we might want to add something akin to distribute_over_set_volume!(du += input_wattage)
total_heating_volume = get_cell_set_total_volume(grid, "heating_areas", config.geo)

input_wattage = 100.0 # W
corrected_volumetric_heating = input_wattage / total_heating_volume

add_region!(
    config, "heating_areas";
    type = Solid(),
    initial_conditions = (
        temp = ustrip(270.0u"°C" |> u"K"),
    ),
    properties = (
        k = 237.0, # k (W/(m*K))
        rho = 2700.0, # rho (kg/m^3)
        cp = 921.0, # cp (J/(kg*K))
    ),
    state_syms = [:temp],
    cache_syms = [:heat],
    region_function=
    function heating_area!(du, u, cell_id, vol)
        #property retrieval

        #internal physics

        #sources
        du.heat[cell_id] += corrected_volumetric_heating * vol

        #boundary conditions

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    mw_avg!(u, idx_a)
    rho_ideal!(u, idx_a)
    molar_concentrations!(u, idx_a)

    mw_avg!(u, idx_b)
    rho_ideal!(u, idx_b)
    molar_concentrations!(u, idx_b)
    
    continuity_and_momentum_darcy(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )

    #=diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )

    all_species_advection!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )

    enthalpy_advection!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#

    mass_fraction_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
end

function solid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)

    #=diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )=#
end

function fluid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    mw_avg!(u, idx_a)
    rho_ideal!(u, idx_a)
    molar_concentrations!(u, idx_a)

    #=diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )=#
end

function solid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    mw_avg!(u, idx_b)
    rho_ideal!(u, idx_b)
    molar_concentrations!(u, idx_b)

    #=diffusion_temp_exchange!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    =#
end

#this is the smallest I could make this function
#I like using <: here because it makes it look nice with syntax highlighting
function connection_map_function(type_a, type_b)
    typeof(type_a) <: Fluid && typeof(type_b) <: Fluid && return fluid_fluid_flux!
    typeof(type_a) <: Solid && typeof(type_b) <: Solid && return solid_solid_flux!
    typeof(type_a) <: Fluid && typeof(type_b) <: Solid && return fluid_solid_flux!
    typeof(type_a) <: Solid && typeof(type_b) <: Fluid && return solid_fluid_flux!
end

n_faces = length(config.geo.cell_neighbor_areas[1])
n_cells = length(config.geo.cell_volumes)
n_reactions = length(config.regions[1].properties.reactions.reforming_reactions)
reaction_names = keys(config.regions[1].properties.reactions.reforming_reactions)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))), 
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
)

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, special_caches)

f_closure_implicit = (du, u, p) -> methanol_reformer_f_steady_state!(
    du, u, p, 

    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,

    system.connection_groups, system.controller_groups, system.region_groups, 
    system.merged_properties,

    system.du_diff_cache_vec, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)
#just remove t from the above closure function and from methanol_reformer_f_test! itself to NonlinearSolve this system

#test = get_tmp(system.u_diff_cache_vec, 1.0)

#named_test = create_views_inline(test, system.u_cache_axes)

#named_test.mass_face[26000]

#named_test.rho

#filter(x -> x == 2700.0, named_test.rho)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

p_guess = 0.0

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess), du0_vec, u0_vec, detector
)

nl_func = NonlinearFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

prob = NonlinearProblem(nl_func, u0_vec, p_guess)

VSCodeServer.@profview sol = solve(prob, NonlinearSolve.NewtonRaphson(concrete_jac = true))

sol.u

record_sol = true

sim_file = @__FILE__

u_proto_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

#u_named = rebuild_u_named(sol.u, u_proto_named)

cell = grid.cellsets["inlet"][1]

inlet_pressure = []

for step in eachindex(sol.u)
    push!(inlet_pressure, u_proto_named[step].pressure[cell])
end

inlet_pressure

sol.u[1] == sol.u[end]

mass_fractions_beginning = u_proto_named[1].mass_fractions.methanol[6389]

mass_fractions_end = u_proto_named[end].mass_fractions.methanol[6389]

mass_fractions_beginning == mass_fractions_end
mass_fractions_beginning - mass_fractions_end

length(u_proto_named[1])

if record_sol == true
    sol_to_vtk(sol, u_proto_named, grid, sim_file)
end