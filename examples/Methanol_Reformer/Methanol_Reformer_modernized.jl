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

mesh_path = joinpath(@__DIR__, "Methanol Reformer V1/output.msh")

grid = togrid(mesh_path)

grid.cellsets

n_cells = length(grid.cells)
u_proto = ComponentVector(
    pressure = zeros(n_cells)u"Pa",
    mass_fractions = (
        methanol = zeros(n_cells)u"kg/kg", 
        water = zeros(n_cells)u"kg/kg", 
        carbon_monoxide = zeros(n_cells)u"kg/kg", 
        hydrogen = zeros(n_cells)u"kg/kg", 
        carbon_dioxide = zeros(n_cells)u"kg/kg"
    ),
    temp = zeros(n_cells)u"K"
)

config = create_fvm_config(grid, u_proto)

#for CH3OH, HCOO, OH
van_t_hoff_A = (CH3O = 1.7e-6, HCOO = 4.74e-13, OH = 3.32e-14)
van_t_hoff_dH = (CH3O = -46800.0u"J/mol", HCOO = -115000.0u"J/mol", OH = -110000.0u"J/mol")

MSR_rxn = (
    heat_of_reaction = 49500.0u"J/mol", 
    ref_delta_G = -3800.0u"J/mol", 
    ref_temp = 298.15u"K",
    kf_A = 1.25e7u"s^-1",
    kf_Ea = 103000.0u"J/mol",
    reactant_stoich_coeffs = (methanol = 1, water = 1),
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3), 
    stoich_coeffs = (methanol = -1, water = -1, carbon_monoxide = 0, hydrogen = 3, carbon_dioxide = 1), 
    van_t_hoff_A = van_t_hoff_A,  
    van_t_hoff_dH = van_t_hoff_dH
)

MD_rxn = (
    heat_of_reaction = 90200.0u"J/mol", 
    ref_delta_G = 24800.0u"J/mol", 
    ref_temp = 298.15u"K",
    kf_A = 1.15e11u"s^-1",
    kf_Ea = 170000.0u"J/mol",
    reactant_stoich_coeffs = (methanol = 1,),
    product_stoich_coeffs = (carbon_monoxide = 1, hydrogen = 2),
    stoich_coeffs = (methanol = -1, water = 0, carbon_monoxide = 1, hydrogen = 2, carbon_dioxide = 0), 
    van_t_hoff_A = van_t_hoff_A,
    van_t_hoff_dH = van_t_hoff_dH
)

WGS_rxn = (
    heat_of_reaction = -41100.0u"J/mol", 
    ref_delta_G = -28600.0u"J/mol", 
    ref_temp = 298.15u"K",
    kf_A = 3.65e7u"s^-1",
    kf_Ea = 87500.0u"J/mol",
    reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1),
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1),
    stoich_coeffs = (methanol = 0, water = -1, carbon_monoxide = -1, hydrogen = 1, carbon_dioxide = 1), 
    van_t_hoff_A = van_t_hoff_A,
    van_t_hoff_dH = van_t_hoff_dH
)

#since each mass fraction is modified, they have to be vectors
initial_mass_fractions = ComponentVector(
    methanol = 1.0u"kg/kg",
    water = 1.3u"kg/kg",
    carbon_monoxide = 0.0001u"kg/kg",
    hydrogen = 0.02u"kg/kg",
    carbon_dioxide = 0.0001u"kg/kg"
)

total_mass_fractions = sum(initial_mass_fractions)
initial_mass_fractions ./= total_mass_fractions

reforming_area_properties = ComponentVector(
    k = 237.0u"W/(m*K)", # k (W/(m*K))
    cp = 4.184u"J/(kg*K)", # cp (J/(kg*K))
    mu = 1e-5u"Pa*s", # mu (Pa*s)
    R_gas = 8.314u"J/(mol*K)",
    permeability = 0.6e-11u"m^2",
    species_ids = (methanol = 1, water = 2, carbon_monoxide = 3, hydrogen = 4, carbon_dioxide = 5), #we could use mass_fractions for species loops, but this is just more consistent
    diffusion_coefficients = (
        methanol = 1e-5u"m^2/s",
        water = 1e-5u"m^2/s",
        carbon_monoxide = 1e-5u"m^2/s",
        hydrogen = 1e-5u"m^2/s",
        carbon_dioxide = 1e-5u"m^2/s"
    ), #diffusion coefficients (m^2/s)
    molecular_weights = (
        methanol = 0.03204u"kg/mol",
        water = 0.01802u"kg/mol",
        carbon_monoxide = 0.02801u"kg/mol",
        hydrogen = 0.00202u"kg/mol",
        carbon_dioxide = 0.04401u"kg/mol"
    ), #species_molecular_weights [kg/mol]
    reactions = (reforming_reactions = (MSR_rxn = MSR_rxn, MD_rxn = MD_rxn, WGS_rxn = WGS_rxn),),
    reactions_kg_cat = (reforming_reactions = (MSR_rxn = 1250.0u"kg/m^3", MD_rxn = 1250.0u"kg/m^3", WGS_rxn = 1250.0u"kg/m^3"),), # cell_kg_cat_per_m3_for_each_reaction
)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

common_cache_syms_and_units = (
    mass_face = u"kg",
    mass = u"kg",
    heat = u"J",
    mw_avg = u"kg/mol",
    rho = u"kg/m^3",
    molar_concentrations = u"mol/m^3",
    net_rates = u"mol/(m^3*s)",
    #species_mass_flows = u"kg"
)

add_region!(
    config, "reforming_area";
    type = Fluid(),
    initial_conditions = ComponentVector(
        pressure = 100000.0u"Pa",
        mass_fractions = initial_mass_fractions,
        temp = 270.0u"°C"
    ),
    properties = reforming_area_properties,
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
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
        #cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    end
)

inlet_total_volume = get_cellset_volume(config, "inlet")u"m^3"

desired_m_dot = 0.2u"g/s"
corrected_m_dot_per_volume = desired_m_dot / inlet_total_volume

#we should probably create methods to automatically copy the parameters from other already defined regions
add_region!(
    config, "inlet",
    type = Fluid(),
    initial_conditions = ComponentVector(
        pressure = 110000.0u"Pa",
        mass_fractions = initial_mass_fractions,
        temp = 270.0u"°C"
    ),
    properties = merge_properties(reforming_area_properties, ComponentVector(corrected_m_dot_per_volume = corrected_m_dot_per_volume,)),
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
    #anything else is assumed to be a fixed variable
    region_function =
    function inlet_area!(du, u, cell_id, vol)
        #property retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #sources
        du.mass[cell_id] += u.corrected_m_dot_per_volume[cell_id] * vol

        #boundary conditions
        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
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
    initial_conditions = ComponentVector(
        pressure = 100000.0u"Pa",
        mass_fractions = initial_mass_fractions,
        temp = 270.0u"°C"
    ),
    properties = reforming_area_properties,
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
    region_function = 
    function outlet_area!(du, u, cell_id, vol)
        #property retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics

        #sources

        #boundary conditions
        du.mass[cell_id] -= u.corrected_m_dot_per_volume[cell_id] * vol

        #variable summations
        sum_mass_flux_face_to_cell!(du, u, cell_id)

        #capacities
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
        #cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "wall";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 270.0u"°C",
    ),
    properties = ComponentVector(
        k = 237.0u"W/(m*K)",
        rho = 2700.0u"kg/m^3",
        cp = 921.0u"J/(kg*K)",
    ),
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
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

total_heating_volume = get_cellset_volume(config, "heating_areas")u"m^3"

input_wattage = 100.0u"W"
corrected_volumetric_heating = input_wattage / total_heating_volume

add_region!(
    config, "heating_areas";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 270.0u"°C",
    ),
    properties = ComponentVector(
        k = 237.0u"W/(m*K)",
        rho = 2700.0u"kg/m^3",
        cp = 921.0u"J/(kg*K)",
        heater_volumetric_heating = corrected_volumetric_heating
    ),
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
    region_function=
    function heating_area!(du, u, cell_id, vol)
        #property retrieval

        #internal physics

        #sources
        du.heat[cell_id] += u.heater_volumetric_heating[cell_id] * vol

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
    
    pressure_driven_mass_flux!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )

    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )

    #=all_species_advection!(
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

    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function fluid_solid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    mw_avg!(u, idx_a)
    rho_ideal!(u, idx_a)
    molar_concentrations!(u, idx_a)

    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
end

function solid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances
)
    mw_avg!(u, idx_b)
    rho_ideal!(u, idx_b)
    molar_concentrations!(u, idx_b)

    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx]
    )
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
special_caches = ComponentVector(
    mass_face = zeros(n_cells, n_faces)u"kg",
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(
            Tuple(zeros(n_cells)u"mol/s" for _ in 1:length(reaction_names))
        ), #don't forget the comma!
    ), 
    molar_concentrations = NamedTuple{species_names}(
        Tuple(zeros(n_cells)u"mol/m^3" for _ in 1:length(species_names))
    ),
    #=species_mass_flows = NamedTuple{species_names}(
        Tuple(zeros(n_cells)u"kg" for _ in 1:length(species_names))
    ),=#
)

du0_vec, u0_vec, state_axes, system, geo = finish_fvm_config(config, connection_map_function, special_caches, check_units = false);

function solve_system!(du, u, p, t, system, geo)
    properties = ComponentVector(system.properties_vec, system.properties_axes)
    u.rho .= properties.rho
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_controller_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end


f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

p_guess = 0.0

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 1000.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
#@time sol = solve(implicit_prob, AutoTsit5(FBDF()), callback = approximate_time_to_finish_cb)

prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 1e-2), p_guess)
VSCodeServer.@profview sol = solve(prob, Tsit5(), callback = approximate_time_to_finish_cb)

u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)]

sim_file = @__FILE__

sol_to_vtk(sol, u_named, grid, sim_file)