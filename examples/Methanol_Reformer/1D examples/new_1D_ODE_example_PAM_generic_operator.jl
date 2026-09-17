using Unitful
using OrdinaryDiffEq
using Ferrite
using SparseConnectivityTracer
using ComponentArrays
import ADTypes

using FVMFramework

total_pipe_length = 16.5u"cm"
stripped_pipe_length = ustrip(total_pipe_length |> u"m")
pipe_width = ustrip(0.5u"cm" |> u"m")

grid_dimensions = (100, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((stripped_pipe_length, pipe_width, pipe_width))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

total_pipe_segments = 100

addcellset!(grid, "inlet", xyz -> xyz[1] <= (1 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "inlet")

addcellset!(grid, "vaporization_area", xyz -> xyz[1] >= (1 * (stripped_pipe_length / total_pipe_segments)) && xyz[1] <= (19 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "vaporization_area")

addcellset!(grid, "reforming_area", xyz -> xyz[1] >= (19 * (stripped_pipe_length / total_pipe_segments)) && xyz[1] <= (99 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "reforming_area")

addcellset!(grid, "outlet", xyz -> xyz[1] >= (99 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "outlet")

n_cells = total_pipe_segments
u_proto = ComponentVector(
    mass_fractions = (
        methanol = zeros(n_cells)u"kg/kg",
        water = zeros(n_cells)u"kg/kg",
        carbon_monoxide = zeros(n_cells)u"kg/kg",
        hydrogen = zeros(n_cells)u"kg/kg",
        carbon_dioxide = zeros(n_cells)u"kg/kg",
        air = zeros(n_cells)u"kg/kg"
    ),
    pressure = zeros(n_cells)u"Pa",
    temp = zeros(n_cells)u"K"
)

config = create_fvm_config(grid, u_proto)

van_t_hoff_A = ComponentVector(CH3O = 1.7e-6u"s^-1", HCOO = 4.74e-13u"s^-1", OH = 3.32e-14u"s^-1")
van_t_hoff_dH = ComponentVector(CH3O = -46800.0u"J/mol", HCOO = -115000.0u"J/mol", OH = -110000.0u"J/mol")

MSR_rxn = ComponentVector(
    heat_of_reaction = 49500.0u"J/mol", 
    ref_delta_G = -3800.0u"J/mol", 
    ref_temp = 298.15u"K", 
    kf_A = 1.59e10u"s^-1", #sources online point to values around 1.25e7 mol / (kg * s * bar)
    kf_Ea = 103000.0u"J/mol",
    reactant_stoich_coeffs = (methanol = 1, water = 1),
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3), 
    stoich_coeffs = (methanol = -1, water = -1, hydrogen = 3, carbon_dioxide = 1), # CH3OH + H2O ⇋ CO2 + 3(H2)
    van_t_hoff_A = van_t_hoff_A, 
    van_t_hoff_dH = van_t_hoff_dH
)

MD_rxn = ComponentVector(
    heat_of_reaction = 90200.0u"J/mol", 
    ref_delta_G = 24800.0u"J/mol", 
    ref_temp = 298.15u"K", 
    kf_A = 1.46e13u"s^-1", #sources online point to values around 1.15e11 mol / (kg * s * bar)
    kf_Ea = 170000.0u"J/mol",
    reactant_stoich_coeffs = (methanol = 1,), 
    product_stoich_coeffs = (carbon_monoxide = 1, hydrogen = 2), 
    stoich_coeffs = (methanol = -1, carbon_monoxide = 1, hydrogen = 2), # CH3OH ⇋ CO + 2(H2)
    van_t_hoff_A = van_t_hoff_A,
    van_t_hoff_dH = van_t_hoff_dH
)

WGS_rxn = ComponentVector(
    heat_of_reaction = -41100.0u"J/mol", 
    ref_delta_G = -28600.0u"J/mol", 
    ref_temp = 298.15u"K", 
    kf_A = 4.63e10u"s^-1", #sources online point to values around 3.65e7 mol / (kg * s * bar)
    kf_Ea = 87500.0u"J/mol",
    reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1), 
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1), 
    stoich_coeffs = (water = -1, carbon_monoxide = -1, hydrogen = 1, carbon_dioxide = 1), # CO + H2O ⇋ CO2 + H2
    van_t_hoff_A = van_t_hoff_A, 
    van_t_hoff_dH = van_t_hoff_dH
)

initial_mass_fractions = ComponentVector(
    methanol = 1.0u"kg/kg",
    water = 1.3u"kg/kg",
    carbon_monoxide = 0.0001u"kg/kg",
    hydrogen = 0.02u"kg/kg",
    carbon_dioxide = 0.0001u"kg/kg",
    air = 0.0001u"kg/kg"
)

empty_mass_fractions = ComponentVector(
    methanol = 1e-20u"kg/kg",
    water = 1e-20u"kg/kg",
    carbon_monoxide = 1e-20u"kg/kg",
    hydrogen = 1e-6u"kg/kg",
    carbon_dioxide = 1e-20u"kg/kg",
    air = 1.0u"kg/kg"
)

initial_total_mass_fractions = sum(initial_mass_fractions)
empty_total_mass_fractions = sum(empty_mass_fractions)

initial_mass_fractions ./= initial_total_mass_fractions
empty_mass_fractions ./= empty_total_mass_fractions

total_pipe_length = 16.5u"cm"
n_segments = length(grid.cells)

pipe_width = 0.5u"cm"
pipe_height = 0.5u"cm"
pipe_area = pipe_width * pipe_height

pipe_mass_flow = 0.0278u"g/s"
pipe_volumetric_flow = pipe_mass_flow / 0.54u"kg/m^3"
velocity = pipe_volumetric_flow / pipe_area
superficial_mass_velocity = pipe_mass_flow / pipe_area #forgot to use bed void fraction here

pipe_hydraulic_diameter = (2 * pipe_width * pipe_height) / (pipe_width + pipe_height)

pipe_thickness = 1.0u"cm"
pipe_k = 237.0u"W/(m*K)"

approximate_residence_time = total_pipe_length / velocity |> u"s"   

reforming_area_properties = ComponentVector(
    k = 0.025u"W/(m*K)", 
    cp = 4184u"J/(kg*K)",
    mu = 1e-5u"Pa*s",
    rho = 791.0u"kg/m^3",
    viscosity = 1e-5u"Pa*s",
    R_gas = 8.314u"J/(mol*K)", 
    
    pipe_mass_flow = pipe_mass_flow,
    pipe_area = pipe_area,
    pipe_length = total_pipe_length / n_segments,
    pipe_hydraulic_diameter = pipe_hydraulic_diameter,
    pipe_thickness = pipe_thickness,
    pipe_inside_diameter = pipe_hydraulic_diameter,
    pipe_outside_diameter = pipe_hydraulic_diameter + 2 * pipe_thickness,
    pipe_k = pipe_k,
    velocity = velocity,
    superficial_mass_velocity = superficial_mass_velocity,
    bed_void_fraction = 0.80,
    catalyst_particle_diameter = 1.0u"mm",
    external_temp = 270.0u"°C",

    diffusion_coefficients = (
        methanol = 1e-5u"m^2/s",
        water = 1e-5u"m^2/s",
        carbon_monoxide = 1e-5u"m^2/s",
        hydrogen = 1e-5u"m^2/s",
        carbon_dioxide = 1e-5u"m^2/s",
        air = 1e-5u"m^2/s"
    ), 
    molecular_weights = (
        methanol = 32.04u"g/mol",
        water = 18.02u"g/mol",
        carbon_monoxide = 28.01u"g/mol",
        hydrogen = 2.02u"g/mol",
        carbon_dioxide = 44.01u"g/mol",
        air = 28.97u"g/mol"
    ), 
    reactions = (reforming_reactions = (MSR_rxn = MSR_rxn, MD_rxn = MD_rxn, WGS_rxn = WGS_rxn),),
    reactions_kg_cat = (reforming_reactions = (MSR_rxn = 1250.0u"kg/m^3", MD_rxn = 1250.0u"kg/m^3", WGS_rxn = 1250.0u"kg/m^3"),), 
)
#remember, for any named tuple with a single field inside it, remember to add a comma to the end 
#good: (reforming_reactions = (MSR_rxn = MSR_rxn,),)
#bad: (reforming_reaction = (MSR_rxn = MSR_rxn))

#these are just for classifying regions to make sure they do the right connection functions
struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "1D physics", "ergun_pressure_drop.jl"))
Revise.includet(joinpath(@__DIR__, "1D physics", "empirical_hx_correlations.jl"))

function update_properties!(u, cell_id)
    mw_avg!(u, cell_id)
    rho_ideal!(u, cell_id)
    molar_concentrations!(u, cell_id)
end

#you could also add any of these functions individually

function sum_and_cap_fluxes!(du, u, cell_id, vol)
    sum_mass_flux_face_to_cell!(du, u, cell_id)

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
end

common_cache_syms_and_units = (
    heat = u"J",
    mw_avg = u"kg/mol",
    rho = u"kg/m^3",
    molar_concentrations = u"mol/m^3",
    species_mass_flows = u"kg/s",
    net_rates = u"mol/s",
    mass = u"kg",
    mass_face = u"kg",
    Pr = u"NoUnits"
)

add_region!(
    config, "inlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = initial_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
    region_function =
    function inlet!(du, u, cell_id, vol)
        update_properties!(u, cell_id)

        du.mass_face[cell_id, 5] += u.pipe_mass_flow[cell_id]

        du.heat[cell_id] *= 0.0

        sum_and_cap_fluxes!(du, u, cell_id, vol)

        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
        end
    end
)

add_region!(
    config, "vaporization_area";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
    region_function =
    function vaporization_area!(du, u, cell_id, vol)
        update_properties!(u, cell_id)

        UA = overall_heat_transfer_coefficient(du, u, cell_id, vol) #W/(m^2 * K)
        surface_area = pi * u.pipe_inside_diameter[cell_id] * u.pipe_length[cell_id]
        du.heat[cell_id] += UA * (u.external_temp[cell_id] - u.temp[cell_id]) * surface_area

        sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "reforming_area";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        update_properties!(u, cell_id)

        PAM_reforming_react_cell!(du, u, cell_id, vol)

        UA = overall_heat_transfer_coefficient(du, u, cell_id, vol) #W/(m^2 * K)
        surface_area = pi * u.pipe_inside_diameter[cell_id] * u.pipe_length[cell_id]
        du.heat[cell_id] += UA * (u.external_temp[cell_id] - u.temp[cell_id]) * surface_area

        sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "outlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
    ),
    properties = reforming_area_properties,
    optimized_syms = (),
    cache_syms_and_units = common_cache_syms_and_units,
    region_function =
    function outlet!(du, u, cell_id, vol)
        update_properties!(u, cell_id)

        du.mass_face[cell_id, 3] -= u.pipe_mass_flow[cell_id]
        #this is simulating the mass flow out of the system

        for_fields!(u.mass_fractions, du.species_mass_flows) do species, u_mass_fractions, du_species_mass_flows
            du_species_mass_flows[species[cell_id]] -= u.pipe_mass_flow[cell_id] * u_mass_fractions[species[cell_id]]
        end
        #this is to prevent the concentration of all species from building up at the outlet

        du.heat[cell_id] -= u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id] 
        #this is to prevent the temperature from building up at the outlet

        UA = overall_heat_transfer_coefficient(du, u, cell_id, vol) #W/(m^2 * K)
        surface_area = pi * u.pipe_inside_diameter[cell_id] * u.pipe_length[cell_id]
        du.heat[cell_id] += UA * (u.external_temp[cell_id] - u.temp[cell_id]) * surface_area

        sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
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
    )

    #=mass_fraction_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#
end

function connection_map_function(type_a, type_b)
    typeof(type_a) <: Fluid && typeof(type_b) <: Fluid && return fluid_fluid_flux!
end

n_faces = length(config.geo.cell_neighbor_areas[1])
n_cells = length(config.geo.cell_volumes)
reaction_names = keys(config.regions[1].properties.reactions.reforming_reactions)
species_names = keys(config.regions[1].properties.molecular_weights)
 
#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = ComponentArray(
    mass_face = zeros(n_cells, n_faces)u"kg",
    net_rates = (
        reforming_reactions = NamedTuple{reaction_names}(
            Tuple(zeros(n_cells)u"mol/s" for _ in 1:length(reaction_names))
        ), #don't forget the comma!
    ), 
    molar_concentrations = NamedTuple{species_names}(
        Tuple(zeros(n_cells)u"mol/m^3" for _ in 1:length(species_names))
    ),
    species_mass_flows = NamedTuple{species_names}(
        Tuple(zeros(n_cells)u"kg" for _ in 1:length(species_names))
    ),
)

#you can check units by setting check_units = true and du0_vec and u0_vec will be returned as unitful ComponentVectors
du0_vec, u0_vec, state_axes, system, geo = finish_fvm_config(config, connection_map_function, special_caches, check_units = false);

function one_dimensional_pipe_f!(du, u, p, t, system, geo)
    #du, u = unpack_fvm_state(du_vec, u_vec, p, t, system) 
    #I put this here just so you know that you can modify this if you want to
    #that would probably only be useful when doing optimization 
    #For example, making u.diffusion_coefficients[:] .= p.diffusion_coefficients
    
    for cell_id in 1:length(geo.cell_volumes)-1
        du.mass_face[cell_id, 3] -= u.pipe_mass_flow[cell_id]
        du.mass_face[cell_id + 1, 5] += u.pipe_mass_flow[cell_id]
    end

    solve_connection_groups!(du, u, p, t, system, geo)
    solve_controller_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)

    #you could also just do: default_order_solve_all_groups!(du, u, p, t, system, geo),
    #but then you can't print stuff between group type loops.
    #For example, you can't do @show u.heat[1:5] after the connection groups are solved to see if some heat flux is screwing everything up.

    #also, if you don't need any custom logic at all, you can directly pass default_order_solve_all_groups!() to the closure instead of defining this function.
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, one_dimensional_pipe_f!, system, geo)

p_guess = 0.0

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 100.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF(), callback = approximate_time_to_finish_cb)
#VSCodeServer.@profview solve(implicit_prob, FBDF(), callback = approximate_time_to_finish_cb)

u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)]

record_sol = true

sim_file = @__FILE__
#u_named[end].mass_fractions.methanol[1]

u_named[end].mass_fractions.methanol[1]
u_named[end].mass_fractions.methanol[99]

conversion = ((u_named[end].mass_fractions.methanol[1] - u_named[end].mass_fractions.methanol[99]) / u_named[end].mass_fractions.methanol[1]
)

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end

using BenchmarkTools

explicit_prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 1e-2), p_guess)
@time sol = solve(explicit_prob, AutoTsit5(FBDF()))#, callback = approximate_time_to_finish_cb)