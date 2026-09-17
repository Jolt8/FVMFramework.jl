using FVMFramework
using Unitful
using OrdinaryDiffEq
using Ferrite
using SparseConnectivityTracer
import ADTypes

total_pipe_length = 50.0u"cm" |> u"m"
stripped_pipe_length = ustrip(total_pipe_length)
pipe_width = ustrip(1.0u"cm" |> u"m")

grid_dimensions = (100, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((stripped_pipe_length, pipe_width, pipe_width))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

total_pipe_segments = 100

addcellset!(grid, "inlet", xyz -> xyz[1] <= (1 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "inlet")

addcellset!(grid, "vaporization_area", xyz -> xyz[1] >= (1 * (stripped_pipe_length / total_pipe_segments)) && xyz[1] <= (20 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "vaporization_area")

addcellset!(grid, "reforming_area", xyz -> xyz[1] >= (20 * (stripped_pipe_length / total_pipe_segments)) && xyz[1] <= (99 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "reforming_area")

addcellset!(grid, "outlet", xyz -> xyz[1] >= (99 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "outlet")

n_cells = total_pipe_segments
u_proto = (
    mass_fractions = (
        methanol = zeros(n_cells),
        water = zeros(n_cells),
        carbon_monoxide = zeros(n_cells),
        hydrogen = zeros(n_cells),
        carbon_dioxide = zeros(n_cells),
        air = zeros(n_cells)
    ),
    pressure = zeros(n_cells),
    temp = zeros(n_cells)
)

config = create_fvm_config(grid, u_proto)

van_t_hoff_A = (CH3O = 1.7e-6, HCOO = 4.74e-13, OH = 3.32e-14)
van_t_hoff_dH = (CH3O = -46800.0, HCOO = -115000.0, OH = -110000.0)

function K_gibbs_free(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
    K_ref = exp(-ΔG_rxn_ref / (R_gas * T_ref))

    ln_K_ratio = (-ΔH_rxn_ref / R_gas) * (1/T_actual - 1/T_ref)
    
    K_T = K_ref * exp(ln_K_ratio)
    
    return ustrip(K_T) # Return a unitless equilibrium constant
end

MSR_rxn = (
    heat_of_reaction = 49500.0, # [J/mol]
    ref_delta_G = -3800.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 1.25e8,#1.25e7, # [s^-1] #sources online point to values around 1.25e7 mol / (kg * s * bar)
    kf_Ea = 103000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs = (methanol = -1, water = -1, hydrogen = 3, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

K_gibbs_free(MSR_rxn.ref_temp, ustrip(270u"°C" |> u"K"), MSR_rxn.ref_delta_G, MSR_rxn.heat_of_reaction)

MD_rxn = (
    heat_of_reaction = 90200.0, # [J/mol]
    ref_delta_G = 24800.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 1.15e12,#1.15e11, # [s^-1] #sources online point to values around 1.15e11 mol / (kg * s * bar)
    kf_Ea = 170000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1,), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_monoxide = 1, hydrogen = 2), # product_stoich_coeffs: 1 CO + 2 H2
    stoich_coeffs = (methanol = -1, carbon_monoxide = 1, hydrogen = 2), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

K_gibbs_free(MD_rxn.ref_temp, ustrip(270u"°C" |> u"K"), MD_rxn.ref_delta_G, MD_rxn.heat_of_reaction)

WGS_rxn = (
    heat_of_reaction = -41100.0, # [J/mol]
    ref_delta_G = -28600.0, # [J/mol]
    ref_temp = 298.15, # [K]
    kf_A = 3.65e8, #3.65e7, # [s^-1] #sources online point to values around 3.65e7 mol / (kg * s * bar)
    kf_Ea = 87500.0, # [J/mol]
    reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1), # product_stoich_coeffs: 1 CO2 + 1 H2
    stoich_coeffs = (water = -1, carbon_monoxide = -1, hydrogen = 1, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A = van_t_hoff_A,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH = van_t_hoff_dH # dH vector (CH3O, HCOO, OH) [J/mol]
)

K_gibbs_free(WGS_rxn.ref_temp, ustrip(270u"°C" |> u"K"), WGS_rxn.ref_delta_G, WGS_rxn.heat_of_reaction)

#since each mass fraction is modified, they have to be vectors
initial_mass_fractions = (
    methanol = [1.0],
    water = [1.3],
    carbon_monoxide = [0.0001],
    hydrogen = [0.02],
    carbon_dioxide = [0.0001],
    air = [0.0001]
)

empty_mass_fractions = (
    methanol = [1e-20],
    water = [1e-20],
    carbon_monoxide = [1e-20],
    hydrogen = [1e-6],
    carbon_dioxide = [1e-20],
    air = [1.0]
)

initial_total_mass_fractions = 0.0
empty_total_mass_fractions = 0.0

for (species_name, mass_fraction) in pairs(initial_mass_fractions)
    initial_total_mass_fractions += initial_mass_fractions[species_name][1]
end

for (species_name, mass_fraction) in pairs(empty_mass_fractions)
    empty_total_mass_fractions += empty_mass_fractions[species_name][1]
end

for (species_name, mass_fraction) in pairs(initial_mass_fractions)
    initial_mass_fractions[species_name][1] /= initial_total_mass_fractions
end

for (species_name, mass_fraction) in pairs(empty_mass_fractions)
    empty_mass_fractions[species_name][1] /= empty_total_mass_fractions
end

initial_mass_fractions = NamedTuple{keys(initial_mass_fractions)}(first.(values(initial_mass_fractions)))
empty_mass_fractions = NamedTuple{keys(empty_mass_fractions)}(first.(values(empty_mass_fractions)))
#this is not ideal, but we can't use scalars inside the initial mass_fractions because the values of NamedTuples can't be modified 

total_pipe_length = 50.1u"cm" |> u"m"
n_segments = length(grid.cells)

pipe_width = 0.35u"cm" |> u"m"
pipe_height = 0.35u"cm" |> u"m"
pipe_area = pipe_width * pipe_height

pipe_mass_flow = 0.0278u"g/s" |> u"kg/s"
pipe_volumetric_flow = pipe_mass_flow / 791u"kg/m^3"
velocity = pipe_volumetric_flow / pipe_area
superficial_mass_velocity = pipe_mass_flow / pipe_area #forgot to use bed void fraction here

pipe_hydraulic_diameter = (2 * pipe_width * pipe_height) / (pipe_width + pipe_height)

pipe_thickness = 1.0u"cm" |> u"m"
pipe_k = 237.0u"W/(m*K)"

approximate_residence_time = total_pipe_length / velocity

reforming_area_properties = (
    k = 237.0, # k (W/(m*K))
    cp = 4.184, # cp (J/(kg*K))
    mu = 1e-5, # mu (Pa*s)
    rho = 791.0, # rho (kg/m^3)
    viscosity = ustrip(1e-5u"Pa*s" |> u"Pa*s"),
    
    pipe_mass_flow = ustrip(pipe_mass_flow |> u"kg/s"),
    pipe_area = ustrip(pipe_area |> u"m^2"),
    pipe_length = ustrip((total_pipe_length / n_segments) |> u"m"),
    pipe_hydraulic_diameter = ustrip(pipe_hydraulic_diameter |> u"m"),
    pipe_thickness = ustrip(pipe_thickness |> u"m"),
    pipe_k = ustrip(pipe_k |> u"W/(m*K)"),
    velocity = ustrip(velocity |> u"m/s"),
    superficial_mass_velocity = ustrip(superficial_mass_velocity |> u"kg/(m^2*s)"),
    bed_void_fraction = 0.80,
    catalyst_particle_diameter = ustrip(1.0u"mm" |> u"m"),
    external_temp = ustrip(270.0u"°C" |> u"K"),
    
    species_ids = (methanol = 1, water = 2, carbon_monoxide = 3, hydrogen = 4, carbon_dioxide = 5, air = 6), #we could use mass_fractions for species loops, but this is just more consistent
    diffusion_coefficients = (
        methanol = 1e-5,
        water = 1e-5,
        carbon_monoxide = 1e-5,
        hydrogen = 1e-5,
        carbon_dioxide = 1e-5,
        air = 1e-5
    ), #diffusion coefficients (m^2/s)
    molecular_weights = (
        methanol = 0.03204,
        water = 0.01802,
        carbon_monoxide = 0.02801,
        hydrogen = 0.00202,
        carbon_dioxide = 0.04401,
        air = 0.02897
    ), #species_molecular_weights [kg/mol]
    reactions = (reforming_reactions = (MSR_rxn = MSR_rxn, MD_rxn = MD_rxn, WGS_rxn = WGS_rxn),),
    reactions_kg_cat = (reforming_reactions = (MSR_rxn = 1250.0, MD_rxn = 1250.0, WGS_rxn = 1250.0),), # cell_kg_cat_per_m3_for_each_reaction
)
#remember, for any named tuple with a single field inside it, remember to add a comma to the end 
#ex. (reforming_reactions = (MSR_rxn = MSR_rxn,),)
#not: (reforming_reaction = (MSR_rxn = MSR_rxn))

#these are just for classifying regions to make sure they do the right connection functions
struct Inlet <: AbstractPhysics end
struct Fluid <: AbstractPhysics end
struct Outlet <: AbstractPhysics end

include(joinpath(@__DIR__, "1D physics", "ergun_pressure_drop.jl"))
include(joinpath(@__DIR__, "1D physics", "empirical_hx_correlations.jl"))

add_region!(
    config, "inlet";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = initial_mass_fractions,
        pressure = ustrip(1.0u"atm" |> u"Pa"),
        temp = ustrip(270.0u"°C" |> u"K"),
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :species_mass_flows, :net_rates, :mass, :mass_face, :Pr],
    region_function =
    function inlet!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics
        #ergun_pressure_drop!(du, u, cell_id, vol)

        #du.mass_face[cell_id][5] += u.pipe_mass_flow[cell_id]
        du.mass[cell_id] += u.pipe_mass_flow[cell_id]

        map(keys(u.species_mass_flows)) do species_name
            #du.mass_fractions[species_name][cell_id] += u.pipe_mass_flow[cell_id] * initial_mass_fractions[species_name][1]
            du.species_mass_flows[species_name][cell_id] *= 0.0
        end
        
        du.heat[cell_id] *= 0.0

        sum_mass_flux_face_to_cell!(du, u, cell_id)

        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
        cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "vaporization_area";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = empty_mass_fractions,
        pressure = ustrip(1.0u"atm" |> u"Pa"),
        temp = ustrip(270.0u"°C" |> u"K"),
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :species_mass_flows, :net_rates, :mass, :mass_face, :Pr],
    region_function =
    function vaporization_area!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics
        #prandtl_number!(du, u, cell_id, vol)
        UA = UA_prime(du, u, cell_id, vol)
        du.heat[cell_id] += UA * (u.external_temp[cell_id] - u.temp[cell_id]) * u.pipe_area[cell_id]

        sum_mass_flux_face_to_cell!(du, u, cell_id)

        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
        cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    end
)

species_names = keys(config.regions[1].properties.species_ids)

add_region!(
    config, "reforming_area";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = empty_mass_fractions,
        pressure = ustrip(1.0u"atm" |> u"Pa"),
        temp = ustrip(270.0u"°C" |> u"K"),
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :species_mass_flows, :net_rates, :mass, :mass_face, :Pr],
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics
        #println(du.mass_fractions.methanol[cell_id])
        PAM_reforming_react_cell!(du, u, cell_id, vol)
        #println(du.mass_fractions.methanol[cell_id])

        #prandtl_number!(du, u, cell_id, vol)
        UA = UA_prime(du, u, cell_id, vol)
        du.heat[cell_id] += UA * (u.external_temp[cell_id] - u.temp[cell_id]) * u.pipe_area[cell_id]

        sum_mass_flux_face_to_cell!(du, u, cell_id)

        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
        cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "outlet";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = empty_mass_fractions,
        pressure = ustrip(1.0u"atm" |> u"Pa"),
        temp = ustrip(270.0u"°C" |> u"K"),
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :species_mass_flows, :net_rates, :mass, :mass_face, :Pr],
    region_function =
    function outlet!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics
        #du.mass_face[cell_id][5] -= u.pipe_mass_flow[cell_id]
        du.mass[cell_id] -= u.pipe_mass_flow[cell_id]
       
        map(keys(u.species_mass_flows)) do species_name
            #du.mass_fractions[species_name][cell_id] += u.pipe_mass_flow[cell_id] * initial_mass_fractions[species_name][1]
            du.species_mass_flows[species_name][cell_id] *= 0.0
        end

        #UA = UA_prime(du, u, cell_id, vol)
        #du.heat[cell_id] -= UA * (u.external_temp[cell_id] - u.temp[cell_id]) * u.pipe_area[cell_id]
        du.heat[cell_id] *= 0.0

        #ergun_pressure_drop!(du, u, cell_id, vol)

        sum_mass_flux_face_to_cell!(du, u, cell_id)

        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
        cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
        cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
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
    #=
    enthalpy_advection!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )
        =#

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
config.regions[1].properties.reactions
n_reactions = length(config.regions[1].properties.reactions.reforming_reactions)
reaction_names = keys(config.regions[1].properties.reactions.reforming_reactions)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = (reforming_reactions = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))), ), 
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
    species_mass_flows = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))),
)

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, special_caches)

system.connection_groups[1].cell_neighbors[end]

f_closure_implicit = (du, u, p, t) -> pipe_f!(
    du, u, p, t, 

    geo.cell_volumes, geo.cell_centroids,
    geo.cell_neighbor_areas, geo.cell_neighbor_normals, geo.cell_neighbor_distances,
    geo.unconnected_cell_face_map, geo.cell_face_areas, geo.cell_face_normals,

    system.connection_groups, system.controller_groups, system.region_groups, system.patch_groups,
    system.merged_properties,

    system.du_diff_cache, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)

p_guess = 0.0
#=
prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 5.0), p_guess)
#@time sol = solve(prob, Tsit5(), callback = approximate_time_to_finish_cb)

sol_u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)
sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

sol_u_named_0.pressure
sol_u_named_end.pressure

sol_u_named_0.mass_fractions.methanol
sol_u_named_end.mass_fractions.methanol

sol_u_named_0.mass_fractions.methanol[80] - sol_u_named_end.mass_fractions.methanol[80]

methanol_conversion = ((sol_u_named_0.mass_fractions.methanol[80] - sol_u_named_end.mass_fractions.methanol[80]) / sol_u_named_0.mass_fractions.methanol[80]
)

sol_u_named_0.mass_fractions.hydrogen
sol_u_named_end.mass_fractions.hydrogen
=#

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 10000.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

desired_steps = 100
save_interval = (tspan[end] / desired_steps)

#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)
sol = solve(implicit_prob, FBDF(), callback = approximate_time_to_finish_cb)

record_sol = true

sim_file = @__FILE__

u_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

u_named[end].mass_fractions.methanol[1]
u_named[end].mass_fractions.methanol[99]

inlet_cell_id = 1
outlet_cell_id = 99

conversion = ((u_named[end].mass_fractions.methanol[inlet_cell_id] - u_named[end].mass_fractions.methanol[outlet_cell_id]) / u_named[end].mass_fractions.methanol[inlet_cell_id]
)

if record_sol == true
    sol_to_vtk(sol, u_named, grid, sim_file)
end

hi = 1


#=
function conversion_from_residence_time(tMax)
    tspan = (t0, tMax)
    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)
    sol = solve(implicit_prob, FBDF())
    
    sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

    methanol_conversion = ((sol_u_named_0.mass_fractions.methanol[80] - sol_u_named_end.mass_fractions.methanol[80]) / sol_u_named_0.mass_fractions.methanol[80])

    #oh wait, I just realized that this doesn't make any sense because the conversion is spatial rather than temporal

    return methanol_conversion
end

using Roots
desired_conversion = 0.95
time = 1.0

#find_zero(time -> conversion_from_residence_time(time) - desired_conversion, time)


using Plots

plot(sol.t, [u.mass_fractions.methanol[80] for u in u_named], label = "Methanol")
plot!(sol.t, [u.mass_fractions.water[80] for u in u_named], label = "Water")
plot!(sol.t, [u.mass_fractions.carbon_monoxide[80] for u in u_named], label = "Carbon Monoxide")
plot!(sol.t, [u.mass_fractions.hydrogen[80] for u in u_named], label = "Hydrogen")
plot!(sol.t, [u.mass_fractions.carbon_dioxide[80] for u in u_named], label = "Carbon Dioxide")

methanol_molar_concentration = [(u.mass_fractions.methanol[80] * 1000 / reforming_area_properties.molecular_weights.methanol) for u in u_named]
water_molar_concentration = [(u.mass_fractions.water[80] * 1000 / reforming_area_properties.molecular_weights.water) for u in u_named]
carbon_monoxide_molar_concentration = [(u.mass_fractions.carbon_monoxide[80] * 1000 / reforming_area_properties.molecular_weights.carbon_monoxide) for u in u_named]
hydrogen_molar_concentration = [(u.mass_fractions.hydrogen[80] * 1000 / reforming_area_properties.molecular_weights.hydrogen) for u in u_named]
carbon_dioxide_molar_concentration = [(u.mass_fractions.carbon_dioxide[80] * 1000 / reforming_area_properties.molecular_weights.carbon_dioxide) for u in u_named]

plot(sol.t, methanol_molar_concentration, label = "Methanol")
plot!(sol.t, water_molar_concentration, label = "Water")
plot!(sol.t, carbon_monoxide_molar_concentration, label = "Carbon Monoxide")
plot!(sol.t, hydrogen_molar_concentration, label = "Hydrogen")
plot!(sol.t, carbon_dioxide_molar_concentration, label = "Carbon Dioxide")

total_molar_concentration = methanol_molar_concentration + water_molar_concentration + carbon_monoxide_molar_concentration + hydrogen_molar_concentration + carbon_dioxide_molar_concentration

methanol_molar_fractions = methanol_molar_concentration ./ total_molar_concentration
water_molar_fractions = water_molar_concentration ./ total_molar_concentration
carbon_monoxide_molar_fractions = carbon_monoxide_molar_concentration ./ total_molar_concentration
hydrogen_molar_fractions = hydrogen_molar_concentration ./ total_molar_concentration
carbon_dioxide_molar_fractions = carbon_dioxide_molar_concentration ./ total_molar_concentration

plot(sol.t, methanol_molar_fractions, label = "Methanol")
plot!(sol.t, water_molar_fractions, label = "Water")
plot!(sol.t, carbon_monoxide_molar_fractions, label = "Carbon Monoxide")
plot!(sol.t, hydrogen_molar_fractions, label = "Hydrogen")
plot!(sol.t, carbon_dioxide_molar_fractions, label = "Carbon Dioxide")


1.0e6u"mol/(kg*s*bar)" |> u"mol/(kg*s*Pa)"

1.0e6u"bar/m" |> u"Pa/m"
1.0e6u"m/bar" |> u"m/Pa"
=#