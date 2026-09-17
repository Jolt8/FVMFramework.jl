using FVMFramework
using Unitful
using OrdinaryDiffEq

#Chemistry functions and definitions 
function arrenhius_equation_pre_exponential_factor(k, Ea, T)
    R_GAS = 8.314u"J/(mol*K)"
    return (k / exp(-Ea / (R_GAS * T)))
end

function arrenhius_equation_rate_constant(A, Ea, T)
    R_GAS = 8.314u"J/(mol*K)"
    return (A * exp(-Ea / (R_GAS * T)))
end

ref_T_SMR = 700.0u"°C" |> u"K" 
kf_ref_smr = 5.0e-1u"mol/(kg*s*bar^2)" 
Ea_f_smr = 240.0u"kJ/mol" 

# CH4 + H2O <-> CO + 3H2
SMR_reaction = (
    heat_of_reaction = 206.0u"kJ/mol" |> u"J/mol", # [J/mol]
    ref_delta_G = 142.0u"kJ/mol" |> u"J/mol", # [J/mol]
    ref_temp = ref_T_SMR, # [K]
    kf_A = arrenhius_equation_pre_exponential_factor(kf_ref_smr, Ea_f_smr, ref_T_SMR), # [s^-1]
    kf_Ea = 103000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs = (methanol = -1, water = -1, carbon_monoxide = 0, hydrogen = 3, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
)


ref_T_WGS = 300.13u"°C" |> u"K"
kf_ref_wgs = 0.2e-1u"mol/(kg*s*bar^2)"
Ea_f_wgs = 60u"kJ/mol"
kf_A_wgs = arrenhius_equation_pre_exponential_factor(kf_ref_wgs, Ea_f_wgs, ref_T_WGS)

# CO + H2O <-> CO2 + H2
WGS_reaction = (
    heat_of_reaction = -41.1u"kJ/mol" |> u"J/mol", # [J/mol]
    ref_delta_G = -28.63u"kJ/mol" |> u"J/mol", # [J/mol]
    ref_temp = ref_T_WGS, # [K]
    kf_A = arrenhius_equation_pre_exponential_factor(kf_ref_wgs, Ea_f_wgs, ref_T_WGS), # [s^-1]
    kf_Ea = 60000.0, # [J/mol]
    reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs = (carbon_monoxide = -1, water = -1, hydrogen = 1, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
)

# --- Reaction 3: Coking (Methane Cracking) ---
# CH4 -> C(s) + 2H2
# Note: Irreversible approximation for this estimation
ref_T_Coking = 700.0u"°C" |> u"K"
kf_ref_coke = 1.0e-4u"mol/(kg*s*bar)" # Usually slow compared to SMR if steam is high
Ea_f_coke = 150.0u"kJ/mol"
kf_A_coke = arrenhius_equation_pre_exponential_factor(kf_ref_coke, Ea_f_coke, ref_T_Coking)

Coking_reaction = (
    heat_of_reaction = 74.8u"kJ/mol" |> u"J/mol", # [J/mol]
    ref_delta_G = 142.0u"kJ/mol" |> u"J/mol", # [J/mol]
    ref_temp = ref_T_Coking, # [K]
    kf_A = arrenhius_equation_pre_exponential_factor(kf_ref_coke, Ea_f_coke, ref_T_Coking), # [s^-1]
    kf_Ea = 150000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 3),     # product_stoich_coeffs: 1 CO2 + 3 H2
    stoich_coeffs = (methanol = -1, water = -1, carbon_monoxide = 0, hydrogen = 3, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
)

using Ferrite

total_pipe_length = 50.0u"cm" |> u"m"
stripped_pipe_length = ustrip(total_pipe_length)

grid_dimensions = (100, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((stripped_pipe_length, stripped_pipe_length, stripped_pipe_length))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

total_pipe_segments = 100

addcellset!(grid, "vaporization_area", xyz -> xyz[1] <= (20 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "vaporization_area")
addcellset!(grid, "reforming_area", xyz -> xyz[1] >= (20 * (stripped_pipe_length / total_pipe_segments)) && xyz[1] <= (100 * (stripped_pipe_length / total_pipe_segments)))
getcellset(grid, "reforming_area")

n_cells = total_pipe_segments
u_proto = (
    mass_fractions = (
        methanol = zeros(n_cells),
        water = zeros(n_cells),
        carbon_monoxide = zeros(n_cells),
        hydrogen = zeros(n_cells),
        carbon_dioxide = zeros(n_cells)
    ),
    pressure = zeros(n_cells),
)

config = create_fvm_config(grid, u_proto)

# CH3OH -> CO + 2H2
ref_T = 300.13u"°C" |> u"K"
kf_ref = 1e-2u"mol/(kg*s*bar)" 
Ea_f = 90u"kJ/mol" 

arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T) |> u"mol/(kg*s*Pa)"


MD_reaction = (
    heat_of_reaction = -49400.0, # [J/mol]
    ref_delta_G = 33200.0, # [J/mol]
    ref_temp = ustrip(ref_T), # [K]
    kf_A = ustrip(arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T) |> u"mol/(kg*s*Pa)"), # [s^-1]
    kf_Ea = 90000.0, # [J/mol]
    reactant_stoich_coeffs = (methanol = 1,), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_monoxide = 1, hydrogen = 2), # product_stoich_coeffs: 1 CO2 + 1 H2
    stoich_coeffs = (methanol = -1, water = 0, carbon_monoxide = 1, hydrogen = 2, carbon_dioxide = 0),
)

K_ref = exp(-MD_reaction.ref_delta_G / (R_gas * MD_reaction.ref_temp)) 

ln_K_ratio = (-MD_reaction.heat_of_reaction / R_gas) * (1 / ustrip(270u"°C" |> u"K") - 1 / MD_reaction.ref_temp)

K_T = K_ref * exp(ln_K_ratio)

kf_A = MD_reaction.kf_A

kr_A = (MD_reaction.kf_A / K_T) * exp(-MD_reaction.heat_of_reaction / (R_gas * ustrip(270u"°C" |> u"K"))) 

#find reverse Ea
kr_Ea = MD_reaction.kf_Ea - MD_reaction.heat_of_reaction 


# CO + H2O ⇋ CO2 + H2
ref_T = 300.13u"°C" |> u"K"
kf_ref = 1e-2u"mol/(kg*s*bar^2)" #around 1e-6 to 1e-5
Ea_f = 60u"kJ/mol" #usually around 55-100 kJ/mol

WGS_reaction = (
    heat_of_reaction = -41100.0, # [J/mol]
    ref_delta_G = -28600.0, # [J/mol]
    ref_temp = ustrip(ref_T), # [K]
    kf_A = ustrip(arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T) |> u"mol/(kg*s*Pa^2)"), # [s^-1]
    kf_Ea = 60000.0, # [J/mol]
    reactant_stoich_coeffs = (carbon_monoxide = 1, water = 1), # reactant_stoich_coeffs
    product_stoich_coeffs = (carbon_dioxide = 1, hydrogen = 1), # product_stoich_coeffs: 1 CO2 + 1 H2
    stoich_coeffs = (methanol = 0, water = -1, carbon_monoxide = -1, hydrogen = 1, carbon_dioxide = 1), # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
)

#since each mass fraction is modified, they have to be vectors
initial_mass_fractions = (
    methanol = [1.0],
    water = [1.3],
    carbon_monoxide = [0.0001],
    hydrogen = [0.02],
    carbon_dioxide = [0.0001]
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

total_pipe_length = 50.1u"cm" |> u"m"
n_segments = length(grid.cells)

pipe_mass_flow = 0.0278u"g/s" |> u"kg/s"
pipe_area = 2.0u"cm^2" |> u"m^2"
superficial_mass_velocity = pipe_mass_flow / pipe_area

reforming_area_properties = (
    k = 237.0, # k (W/(m*K))
    cp = 4.184, # cp (J/(kg*K))
    mu = 1e-5, # mu (Pa*s)
    rho = 791.0, # rho (kg/m^3)
    temp = ustrip(270.0u"°C" |> u"K"),
    viscosity = ustrip(1e-5u"Pa*s" |> u"Pa*s"),
    
    pipe_mass_flow = ustrip(pipe_mass_flow |> u"kg/s"),
    pipe_area = ustrip(pipe_area |> u"m^2"),
    pipe_length = ustrip((total_pipe_length / n_segments) |> u"m"),
    superficial_mass_velocity = ustrip(superficial_mass_velocity |> u"kg/(m^2*s)"),
    bed_void_fraction = 0.80,
    catalyst_particle_diameter = ustrip(1.0u"mm" |> u"m"),
    
    species_ids = (methanol = 1, water = 2, carbon_monoxide = 3, hydrogen = 4, carbon_dioxide = 5), #we could use mass_fractions for species loops, but this is just more consistent
    molecular_weights = (
        methanol = 0.03204,
        water = 0.01802,
        carbon_monoxide = 0.02801,
        hydrogen = 0.00202,
        carbon_dioxide = 0.04401
    ), #species_molecular_weights [kg/mol]
    reactions = (MD_reaction = MD_reaction, WGS_reaction = WGS_reaction),
    reactions_kg_cat = (MD_reaction = 1250.0, WGS_reaction = 1250.0), # cell_kg_cat_per_m3_for_each_reaction
)
#remember, for any named tuple with a single field inside it, remember to add a comma to the end 
#ex. (reforming_reactions = (MSR_rxn = MSR_rxn, ),)
#not: (reforming_reaction = (MSR_rxn = MSR_rxn))

#these are just for classifying regions to make sure they do the right connection functions
struct Fluid <: AbstractPhysics end

include(joinpath(@__DIR__, "1D physics", "ergun_pressure_drop.jl"))

add_region!(
    config, "vaporization_area";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = initial_mass_fractions,
        pressure = ustrip(1.0u"atm" |> u"Pa"),
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :net_rates, :mass, :mass_face],
    region_function =
    function vaporization_area!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics
        #ergun_pressure_drop!(du, u, cell_id, vol)

        #sum_mass_flux_face_to_cell!(du, u, cell_id)
    end
)

add_region!(
    config, "reforming_area";
    type = Fluid(),
    initial_conditions = (
        mass_fractions = initial_mass_fractions,
        pressure = ustrip(1.0u"atm" |> u"Pa"),
    ),
    properties = reforming_area_properties,
    optimized_syms = [],
    cache_syms = [:heat, :mw_avg, :rho, :molar_concentrations, :net_rates, :mass, :mass_face],
    region_function =
    function reforming_area!(du, u, cell_id, vol)
        #property updating/retrieval
        mw_avg!(u, cell_id)
        rho_ideal!(u, cell_id)
        molar_concentrations!(u, cell_id)

        #internal physics
        power_law_react_cell!(du, u, cell_id, u.reactions.MD_reaction, :MD_reaction, vol)
        power_law_react_cell!(du, u, cell_id, u.reactions.WGS_reaction, :WGS_reaction, vol)

        #ergun_pressure_drop!(du, u, cell_id, vol)

        #sum_mass_flux_face_to_cell!(du, u, cell_id)
    end
)


#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    #=all_species_advection!(
        du, u, 
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#

    #=enthalpy_advection!(
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
n_reactions = length(config.regions[1].properties.reactions)
reaction_names = keys(config.regions[1].properties.reactions)
species_names = keys(config.regions[1].properties.species_ids)

#species caches are for things like mass_face, which has an entry for every face of every cell rather than entries for each cell
special_caches = (
    mass_face = fill(zeros(n_faces), n_cells),
    net_rates = NamedTuple{reaction_names}(fill(0.0, length(reaction_names))), 
    molar_concentrations = NamedTuple{species_names}(fill(zeros(n_cells), length(species_names))), #I'm starting to really enjoy these NamedTuple constructors
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

    system.du_diff_cache_vec, system.u_diff_cache_vec,
    system.du_proto_axes, system.u_proto_axes,
    system.du_cache_axes, system.u_cache_axes
)

p_guess = 0.0

prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 10000000.0), p_guess)
@time sol = solve(prob, Tsit5(), callback = approximate_time_to_finish_cb)

sol_u_named_0 = create_views_inline(sol.u[1], system.u_proto_axes)
sol_u_named_end = create_views_inline(sol.u[end], system.u_proto_axes)

sol_u_named_0.pressure
sol_u_named_end.pressure

sol_u_named_0.mass_fractions.methanol
sol_u_named_end.mass_fractions.methanol

sol_u_named_0.mass_fractions.methanol[80] - sol_u_named_end.mass_fractions.methanol[80]

sol_u_named_0.mass_fractions.hydrogen
sol_u_named_end.mass_fractions.hydrogen

record_sol = false

sim_file = @__FILE__

u_proto_named = [create_views_inline(sol.u[i], system.u_proto_axes) for i in eachindex(sol.u)]

if record_sol == true
    sol_to_vtk(sol, u_proto_named, grid, sim_file)
end


using Plots

plot(sol.t, [u.mass_fractions.methanol[80] for u in u_proto_named])
plot(sol.t, [u.mass_fractions.carbon_dioxide[80] for u in u_proto_named])

