using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ComponentArrays
import ADTypes
using NonlinearSolve

using FVMFramework

pipe_inside_diameter = 0.5u"inch" |> u"m"
pipe_length = 12.1u"inch" |> u"m"

stripped_pipe_length = ustrip(pipe_length |> u"m")
pipe_width = ustrip(pipe_inside_diameter |> u"m")

n_cells = 100

grid_dimensions = (1, 1, n_cells)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((pipe_width, pipe_width, stripped_pipe_length))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "inlet", xyz -> xyz[3] <= (1 * (stripped_pipe_length / n_cells)))
getcellset(grid, "inlet")

addcellset!(grid, "evaporator", xyz -> xyz[3] >= (1 * (stripped_pipe_length / n_cells)) && xyz[3] <= (19 * (stripped_pipe_length / n_cells)))
getcellset(grid, "evaporator")

addcellset!(grid, "reactor", xyz -> xyz[3] >= (19 * (stripped_pipe_length / n_cells)) && xyz[3] <= (99 * (stripped_pipe_length / n_cells)))
getcellset(grid, "reactor")

addcellset!(grid, "outlet", xyz -> xyz[3] >= (99 * (stripped_pipe_length / n_cells)))
getcellset(grid, "outlet")

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
    temp = zeros(n_cells)u"K",
    liquid_holdup = zeros(n_cells)u"m^3/m^3",
    gas_holdup = zeros(n_cells)u"m^3/m^3"
)

config = create_fvm_config(grid, u_proto);

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
    stoich_coeffs = (methanol = -1, water = -1, carbon_dioxide = 1, hydrogen = 3), # CH3OH + H2O ⇋ CO2 + 3(H2)
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
    stoich_coeffs = (carbon_monoxide = -1, water = -1, carbon_dioxide = 1, hydrogen = 1), # CO + H2O ⇋ CO2 + H2
    van_t_hoff_A = van_t_hoff_A, 
    van_t_hoff_dH = van_t_hoff_dH
)

inlet_mass_fractions = ComponentVector(
    methanol = 1.0u"kg/kg",
    water = 1.3u"kg/kg",
    carbon_monoxide = 0.0001u"kg/kg",
    hydrogen = 0.001u"kg/kg",
    carbon_dioxide = 0.0001u"kg/kg",
    air = 0.0001u"kg/kg"
)

total_inlet_mass_fractions = sum(inlet_mass_fractions)
inlet_mass_fractions ./= total_inlet_mass_fractions

empty_mass_fractions = ComponentVector(
    methanol = 1e-20u"kg/kg",
    water = 1e-20u"kg/kg",
    carbon_monoxide = 1e-20u"kg/kg",
    hydrogen = 1e-6u"kg/kg",
    carbon_dioxide = 1e-20u"kg/kg",
    air = 1.0u"kg/kg"
)

total_empty_mass_fractions = sum(empty_mass_fractions)
empty_mass_fractions ./= total_empty_mass_fractions

pipe_area = pi * (pipe_inside_diameter / 2)^2

cell_lengths_along_pipe = [config.geo.cell_centroids[i][3]u"m" for i in 1:length(config.geo.cell_centroids)] 
#this only works because the base of the pipe is at z = 0.0

reforming_area_properties = ComponentVector(
    k = 0.025u"W/(m*K)", 
    cp = 4184u"J/(kg*K)",
    mu = 1e-5u"Pa*s",
    rho = 791.0u"kg/m^3",
    viscosity = 1e-5u"Pa*s",
    R_gas = 8.314u"J/(mol*K)",

    pipe_mass_flow = 1.0u"g/minute",

    pipe_inside_diameter = pipe_inside_diameter,
    pipe_area = pipe_area,
    per_cell_pipe_length = pipe_length / n_cells,
    cell_lengths_along_pipe = cell_lengths_along_pipe,

    bed_void_fraction = 0.4,
    packing_surface_area = 100.0u"m^2/m^3",
    particle_diameter = 5u"mm",

    overall_heat_transfer_coefficient = 1000.0u"W/(m^2*K)",

    external_temp = 300.0u"°C",
    saturation_temp = 72.4u"°C",
    liquid_rho = 791.0u"kg/m^3",
    gas_rho = 1.225u"kg/m^3",
    mass_transfer_coeff_vap = 0.001u"kg/(m^2*s*K)",
    heat_of_vaporization = 1.5u"kJ/g",
    liquid_feed_mass_fractions = inlet_mass_fractions,

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

permeability = (reforming_area_properties.particle_diameter^2 * reforming_area_properties.bed_void_fraction^3) / (150.0 * (1.0 - reforming_area_properties.bed_void_fraction)^2)

reforming_area_properties = merge_properties(reforming_area_properties, ComponentVector(permeability = permeability))

struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "physics", "multiphase.jl"))
Revise.includet(joinpath(@__DIR__, "physics", "energy.jl"))
Revise.includet(joinpath(@__DIR__, "physics", "momentum.jl"))

function update_properties!(du, u, cell_id, vol)
    mw_avg!(u, cell_id)
    rho_ideal!(u, cell_id)
    #rho_multiphase!(du, u, cell_id, vol)
    molar_concentrations!(u, cell_id)
    update_velocity!(du, u, cell_id, vol)
end

#you could also add any of these functions individually

function sum_and_cap_fluxes!(du, u, cell_id, vol)
    sum_mass_flux_face_to_cell!(du, u, cell_id, vol) #this always has to go before cap_mass_flux_to_pressure_change!

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)

    cap_evaporation_rate_to_phase_holdup!(du, u, cell_id, vol)
end

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])
reaction_names = keys(reforming_area_properties.reactions.reforming_reactions)
species_names = keys(reforming_area_properties.molecular_weights)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        rho = u"kg/m^3",
        molar_concentrations = u"mol/m^3",
        species_masses = u"kg",
        net_rates = u"mol/s",
        mass = u"kg",
        mass_face = u"kg",
        mass_evaporated = u"kg",
        superficial_velocity = u"m/s"
    ),
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
        species_masses = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"kg" for _ in 1:length(species_names))
        )
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector()
)

function common_physics_functions!(du, u, cell_id, vol)
    #vaporization_model!(du, u, cell_id, vol)
    #ergun_momentum_friction!(du, u, cell_id, vol)
end

add_region!(
    config, "inlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = inlet_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
        liquid_holdup = 1.0,
        gas_holdup = 0.0
    ),
    properties = reforming_area_properties,
    region_function =
    function inlet!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 1] += u.pipe_mass_flow[cell_id]

        du.heat[cell_id] *= 0.0

        sum_and_cap_fluxes!(du, u, cell_id, vol)

        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
        end
    end
)

add_region!(
    config, "evaporator";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
        liquid_holdup = 1.0,
        gas_holdup = 0.0
    ),
    properties = reforming_area_properties,
    region_function =
    function evaporator!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        #UA = overall_heat_transfer_coefficient(du, u, cell_id, vol) #W/(m^2 * K)
        surface_area = pi * u.pipe_inside_diameter[cell_id] * u.per_cell_pipe_length[cell_id]
        du.heat[cell_id] += u.overall_heat_transfer_coefficient[cell_id] * (u.external_temp[cell_id] - u.temp[cell_id]) * surface_area
        #wall_heat_flux!(du, u, cell_id, vol)

        sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

add_region!(
    config, "reactor";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = empty_mass_fractions,
        pressure = 1.0u"atm",
        temp = 21.0u"°C",
        liquid_holdup = 1.0,
        gas_holdup = 0.0
    ),
    properties = reforming_area_properties,
    region_function =
    function reactor!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        PAM_reforming_react_cell!(du, u, cell_id, vol)

        #UA = overall_heat_transfer_coefficient(du, u, cell_id, vol) #W/(m^2 * K)
        surface_area = pi * u.pipe_inside_diameter[cell_id] * u.per_cell_pipe_length[cell_id]
        du.heat[cell_id] += u.overall_heat_transfer_coefficient[cell_id] * (u.external_temp[cell_id] - u.temp[cell_id]) * surface_area
        #wall_heat_flux!(du, u, cell_id, vol)

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
        liquid_holdup = 1.0,
        gas_holdup = 0.0
    ),
    properties = reforming_area_properties,
    region_function =
    function outlet!(du, u, cell_id, vol)
        common_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 6] -= u.pipe_mass_flow[cell_id]
        #this is simulating the mass flow out of the system

        for_fields!(u.mass_fractions, du.species_masses) do species, u_mass_fractions, du_species_mass
            du_species_mass[species[cell_id]] -= u.pipe_mass_flow[cell_id] * u_mass_fractions[species[cell_id]]
        end
        #this is to prevent the concentration of all species from building up at the outlet

        du.heat[cell_id] -= u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id] 
        #this is to prevent the temperature from building up at the outlet

        #UA = overall_heat_transfer_coefficient(du, u, cell_id, vol) #W/(m^2 * K)
        #surface_area = pi * u.pipe_inside_diameter[cell_id] * u.per_cell_pipe_length[cell_id]
        #du.heat[cell_id] += u.overall_heat_transfer_coefficient[cell_id] * (u.external_temp[cell_id] - u.temp[cell_id]) * surface_area

        sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
)
    new_pressure_driven_mass_flux!(
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
    )

    #=mass_fraction_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
end

#you can check units by setting check_units = true and du0_vec and u0_vec will be returned as unitful ComponentVectors
du0_vec, u0_vec, state_axes, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

function solve_system!(du, u, p, t, geo, system)
    #sus_cell_id = 5162
    #VERY IMPORTANT: since most software uses 0-based indexing, you need to adjust the cell id by +1
    #for example, if you mouse over cell_id 5161 in paraview, you need to use 5162 in the code because julia uses 1-based indexing 

    for cell_id in eachindex(geo.cell_volumes)
        update_properties!(du, u, cell_id, geo.cell_volumes[cell_id])
    end

    #=for cell_id in 1:length(geo.cell_volumes)-1
        du.mass_face[cell_id, 6] -= u.pipe_mass_flow[cell_id]
        du.mass_face[cell_id + 1, 1] += u.pipe_mass_flow[cell_id]
    end=#

    solve_connection_groups!(du, u, geo, system)
    solve_patch_groups!(du, u, geo, system)
    solve_region_groups!(du, u, geo, system)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

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
#KrylovJL_BICGSTAB is another option, so is using algebraicmultigrid
sol.destats

prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 1e-5), p_guess)
#@time sol = solve(prob, Tsit5(), callback = approximate_time_to_finish_cb)

u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)]

sim_file = @__FILE__

sol_to_vtk(sol, u_named, grid, sim_file)