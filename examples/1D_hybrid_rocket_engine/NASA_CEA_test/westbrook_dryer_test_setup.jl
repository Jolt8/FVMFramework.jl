using Revise

using FVMFramework

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
using Clapeyron
using PreallocationTools
using DataInterpolations
using JSON

grid_dimensions = (1, 1, 1)

left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0)) 

grid = generate_grid(Hexahedron, grid_dimensions, left, right)
addcellset!(grid, "fluid_domain", x -> true)
getcellset(grid, "fluid_domain")

n_cells = length(grid.cells)

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3",
    mass_fractions = (
        nitrous_oxide = zeros(n_cells)u"kg/kg",
        ethylene = zeros(n_cells)u"kg/kg",
        nitrogen = zeros(n_cells)u"kg/kg",
        oxygen = zeros(n_cells)u"kg/kg",
        carbon_dioxide = zeros(n_cells)u"kg/kg",
        carbon_monoxide = zeros(n_cells)u"kg/kg",
        water = zeros(n_cells)u"kg/kg",
        hydrogen = zeros(n_cells)u"kg/kg"
    ),
    internal_energy_density = zeros(n_cells)u"J/m^3"
)

config = create_fvm_config(grid, u_proto);

Revise.includet(joinpath(@__DIR__, "CEA_retrieval", "CEA_lookup_table.jl"));

config.data = (
    chamber_cp_interp = chamber_cp_interp,
    chamber_mu_interp = chamber_mu_interp,
    chamber_k_interp = chamber_k_interp,
    chamber_prandtl_interp = chamber_prandtl_interp,
    chamber_isp_vac_interp = chamber_isp_vac_interp,
    chamber_cstar_interp =  chamber_cstar_interp,
    chamber_T_comb_interp = chamber_T_comb_interp,
    chamber_mw_avg_interp = chamber_mw_avg_interp,
    chamber_gamma_interp = chamber_gamma_interp,
    mass_fractions_interps_vec = mass_fractions_interps_vec,
    use_smooth_reaction_ramp = false
);

add_setup_syms!(
    config;
    cache_syms_and_units = (
        temperature = u"K",
        #updated based on temperature, density, and composition:
        pressure = u"Pa",

        #might be updated based on temp in the future, but for now they're fixed
        #k = u"W/(m*K)",
        #mu = u"Pa*s",
        
        #updated based on temperature, density, and composition:
        cp = u"J/(kg*K)",
        cv = u"J/(kg*K)",

        molar_concentrations = (
            nitrous_oxide = u"mol/m^3",
            ethylene = u"mol/m^3",
            nitrogen = u"mol/m^3",
            oxygen = u"mol/m^3",
            carbon_dioxide = u"mol/m^3",
            carbon_monoxide = u"mol/m^3",
            water = u"mol/m^3",
            hydrogen = u"mol/m^3"
        ),
        
        #updated based on composition:
        mw_avg = u"g/mol",
    ),
    special_caches = ComponentVector(
        
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(
        # Nitrous Oxide Decomposition kinetics:
        nitrous_oxide_decomposition_pre_exponential_factor = 4.2e9u"s^-1",
        nitrous_oxide_decomposition_activation_energy = 222_000u"J/mol",
        nitrous_oxide_decomposition_heat_of_reaction = -82_000u"J/mol",

        # Westbrook kinetics:
        westbrook_dryer_step_1_pre_exponential_factor = 7.5895e7u"m^2.25 * mol^-0.75 * s^-1",
        westbrook_dryer_step_1_activation_energy = 125_520u"J/mol",
        westbrook_dryer_step_1_heat_of_reaction = -757_120u"J/mol",

        westbrook_dryer_step_2_forward_pre_exponential_factor = 1.2589e10u"m^2.25 * mol^-0.75 * s^-1",
        westbrook_dryer_step_2_forward_activation_energy = 167_360u"J/mol",
        westbrook_dryer_step_2_forward_heat_of_reaction = -282_980u"J/mol",

        westbrook_dryer_step_2_reverse_pre_exponential_factor = 5.0000e8u"s^-1",
        westbrook_dryer_step_2_reverse_activation_energy = 167_360u"J/mol",
        westbrook_dryer_step_2_reverse_heat_of_reaction = 282_980u"J/mol",

        oxidizer_to_fuel_ratio = 7.0,
        chamber_pressure = 20.0u"bar",
        expansion_ratio = 4.0,

        species_min_mass_fraction_ramp_activation = 1e-6, #typically 1e-8
        species_max_mass_fraction_ramp_max = 1e-4 #typically 1e-6
    ),
)

struct Fluid <: AbstractPhysics end

function no_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    #nothing here since we handle these BCs at the patches
end

function connection_map_function(phys_a, phys_b)
    return no_flux!
end

u_proto = (
    density = zeros(n_cells)u"kg/m^3",
    mass_fractions = (
        nitrous_oxide = zeros(n_cells)u"kg/kg",
        ethylene = zeros(n_cells)u"kg/kg",
        nitrogen = zeros(n_cells)u"kg/kg",
        oxygen = zeros(n_cells)u"kg/kg",
        carbon_dioxide = zeros(n_cells)u"kg/kg",
        carbon_monoxide = zeros(n_cells)u"kg/kg",
        water = zeros(n_cells)u"kg/kg",
        hydrogen = zeros(n_cells)u"kg/kg"
    ),
    internal_energy_density = zeros(n_cells)u"J/m^3"
)

Revise.includet(joinpath(@__DIR__, "fuel_oxidizer_reaction_physics", "combustion_reaction_model", "simple_arrhenius_reaction_model.jl"))
Revise.includet(joinpath(@__DIR__, "fuel_oxidizer_reaction_physics", "combustion_reaction_model", "westbrook_dryer_reaction_mechanism.jl"))
Revise.includet(joinpath(@__DIR__, "fuel_oxidizer_reaction_physics", "fuel_oxidizer_reaction_physics.jl"))

function molar_concentrations!(du, u, p, t, system, geo, cell_id)
    for_fields!(u.mass_fractions, u.molecular_weights, u.molar_concentrations) do species, u_mass_fractions, u_molecular_weights, u_molar_concentrations
        u_molar_concentrations[species[cell_id]] = (u.density[cell_id] * u_mass_fractions[species[cell_id]]) / u_molecular_weights[species[cell_id]]
    end
end

function smooth_ramp(value, min_value, max_value)
    return clamp((value - min_value) / (max_value - min_value), 0.0, 1.0)
end

function overall_combustion_model!(du, u, p, t, system, geo, cell_id)
    decompose_nitrous_oxide!(du, u, p, t, system, geo, cell_id)
    combust_ethylene_and_oxygen_westbrook_dryer!(du, u, p, t, system, geo, cell_id) 
    #combust_ethylene_and_oxygen_simple!(du, u, p, t, system, geo, cell_id) #is also an option
end

function normalize!(vec)
    total = sum(vec)
    vec ./= total
end

add_region!(
    config, "fluid_domain";
    type = Fluid(),
    initial_conditions = ComponentVector(
        density = 2.7u"kg/m^3",
        mass_fractions = normalize!(ComponentVector(
            nitrous_oxide = 1.0u"kg/kg" * config.optimized_parameters.oxidizer_to_fuel_ratio,
            ethylene = 1.0u"kg/kg",
            nitrogen = 0.0u"kg/kg",
            oxygen = 0.0u"kg/kg",
            carbon_dioxide = 0.0u"kg/kg",
            carbon_monoxide = 0.0u"kg/kg",
            water = 0.0u"kg/kg",
            hydrogen = 0.0u"kg/kg"
        )),
        internal_energy_density = 0.0u"J/m^3"
    ),
    properties = ComponentVector(
        R_gas = 8.314u"J/(mol*K)",
        molecular_weights = (
            nitrous_oxide = 44.01u"g/mol", 
            ethylene = 28.05u"g/mol", 
            nitrogen = 28.01u"g/mol", 
            oxygen = 32.00u"g/mol", 
            carbon_dioxide = 44.01u"g/mol",
            carbon_monoxide = 28.01u"g/mol",
            water = 18.02u"g/mol",
            hydrogen = 2.02u"g/mol",
        ),
    ),
    property_update_function = 
    function update_pre_chamber_region!(du, u, p, t, system, geo, cell_id)
        #u.cp[cell_id] = system.chamber_cp_interp(p.chamber_pressure, p.oxidizer_to_fuel_ratio, p.expansion_ratio)
        #u.cv[cell_id] = u.cp[cell_id] / system.chamber_gamma_interp(p.chamber_pressure, p.oxidizer_to_fuel_ratio, p.expansion_ratio)
        #u.temperature[cell_id] = 
        u.temperature[cell_id] = 1300.0
        molar_concentrations!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function pre_chamber_region!(du, u, p, t, system, geo, cell_id)
        overall_combustion_model!(du, u, p, t, system, geo, cell_id)
    end
)

function solve_system!(du, u, p, t, system, geo)
    update_region_groups!(du, u, p, t, system, geo)
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end

system_data = (;
    config.data...,
);

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, system_data, check_units = false);

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

p_guess = ustrip.(upreferred.(deepcopy(config.optimized_parameters)))

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

t0 = 0.0
tMax = 10.0

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
#ode_func = ODEFunction(f_closure_implicit)

westbrook_dryer_prob = ODEProblem(ode_func, u0_vec, (t0, tMax), p_guess)

@time sol = solve(
    westbrook_dryer_prob,
    FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff())
    #callback = approximate_time_to_finish_cb,
)

@time sol = solve(
    westbrook_dryer_prob,
    #Tsit5(),
    #FBDF(linsolve = SparSpakFactorization())
    #FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff())
    FBDF(linsolve = KLUFactorization()),
    #callback = approximate_time_to_finish_cb,
)

@time sol = solve(
    westbrook_dryer_prob,
    Tsit5(),
    callback = approximate_time_to_finish_cb,
)

desired_shown_steps = 250
save_fraction = desired_shown_steps / length(sol.t)

reduced_westbrook_dryer_sol = (
    u = [[sol.u[i] for i in 1:(Int(ceil(1/save_fraction))):length(sol.u)]..., sol.u[end]],
    t = [sol.t[1:(Int(ceil(1/save_fraction))):length(sol.t)]..., sol.t[end]]
)

du_named, u_named = regenerate_fvm_state(reduced_westbrook_dryer_sol, system, solve_system!, config.geo, p_guess, save_properties_and_caches = true, track_progress = true);

using Plots

plot(reduced_westbrook_dryer_sol.t, [u_named[i].internal_energy_density[1] for i in eachindex(reduced_westbrook_dryer_sol.t)])

plot(reduced_westbrook_dryer_sol.t, [u_named[i].mass_fractions.nitrous_oxide[1] for i in eachindex(reduced_westbrook_dryer_sol.t)], label = "nitrous_oxide")
plot!(reduced_westbrook_dryer_sol.t, [u_named[i].mass_fractions.ethylene[1] for i in eachindex(reduced_westbrook_dryer_sol.t)], label = "ethylene")
plot!(reduced_westbrook_dryer_sol.t, [u_named[i].mass_fractions.nitrogen[1] for i in eachindex(reduced_westbrook_dryer_sol.t)], label = "nitrogen")
plot!(reduced_westbrook_dryer_sol.t, [u_named[i].mass_fractions.oxygen[1] for i in eachindex(reduced_westbrook_dryer_sol.t)], label = "oxygen")
plot!(reduced_westbrook_dryer_sol.t, [u_named[i].mass_fractions.carbon_dioxide[1] for i in eachindex(reduced_westbrook_dryer_sol.t)], label = "carbon_dioxide")
plot!(reduced_westbrook_dryer_sol.t, [u_named[i].mass_fractions.carbon_monoxide[1] for i in eachindex(reduced_westbrook_dryer_sol.t)], label = "carbon_monoxide")
plot!(reduced_westbrook_dryer_sol.t, [u_named[i].mass_fractions.water[1] for i in eachindex(reduced_westbrook_dryer_sol.t)], label = "water")
plot!(reduced_westbrook_dryer_sol.t, [u_named[i].mass_fractions.hydrogen[1] for i in eachindex(reduced_westbrook_dryer_sol.t)], label = "hydrogen")

