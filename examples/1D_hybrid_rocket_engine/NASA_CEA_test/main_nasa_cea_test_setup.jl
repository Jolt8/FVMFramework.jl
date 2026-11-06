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

Revise.includet(joinpath(@__DIR__, "CEA_retrieval", "CEA_lookup_table.jl"))

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
    mass_fractions_interps_vec = mass_fractions_interps_vec
);

add_setup_syms!(
    config;
    cache_syms_and_units = (
        #updated based on temperature, density, and composition:
        pressure = u"Pa",

        #might be updated based on temp in the future, but for now they're fixed
        #k = u"W/(m*K)",
        #mu = u"Pa*s",
        
        #updated based on temperature, density, and composition:
        cp = u"J/(kg*K)",
        cv = u"J/(kg*K)",
        
        #updated based on composition:
        mw_avg = u"g/mol",
    ),
    special_caches = ComponentVector(
        
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(
        relaxation_factor = 1e-4u"s",

        oxidizer_to_fuel_ratio = 7.0,
        chamber_pressure = 20.0u"bar",
        expansion_ratio = 4.0,
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

add_region!(
    config, "fluid_domain";
    type = Fluid(),
    initial_conditions = ComponentVector(
        density = 2.7u"kg/m^3",
        mass_fractions = (
            nitrous_oxide = 0.6,
            ethylene = 0.4,
            nitrogen = 0.0u"kg/kg",
            oxygen = 0.0u"kg/kg",
            carbon_dioxide = 0.0u"kg/kg",
            carbon_monoxide = 0.0u"kg/kg",
            water = 0.0u"kg/kg",
            hydrogen = 0.0u"kg/kg"
        ),
        internal_energy_density = 0.0u"J/m^3"
    ),
    properties = ComponentVector(
        cp = 1.0u"J/(kg*K)"
    ),
    property_update_function = 
    function update_pre_chamber_region!(du, u, p, t, system, geo, cell_id)
        #update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function pre_chamber_region!(du, u, p, t, system, geo, cell_id)
        #oxidizer_to_fuel_ratio = u.mass_fractions.nitrous_oxide[cell_id] / u.mass_fractions.ethylene[cell_id]

        interp_idx = 1
        
        for_fields!(du.mass_fractions, u.mass_fractions) do species, du_mass_fractions, u_mass_fractions
            du_mass_fractions[species[cell_id]] -= (1 / p.relaxation_factor) * 
                (u_mass_fractions[species[cell_id]] - system.mass_fractions_interps_vec[interp_idx](p.chamber_pressure, p.oxidizer_to_fuel_ratio, p.expansion_ratio))

            interp_idx += 1
        end
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