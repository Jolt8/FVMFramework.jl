using Revise
using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using NonlinearSolve
using Sparspak
using Dates
using DataInterpolations
using Clapeyron

Revise.includet(joinpath(@__DIR__, "overloads", "clapeyron_tracer_overloads.jl"))

#=
using XLSX
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using ForwardDiff
using DataFrames
using CSV
=#

using FVMFramework

grid_dimensions = (1, 1, 100)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 20.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)
#=
addcellset!(grid, "wall", xyz -> 
    xyz[1] <= 1.0 &&
    xyz[1] >= 99.0 &&
    xyz[2] <= 1.0 &&
    xyz[2] >= 99.0 &&
    xyz[3] <= 1.0 &&
    xyz[3] >= 99.0
)
getcellset(grid, "wall")

addcellset!(grid, "fluid", xyz -> 
    xyz[1] < 99.0 &&
    xyz[1] > 1.0 &&
    xyz[2] < 99.0 &&
    xyz[2] > 1.0 &&
    xyz[3] < 99.0 &&
    xyz[3] > 1.0
)
getcellset(grid, "fluid")
=#

addcellset!(grid, "fluid", xyz -> true)
getcellset(grid, "fluid")

n_cells = length(grid.cells)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

#properties
Revise.includet(joinpath(@__DIR__, "properties", "water_methanol_properties.jl"))
water_methanol_properties = get_water_methanol_properties()

#physics
Revise.includet(joinpath(@__DIR__, "physics", "drift_flux_vaporization.jl"))
Revise.includet(joinpath(@__DIR__, "physics", "eos_stuff.jl")) #for update_eos_densities! and update_K_vle!
#NOTE: if the solver ever crashes, it's most likely due to the rate in which the liquid boils in the drift_flux_vaporization physics functions
#this is on lines 90-92 
#the culprit is usually: kinetic_constant = u.phase_change_mass_transfer_coefficient[cell_id] * effective_bubble_area
#if the phase_change_mass_transfer_coefficient is too high, the solver will be hit with extreme stiffness and never solve

function update_solid_properties!(du, u, p, t, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.solid_rho[cell_id]
end

#clapeyron_model = CPA(["methanol", "water"])
clapeyron_model = PR(["methanol", "water"])

Clapeyron.volume(clapeyron_model, 100000u"Pa", 273.13u"K", [0.5, 0.5], phase = :liquid)
Clapeyron.volume(clapeyron_model, 100000u"Pa", 273.13u"K", [1.0, 1.0], phase = :liquid)

Clapeyron.isothermal_compressibility(clapeyron_model, 100000u"Pa", 273.13u"K", [0.5, 0.5], phase = :liquid)
Clapeyron.isothermal_compressibility(clapeyron_model, 100000u"Pa", 273.13u"K", [0.5, 0.5], phase = :gas)
isobaric_heat_capacity(clapeyron_model, 100000u"Pa", 273.13u"K", [0.5, 0.5], phase = :gas)
isobaric_heat_capacity(clapeyron_model, 100000u"Pa", 273.13u"K", [0.5, 0.5])

function update_fluid_properties!(du, u, p, t, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    multi_species_overall_drift_flux_state_update!(du, u, cell_id, vol, clapeyron_model)

    gas_mole_fractions_vec, liquid_mole_fractions_vec = get_mole_fractions_vec(du, u, cell_id, vol)
    
    u.k[cell_id] = properties.k[cell_id]

    u.rho[cell_id] = u.bed_porosity[cell_id] * u.fluid_rho[cell_id] + (1.0 - u.bed_porosity[cell_id]) * u.solid_rho[cell_id]

    #=
    gas_cp_mass = isobaric_heat_capacity(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], gas_mole_fractions_vec, phase = :gas) / u.gas_mw_avg[cell_id]
    liquid_cp_mass = isobaric_heat_capacity(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], liquid_mole_fractions_vec, phase = :liquid) / u.liquid_mw_avg[cell_id]
    
    fluid_rho_cp = u.gas_holdup[cell_id] * u.eos_gas_density[cell_id] * gas_cp_mass + u.liquid_holdup[cell_id] * u.eos_liquid_density[cell_id] * liquid_cp_mass

    rho_cp_total = u.bed_porosity[cell_id] * fluid_rho_cp + (1.0 - u.bed_porosity[cell_id]) * u.solid_rho[cell_id] * u.solid_cp[cell_id]
    u.cp[cell_id] = rho_cp_total / u.rho[cell_id]
    =#
    u.cp[cell_id] = u.solid_cp[cell_id]

    gas_compressibility = isothermal_compressibility(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], gas_mole_fractions_vec, phase = :gas)
    liquid_compressibility = isothermal_compressibility(clapeyron_model, u.pressure[cell_id], u.temp[cell_id], liquid_mole_fractions_vec, phase = :liquid)
    u.compressibility_effective[cell_id] = gas_compressibility * u.gas_holdup[cell_id] + liquid_compressibility * u.liquid_holdup[cell_id]
end

function fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    #it's definitely one of these
    sum_mass_flux_face_to_cell!(du, u, cell_id, vol) #this always has to go before cap_mass_flux_to_pressure_change!

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change_with_compressibility!(du, u, cell_id, vol) 
    #cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    #compared to cap_mass_flux_to_pressure_change! that uses the ideal gas law, the version that uses 
    #compressibility makes the solver take 3x longer to solve
end

function solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
end

u_proto = ComponentVector(
    gas_densities = (
        methanol = zeros(n_cells)u"kg/m^3",
        water = zeros(n_cells)u"kg/m^3",
    ),
    liquid_densities = (
        methanol = zeros(n_cells)u"kg/m^3",
        water = zeros(n_cells)u"kg/m^3",
    ),
    #I think overall mass fractions would be a cached variable then
    pressure = zeros(n_cells)u"Pa",
    temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto);

n_faces = length(config.geo.cell_neighbor_areas[1])
#reaction_names = keys(water_methanol_properties.reactions.reforming_reactions)
species_names = keys(u_proto.gas_densities)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        k = u"W/(m*K)",
        cp = u"J/(kg*K)",
        rho = u"kg/m^3",
        mass = u"kg",
        gas_holdup = u"1",
        liquid_holdup = u"1",
        gas_enthalpy = u"J/kg",
        liquid_enthalpy = u"J/kg",
        overall_gas_density = u"kg/m^3",
        overall_liquid_density = u"kg/m^3",
        fluid_rho = u"kg/m^3",

        gas_mw_avg = u"kg/mol",
        liquid_mw_avg = u"kg/mol",

        eos_gas_density = u"kg/m^3", #these are the intrinsic density of the current phase found with the eos model
        eos_liquid_density = u"kg/m^3",
        
        mass_face = u"kg",
        
        viscous_drag = u"kg/(m*s)",
        inertial_drag = u"kg/(m^3)",
        driving_force = u"N/m^3",
        mixture_superficial_velocity = u"m/s",
        mixture_pore_velocity = u"m/s",

        gas_velocity_face = u"m/s",
        liquid_velocity_face = u"m/s",

        gas_mass_fractions = u"1",
        liquid_mass_fractions = u"1",

        gas_mole_fractions = u"1",
        liquid_mole_fractions = u"1",

        gas_generation = u"kg/(m^3*s)",

        K_vle = u"1",

        compressibility_effective = u"1/Pa"
    ),
    special_caches = ComponentArray(
        mass_face = zeros(n_cells, n_faces)u"kg",

        viscous_drag = zeros(n_cells, n_faces)u"kg/(m*s)",
        inertial_drag = zeros(n_cells, n_faces)u"kg/(m^3)",
        driving_force = zeros(n_cells, n_faces)u"N/m^3",
        mixture_superficial_velocity = zeros(n_cells, n_faces)u"m/s",
        mixture_pore_velocity = zeros(n_cells, n_faces)u"m/s",

        gas_velocity_face = zeros(n_cells, n_faces)u"m/s",
        liquid_velocity_face = zeros(n_cells, n_faces)u"m/s",

        gas_mass_fractions = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),
        liquid_mass_fractions = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),

        gas_mole_fractions = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),
        liquid_mole_fractions = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),

        gas_generation = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"kg/(m^3*s)" for _ in 1:length(species_names))
        ),

        K_vle = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"1" for _ in 1:length(species_names))
        ),
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector()
)

add_region!(
    config, "fluid";
    type = Fluid(),
    initial_conditions = ComponentVector(
        gas_densities = ComponentVector(
            methanol = 0.1u"kg/m^3",
            water = 0.2u"kg/m^3",
        ),
        liquid_densities = ComponentVector(
            methanol = 80.0u"kg/m^3",
            water = 100.0u"kg/m^3",
        ),
        pressure = 1.0u"atm",
        temp = 60.0u"°C",
    ),
    properties = water_methanol_properties,
    property_update_function = 
    function update_inlet!(du, u, p, t, cell_id, vol, system)
        update_fluid_properties!(du, u, p, t, cell_id, vol, system)
    end,
    region_function =
    function inlet!(du, u, p, t, cell_id, vol)
        #update_fluid_properties!(du, u, cell_id, vol, system)

        if t >= 30000.0
            du.heat[cell_id] += 1000.0
        end

        fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    end
)
#Connection functions
function fluid_fluid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    overall_drift_flux_model!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )

    #this is the only equation that's still needed
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function fluid_solid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function solid_solid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Solid && return fluid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Fluid && return fluid_solid_flux!

    typeof(phys_a) <: Solid && typeof(phys_b) <: Solid && return solid_solid_flux!
end

function solve_system!(du, u, p_vec, t, geo, system)
    p = ComponentVector(p_vec, system.p_axes)

    update_region_groups!(du, u, p, t, geo, system)

    solve_connection_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

p_guess = [0.0]

@time f_closure(du0_vec, u0_vec, p_guess, 0.0)

#test_prob = ODEProblem(f_closure, u0_vec, (0.0, 0.01), p_guess)
#@time sol = solve(test_prob, Tsit5(), callback = approximate_time_to_finish_cb)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 30000.0
tspan = (t0, tMax)

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

#@time sol = solve(implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)
@time sol = solve(implicit_prob, FBDF(), callback = approximate_time_to_finish_cb)
#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)

du_complete, u_complete = regenerate_fvm_state(sol, system, solve_system!, geo, p_guess, u_additional_information = ComponentVector());

root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files"

sol_to_vtk(sol, du_complete, u_complete, grid, geo, @__FILE__, root_dir; include_zeros_fields = false)