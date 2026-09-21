using Revise
using Logging
using Unitful
using OrdinaryDiffEq
using NonlinearSolve
using Sparspak
using Ferrite
using SparseConnectivityTracer
using ForwardDiff
using ComponentArrays
import ADTypes
using ILUZero
using LinearAlgebra
using StaticArrays
using FVMFramework
using Clapeyron

pre_chamber_length = ustrip(5.0u"cm" |> u"m")
grain_length = ustrip(30.0u"cm" |> u"m")
post_combustion_chamber_length = ustrip(5.0u"cm" |> u"m")
nozzle_length = ustrip(10.0u"cm" |> u"m")

grid_x_length = ustrip(pre_chamber_length + grain_length + post_combustion_chamber_length + nozzle_length)
grid_y_length = 0.4
grid_z_length = 0.1

grid_dimensions = (100, 4, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((grid_x_length, grid_y_length, grid_z_length))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

n_cells_axial = grid_dimensions[1]
n_layers = grid_dimensions[2]
n_cells = length(grid.cells)
n_faces = nfacets(grid.cells[1])

# Grid sets
delta_x = grid_x_length / n_cells_axial
delta_y = grid_y_length / n_layers

# Axial boundary definitions
x_0 = 0.0
x_1 = 1 * delta_x
x_2 = 2 * delta_x
x_3 = pre_chamber_length
x_4 = pre_chamber_length + grain_length
x_5 = pre_chamber_length + grain_length + post_combustion_chamber_length
x_6 = grid_x_length

# Radial boundary definitions
y_0 = 0.0
y_1 = 1 * delta_y
y_2 = 2 * delta_y
y_3 = 3 * delta_y
y_4 = 4 * delta_y

"""
Potential Axial Layers:
    - oxidizer_tank (represents N2O tank) (no radial layers)
    - Post valve (represents valve after N2O tank) (no radial layers)
    - Pre chamber (surrounded by ablative layer, phenolic liner, and structural liner)
    - Fuel grain void (empty space in grain port) (surrounded by fuel grain, phenolic liner, and structural liner)
    - Fuel grain (this would be another cell set 1 cell out radially than the fuel grain void) (surrounded by phenolic liner and structural liner)
    - Post combustion chamber (right before nozzle) (surrounded by ablative layer, phenolic liner, and structural liner)
    - Nozzle (surrounded by a half or full thickness phenolic liner and structural liner)

Potential Radial Layers
    - Ablative layer (This would be next to the pre chamber and post combustion champer)
    - Nozzle material
    - Phenolic liner (this would encompass all other radial layers)
        - This will also encompass the nozzle material
        - This isn't exactly accurate though because usually the nozzle has a step cut into it that's about half the thickness of the phenolic liner
    - Structural liner (this would encompass the phenolic liner)
        - Not sure if an aluminum structural liner would actually be required, it looks like some hybrids ditch this and just have the carbon composite right next to the phenolic liner

Max amount of radial layers is 4
"""

# Helper to avoid floating point issues with boundaries
tol = 1e-6
in_bounds(v, lower, upper) = (v >= lower - tol) && (v <= upper + tol)

# ================================
# Axial Layer 1 (Inner Core Flow)
# ================================
addcellset!(grid, "oxidizer_tank", xyz -> 
    in_bounds(xyz[1], x_0, x_1) && in_bounds(xyz[2], y_0, y_1)
)

addcellset!(grid, "post_valve", xyz -> 
    in_bounds(xyz[1], x_1, x_2) && in_bounds(xyz[2], y_0, y_1)
)

addcellset!(grid, "pre_chamber", xyz -> 
    in_bounds(xyz[1], x_2, x_3) && in_bounds(xyz[2], y_0, y_1)
)

addcellset!(grid, "fuel_grain_void", xyz -> 
    in_bounds(xyz[1], x_3, x_4) && in_bounds(xyz[2], y_0, y_1)
)

addcellset!(grid, "post_combustion_chamber", xyz -> 
    in_bounds(xyz[1], x_4, x_5) && in_bounds(xyz[2], y_0, y_1)
)

addcellset!(grid, "nozzle", xyz -> 
    in_bounds(xyz[1], x_5, x_6) && in_bounds(xyz[2], y_0, y_1)
)

# ================================
# Axial Layer 2 (Inner Walls)
# ================================
addcellset!(grid, "ablative_layer", xyz -> 
    (in_bounds(xyz[1], x_2, x_3) || in_bounds(xyz[1], x_4, x_5)) && 
    in_bounds(xyz[2], y_1, y_2)
)

addcellset!(grid, "fuel_grain", xyz -> 
    in_bounds(xyz[1], x_3, x_4) && in_bounds(xyz[2], y_1, y_2)
)

addcellset!(grid, "nozzle_material", xyz -> 
    in_bounds(xyz[1], x_5, x_6) && in_bounds(xyz[2], y_1, y_2)
)

# ================================
# Axial Layer 3 (Phenolic Liner)
# ================================
addcellset!(grid, "phenolic_liner", xyz -> 
    in_bounds(xyz[1], x_2, x_6) && in_bounds(xyz[2], y_2, y_3)
)

# ================================
# Axial Layer 4 (Structural Liner)
# ================================
addcellset!(grid, "structural_liner", xyz -> 
    in_bounds(xyz[1], x_2, x_6) && in_bounds(xyz[2], y_3, y_4)
)

# ================================
# Facet Sets
# ================================
addfacetset!(grid, "adjustable_valve", xyz ->
    abs(xyz[1] - x_1) < tol && in_bounds(xyz[2], y_0, y_1)
)

addfacetset!(grid, "injector_plate", xyz ->
    abs(xyz[1] - x_2) < tol && in_bounds(xyz[2], y_0, y_1)
)

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3",
    momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
    volumetric_energy = zeros(n_cells)u"J/m^3",
    mass_fractions = (
        air = zeros(n_cells)u"kg/kg",
        nitrous_oxide = zeros(n_cells)u"kg/kg",
        hdpe_vapor = zeros(n_cells)u"kg/kg",
        nitrogen = zeros(n_cells)u"kg/kg",
        oxygen = zeros(n_cells)u"kg/kg",
        carbon_monoxide = zeros(n_cells)u"kg/kg",
        carbon_dioxide = zeros(n_cells)u"kg/kg",
        hydrogen = zeros(n_cells)u"kg/kg",
        water = zeros(n_cells)u"kg/kg",
    )
)

model = PR(
    [
        "air", 
        "nitrous oxide", 
        "hdpe vapor", 
        "nitrogen", 
        "oxygen", 
        "carbon monoxide", 
        "carbon dioxide", 
        "hydrogen", 
        "water",
    ], 
    userlocations = joinpath(@__DIR__, "physics/clapeyron_additional_data.csv")
)
# Note: 'hdpe vapor' uses ethylene properties as a surrogate gas for HDPE vapor

config = create_fvm_config(grid, u_proto);

add_setup_syms!(
    config;
    cache_syms_and_units = (
        vel_u = u"m/s",
        pressure = u"Pa",
        temperature = u"K",
        speed_of_sound = u"m/s",

        #flows for later integrating fluxes
        density_flow = u"kg/s",
        momentum_density_u_flow = u"N",
        volumetric_energy_flow = u"W",
        k = u"W/(m*K)",
        
        #specific_internal_energy = u"J/kg", #not cached right now, just a property
    ),
    special_caches = ComponentVector(
        grad_density = zeros(n_cells, 3)u"kg/m^4",
        grad_vel_u = zeros(n_cells, 3)u"m/(s*m)",
        grad_temperature = zeros(n_cells, 3)u"K/m",

        
        vel_u_face = zeros(n_cells, n_faces)u"m/s",
        temperature_face = zeros(n_cells, n_faces)u"J/m^3",
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(),
)

struct Fluid <: AbstractPhysics end
#=
Revise.includet(joinpath(@__DIR__, "face_reconstructors/first_order_face_reconstruction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC_low_mach_correction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC.jl"))
Revise.includet(joinpath(@__DIR__, "weighted_least_squares/weighted_least_squares.jl"))
Revise.includet(joinpath(@__DIR__, "face_reconstructors/MUSCL_face_reconstruction.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/fluid_viscous_and_diffusive_fluxes.jl"))
Revise.includet(joinpath(@__DIR__, "viscous_and_diffusive_terms/wall_viscous_and_diffusive_fluxes.jl"))
=#

function fluid_fluid_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    HLLC!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b,
        #MUSCL_face_reconstruction!,
        first_order_face_reconstruction!,
        #thornber_low_mach_correction,
    )
end

nitrous_oxide_density = mass_density(
    model, ustrip(upreferred(tank_initial_pressure)), ustrip(upreferred(tank_initial_temperature)), Vector(ComponentVector(
        air = 0.0,
        nitrous_oxide = 1.0,
        hdpe_vapor = 0.0,
        nitrogen = 0.0,
        oxygen = 0.0,
        carbon_monoxide = 0.0,
        carbon_dioxide = 0.0,
        hydrogen = 0.0,
        water = 0.0,
    ))
) * u"kg/m^3"

oxidizer_volumetric_flow = oxidizer_mass_flow / nitrous_oxide_density 

Revise.includet(joinpath(@__DIR__, "physics/valve_flow_coefficient.jl"))

min_valve_openness = 0.1
#flow_coefficient = valve_flow_coefficient(min_valve_openness)
max_pressure_drop = 10u"bar"
min_flow_cofficient = oxidizer_volumetric_flow * sqrt((0.5 * tank_initial_pressure) / max_pressure_drop)

#pressure_drop = 0.5 * tank_initial_pressure * (oxidizer_volumetric_flow / flow_coefficient)^2 |> u"bar"

function tank_post_valve_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    oxidizer_volumetric_flow = ustrip(oxidizer_mass_flow |> u"kg/s") / u.rho[idx_a]

    valve_openness = system.experimental_valve_opening_percentage(t) #get the valve's openness based on experimental data through time
    flow_coefficient = system.adjustable_valve_flow_coefficient(valve_openness) #use the valve's openness to calculate its flow coefficient

    pressure_drop = 0.5 * u.pressure[idx_a] * (oxidizer_volumetric_flow / flow_coefficient)^2

    u.pressure[idx_b] = u.pressure[idx_a] - pressure_drop

    flash_result = ph_flash(model, u.pressure[idx_b], u.temperature[idx_a], system.mole_fractions_vec[idx_a]).data

    u.temperature[idx_b] = flash_result.T
end

discharge_coefficient = 0.7

orifice_diameter = 0.5u"cm"
valve_cross_sectional_area = pi * orifice_diameter^2 / 4

function post_valve_pre_chamber_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, 
    idx_b, face_b,
)
    oxidizer_volumetric_flow = ustrip(oxidizer_mass_flow |> u"kg/s") / u.rho[idx_a]

    pressure_drop = 0.5 * u.pressure[idx_a] * (oxidizer_volumetric_flow / (discharge_coefficient * ustrip(upreferred(valve_cross_sectional_area))))^2
end




#=
# Mapping connection functions
#IMPORTANT: when we switch to evaluating only one flux update per face rather than two, this will determine which idx_a will be and idx_b will be
#This makes sense because if we have two fluids that have the same type (flux functions should be the same) they will both be updated randomly
function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    typeof(phys_a) <: Tank && typeof(phys_b) <: PostValve && return tank_post_valve_flux!
end

# Add region and patches
add_region!(
    config, "oxidizer_tank";
    type = Tank(),
    initial_conditions = (
        density = 0.0u"kg/m^3",
        momentum_density_u = 0.0u"kg/(m^2*s)",
        volumetric_energy = 0.0u"J/m^3",
        mass_fractions = (
            air = 0.0u"kg/kg",
            nitrous_oxide = 1.0u"kg/kg",
            hdpe_vapor = 0.0u"kg/kg",
            nitrogen = 0.0u"kg/kg",
            oxygen = 0.0u"kg/kg",
            carbon_monoxide = 0.0u"kg/kg",
            carbon_dioxide = 0.0u"kg/kg",
            hydrogen = 0.0u"kg/kg",
            water = 0.0u"kg/kg",
        )
    ),
    properties = ComponentVector(
        nitrous_oxide_mass = 0.5u"kg",
        tank_volume = 0.1u"m^3",
        temperature = 21.0u"°C"
    ),
    property_update_function = 
    function update_oxidizer_tank_properties!(du, u, p, t, system, geo, cell_id)
        mw_avg!(du, u, p, t, system, geo, cell_id)
        
        u.rho[cell_id] = p.nitrous_oxide_mass / p.tank_volume

        molar_density = u.rho[cell_id] / u.mw_avg[cell_id]

        molar_volume = 1 / molar_density

        u.pressure[cell_id] = pressure(system.model, molar_volume, u.temperature[cell_id], system.mole_fractions_vec[cell_id])
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, system, geo, cell_id)
        
    end,
)

add_region!(
    config, "post_valve";
    type = PostValve(),
    initial_conditions = (
        density = 0.0u"kg/m^3",
        momentum_density_u = 0.0u"kg/(m^2*s)",
        volumetric_energy = 0.0u"J/m^3",
        mass_fractions = (
            air = 0.0u"kg/kg",
            nitrous_oxide = 1.0u"kg/kg",
            hdpe_vapor = 0.0u"kg/kg",
            nitrogen = 0.0u"kg/kg",
            oxygen = 0.0u"kg/kg",
            carbon_monoxide = 0.0u"kg/kg",
            carbon_dioxide = 0.0u"kg/kg",
            hydrogen = 0.0u"kg/kg",
            water = 0.0u"kg/kg",
        )
    ),
    properties = ComponentVector(
        nitrous_oxide_mass = 0.5u"kg",
        tank_volume = 0.1u"m^3",
        temperature = 21.0u"°C"
    ),
    property_update_function = 
    function update_oxidizer_tank_properties!(du, u, p, t, system, geo, cell_id)
        mw_avg!(du, u, p, t, system, geo, cell_id)
        
        u.rho[cell_id] = p.nitrous_oxide_mass / p.tank_volume

        molar_density = u.rho[cell_id] / u.mw_avg[cell_id]

        molar_volume = 1 / molar_density

        u.pressure[cell_id] = pressure(system.model, molar_volume, u.temperature[cell_id], system.mole_fractions_vec[cell_id])
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, system, geo, cell_id)
        
    end,
)

#I think we're going to add an additional_data field of the system so that arbitrary data can be passed into any function without allocations
additional_data = (
    mole_fractions = [zeros(length(u_proto.mass_fractions[1])) for _ in 1:n_cells],
)

du0_vec, u0_vec, system, geo = finish_fvm_config(config, connection_map_function, additional_data, check_units = false);
WLS_STENCIL = build_weighted_least_squares_stencil(geo)

# System solver function
function solve_system!(du, u, p, t, system, geo)
    populate_system_mole_fractions_vec!(du, u, p, t, system, geo)
    update_region_groups!(du, u, p, t, system, geo)
    update_weighted_least_squares_gradients!(u, WLS_STENCIL)
    update_MUSCL_gradients!(u, WLS_STENCIL)
    solve_connection_groups!(du, u, p, t, system, geo)
    solve_patch_groups!(du, u, p, t, system, geo)
    solve_region_groups!(du, u, p, t, system, geo)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_system!)

#=
function f_closure_implicit(du, u, p, t)
    fvm_operator!(du, u, p, t, system, geo, solve_system!)
end
=#

p_guess = 0.0

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

# Detect Jacobian sparsity for NonlinearProblem
jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
) 
#this scales absolutely abysmally as the number of cells goes up
#for example, a 10x increase in the amount of cells made this take around 100x longer! 


function state_is_invalid(u, p, t, system)
    u_named = ComponentVector(u, system.state_axes)

    for i in 1:n_cells
        kinetic_energy_density = 0.5 * (
            u_named.momentum_density_u[i]^2 +
            u_named.momentum_density_v[i]^2 +
            u_named.momentum_density_w[i]^2
        ) / u_named.density[i]
        
        internal_energy_density = u_named.volumetric_energy[i] - kinetic_energy_density

        if !all(
            isfinite, (
                u_named.density[i],
                u_named.momentum_density_u[i],
                u_named.momentum_density_v[i],
                u_named.momentum_density_w[i],
                u_named.volumetric_energy[i],
                internal_energy_density,
            )) || u_named.density[i] <= 0.0 || internal_energy_density <= 0.0
            return true
        end
    end

    return false
end

state_is_invalid_closure = (u, p, t) -> state_is_invalid(u, p, t, system);

#transient
t0 = 0.0
tMax = 1000.0
tspan = (t0, tMax)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))
implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

early_dtmax = 100.0
late_dtmax  = 1000.0
switch_time = 400.0

function increase_dtmax!(integrator)
    integrator.opts.dtmax = late_dtmax

    # We changed solver settings, not the state vector.
    u_modified!(integrator, false)

    println("Raised dtmax to $late_dtmax at t = $(integrator.t)")
end

increase_dtmax_cb = DiscreteCallback(
    (u, t, integrator) -> t >= switch_time && integrator.opts.dtmax < late_dtmax,
    increase_dtmax!;
    save_positions = (false, false),
)

callbacks = CallbackSet(
    approximate_time_to_finish_cb,
    #increase_dtmax_cb,
)

@time sol = solve(
    implicit_prob,
    #Tsit5(),
    #AutoTsit5(FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true)),
    #FBDF(linsolve = SparspakFactorization()),
    #FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true),
    #FBDF(linsolve = KrylovJL_GMRES(), nlsolve = NLNewton(relax = 0.7), precs = iluzero, concrete_jac = true),
    callback = callbacks,
    #isoutofdomain = state_is_invalid_closure,
    #saveat = (tMax / 300),
    #dtmax = 100
    #dtmax = early_dtmax
    #dt = 1e-5
)

sol.alg

f_closure_steady = (du, u, p) -> f_closure_implicit(du, u, p, 0.0)

nl_func = NonlinearFunction(f_closure_steady, jac_prototype = float.(jac_sparsity))

prob = NonlinearProblem(nl_func, u0_vec, p_guess)

@time sol_steady = solve(prob, NonlinearSolve.NewtonRaphson(concrete_jac = true))

du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, geo, p_guess, track_progress = false);

root_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\julia_cfd_output_files"

sol_to_vtk(sol, du_named, u_named, grid, geo, @__FILE__, root_dir, track_progress = false)
