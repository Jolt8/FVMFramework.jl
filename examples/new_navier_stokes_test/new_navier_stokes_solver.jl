using Revise
using Logging
using Unitful
using OrdinaryDiffEq
using NonlinearSolve
using Sparspak
using Ferrite
using SparseConnectivityTracer
using ComponentArrays
import ADTypes
using ILUZero
using LinearAlgebra
using StaticArrays
using FVMFramework


grid_dimensions = (100, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

n_cells = length(grid.cells)
n_faces = nfacets(grid.cells[1])

# Grid sets
addcellset!(grid, "fluid", xyz -> xyz[1] >= (1.0 / n_cells) && xyz[1] <= (1.0 - 1.0 / n_cells))
getcellset(grid, "fluid")

addcellset!(grid, "supersonic_inlet", xyz -> xyz[1] <= (1.0 / n_cells))
getcellset(grid, "supersonic_inlet")

addcellset!(grid, "supersonic_outlet", xyz -> xyz[1] >= (1.0 - 1.0 / n_cells))
getcellset(grid, "supersonic_outlet")

u_proto = ComponentVector(
    density = zeros(n_cells)u"kg/m^3",
    momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_v = zeros(n_cells)u"kg/(m^2*s)",
    momentum_density_w = zeros(n_cells)u"kg/(m^2*s)",
    volumetric_energy = zeros(n_cells)u"J/m^3"
)

config = create_fvm_config(grid, u_proto);

add_setup_syms!(
    config;
    cache_syms_and_units = (
        vel_u = u"m/s",
        vel_v = u"m/s",
        vel_w = u"m/s",
        pressure = u"Pa",
        temperature = u"K",
        speed_of_sound = u"m/s",

        #flows for later integrating fluxes
        density_flow = u"kg/s",
        momentum_density_u_flow = u"N",
        momentum_density_v_flow = u"N",
        momentum_density_w_flow = u"N",
        volumetric_energy_flow = u"W",
        
        #specific_internal_energy = u"J/kg", #not cached right now, just a property
    ),
    special_caches = ComponentVector(
    ),
    second_order_syms = [],
    optimized_parameters = ComponentVector(),
)


struct Fluid <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "face_reconstructors/first_order_face_reconstruction.jl"))
Revise.includet(joinpath(@__DIR__, "riemann_solvers/HLLC.jl"))

HLLC_closure! = (
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals, cell_face_distances,
    cell_neighbor_normals, cell_neighbor_distances, 
    cell_volumes
) -> HLLC!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals, cell_face_distances,
    cell_neighbor_normals, cell_neighbor_distances, 
    cell_volumes,
    first_order_face_reconstruction!
)

function fluid_fluid_flux!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    cell_face_areas, cell_face_normals, cell_face_distances,
    cell_neighbor_normals, cell_neighbor_distances, 
    cell_volumes
)
    HLLC_closure!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_face_areas[idx_a][face_idx], cell_face_normals[idx_a][face_idx], cell_face_distances[idx_a][face_idx],
        cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a],
    )
end

# Mapping connection functions
function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
end

Revise.includet(joinpath(@__DIR__, "navier_stokes_fluid_property_update_functions/fluid_property_update_functions.jl"))
Revise.includet(joinpath(@__DIR__, "sum_and_cap_functions/cap_functions.jl"))

function construct_initial_conditions_from_intuitive_inputs(u)
    momentum_density_u = u.density * u.vel_u
    momentum_density_v = u.density * u.vel_v
    momentum_density_w = u.density * u.vel_w

    internal_energy = u.cv * u.temperature
    kinetic_energy = 0.5 * (u.vel_u^2 + u.vel_v^2 + u.vel_w^2)
    
    volumetric_energy = u.density * (internal_energy + kinetic_energy)

    pressure = (u.R_gas / u.mw) * u.density * u.temperature

    @show speed_of_sound = sqrt((u.cp / u.cv) * pressure / u.density) |> u"m/s"

    return ComponentVector(
        density = u.density,
        momentum_density_u = momentum_density_u,
        momentum_density_v = momentum_density_v,
        momentum_density_w = momentum_density_w,
        volumetric_energy = volumetric_energy
    ), ComponentVector(
        cp = u.cp,
        cv = u.cv,
        R_gas = u.R_gas,
        mw = u.mw
    )
end

inital_conditions, properties = construct_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 0.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = 0.0u"m/s",
        density = 1.18u"kg/m^3",
        temperature = 300.0u"K",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol"
    )
)

# Add region and patches
add_region!(
    config, "fluid";
    type = Fluid(),
    initial_conditions = initial_conditions,
    properties = properties,
    property_update_function = 
    function update_fluid_properties!(du, u, p, t, cell_id, vol)
        overall_navier_stokes_property_update!(du, u, p, t, cell_id, vol)
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, cell_id, vol)
        cap_navier_stokes_flow!(du, u, p, t, cell_id, vol)
    end,
)

inital_conditions, properties = construct_initial_conditions_from_intuitive_inputs(
    ComponentVector(
        vel_u = 600.0u"m/s",
        vel_v = 0.0u"m/s",
        vel_w = 0.0u"m/s",
        density = 1.18u"kg/m^3",
        temperature = 300.0u"K",
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",
        R_gas = 8.314u"J/(mol*K)",
        mw = 28.97u"g/mol"
    )
)

add_region!(
    config, "supersonic_inlet";
    type = Fluid(),
    initial_conditions = initial_conditions,
    properties = properties,
    property_update_function = 
    function update_fluid_properties!(du, u, p, t, cell_id, vol)
        overall_navier_stokes_property_update!(du, u, p, t, cell_id, vol)
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, cell_id, vol)
        cap_navier_stokes_flow!(du, u, p, t, cell_id, vol)
    end,
)

add_region!(
    config, "supersonic_outlet";
    type = Fluid(),
    initial_conditions = initial_conditions,
    properties = properties,
    property_update_function = 
    function update_fluid_properties!(du, u, p, t, cell_id, vol)
        overall_navier_stokes_property_update!(du, u, p, t, cell_id, vol)
    end,
    region_function = 
    function fluid_physics!(du, u, p, t, cell_id, vol)
        du.density *= 0.0
        du.
        cap_navier_stokes_flow!(du, u, p, t, cell_id, vol)
    end,
)

function slip_wall_flux!(
    du, u, p, t,
    idx_a,
    area,
    normal
)
    gamma = u.cp[idx_a] / u.cv[idx_a]

    _, _, _, pressure, _ = primitive_from_conservative(
        u.density[idx_a],
        u.momentum_density_u[idx_a],
        u.momentum_density_v[idx_a],
        u.momentum_density_w[idx_a],
        u.volumetric_energy[idx_a],
        gamma,
    )

    du.momentum_density_u_flow[idx_a] -= area * pressure * normal[1]
    du.momentum_density_v_flow[idx_a] -= area * pressure * normal[2]
    du.momentum_density_w_flow[idx_a] -= area * pressure * normal[3]

    # mass flux = 0
    # energy flux = 0
end

# Finish FVM configuration
du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

# System solver function
function solve_system!(du, u, p, t, geo, system)
    solve_connection_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end

f_closure_implicit = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)
p_guess = 0.0
detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

# 1. Direct steady-state solve using NonlinearSolve.jl
f_closure_steady = (du, u, p) -> f_closure_implicit(du, u, p, 0.0)

# Detect Jacobian sparsity for NonlinearProblem
nl_jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_steady(du, u, p_guess), du0_vec, u0_vec, detector
) 
#this scales absolutely abysmally as the number of cells goes up
#for example, a 10x increase in the amount of cells made this take around 100x longer! 

t0 = 0.0
tMax = 100000.0
tspan = (t0, tMax)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(nl_jac_sparsity))
implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

println("Solving the Navier-Stokes ODE system...")
@time sol = solve(
    implicit_prob,
    FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true),
    callback = approximate_time_to_finish_cb,
)
