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

using SparseArrays
using JLD2

Revise.includet(joinpath(@__DIR__, "..", "multiphase_box", "overloads", "clapeyron_tracer_overloads.jl"))

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

addcellset!(grid, "fluid", xyz -> true)
getcellset(grid, "fluid")

n_cells = length(grid.cells)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end

#properties
Revise.includet(joinpath(@__DIR__, "properties", "water_properties.jl"))
water_properties = get_water_properties()

#physics
Revise.includet(joinpath(@__DIR__, "..", "multiphase_box", "physics", "drift_flux_vaporization.jl"))
Revise.includet(joinpath(@__DIR__, "..", "multiphase_box", "physics", "eos_stuff.jl")) #for update_eos_densities! and update_K_vle!
#NOTE: if the solver ever crashes, it's most likely due to the rate in which the liquid boils in the drift_flux_vaporization physics functions
#this is on lines 90-92 
#the culprit is usually: kinetic_constant = u.phase_change_mass_transfer_coefficient[cell_id] * effective_bubble_area
#if the phase_change_mass_transfer_coefficient is too high, the solver will be hit with extreme stiffness and never solve

function update_solid_properties!(du, u, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id] * u.steel_thermal_mass_multiplier[cell_id]
    u.rho[cell_id] = properties.solid_rho[cell_id]
end

#clapeyron_model = CPA(["methanol", "water"])
clapeyron_model = IAPWS95(["water"])
n_species = length(clapeyron_model.components)

if n_species == 1
    overall_drift_flux_state_update!(du, u, cell_id, vol, clapeyron_model) = single_species_overall_drift_flux_state_update!(du, u, cell_id, vol, clapeyron_model)
else
    overall_drift_flux_state_update!(du, u, cell_id, vol, clapeyron_model) = multi_species_overall_drift_flux_state_update!(du, u, cell_id, vol, clapeyron_model)
end

function update_fluid_properties!(du, u, p, t, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    overall_drift_flux_state_update!(du, u, cell_id, vol, clapeyron_model)
    
    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = u.bed_porosity[cell_id] * u.fluid_rho[cell_id] + (1.0 - u.bed_porosity[cell_id]) * properties.solid_rho[cell_id]

    compressibility_effective = u.gas_holdup[cell_id] * properties.gas_compressibility[cell_id] + u.liquid_holdup[cell_id] * properties.liquid_compressibility[cell_id]
    #isothermal_compressibility(clapeyron_model, u.pressure[cell_id], u.temp[cell_id])
    
    u.speed_of_sound[cell_id] = 1 / sqrt(u.fluid_rho[cell_id] * compressibility_effective)
end

function fluid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    #it's definitely one of these
    sum_mass_flux_face_to_cell!(du, u, cell_id, vol) #this always has to go before cap_mass_flux_to_pressure_change!

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
end

function solid_sum_and_cap_fluxes!(du, u, cell_id, vol)
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
end

u_proto = ComponentVector(
    gas_densities = (
        #methanol = zeros(n_cells)u"kg/m^3",
        water = zeros(n_cells)u"kg/m^3",
    ),
    liquid_densities = (
        #methanol = zeros(n_cells)u"kg/m^3",
        water = zeros(n_cells)u"kg/m^3",
    ),
    #I think overall mass fractions would be a cached variable then
    pressure = zeros(n_cells)u"Pa",
    temp = zeros(n_cells)u"K",
)

config = create_fvm_config(grid, u_proto);

n_faces = length(config.geo.cell_neighbor_areas[1])
species_names = keys(u_proto.gas_densities)

transducer_node_ids = [1, 10, 50, 100, 150, 200, 250, 300, 350, 400]
transducer_ids = collect(1:length(transducer_node_ids))

#Ray mapping
Revise.includet(joinpath(@__DIR__, "transducer_profiles", "ray_mapping", "ray_tracers", "default_ray_tracer.jl"))
Revise.includet(joinpath(@__DIR__, "transducer_profiles", "ray_mapping", "generate_ray_map.jl"))
#generate_ray_map(transducer_ids, transducer_node_ids, grid, config.geo) #this takes a while

ray_map_intersected_cells, ray_map_distances_through_cells, ray_map_ray_lengths = load_object(
    joinpath(
        @__DIR__,
        "transducer_profiles",
        "ray_mapping",
        "saved_ray_maps", 
        "ray_map_$(length(transducer_ids))_transducers_$(length(grid.cells))_cells.jls"
    )
)

#Projections 
#remember, although this seems pointless because projections are just for 1 transducer, 
#we use the vector between the two transducers to determine the direction of the projection
transducer_opposing_pairs_ids = [
    (1, 6),
    (2, 7),
    (3, 8),
    (4, 9),
    (5, 10)
]

#travel_time = experimental_travel_time_interp_matrix[opposing_transducer_node_id, source_transducer_node_id](most_consistent_timestamp_for_speed_of_sound)
travel_time = 100u"μs"
#speed_of_sound = 1500u"m/s"
projection_distance = 100.0u"mm"
transducer_frequency = 2.0u"MHz"
transducer_diameter = 10.0u"mm"

speed_of_sound = projection_distance / travel_time
wavelength = speed_of_sound / transducer_frequency
arg = ustrip(upreferred(1.22 * wavelength / transducer_diameter))
transducer_angle_of_projection = asind(min(arg, 1.0))
near_field_distance = ustrip(upreferred(transducer_diameter^2 / (4 * wavelength)))
transducer_radius_val = ustrip(upreferred(transducer_diameter / 2))

Revise.includet(joinpath(@__DIR__, "transducer_profiles", "transducer_projecting", "projectors", "default_projector.jl"))
Revise.includet(joinpath(@__DIR__, "transducer_profiles", "transducer_projecting", "generate_transducer_projections.jl"))
beam_profile = BeamProfile(transducer_angle_of_projection, near_field_distance, transducer_radius_val)
#generate_transducer_projections(transducer_opposing_pairs_ids, transducer_node_ids, grid, config.geo, beam_profile)

transducer_projection_intersected_cells, 
transducer_projection_volume_in_cells, 
transducer_projection_cell_distances, 
transducer_projection_center_ray_intersected_cells,
transducer_projection_center_ray_distances_through_cells,
transducer_projection_center_ray_distances, 
transducer_projection_slant_distances, 
cell_projection_counts, 
projection_unit_vectors_to_cell_centers = load_object(
    joinpath(
        @__DIR__, 
        "transducer_profiles",
        "transducer_projecting",
        "saved_transducer_projections", 
        "transducer_projection_$(length(transducer_opposing_pairs_ids))_pairs_$(length(grid.cells))_cells.jls"
    )
);

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        k = u"W/(m*K)",
        cp = u"J/(kg*K)",
        rho = u"kg/m^3",
        compressibility = u"Pa^-1",
        speed_of_sound = u"m/s",
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

        # Feature 1: Backscatter amplitude (independent constraint on void fraction)
        backscatter_intensity = u"1",

        # Feature 2: Turbulence intensity (velocity variance from xcorr peak width)
        turbulence_intensity = u"m/s",

        # Feature 3: Bubble size distribution spread (log-normal σ of ln(radius))
        bubble_radii_sigma = u"1",

        # Feature 4: Bubble number density and passage frequency
        bubble_number_density = u"m^-3",
        bubble_passage_frequency = u"Hz",
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
            #methanol = 0.1u"kg/m^3",
            water = 0.4u"kg/m^3",
        ),
        liquid_densities = ComponentVector(
            #methanol = 80.0u"kg/m^3",
            water = 200.0u"kg/m^3",
        ),
        pressure = 1.0u"atm",
        temp = 60.0u"°C",
    ),
    properties = water_properties,
    property_update_function = 
    function update_inlet!(du, u, p, t, cell_id, vol, system)
        update_fluid_properties!(du, u, p, t, cell_id, vol, system)
    end,
    region_function =
    function inlet!(du, u, p, t, cell_id, vol)
        #update_fluid_properties!(du, u, cell_id, vol, system)

        #du.heat[cell_id] += 10.0

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
    #no=slip condition at walls to prevent liquid and gas from building up at the top and bottom cells
    u.gravity[1] = 0.0 
    u.local_gas_drift_velocity[end] = 0.0

    p = ComponentVector(p_vec, system.p_axes)

    update_region_groups!(du, u, p, t, geo, system)
    solve_connection_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end

du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);

f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)

p_guess = 0.0

#test_prob = ODEProblem(f_closure, u0_vec, (0.0, 0.01), p_guess)
#@time sol = solve(test_prob, Tsit5(), callback = approximate_time_to_finish_cb)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)

ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

t0 = 0.0
tMax = 2300.0

tspan = (t0, tMax)

saveat = 10.0

implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

@time sol = solve(implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)
#@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)

#u_named = [ComponentVector(sol.u[i], state_axes) for i in 1:length(sol.u)];

#root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files"

#sol_to_vtk(sol, u_named, grid, @__FILE__, root_dir)

Revise.includet(joinpath(@__DIR__, "sonic_helper_functions.jl"))

Revise.includet(joinpath(@__DIR__, "precompute_experimental_velocities.jl"))
experimental_vel_x_interps, experimental_vel_y_interps, experimental_vel_z_interps, 
cell_vel_x_uncertainty_interps, cell_vel_y_uncertainty_interps, cell_vel_z_uncertainty_interps, 
speeds_of_sound_interps, speeds_of_sound_uncertainty_interps,
backscatter_intensity_interps, turbulence_intensity_interps, bubble_passage_frequency_interps = precompute_experimental_velocities_over_time(
    grid, geo,
    #cell_speeds_of_sound, 
    transducer_ids,
    transducer_opposing_pairs_ids,
    transducer_projection_intersected_cells,
    transducer_projection_cell_distances,
    transducer_projection_center_ray_distances,
    transducer_projection_center_ray_intersected_cells,
    transducer_projection_center_ray_distances_through_cells,
    projection_unit_vectors_to_cell_centers,
    cell_projection_counts,

    experimental_sonic_data, #TODO: feed in the actual data

    window_size_indices,
    reltol = 1.5
)

# Feature 3: Precompute spectral attenuation from transducer-to-transducer data
transducer_center_frequency_hz = ustrip(upreferred(transducer_frequency)) # 2.0 MHz → Hz
spectral_frequencies, spectral_attenuation_interps = precompute_experimental_spectral_attenuation(
    transducer_ids,
    transducer_opposing_pairs_ids,
    experimental_sonic_data,
    transducer_center_frequency_hz;
    n_frequency_bins = 10,
    bandwidth_factor = 0.25
)

# Measurement cross-section for bubble passage frequency (beam area at near-field boundary)
measurement_cross_section = pi * transducer_radius_val^2

function build_prob(f_closure, du0_vec, u0_vec, tMax, p_guess)
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
    )

    ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

    t0 = 0.0
    tspan = (t0, tMax)

    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

    return implicit_prob
end

test_interp = LinearInterpolation(
    [0.0, 1.0, 2.0, 3.0],
    [0.0, 1500.0, 4000.0, 8000.0]
)

experimental_travel_time_interp_matrix = [deepcopy(test_interp) for i in 1:length(transducer_ids), j in 1:length(transducer_ids)]

function loss(θ)
    prob, system_copy, geo_copy = get!(task_local_storage(), :implicit_prob_setup) do
        system_copy = deepcopy(system)
        geo_copy = deepcopy(geo)
        
        f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo_copy, system_copy)
        
        prob = build_prob(f_closure, du0_vec, u0_vec, experimental_tMax, θ)
        return (prob, system_copy, geo_copy)
    end

    loss_prob = remake(prob, p = θ)

    sol = solve(
        loss_prob,
        FBDF(linsolve = SparspakFactorization()), 
        sensealg = ForwardSensitivity(),
        saveat = saveat,
        callback = approximate_time_to_finish_cb
    )

    n_saves = Int(experimental_tMax / saveat)

    if length(sol.t) < n_saves + 1
        return nothing, nothing, 1e10
    end

    u_named = [ComponentVector(sol.u[i], state_axes) for i in eachindex(sol.u)]

    du_temporary = similar(sol.u[1])

    mean_squared_error = 0.0

    p_named = ComponentVector(θ, system_copy.p_axes)

    for (i, t) in enumerate(sol.t)
        du, u = unpack_fvm_state(du_temporary, sol.u[i], θ, t, system_copy)

        update_region_groups!(du, u, p_named, t, geo_copy, system_copy)

        reconstructed_speeds_of_sound = reconstruct_simulated_speeds_of_sound(
            du, u, geo_copy, system_copy,
            transducer_opposing_pairs_ids, experimental_sonic_data,
            center_ray_intersected_cells, center_ray_distances_through_cells;
            regularization_lambda = 1e-4
        )

        #getting the scaling right between these three so that one of the them doesn't overpower the other is going to be tricky
        #we also have to get the scaling right for the two simulated transducer paths
        for cell_id in 1:length(geo.cell_volumes)
            mean_squared_error += abs2(reconstructed_speeds_of_sound[cell_id] - experimental_speeds_of_sound_interps[cell_id](t)) / speeds_of_sound_uncertainty_interps[cell_id](t)
        end

        for cell_id in 1:length(geo.cell_volumes)
            mean_squared_error += abs2(u.local_gas_drift_velocity[cell_id] - experimental_vel_z_interps[cell_id](t)) / cell_vel_z_uncertainty_interps[cell_id](t)
        end

        for cell_id in 1:length(geo.cell_volumes)
            mean_squared_error += abs2(u.bubble_acceleration[cell_id] - DataInterpolations.derivative(experimental_vel_z_interps[cell_id], t)) / cell_vel_z_uncertainty_interps[cell_id](t) 
        end

        # ── Feature 4: Bubble passage frequency constraint ────────────────────
        # Constrains the relationship between gas_holdup, bubble_radii, and gas velocity.
        # passage_freq = n_density × |v_gas| × A_measurement
        for cell_id in 1:length(geo.cell_volumes)
            sim_passage_freq, _ = calculate_simulated_bubble_passage_frequency(u, cell_id, measurement_cross_section)
            exp_passage_freq = bubble_passage_frequency_interps[cell_id](t)
            if exp_passage_freq > 0.0
                mean_squared_error += abs2(sim_passage_freq - exp_passage_freq) / exp_passage_freq
            end
        end

        # ── Feature 1: Backscatter amplitude constraint ────────────────────────
        # Independent per-cell constraint on void fraction from self-to-self echo path.
        # Does not go through the lossy tomographic inversion.
        for transducer_id in transducer_ids
            ray_cells, sim_backscatter = calculate_simulated_backscatter_amplitudes(
                u, transducer_id,
                transducer_projection_center_ray_intersected_cells,
                transducer_projection_center_ray_distances_through_cells
            )
            for (idx, cell_id) in enumerate(ray_cells)
                exp_backscatter = backscatter_intensity_interps[cell_id](t)
                if exp_backscatter > 0.0
                    mean_squared_error += abs2(sim_backscatter[idx] - exp_backscatter) / exp_backscatter
                end
            end
        end

        for origin_transducer_idx in transducer_ids
            for destination_transducer_idx in transducer_ids
                if origin_transducer_idx != destination_transducer_idx
                    #=
                    this part of the loss function will force the optimizer to get the density, bulk modulus, and gas_holdup/void_fraction right
                    bulk modulus is probably going to come from an EOS model though, so that only density and gas_holdup are variables that must be satisfied through other conditions
                    however, we also have access to experimental speeds of sound from our precompute_experimental_velocities_over_time function, so we can still constrain the solver by making it satisfy 
                    that condition as well
                    don't know if it's going to be better to constrain it by simulating the sound waves and comparing it to the experimental travel times, or compare each cell's speed of sound to 
                    the approximate_experimental_speeds_of_sound found in precompute_experimental_velocities_over_time

                    if we do the precompute_experimental_velocities_over_time, I wonder if it would be a good idea to plug the simulated travel times into the same function that finds 
                    the approximate_experimental_speeds_of_sound and see if that matches the experimental one or if it would be better to compare the speeds of sound f each cell direclty
                    maybe then we could ignore cells that have not been observed by a transducer instead of filling them in
                    this is going to be challenging to figure out the best solution
                    actually, would it make sense to do this by also simulating sound waves in the simulation to find simulated_velocities and then compare that to experimental_velocities
                        - notes on this:

                        - I think this would actually be a good idea as it would impose the same lossiness and limitations of the velocity reconstruction
                        - However, the only issue is that running it through this pipeline might be too performance intensive and thus comparing the raw simulated values 
                        to the experimentally-derived values would be better
                        - There's no way to predict this, so we'll just see how slow precompute_experimental_velocities_over_time is
                        - TODO: check how slow precompute_experimental_velocities_over_time is
                        - Wait a minute, we can't do this because precompute_experimental_velocities_over_time relies on actual experimental data
                        of the frequency response from the transducers and we can't simulate that
                        - I guess we could maybe built an uncertainty quantification map for each cell's experimental velocity 
                            - for example, this would be the map that we would probably go for
                                uncertainty = 0.1^n (cell has been observed by n tansducers)
                                uncertainty = 0.1^2 (cell has been observed by 2 transducer)
                                uncertianty = 1 (cell has been observed by only 1 transducer)
                                uncertainty = 10 (cell has had no transducers observe it and requied a neighbor fill-in)
                                NOTE: the uncertainty would also go up if a cell was only observed by two transducers that were roughly parallel to each other,
                                beacause then the measured velocity of the cell can "notch" both transducers 
                                for example, if the difference  between the angle formed by the ray projected by transducer 1 to the cell 
                                and the ray projected by transducer 2 to the cell is small, the uncertainty would be greater
                                - I'm sure there's a smart way to do this, I'll look into it later
                                - Smart way: You could use 1 / cond(A[1:n, :]) or the smallest singular value as your uncertainty metric — 
                                it's cheap to compute and directly measures how well-constrained the velocity reconstruction is for that cell.
                                 
                            - We could also use a greater weighting for simulated cell speeds of sound that are observed by two transducers that are closer together
                            since if the transducers are closer together, there are less cells along the path between them, meaning that that speed of sound is more accurate for 
                            that specific subset of cells and a cell within that should be weighted more by that speed of sound calculation than when it's 
                            calculated from transducers that are further away from each other
                            
                        - Apparently this is an "Inverse Crime" as I would normally be comparing raw simulated cell values against reconstructed experimental values

                        - NOTE: we can still use precompute_experimental_speeds_of_sound though because it only relies on travel times
                        - The only uncertainty would be in finding velocities


                                
                    =#
                    experimental_travel_time = experimental_travel_time_interp_matrix[origin_transducer_idx, destination_transducer_idx](t)
                    sim_travel_time = calculate_simulated_travel_time(
                        u, origin_transducer_idx, destination_transducer_idx, 
                        ray_map_intersected_cells, ray_map_distances_through_cells
                    )
                    mean_squared_error += abs2(sim_travel_time - experimental_travel_time)
                    
                    ping_volume = experimental_ping_volumes_interp_matrix[origin_transducer_idx, destination_transducer_idx](t)
                    experimental_returned_volume = experimental_return_volumes_interp_matrix[origin_transducer_idx, destination_transducer_idx](t)
                    sim_returned_volume = calculate_simulated_received_returned_volume(
                        u, t, origin_transducer_idx, destination_transducer_idx, 
                        ray_map_intersected_cells, ray_map_distances_through_cells, ping_volume
                    )
                    mean_squared_error += abs2(sim_returned_volume - experimental_returned_volume)

                    # ── Feature 3: Spectral attenuation constraint ────────────────
                    # Compares the frequency-dependent attenuation along the ray path
                    # to constrain the bubble size distribution (log-normal: median + σ).
                    # Both sides are normalized by the center frequency, so the
                    # transmitted spectrum shape cancels exactly.
                    if length(spectral_frequencies) > 0
                        sim_total_atten = calculate_simulated_total_attenuation_per_frequency(
                            u, origin_transducer_idx, destination_transducer_idx,
                            ray_map_intersected_cells, ray_map_distances_through_cells,
                            spectral_frequencies
                        )
                        center_f_idx = argmin(abs.(spectral_frequencies .- transducer_center_frequency_hz))
                        sim_center_atten = sim_total_atten[center_f_idx]

                        for f_idx in eachindex(spectral_frequencies)
                            sim_normalized = -(sim_total_atten[f_idx] - sim_center_atten)
                            exp_normalized = spectral_attenuation_interps[origin_transducer_idx, destination_transducer_idx][f_idx](t)
                            mean_squared_error += abs2(sim_normalized - exp_normalized)
                        end
                    end
                end
            end
        end
    end

    return sol.t, u_named, mean_squared_error
end

#=
u_additional_information = ComponentVector(
    speeds_of_sound_uncertainties = speeds_of_sound_uncertainty_interps[i](t),
    experimental_speeds_of_sound = experimental_speeds_of_sound_interps[i](t),
    experimental_vel_vec = [
        SVector{3, Float64}(
            experimental_vel_x_interps[i](t)[cell_id], 
            experimental_vel_y_interps[i](t)[cell_id], 
            experimental_vel_z_interps[i](t)[cell_id]
        )
        for cell_id in 1:n_cells
    ],
    experimental_vel_vec_uncertainties = [
        SVector{3, Float64}(
            cell_vel_x_uncertainty_interps[i](t)[cell_id], 
            cell_vel_y_uncertainty_interps[i](t)[cell_id], 
            cell_vel_z_uncertainty_interps[i](t)[cell_id]
        )
        for cell_id in 1:n_cells
    ]
)
    =#

du_complete, u_complete = regenerate_fvm_state(sol, system, solve_system!, geo, p_guess; u_additional_information = ComponentVector());

root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files"

sol_to_vtk(sol, du_complete, u_complete, grid, geo, @__FILE__, root_dir; include_zeros_fields = false) 

loss(1.0)