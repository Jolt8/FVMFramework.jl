include(joinpath(
    @__DIR__, "..", "..", "navier_stokes_testing_grounds",
    "face_reconstructors", "first_order_face_reconstruction.jl",
))
include(joinpath(
    @__DIR__, "..", "..", "navier_stokes_testing_grounds",
    "riemann_solvers", "HLLC_low_mach_correction.jl",
))
include(joinpath(
    @__DIR__, "..", "..", "navier_stokes_testing_grounds",
    "riemann_solvers", "HLLC.jl",
))
include(joinpath(
    @__DIR__, "..", "..", "navier_stokes_testing_grounds",
    "weighted_least_squares", "weighted_least_squares.jl",
))
include(joinpath(
    @__DIR__, "..", "..", "navier_stokes_testing_grounds",
    "face_reconstructors", "MUSCL_face_reconstruction.jl",
))
include(joinpath(
    @__DIR__, "..", "..", "navier_stokes_testing_grounds",
    "navier_stokes_fluid_property_update_functions", "fluid_property_update_functions.jl",
))
include(joinpath(
    @__DIR__, "..", "..", "navier_stokes_testing_grounds",
    "sum_and_cap_functions", "cap_functions.jl",
))

struct VerificationFluid <: AbstractPhysics end

struct CompressibleCase{S, G, F}
    name::String
    system::S
    geo::G
    rhs!::F
    u0::Vector{Float64}
    p::Float64
    gamma::Float64
    cv::Float64
    boundary_condition::Symbol
    spatial_method::Symbol
    low_mach_method::Symbol
    species_diffusion::Bool
end

struct CompressibleFlowVerification
    random_seed::Int
    unit_cells::Int
    integration_cells::Int
    full_cells::Int
end

function CompressibleFlowVerification(;
    random_seed = 0x46564d,
    unit_cells = 5,
    integration_cells = 32,
    full_cells = 64,
)
    return CompressibleFlowVerification(
        random_seed,
        unit_cells,
        integration_cells,
        full_cells,
    )
end

function _conservative_values(density, vel_u, vel_v, vel_w, pressure, gamma)
    momentum_density_u = density * vel_u
    momentum_density_v = density * vel_v
    momentum_density_w = density * vel_w
    volumetric_energy = pressure / (gamma - 1.0) + 0.5 * density * (
        vel_u^2 + vel_v^2 + vel_w^2
    )
    return (
        density,
        momentum_density_u,
        momentum_density_v,
        momentum_density_w,
        volumetric_energy,
    )
end

function _set_cell_state!(
    state,
    cell_id,
    density,
    vel_u,
    vel_v,
    vel_w,
    pressure,
    gamma,
    species_a_mass_fraction,
)
    (
        state.density[cell_id],
        state.momentum_density_u[cell_id],
        state.momentum_density_v[cell_id],
        state.momentum_density_w[cell_id],
        state.volumetric_energy[cell_id],
    ) = _conservative_values(density, vel_u, vel_v, vel_w, pressure, gamma)
    state.species_densities.species_a[cell_id] =
        density * species_a_mass_fraction
    state.species_densities.species_b[cell_id] =
        density * (1.0 - species_a_mass_fraction)
    return nothing
end

function _populate_profile!(state_vector, state_axes, geo, profile, gamma)
    state = ComponentVector(state_vector, state_axes)
    for cell_id in eachindex(geo.cell_centroids)
        x = geo.cell_centroids[cell_id][1]
        if profile == :uniform
            density = 1.18
            vel_u = 120.0
            vel_v = 3.0
            vel_w = -2.0
            pressure = 101325.0
            species_a_mass_fraction = 0.25
        elseif profile == :stationary_smooth
            density = 1.0 + 0.08 * sinpi(2.0 * x)
            vel_u = 0.0
            vel_v = 0.0
            vel_w = 0.0
            pressure = 101325.0
            species_a_mass_fraction = 0.35 + 0.1 * sinpi(2.0 * x)
        elseif profile == :jacobian
            density = 1.0 + 0.04 * x + 0.001 * sinpi(3.7 * x + 0.2)
            vel_u = 35.0 + 2.0 * x + 0.03 * sinpi(2.3 * x + 0.1)
            vel_v = 0.2 + 0.4 * x + 0.01 * sinpi(2.9 * x + 0.4)
            vel_w = -0.3 + 0.2 * x + 0.005 * sinpi(3.1 * x + 0.3)
            pressure = 100000.0 + 1200.0 * x + 20.0 * sinpi(2.7 * x + 0.15)
            species_a_mass_fraction = 0.3 + 0.08 * x +
                0.002 * sinpi(2.5 * x + 0.35)
        elseif profile == :stationary_contact
            if x < 0.5
                density = 1.0
            else
                density = 0.125
            end
            vel_u = 0.0
            vel_v = 0.0
            vel_w = 0.0
            pressure = 1.0
            if x < 0.5
                species_a_mass_fraction = 0.2
            else
                species_a_mass_fraction = 0.8
            end
        elseif profile == :sod
            if x < 0.5
                density = 1.0
                pressure = 1.0
            else
                density = 0.125
                pressure = 0.1
            end
            vel_u = 0.0
            vel_v = 0.0
            vel_w = 0.0
            species_a_mass_fraction = 0.4
        elseif profile == :species_advection
            density = 1.0
            vel_u = 0.5
            vel_v = 0.0
            vel_w = 0.0
            pressure = 1.0
            if x < 0.5
                species_a_mass_fraction = 0.8
            else
                species_a_mass_fraction = 0.2
            end
        else
            throw(ArgumentError("unknown compressible verification profile: $profile"))
        end
        _set_cell_state!(
            state,
            cell_id,
            density,
            vel_u,
            vel_v,
            vel_w,
            pressure,
            gamma,
            species_a_mass_fraction,
        )
    end
    return state_vector
end

function _verification_internal_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b,
)
    HLLC!(
        du, u, p, t, system, geo,
        idx_a, face_a,
        idx_b, face_b,
        system.additional_data.face_reconstructor,
        system.additional_data.low_mach_correction,
    )
    if system.additional_data.species_diffusion
        (
            distance,
            face_area_a, face_normal_a, face_distance_a, volume_a,
            face_area_b, face_normal_b, face_distance_b, volume_b,
        ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)
        add_conservative_species_diffusion_flux!(
            du,
            u,
            idx_a,
            idx_b,
            face_area_a,
            distance,
        )
    end
    return nothing
end

function _verification_boundary_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b,
)
    face_area, face_normal, face_distance, cell_volume = boundary_geometry(
        geo,
        idx_a,
        face_a,
    )
    boundary_condition = system.additional_data.boundary_condition

    if boundary_condition == :transmissive
        (
            density_flux,
            momentum_density_u_flux,
            momentum_density_v_flux,
            momentum_density_w_flux,
            volumetric_energy_flux,
        ) = physical_flux(
            u.density[idx_a],
            u.momentum_density_u[idx_a],
            u.momentum_density_v[idx_a],
            u.momentum_density_w[idx_a],
            u.volumetric_energy[idx_a],
            u.pressure[idx_a],
            face_normal,
        )

        du.density_flow[idx_a] -= face_area * density_flux
        du.momentum_density_u_flow[idx_a] -= face_area * momentum_density_u_flux
        du.momentum_density_v_flow[idx_a] -= face_area * momentum_density_v_flux
        du.momentum_density_w_flow[idx_a] -= face_area * momentum_density_w_flux
        du.volumetric_energy_flow[idx_a] -= face_area * volumetric_energy_flux
        add_boundary_species_advection_flux!(
            du,
            u,
            idx_a,
            face_area,
            density_flux,
        )
    elseif boundary_condition == :slip_wall
        pressure = u.pressure[idx_a]
        du.momentum_density_u_flow[idx_a] -= face_area * pressure * face_normal[1]
        du.momentum_density_v_flow[idx_a] -= face_area * pressure * face_normal[2]
        du.momentum_density_w_flow[idx_a] -= face_area * pressure * face_normal[3]
    else
        throw(ArgumentError("unsupported verification boundary condition: $boundary_condition"))
    end
    return nothing
end

function _verification_connection_map(physics_a, physics_b)
    if physics_a isa VerificationFluid && physics_b isa VerificationFluid
        return _verification_internal_flux!
    end
    throw(ArgumentError("no verification flux is defined for the requested region pair"))
end

function _add_boundary_sets!(grid)
    addfacetset!(grid, "x_min", xyz -> abs(xyz[1]) < 1e-12)
    addfacetset!(grid, "x_max", xyz -> abs(xyz[1] - 1.0) < 1e-12)
    addfacetset!(grid, "y_min", xyz -> abs(xyz[2]) < 1e-12)
    addfacetset!(grid, "y_max", xyz -> abs(xyz[2] - 1.0) < 1e-12)
    addfacetset!(grid, "z_min", xyz -> abs(xyz[3]) < 1e-12)
    addfacetset!(grid, "z_max", xyz -> abs(xyz[3] - 1.0) < 1e-12)
    return nothing
end

function build_compressible_case(
    n_cells;
    profile = :uniform,
    boundary_condition = :transmissive,
    spatial_method = :first_order,
    low_mach_method = :none,
    species_diffusion = true,
)
    if n_cells < 2
        throw(ArgumentError("a compressible verification case needs at least two cells"))
    end
    if !(boundary_condition in (:transmissive, :slip_wall))
        throw(ArgumentError("boundary_condition must be :transmissive or :slip_wall"))
    end
    if !(spatial_method in (:first_order, :muscl))
        throw(ArgumentError("spatial_method must be :first_order or :muscl"))
    end
    if !(low_mach_method in (:none, :thornber))
        throw(ArgumentError("low_mach_method must be :none or :thornber"))
    end

    gamma = 1.4
    cv = 717.5
    cp = gamma * cv
    grid = generate_grid(
        Hexahedron,
        (n_cells, 1, 1),
        Ferrite.Vec{3}((0.0, 0.0, 0.0)),
        Ferrite.Vec{3}((1.0, 1.0, 1.0)),
    )
    addcellset!(grid, "fluid", xyz -> true)
    _add_boundary_sets!(grid)

    state_prototype = ComponentVector(
        density = zeros(n_cells)u"kg/m^3",
        momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
        momentum_density_v = zeros(n_cells)u"kg/(m^2*s)",
        momentum_density_w = zeros(n_cells)u"kg/(m^2*s)",
        volumetric_energy = zeros(n_cells)u"J/m^3",
        species_densities = ComponentVector(
            species_a = zeros(n_cells)u"kg/m^3",
            species_b = zeros(n_cells)u"kg/m^3",
        ),
    )
    config = create_fvm_config(grid, state_prototype)
    add_setup_syms!(
        config;
        cache_syms_and_units = (
            vel_u = u"m/s",
            vel_v = u"m/s",
            vel_w = u"m/s",
            pressure = u"Pa",
            temperature = u"K",
            speed_of_sound = u"m/s",
            density_flow = u"kg/s",
            momentum_density_u_flow = u"N",
            momentum_density_v_flow = u"N",
            momentum_density_w_flow = u"N",
            volumetric_energy_flow = u"W",
            mass_fractions = (
                species_a = u"kg/kg",
                species_b = u"kg/kg",
            ),
            species_density_flow = (
                species_a = u"kg/s",
                species_b = u"kg/s",
            ),
        ),
        special_caches = ComponentVector(
            grad_density = zeros(n_cells, 3)u"kg/m^4",
            grad_pressure = zeros(n_cells, 3)u"Pa/m",
            grad_vel_u = zeros(n_cells, 3)u"1/s",
            grad_vel_v = zeros(n_cells, 3)u"1/s",
            grad_vel_w = zeros(n_cells, 3)u"1/s",
            grad_temperature = zeros(n_cells, 3)u"K/m",
        ),
        second_order_syms = [],
        optimized_parameters = ComponentVector(),
    )

    initial_density = 1.0u"kg/m^3"
    initial_pressure = 1.0u"Pa"
    initial_energy = initial_pressure / (gamma - 1.0)
    add_region!(
        config,
        "fluid";
        type = VerificationFluid(),
        initial_conditions = ComponentVector(
            density = initial_density,
            momentum_density_u = 0.0u"kg/(m^2*s)",
            momentum_density_v = 0.0u"kg/(m^2*s)",
            momentum_density_w = 0.0u"kg/(m^2*s)",
            volumetric_energy = initial_energy,
            species_densities = ComponentVector(
                species_a = 0.4 * initial_density,
                species_b = 0.6 * initial_density,
            ),
        ),
        properties = ComponentVector(
            cp = cp * u"J/(kg*K)",
            cv = cv * u"J/(kg*K)",
            diffusion_coefficients = ComponentVector(
                species_a = 1e-5u"m^2/s",
                species_b = 1e-5u"m^2/s",
            ),
        ),
        property_update_function = overall_navier_stokes_property_update!,
        region_function = cap_navier_stokes_flow!,
    )

    for patch_name in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
        add_patch!(
            config,
            patch_name;
            properties = ComponentVector(),
            patch_function = _verification_boundary_flux!,
        )
    end

    if spatial_method == :muscl
        face_reconstructor = MUSCL_face_reconstruction!
    else
        face_reconstructor = first_order_face_reconstruction!
    end
    if low_mach_method == :thornber
        low_mach_correction = thornber_low_mach_correction
    else
        low_mach_correction = no_low_mach_correction
    end
    additional_data = (
        boundary_condition = boundary_condition,
        spatial_method = spatial_method,
        low_mach_method = low_mach_method,
        face_reconstructor = face_reconstructor,
        low_mach_correction = low_mach_correction,
        species_diffusion = species_diffusion,
    )
    du0, u0, system, geo = finish_fvm_config(
        config,
        _verification_connection_map,
        additional_data;
        check_units = false,
    )
    _populate_profile!(u0, system.state_axes, geo, profile, gamma)
    weighted_least_squares_stencil = build_weighted_least_squares_stencil(geo)

    function solve_groups!(du, u, p, t, system, geo)
        update_region_groups!(du, u, p, t, system, geo)
        if spatial_method == :muscl
            update_weighted_least_squares_gradients!(u, weighted_least_squares_stencil)
            update_MUSCL_gradients!(u, weighted_least_squares_stencil)
        end
        solve_connection_groups!(du, u, p, t, system, geo)
        solve_patch_groups!(du, u, p, t, system, geo)
        solve_region_groups!(du, u, p, t, system, geo)
        return nothing
    end

    rhs! = (du, u, p, t) -> fvm_operator!(du, u, p, t, system, geo, solve_groups!)
    case_name = "$(profile)_$(boundary_condition)_$(spatial_method)_$(low_mach_method)_species_$(n_cells)_cells"
    return CompressibleCase(
        case_name,
        system,
        geo,
        rhs!,
        u0,
        0.0,
        gamma,
        cv,
        boundary_condition,
        spatial_method,
        low_mach_method,
        species_diffusion,
    )
end

function residual(case::CompressibleCase, state = case.u0, time = 0.0)
    derivative = similar(state)
    case.rhs!(derivative, state, case.p, time)
    return derivative
end

function _append_state_index_metadata!(metadata, indices, variable, species = nothing)
    if indices isa ComponentArray
        for field_name in propertynames(indices)
            nested_variable = species === nothing ? variable : species
            _append_state_index_metadata!(
                metadata,
                getproperty(indices, field_name),
                nested_variable,
                field_name,
            )
        end
        return nothing
    end

    for (cell_id, state_index) in enumerate(indices)
        metadata[state_index] = (
            variable = variable,
            species = species,
            cell = cell_id,
        )
    end
    return nothing
end

function state_index_metadata(case::CompressibleCase)
    index_state = ComponentVector(
        collect(1:length(case.u0)),
        case.system.state_axes,
    )
    metadata = Vector{NamedTuple}(undef, length(case.u0))
    for variable in propertynames(index_state)
        indices = getproperty(index_state, variable)
        _append_state_index_metadata!(metadata, indices, variable)
    end
    return metadata
end
