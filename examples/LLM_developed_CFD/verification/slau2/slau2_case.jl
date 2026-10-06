struct SLAU2VerificationFluid <: AbstractPhysics end

struct SLAU2Case{S, G, F}
    name::String
    system::S
    geo::G
    rhs!::F
    u0::Vector{Float64}
    p::Float64
    gamma::Float64
    cv::Float64
    boundary_condition::Symbol
end

function _slau2_conservative_values(
    density,
    vel_u,
    vel_v,
    vel_w,
    pressure,
    gamma,
)
    return (
        density,
        density * vel_u,
        density * vel_v,
        density * vel_w,
        pressure / (gamma - 1.0) +
            0.5 * density * (vel_u^2 + vel_v^2 + vel_w^2),
    )
end

function _set_slau2_cell_state!(
    state,
    cell_id,
    density,
    vel_u,
    vel_v,
    vel_w,
    pressure,
    gamma,
    species_a_mass_fraction,
    turbulent_kinetic_energy,
    specific_dissipation_rate,
)
    (
        state.density[cell_id],
        state.momentum_density_u[cell_id],
        state.momentum_density_v[cell_id],
        state.momentum_density_w[cell_id],
        state.volumetric_energy[cell_id],
    ) = _slau2_conservative_values(
        density,
        vel_u,
        vel_v,
        vel_w,
        pressure,
        gamma,
    )
    state.species_densities.species_a[cell_id] =
        density * species_a_mass_fraction
    state.species_densities.species_b[cell_id] =
        density * (1.0 - species_a_mass_fraction)
    state.turbulent_kinetic_energy_density[cell_id] =
        density * turbulent_kinetic_energy
    state.specific_dissipation_rate_density[cell_id] =
        density * specific_dissipation_rate
    return nothing
end

function _populate_slau2_profile!(state_vector, state_axes, geo, profile, gamma)
    state = ComponentVector(state_vector, state_axes)
    for cell_id in eachindex(geo.cell_centroids)
        x = geo.cell_centroids[cell_id][1]
        if profile == :uniform
            density = 1.18
            vel_u = 120.0
            vel_v = 3.0
            vel_w = -2.0
            pressure = 101325.0
            species_a_mass_fraction = 0.3
            turbulent_kinetic_energy = 0.25
            specific_dissipation_rate = 4.0
        elseif profile == :smooth
            density = 1.0 + 0.08 * sinpi(2.0 * x)
            vel_u = 45.0 + 4.0 * sinpi(2.0 * x + 0.2)
            vel_v = 2.0 + 0.3 * cospi(2.0 * x)
            vel_w = -1.0 + 0.2 * sinpi(2.0 * x)
            pressure = 100000.0 + 1500.0 * cospi(2.0 * x + 0.1)
            species_a_mass_fraction = 0.3 + 0.1 * sinpi(2.0 * x)
            turbulent_kinetic_energy = 0.25 + 0.05 * cospi(2.0 * x)
            specific_dissipation_rate = 4.0 + 0.4 * sinpi(2.0 * x)
        elseif profile == :stationary_contact
            if x < 0.5
                density = 1.0
                species_a_mass_fraction = 0.2
                turbulent_kinetic_energy = 0.2
                specific_dissipation_rate = 3.0
            else
                density = 0.125
                species_a_mass_fraction = 0.8
                turbulent_kinetic_energy = 0.4
                specific_dissipation_rate = 6.0
            end
            vel_u = 0.0
            vel_v = 0.0
            vel_w = 0.0
            pressure = 1.0
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
            turbulent_kinetic_energy = 0.2
            specific_dissipation_rate = 3.0
        else
            throw(ArgumentError("unknown SLAU2 verification profile: $profile"))
        end
        _set_slau2_cell_state!(
            state,
            cell_id,
            density,
            vel_u,
            vel_v,
            vel_w,
            pressure,
            gamma,
            species_a_mass_fraction,
            turbulent_kinetic_energy,
            specific_dissipation_rate,
        )
    end
    return state_vector
end

function _slau2_internal_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b,
)
    SLAU2!(
        du, u, p, t, system, geo,
        idx_a, face_a,
        idx_b, face_b,
        system.additional_data.face_reconstructor,
    )
    return nothing
end

function _slau2_boundary_physical_flux(
    density,
    momentum_density_u,
    momentum_density_v,
    momentum_density_w,
    volumetric_energy,
    pressure,
    normal,
)
    vel_u = momentum_density_u / density
    vel_v = momentum_density_v / density
    vel_w = momentum_density_w / density
    normal_velocity =
        vel_u * normal[1] +
        vel_v * normal[2] +
        vel_w * normal[3]
    return (
        density * normal_velocity,
        momentum_density_u * normal_velocity + pressure * normal[1],
        momentum_density_v * normal_velocity + pressure * normal[2],
        momentum_density_w * normal_velocity + pressure * normal[3],
        (volumetric_energy + pressure) * normal_velocity,
    )
end

function _slau2_boundary_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b,
)
    face_area, face_normal, face_distance, cell_volume = boundary_geometry(
        geo,
        idx_a,
        face_a,
    )

    if system.additional_data.boundary_condition == :transmissive
        (
            density_flux,
            momentum_density_u_flux,
            momentum_density_v_flux,
            momentum_density_w_flux,
            volumetric_energy_flux,
        ) = _slau2_boundary_physical_flux(
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
        add_boundary_sst_advection_flux!(
            du,
            u,
            idx_a,
            face_area,
            density_flux,
        )
    elseif system.additional_data.boundary_condition == :slip_wall
        pressure = u.pressure[idx_a]
        du.momentum_density_u_flow[idx_a] -= face_area * pressure * face_normal[1]
        du.momentum_density_v_flow[idx_a] -= face_area * pressure * face_normal[2]
        du.momentum_density_w_flow[idx_a] -= face_area * pressure * face_normal[3]
    else
        throw(ArgumentError("unsupported SLAU2 boundary condition"))
    end
    return nothing
end

function _cap_slau2_sst_advection!(du, u, geo, cell_id)
    cell_volume = cell_geometry(geo, cell_id)
    du.turbulent_kinetic_energy_density[cell_id] +=
        du.turbulent_kinetic_energy_density_flow[cell_id] / cell_volume
    du.specific_dissipation_rate_density[cell_id] +=
        du.specific_dissipation_rate_density_flow[cell_id] / cell_volume
    return nothing
end

function _slau2_connection_map(physics_a, physics_b)
    if physics_a isa SLAU2VerificationFluid &&
        physics_b isa SLAU2VerificationFluid
        return _slau2_internal_flux!
    end
    throw(ArgumentError("no SLAU2 flux is defined for the requested region pair"))
end

function _add_slau2_boundary_sets!(grid)
    addfacetset!(grid, "x_min", xyz -> abs(xyz[1]) < 1e-12)
    addfacetset!(grid, "x_max", xyz -> abs(xyz[1] - 1.0) < 1e-12)
    addfacetset!(grid, "y_min", xyz -> abs(xyz[2]) < 1e-12)
    addfacetset!(grid, "y_max", xyz -> abs(xyz[2] - 1.0) < 1e-12)
    addfacetset!(grid, "z_min", xyz -> abs(xyz[3]) < 1e-12)
    addfacetset!(grid, "z_max", xyz -> abs(xyz[3] - 1.0) < 1e-12)
    return nothing
end

function build_slau2_case(
    n_cells;
    profile = :uniform,
    boundary_condition = :transmissive,
)
    if n_cells < 2
        throw(ArgumentError("an SLAU2 verification case needs at least two cells"))
    end
    if !(boundary_condition in (:transmissive, :slip_wall))
        throw(ArgumentError("boundary_condition must be :transmissive or :slip_wall"))
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
    _add_slau2_boundary_sets!(grid)

    state_prototype = ComponentVector(
        density = zeros(n_cells)u"kg/m^3",
        momentum_density_u = zeros(n_cells)u"kg/(m^2*s)",
        momentum_density_v = zeros(n_cells)u"kg/(m^2*s)",
        momentum_density_w = zeros(n_cells)u"kg/(m^2*s)",
        volumetric_energy = zeros(n_cells)u"J/m^3",
        turbulent_kinetic_energy_density = zeros(n_cells)u"J/m^3",
        specific_dissipation_rate_density = zeros(n_cells)u"kg/(m^3*s)",
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
            turbulent_kinetic_energy = u"m^2/s^2",
            specific_dissipation_rate = u"1/s",
            turbulent_kinetic_energy_density_flow = u"W",
            specific_dissipation_rate_density_flow = u"kg/s^2",
            mass_fractions = (
                species_a = u"kg/kg",
                species_b = u"kg/kg",
            ),
            species_density_flow = (
                species_a = u"kg/s",
                species_b = u"kg/s",
            ),
        ),
        special_caches = ComponentVector(),
        second_order_syms = [],
        optimized_parameters = ComponentVector(),
    )

    initial_density = 1.0u"kg/m^3"
    initial_pressure = 1.0u"Pa"
    add_region!(
        config,
        "fluid";
        type = SLAU2VerificationFluid(),
        initial_conditions = ComponentVector(
            density = initial_density,
            momentum_density_u = 0.0u"kg/(m^2*s)",
            momentum_density_v = 0.0u"kg/(m^2*s)",
            momentum_density_w = 0.0u"kg/(m^2*s)",
            volumetric_energy = initial_pressure / (gamma - 1.0),
            turbulent_kinetic_energy_density = 0.2u"J/m^3",
            specific_dissipation_rate_density = 3.0u"kg/(m^3*s)",
            species_densities = ComponentVector(
                species_a = 0.4u"kg/m^3",
                species_b = 0.6u"kg/m^3",
            ),
        ),
        properties = ComponentVector(
            cp = cp * u"J/(kg*K)",
            cv = cv * u"J/(kg*K)",
        ),
        property_update_function = function (du, u, p, t, system, geo, cell_id)
            overall_navier_stokes_property_update!(
                du,
                u,
                p,
                t,
                system,
                geo,
                cell_id,
            )
            update_sst_primitives!(du, u, p, t, system, geo, cell_id)
        end,
        region_function = function (du, u, p, t, system, geo, cell_id)
            cap_navier_stokes_flow!(du, u, p, t, system, geo, cell_id)
            _cap_slau2_sst_advection!(du, u, geo, cell_id)
        end,
    )

    for patch_name in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")
        add_patch!(
            config,
            patch_name;
            properties = ComponentVector(),
            patch_function = _slau2_boundary_flux!,
        )
    end

    additional_data = (
        boundary_condition = boundary_condition,
        face_reconstructor = first_order_face_reconstruction!,
    )
    du0, u0, system, geo = finish_fvm_config(
        config,
        _slau2_connection_map,
        additional_data;
        check_units = false,
    )
    _populate_slau2_profile!(u0, system.state_axes, geo, profile, gamma)

    function solve_groups!(du, u, p, t, system, geo)
        update_region_groups!(du, u, p, t, system, geo)
        solve_connection_groups!(du, u, p, t, system, geo)
        solve_patch_groups!(du, u, p, t, system, geo)
        solve_region_groups!(du, u, p, t, system, geo)
        return nothing
    end

    rhs! = (du, u, p, t) -> fvm_operator!(
        du,
        u,
        p,
        t,
        system,
        geo,
        solve_groups!,
    )
    return SLAU2Case(
        "$(profile)_$(boundary_condition)_slau2_$(n_cells)_cells",
        system,
        geo,
        rhs!,
        u0,
        0.0,
        gamma,
        cv,
        boundary_condition,
    )
end

function slau2_residual(case::SLAU2Case, state = case.u0, time = 0.0)
    derivative = similar(state)
    case.rhs!(derivative, state, case.p, time)
    return derivative
end
