function get_igniter_total_volume(geo, grid)
    igniter_volume = 0.0
    
    for igniter_cell_id in getcellset(grid, "igniter_puck_volume")
        igniter_volume += geo.cell_volumes[igniter_cell_id]
    end

    return igniter_volume
end

function get_igniter_total_surface_area(geo, grid)
    igniter_total_surface_area = 0.0
    igniter_cell_ids = getcellset(grid, "igniter_puck_volume")
    
    for face in getfacetset(grid, "fluid_domain_to_igniter_puck_interface")
        cell_id = face[1]
        face_idx = face[2]

        if cell_id in igniter_cell_ids
            igniter_total_surface_area += geo.cell_face_areas[cell_id, face_idx]
        end
    end

    return igniter_total_surface_area
end

function setup_igniter_properties!(fluid_initial_conditions, fluid_properties, config, grid, input_igniter_properties)
    igniter_puck_initial_conditions = deepcopy(fluid_initial_conditions)
    igniter_puck_properties = deepcopy(fluid_properties)

    initial_gas_density = fluid_initial_conditions.density
    initial_igniter_bulk_density = ustrip(upreferred(input_igniter_properties.density))
    density_scaling_ratio = initial_igniter_bulk_density / initial_gas_density

    igniter_puck_initial_conditions.density = initial_igniter_bulk_density

    igniter_puck_initial_conditions.momentum_density_u = 0.0
    igniter_puck_initial_conditions.momentum_density_v = 0.0
    igniter_puck_initial_conditions.momentum_density_w = 0.0

    # Preserve the copied gas temperature while increasing the initial bulk density
    # to represent the unburned solid puck.
    igniter_puck_initial_conditions.volumetric_energy *= density_scaling_ratio

    # The current chemistry model treats the igniter products as a nitrogen-like
    # surrogate gas and adds the measured heat of combustion separately at emission.
    igniter_puck_initial_conditions.species_densities .= 0.0
    igniter_puck_initial_conditions.species_densities.nitrogen = initial_igniter_bulk_density

    if hasproperty(igniter_puck_initial_conditions, :turbulent_kinetic_energy_density)
        igniter_puck_initial_conditions.turbulent_kinetic_energy_density *= density_scaling_ratio
        igniter_puck_initial_conditions.specific_dissipation_rate_density *= density_scaling_ratio
    end

    igniter_volume = get_igniter_total_volume(config.geo, grid) * u"m^3"
    igniter_total_surface_area = get_igniter_total_surface_area(config.geo, grid) * u"m^2"

    initial_igniter_bulk_mass = (initial_igniter_bulk_density * u"kg/m^3") * igniter_volume
    burnout_gas_mass = (initial_gas_density * u"kg/m^3") * igniter_volume
    initial_burnable_mass = initial_igniter_bulk_mass - burnout_gas_mass

    config.u_proto.igniter_mass[1] = initial_burnable_mass

    igniter_system_properties = ComponentVector(
        volume = igniter_volume,
        initial_bulk_density = initial_igniter_bulk_density * u"kg/m^3",
        burnout_gas_density = initial_gas_density * u"kg/m^3",
        initial_burnable_mass = initial_burnable_mass,
        total_surface_area = igniter_total_surface_area,
    )

    igniter_system_properties = merge_properties(igniter_system_properties, input_igniter_properties)
    igniter_system_properties = ustrip.(upreferred.(igniter_system_properties))

    return igniter_puck_initial_conditions, igniter_puck_properties, igniter_system_properties
end

function igniter_mass_burn_rate_at_t(du, u, p, t, system, geo)
    initial_burnable_mass = system.igniter.initial_burnable_mass
    burn_time = system.igniter.burn_time

    if initial_burnable_mass <= zero(initial_burnable_mass) || burn_time <= zero(burn_time) ||
       t < zero(t) || t >= burn_time ||
       u.igniter_mass[1] <= zero(u.igniter_mass[1])
        return 0.0
    end

    nominal_mass_burn_rate = initial_burnable_mass / burn_time
    return nominal_mass_burn_rate
end

function update_igniter_combustion!(du, u, p, t, system, geo)
    du.igniter_mass[1] = -igniter_mass_burn_rate_at_t(du, u, p, t, system, geo)
end

function update_igniter_puck_combustion_source!(du, u, p, t, system, geo, igniter_cell_id)
    if igniter_burnt(du, u, p, t, system, geo)
        return nothing
    end

    igniter_bulk_density = max(u.density[igniter_cell_id], 1e-10)
    igniter_specific_sensible_energy = u.volumetric_energy[igniter_cell_id] / igniter_bulk_density
    igniter_bulk_density_burn_rate = du.igniter_mass[1] / system.igniter.volume

    # The burnable material is distributed uniformly through the puck volume.
    # Removing the same density rate from every puck cell makes the integrated
    # puck mass loss equal the scalar igniter mass loss.
    du.density[igniter_cell_id] += igniter_bulk_density_burn_rate
    du.species_densities.nitrogen[igniter_cell_id] += igniter_bulk_density_burn_rate
    du.volumetric_energy[igniter_cell_id] += igniter_bulk_density_burn_rate * igniter_specific_sensible_energy

    if system.compressible_solver_options.use_sst_turbulence == true
        igniter_specific_turbulent_kinetic_energy = u.turbulent_kinetic_energy_density[igniter_cell_id] / igniter_bulk_density
        igniter_specific_dissipation_rate = u.specific_dissipation_rate_density[igniter_cell_id] / igniter_bulk_density

        du.turbulent_kinetic_energy_density[igniter_cell_id] += igniter_bulk_density_burn_rate * igniter_specific_turbulent_kinetic_energy
        du.specific_dissipation_rate_density[igniter_cell_id] += igniter_bulk_density_burn_rate * igniter_specific_dissipation_rate
    end

    return nothing
end

# idx_a is the surrounding fluid cell and idx_b is the igniter puck cell.
function fluid_to_igniter_puck_interface_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a, idx_b, face_b
)
    (
        dist,
        face_area_a, face_normal_a, face_distance_a, vol_a,
        face_area_b, face_normal_b, face_distance_b, vol_b
    ) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)

    total_igniter_product_mass_flow_rate = max(-du.igniter_mass[1], zero(du.igniter_mass[1]))
    face_area_fraction = face_area_b / system.igniter.total_surface_area
    igniter_product_mass_flow_rate = total_igniter_product_mass_flow_rate * face_area_fraction
    igniter_product_mass_flux = igniter_product_mass_flow_rate / face_area_b

    # Assume the products enter at the adjacent fluid density. The normal velocity
    # then follows directly from mass flux = density * normal velocity.
    injected_gas_density = max(u.density[idx_a], 1e-10)
    injected_gas_normal_velocity = igniter_product_mass_flux / injected_gas_density
    injected_gas_velocity_u = face_normal_b[1] * injected_gas_normal_velocity
    injected_gas_velocity_v = face_normal_b[2] * injected_gas_normal_velocity
    injected_gas_velocity_w = face_normal_b[3] * injected_gas_normal_velocity

    igniter_bulk_density = max(u.density[idx_b], 1e-10)
    igniter_specific_sensible_energy = u.volumetric_energy[idx_b] / igniter_bulk_density
    injected_gas_specific_kinetic_energy = 0.5 * (injected_gas_velocity_u^2 + injected_gas_velocity_v^2 + injected_gas_velocity_w^2)
    injected_gas_specific_total_energy = igniter_specific_sensible_energy + system.igniter.heat_of_combustion + injected_gas_specific_kinetic_energy

    du.density[idx_a] += igniter_product_mass_flow_rate / vol_a
    du.momentum_density_u[idx_a] += igniter_product_mass_flow_rate * injected_gas_velocity_u / vol_a
    du.momentum_density_v[idx_a] += igniter_product_mass_flow_rate * injected_gas_velocity_v / vol_a
    du.momentum_density_w[idx_a] += igniter_product_mass_flow_rate * injected_gas_velocity_w / vol_a
    du.species_densities.nitrogen[idx_a] += igniter_product_mass_flow_rate / vol_a
    du.volumetric_energy[idx_a] += igniter_product_mass_flow_rate * injected_gas_specific_total_energy / vol_a

    return nothing
end

function igniter_burnt(du, u, p, t, system, geo)
    burnout_mass_tolerance = 1e-9 * system.igniter.initial_burnable_mass
    return t >= system.igniter.burn_time || u.igniter_mass[1] <= burnout_mass_tolerance
end

function igniter_puck_internal_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    if igniter_burnt(du, u, p, t, system, geo)
        fluid_fluid_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    else
        solid_solid_flux!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b)
    end

    return nothing
end
