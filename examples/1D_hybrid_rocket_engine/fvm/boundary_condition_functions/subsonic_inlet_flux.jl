"""
    subsonic_inlet_flux!(
        du, u, p, t, system, geo,
        idx_a, face_a, idx_b, face_b;
        inlet_conditions,
    )

Apply the one-way, subsonic mass-flow boundary returned by
`update_injector_train!`.

`inlet_conditions.mass_flow` is the total injector mass flow and
`system.injector_train_properties.injector_orifice_area` is the corresponding
physical flow area. The boundary-face areas may partition that orifice area,
so every face receives the same mass flux per unit area. A positive injector
mass flow is directed opposite the outward face normal.

The energy flux uses the static specific internal energy supplied by the feed
model plus kinetic energy and flow work. Injected species mass is pure nitrous
oxide. Until the feed model supplies inlet turbulence values, SST scalars are
advected from the adjacent fluid cell.
"""
function subsonic_inlet_flux!(
    du, u, p, t, system, geo,
    idx_a, face_a,
    idx_b, face_b;
    inlet_conditions,
)
    face_area, face_normal, _, _ = boundary_geometry(geo, idx_a, face_a)

    density = inlet_conditions.mid_section_density
    pressure = inlet_conditions.mid_section_pressure
    specific_internal_energy =
        inlet_conditions.mid_section_specific_internal_energy
    mass_flow = inlet_conditions.mass_flow
    injector_velocity = inlet_conditions.injector_velocity
    injector_area = inlet_conditions.injector_area

    density > zero(density) || throw(DomainError(
        density,
        "the injector inlet density must be positive",
    ))
    pressure > zero(pressure) || throw(DomainError(
        pressure,
        "the injector inlet pressure must be positive",
    ))
    specific_internal_energy > zero(specific_internal_energy) ||
        throw(DomainError(
            specific_internal_energy,
            "the injector specific internal energy must be positive",
        ))
    mass_flow >= zero(mass_flow) || throw(DomainError(
        mass_flow,
        "the injector inlet is one-way and requires nonnegative mass flow",
    ))
    injector_area > zero(injector_area) || throw(DomainError(
        injector_area,
        "the injector orifice area must be positive",
    ))

    # Fluxes are positive in the outward-normal direction. The prescribed
    # injector mass flow is positive into the CFD domain, hence the minus sign.
    density_flux = -mass_flow / injector_area

    vel_u = -injector_velocity * face_normal[1]
    vel_v = -injector_velocity * face_normal[2]
    vel_w = -injector_velocity * face_normal[3]

    momentum_density_u_flux = density_flux * vel_u + pressure * face_normal[1]
    momentum_density_v_flux = density_flux * vel_v + pressure * face_normal[2]
    momentum_density_w_flux = density_flux * vel_w + pressure * face_normal[3]

    specific_total_enthalpy =
        specific_internal_energy +
        0.5 * injector_velocity^2 +
        pressure / density
    volumetric_energy_flux = density_flux * specific_total_enthalpy

    du.density_flow[idx_a] -= face_area * density_flux
    du.momentum_density_u_flow[idx_a] -=
        face_area * momentum_density_u_flux
    du.momentum_density_v_flow[idx_a] -=
        face_area * momentum_density_v_flux
    du.momentum_density_w_flow[idx_a] -=
        face_area * momentum_density_w_flux
    du.volumetric_energy_flow[idx_a] -= face_area * volumetric_energy_flux

    if hasproperty(du, :species_density_flow)
        du.species_density_flow.nitrous_oxide[idx_a] -= face_area * density_flux
    end

    if system.compressible_solver_options.use_sst_turbulence
        add_boundary_sst_advection_flux!(du, u, idx_a, face_area, density_flux)
    end

    return nothing
end
