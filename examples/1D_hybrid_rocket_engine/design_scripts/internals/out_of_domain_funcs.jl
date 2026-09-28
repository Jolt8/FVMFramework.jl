function isoutofdomain_feedsystem_expanded(u, p, t, u_axes)
    u_named = ComponentVector(u, u_axes)

    #=
    if u_named.tank_oxidizer_mass < 0.0 || u_named.mid_section_mass < 0.0 || u_named.chamber_gas_mass < 0.0
        @show "caught out of domain"
    end
    =#

    return (
        u_named.tank_oxidizer_mass < 0.0 ||
        u_named.mid_section_mass < 0.0 ||
        u_named.chamber_gas_mass < 0.0
    )
end

isoutofdomain_feedsystem = let
    u_axes_local = u_axes
    (u, p, t) -> isoutofdomain_feedsystem_expanded(u, p, t, u_axes_local)
end