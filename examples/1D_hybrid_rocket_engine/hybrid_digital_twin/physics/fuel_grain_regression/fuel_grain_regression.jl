# Moving fuel-port geometry and fuel regression.
#
# Each axial station pairs one fluid_domain cell (the port gas) with the fuel_grain cell around it.
# The port radius r_p is a state stored in the fuel grain cell (u.fuel_port_radius[fuel_grain_cell_id]).
# Every geometric quantity that depends on r_p is recomputed from it at the start of each RHS call,
# so nothing in geo is mutated and the geometry stays differentiable for ForwardDiff and the tracers.
#
# Geometry of one station (coaxial cylinders, axial length L):
#     port gas volume          V_fluid = π r_p^2 L
#     fuel grain volume        V_grain = π (R_grain^2 - r_p^2) L
#     burning surface area     A_wall  = 2 π r_p L
#     port cross section       A_port  = π r_p^2
#
# r_p: port radius (m); R_grain: outer radius of the fuel grain, the liner's inner radius (m); L: station length (m).


"""
Collect the fixed information about each axial station once, at setup, from the ROM export:
which fluid and fuel grain cells belong together, the station length, the grain's outer radius,
and the faces used by the axial and radial fluxes between moving cells.

Stations are ordered by axial position, so neighbouring stations in the returned vector are
neighbouring cells along the port.
"""
function build_fuel_port_geometry(geo, rom_data)
    fluid_domain_cells = get_cellset(geo, "fluid_domain")
    fuel_grain_cells = get_cellset(geo, "fuel_grain")
    phenolic_liner_cells = get_cellset(geo, "phenolic_liner")

    pre_chamber_cell_id = get_cellset(geo, "pre_chamber")[1]

    # Order the port cells from one end of the grain to the other by their axial position
    sorted_fluid_domain_cells = sort(fluid_domain_cells, by = cell_id -> rom_data["cells"][cell_id]["z_min"])

    stations = []

    for fluid_cell_id in sorted_fluid_domain_cells
        rom_fluid_cell = rom_data["cells"][fluid_cell_id]

        # Float64 and Int conversions below: rom_data comes from JSON, so a value written as 0 instead of 0.0 parses as an Int.
        # One mismatched station would give the stations vector a non-concrete element type, and every station read
        # in the RHS would then be dynamically dispatched.
        axial_length = Float64(rom_fluid_cell["z_max"] - rom_fluid_cell["z_min"])

        # The fuel grain cell for this station is the fuel_grain neighbour across the fluid cell's cylindrical wall
        fuel_grain_cell_id = 0
        fluid_wall_face = 0
        grain_inner_face = 0

        for (neighbor_cell_id, local_face, neighbor_local_face) in geo.cell_neighbors[fluid_cell_id][2]
            if neighbor_cell_id in fuel_grain_cells
                fuel_grain_cell_id = neighbor_cell_id
                fluid_wall_face = local_face
                grain_inner_face = neighbor_local_face
            end
        end

        if fuel_grain_cell_id == 0
            error("fluid_domain cell $fluid_cell_id has no fuel_grain neighbour")
        end

        rom_fuel_grain_cell = rom_data["cells"][fuel_grain_cell_id]

        grain_outer_radius = Float64(rom_fuel_grain_cell["r_outer"])

        # The liner cell across the grain's outer wall. Its radius never changes, but the grain's
        # half-cell distance to it does, so the conduction across that face is recomputed every call.
        phenolic_liner_cell_id = 0
        grain_outer_face = 0
        liner_inner_face = 0

        for (neighbor_cell_id, local_face, neighbor_local_face) in geo.cell_neighbors[fuel_grain_cell_id][2]
            if neighbor_cell_id in phenolic_liner_cells
                phenolic_liner_cell_id = neighbor_cell_id
                grain_outer_face = local_face
                liner_inner_face = neighbor_local_face
            end
        end

        push!(stations, (
            fluid_cell_id = fluid_cell_id,
            fuel_grain_cell_id = fuel_grain_cell_id,
            phenolic_liner_cell_id = phenolic_liner_cell_id,
            axial_length = axial_length,
            grain_outer_radius = grain_outer_radius,
            fluid_wall_face = fluid_wall_face,
            grain_inner_face = grain_inner_face,
            grain_outer_face = grain_outer_face,
            liner_inner_face = liner_inner_face,
        ))
    end

    # Faces between neighbouring stations, for the axial gas flux and the axial conduction in the grain
    axial_faces = []

    for station_index in 1:(length(stations) - 1)
        lower_station = stations[station_index]
        upper_station = stations[station_index + 1]

        fluid_face_lower, fluid_face_upper, fluid_centroid_distance = find_shared_face(geo, lower_station.fluid_cell_id, upper_station.fluid_cell_id)
        grain_face_lower, grain_face_upper, grain_centroid_distance = find_shared_face(geo, lower_station.fuel_grain_cell_id, upper_station.fuel_grain_cell_id)

        push!(axial_faces, (
            lower_station_index = station_index,
            upper_station_index = station_index + 1,
            fluid_face_lower = fluid_face_lower,
            fluid_face_upper = fluid_face_upper,
            fluid_centroid_distance = Float64(fluid_centroid_distance),
            grain_face_lower = grain_face_lower,
            grain_face_upper = grain_face_upper,
            grain_centroid_distance = Float64(grain_centroid_distance),
        ))
    end

    # Lookup from a cell id to its station, so patch functions can find the station of the cell they were given (0 = not a port cell)
    station_index_by_cell = zeros(Int, geo.n_cells)

    for (station_index, station) in enumerate(stations)
        station_index_by_cell[station.fluid_cell_id] = station_index
        station_index_by_cell[station.fuel_grain_cell_id] = station_index
    end

    # The pre-chamber is a full-diameter gas cell. Its radius sets the contraction loss into the port.
    pre_chamber_radius = Float64(rom_data["cells"][pre_chamber_cell_id]["r_outer"])

    return (
        stations = [station for station in stations],
        axial_faces = [axial_face for axial_face in axial_faces],
        station_index_by_cell = station_index_by_cell,
        pre_chamber_cell_id = pre_chamber_cell_id,
        pre_chamber_radius = pre_chamber_radius,
    )
end

function find_shared_face(geo, cell_id_a, cell_id_b)
    for (neighbor_cell_id, local_face, neighbor_local_face) in geo.cell_neighbors[cell_id_a][2]
        if neighbor_cell_id == cell_id_b
            centroid_distance = geo.cell_neighbor_distances[cell_id_a, local_face]

            return local_face, neighbor_local_face, centroid_distance
        end
    end

    error("cells $cell_id_a and $cell_id_b do not share a face")
end


"""
Write the current volume of every cell into u.cell_volume. Fixed cells copy their CAD volume; the port
gas and fuel grain cells are recomputed from the port radius:

    V_fluid = π r_p^2 L        V_grain = π (R_grain^2 - r_p^2) L

Runs first in solve_system!, before any property update or flux, because everything that turns a mass
or energy flow into a change of density or temperature divides by these volumes.
"""
function update_fuel_port_geometry!(du, u, p, t, system, geo)
    for cell_id in 1:geo.n_cells
        u.cell_volume[cell_id] = geo.cell_volumes[cell_id]
    end

    for station in system.fuel_port_geometry.stations
        fuel_port_radius = u.fuel_port_radius[station.fuel_grain_cell_id]

        port_gas_volume = pi * fuel_port_radius^2 * station.axial_length

        fuel_grain_volume = pi * (station.grain_outer_radius^2 - fuel_port_radius^2) * station.axial_length

        u.cell_volume[station.fluid_cell_id] = port_gas_volume
        u.cell_volume[station.fuel_grain_cell_id] = fuel_grain_volume
    end
end


"""
How strongly the fuel is regressing at a station, from 0 (cold) to 1 (burning).

Regression is keyed on the gas temperature next to the fuel, not the grain cell's temperature: the grain
cell is the whole grain thickness, and its average temperature barely rises during a burn even though the
surface itself pyrolyzes within a fraction of a second. A smoothstep from min_fuel_regression_temperature to
100 K above it avoids a hard switch, which would make the implicit solver reject steps at the crossing.
Cold flow (or any run that isn't a hot fire) never regresses.
"""
function fuel_regression_switch(u, p, system, fluid_cell_id)
    if system.cold_flow || !system.hot_fire
        return zero(u.temperature[fluid_cell_id])
    end

    regression_ramp_width = 100.0

    ramp_position = clamp((u.temperature[fluid_cell_id] - p.min_fuel_regression_temperature) / regression_ramp_width, 0.0, 1.0)

    return ramp_position^2 * (3.0 - 2.0 * ramp_position)
end


"""
Regress the fuel at every station and inject the pyrolyzed fuel into the port gas.

Regression law, driven by the injected oxidizer flow over the station's own port area:

    G_ox = ṁ_ox / (π r_p^2)        ṙ = a * G_ox^n

G_ox: oxidizer mass flux (kg/(m^2*s)); ṁ_ox: injector oxidizer mass flow (kg/s); ṙ: regression rate (m/s);
a, n: regression_coefficient_a and regression_exponent_n. Stations with a larger port see a lower flux,
so the port shape evolves along the grain.

Pyrolyzed fuel, entering the port gas as ethylene through the burning surface:

    ṁ_fuel = ρ_fuel * ṙ * A_wall        dV_fluid/dt = A_wall * ṙ

ρ_fuel: solid fuel density (kg/m^3); A_wall = 2 π r_p L, the burning surface area (m^2).

The gas cell's balances, with M = ρ V its gas mass:
    mass:     dρ/dt   = (ṁ_fuel - ρ dV/dt) / V           (the -ρ dV/dt term is the cell growing around the same gas)
    species:  dY_i/dt = (δ_i,ethylene - Y_i) ṁ_fuel / M
    energy:   M de/dt = ṁ_fuel (h_fuel - e) - ṁ_fuel h_gasification - P dV/dt

h_fuel = cp_ethylene(T_surface) * T_surface: the fuel gas enthalpy on the chamber's sensible basis (J/kg);
e = cv T: the cell's specific internal energy (J/kg); h_gasification: fuel_heat_of_pyrolysis, the heat the
gas supplies to turn solid fuel at ambient temperature into gas at the surface (J/kg); P dV/dt: expansion
work as the cell grows (W), about 0.1% of the other terms but kept so energy balances exactly.

The fuel grain cell only shrinks: its density stays constant and its mass falls with its volume.
"""
function regress_fuel_grain!(du, u, p, t, system, geo)
    fuel_port_geometry = system.fuel_port_geometry

    oxidizer_mass_flow = u.injector_oxidizer_mass_flow[fuel_port_geometry.pre_chamber_cell_id]

    # Ethylene is the second species in system.model (ReidIdeal); ideal-gas cp doesn't depend on the 1e5 Pa placeholder
    pure_ethylene_composition = SVector(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    ethylene_cp_at_fuel_surface = mass_isobaric_heat_capacity(system.model, 1e5, p.fuel_surface_temperature, pure_ethylene_composition)

    fuel_gas_specific_enthalpy = ethylene_cp_at_fuel_surface * p.fuel_surface_temperature

    for station in fuel_port_geometry.stations
        fluid_cell_id = station.fluid_cell_id
        fuel_grain_cell_id = station.fuel_grain_cell_id

        # Current geometry of this station
        fuel_port_radius = u.fuel_port_radius[fuel_grain_cell_id]

        port_cross_section_area = pi * fuel_port_radius^2

        burning_surface_area = 2 * pi * fuel_port_radius * station.axial_length

        port_gas_volume = u.cell_volume[fluid_cell_id]


        # Regression rate. G_ox^n has an infinite slope at G_ox = 0 (n < 1), which turns into NaN in the
        # Jacobian before the valve opens, so the flux is floored at a negligible value. The switch is 0 before
        # ignition anyway, so the floor never produces any regression.
        minimum_oxidizer_mass_flux_for_derivatives = 1e-6

        oxidizer_mass_flux = max(oxidizer_mass_flow / port_cross_section_area, minimum_oxidizer_mass_flux_for_derivatives)

        regression_switch = fuel_regression_switch(u, p, system, fluid_cell_id)

        regression_rate = regression_switch * p.regression_coefficient_a * oxidizer_mass_flux^p.regression_exponent_n

        du.fuel_port_radius[fuel_grain_cell_id] += regression_rate


        # Fuel entering the port gas, and the port gas cell growing as the wall recedes
        solid_fuel_density = u.density[fuel_grain_cell_id]

        pyrolyzed_fuel_mass_flow = solid_fuel_density * regression_rate * burning_surface_area

        port_gas_volume_growth_rate = burning_surface_area * regression_rate


        # Mass: the cell gains the fuel, and the same gas also spreads over a larger volume
        port_gas_density = u.density[fluid_cell_id]

        port_gas_mass = port_gas_density * port_gas_volume

        du.density[fluid_cell_id] += (pyrolyzed_fuel_mass_flow - port_gas_density * port_gas_volume_growth_rate) / port_gas_volume


        # Species: the gas mixes toward pure ethylene at the fuel inflow rate.
        # Applied as -Y_i * ṁ / M for every species, plus ṁ / M for ethylene.
        for_fields!(du.mass_fractions, u.mass_fractions) do species, du_mass_fractions, u_mass_fractions
            du_mass_fractions[species[fluid_cell_id]] -= u_mass_fractions[species[fluid_cell_id]] * pyrolyzed_fuel_mass_flow / port_gas_mass
        end

        du.mass_fractions.ethylene[fluid_cell_id] += pyrolyzed_fuel_mass_flow / port_gas_mass


        # Energy, added to du.heat (W) and converted to temperature by cap_heat_flux_to_temp_change!
        port_gas_specific_internal_energy = u.cv[fluid_cell_id] * u.temperature[fluid_cell_id]

        fuel_inflow_power = pyrolyzed_fuel_mass_flow * (fuel_gas_specific_enthalpy - port_gas_specific_internal_energy)

        gasification_heat_power = pyrolyzed_fuel_mass_flow * p.fuel_heat_of_pyrolysis

        expansion_work_power = u.pressure[fluid_cell_id] * port_gas_volume_growth_rate

        du.heat[fluid_cell_id] += fuel_inflow_power - gasification_heat_power - expansion_work_power


        # Diagnostics only
        u.fuel_regression_rate[fuel_grain_cell_id] = regression_rate
        u.pyrolyzed_fuel_mass_flow[fluid_cell_id] = pyrolyzed_fuel_mass_flow
    end
end


"""
Convective heating of the fuel grain by the port gas, through the current burning surface area A_wall = 2 π r_p L:

    Q = h A_wall (T_gas - T_grain) (1 - S)

h: chamber_convective_heat_transfer_coefficient (W/(m^2*K)); S: the regression switch.

Once the fuel is regressing (S = 1), the heat reaching the surface is already accounted for by the gasification
term in regress_fuel_grain! and leaves with the pyrolyzed fuel, so convection into the grain is faded out.
Counting both would pass the same heat twice and heat the grain bulk unphysically.
"""
function fuel_port_wall_convection!(du, u, p, t, system, geo, fluid_cell_id, fuel_grain_cell_id)
    station = system.fuel_port_geometry.stations[system.fuel_port_geometry.station_index_by_cell[fluid_cell_id]]

    fuel_port_radius = u.fuel_port_radius[fuel_grain_cell_id]

    burning_surface_area = 2 * pi * fuel_port_radius * station.axial_length

    regression_switch = fuel_regression_switch(u, p, system, fluid_cell_id)

    convective_heat_flow = p.chamber_convective_heat_transfer_coefficient * burning_surface_area * (u.temperature[fluid_cell_id] - u.temperature[fuel_grain_cell_id]) * (1 - regression_switch)

    du.heat[fluid_cell_id] -= convective_heat_flow
    du.heat[fuel_grain_cell_id] += convective_heat_flow
end


"""
Conduction from the fuel grain into the phenolic liner. The liner side is fixed, but the grain's half-cell
distance grows as the grain thins. With the ROM export's geometric-mean radial node, r_node = sqrt(r_p R_grain),
the equivalent flat-wall distance from the node to the outer face is

    d_grain = R_grain ln(R_grain / r_node) = (R_grain / 2) ln(R_grain / r_p)

which reproduces the export's 4.54 mm at the initial 9.66 mm port radius and goes to zero as r_p reaches R_grain.
"""
function fuel_grain_to_liner_conduction!(du, u, p, t, system, geo, fuel_grain_cell_id, phenolic_liner_cell_id)
    station = system.fuel_port_geometry.stations[system.fuel_port_geometry.station_index_by_cell[fuel_grain_cell_id]]

    fuel_port_radius = u.fuel_port_radius[fuel_grain_cell_id]

    # A rejected trial state can carry a non-positive port radius, where log() would throw and end the run.
    # NaN makes the solver reject it instead (!(r > 0) also catches NaN).
    if !(fuel_port_radius > 0)
        du.heat[fuel_grain_cell_id] += oftype(fuel_port_radius, NaN)
        return nothing
    end

    grain_half_cell_distance = (station.grain_outer_radius / 2) * log(station.grain_outer_radius / fuel_port_radius)

    grain_face_area = geo.cell_face_areas[fuel_grain_cell_id, station.grain_outer_face]

    liner_face_area = geo.cell_face_areas[phenolic_liner_cell_id, station.liner_inner_face]
    liner_half_cell_distance = geo.cell_face_distances[phenolic_liner_cell_id, station.liner_inner_face]

    heat_diffusion_series_resistance_custom_geometry!(
        du, u, p, t, system, geo,
        fuel_grain_cell_id, station.grain_outer_face,
        phenolic_liner_cell_id, station.liner_inner_face,
        grain_face_area, grain_half_cell_distance,
        liner_face_area, liner_half_cell_distance,
    )
end
