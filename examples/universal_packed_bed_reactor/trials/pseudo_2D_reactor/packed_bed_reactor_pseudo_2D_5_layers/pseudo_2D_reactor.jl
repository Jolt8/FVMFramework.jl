using Unitful
using OrdinaryDiffEq
using Ferrite
using FerriteGmsh
using SparseConnectivityTracer
using ComponentArrays
import ADTypes
using NonlinearSolve
using Sparspak

using XLSX
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using ForwardDiff
using DataFrames
using CSV
using DataInterpolations
using Dates

using FVMFramework
using Plots

#to see an explanation of the physical build of this reactor check out:
#joinpath(@__DIR__, "..", "..", "..", "notes", "reactor_design.md")

#plan for mapping layering to FVM groups:
    #cell group 1: fluid 
    #cell group 2: pipe wall (steel)
    #cell group 3: thermocouple with the two layers of mica tape and fiberglass insulation 
    #cell group 4: heating wire (this one should probably be lumped with the thermocouple group)
    #cell group 5: ceramic blanket + aluminum foil insulation

pipe_inside_diameter = 16.0u"mm"
pipe_length = 12.1u"inch" |> u"m"

stripped_pipe_length = ustrip(pipe_length |> u"m")
pipe_width = ustrip(pipe_inside_diameter |> u"m")

n_cells_axial = 100
n_layers = 5

grid_dimensions = (1, n_layers, n_cells_axial)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((pipe_width, pipe_width * n_layers, stripped_pipe_length))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

evaporator_endpoint = 3u"inch" |> u"m"
evaporator_endpoint_cell = Int(round(evaporator_endpoint / (pipe_length / n_cells_axial)))

addcellset!(grid, "pipe_inlet", xyz -> xyz[3] <= (1 * (stripped_pipe_length / n_cells_axial)) && xyz[2] <= pipe_width)
getcellset(grid, "pipe_inlet")

addcellset!(grid, "silicon_carbide_preheater", xyz -> xyz[3] >= (1 * (stripped_pipe_length / n_cells_axial)) && xyz[3] <= (evaporator_endpoint_cell * (stripped_pipe_length / n_cells_axial)) && xyz[2] <= pipe_width)
getcellset(grid, "silicon_carbide_preheater")

addcellset!(grid, "copper_mesh_reformer", xyz -> xyz[3] >= (evaporator_endpoint_cell * (stripped_pipe_length / n_cells_axial)) && xyz[3] <= ((n_cells_axial - 1) * (stripped_pipe_length / n_cells_axial)) && xyz[2] <= pipe_width)
getcellset(grid, "copper_mesh_reformer")

addcellset!(grid, "pipe_outlet", xyz -> xyz[3] >= ((n_cells_axial - 1) * (stripped_pipe_length / n_cells_axial)) && xyz[2] <= pipe_width)
getcellset(grid, "pipe_outlet")

addcellset!(grid, "steel_pipe_wall", xyz -> xyz[2] >= pipe_width && xyz[2] <= 2 * pipe_width)
getcellset(grid, "steel_pipe_wall")

addcellset!(grid, "thermocouple_and_jacket", xyz -> xyz[2] >= 2 * pipe_width && xyz[2] <= 3 * pipe_width)
getcellset(grid, "thermocouple_and_jacket")

addfacetset!(grid, "thermocouple_to_heating_wire", xyz -> xyz[2] == 3 * pipe_width)
getfacetset(grid, "thermocouple_to_heating_wire")

addcellset!(grid, "heating_wire", xyz -> xyz[2] >= 3 * pipe_width && xyz[2] <= 4 * pipe_width)
getcellset(grid, "heating_wire")

addcellset!(grid, "insulation", xyz -> xyz[2] >= 4 * pipe_width)
getcellset(grid, "insulation")

addfacetset!(grid, "insulation_to_air", xyz -> xyz[2] >= 5 * pipe_width)
getfacetset(grid, "insulation_to_air")

addfacetset!(grid, "pipe_endcaps_to_air", xyz -> xyz[2] >= pipe_width && xyz[2] <= 2 * pipe_width && (xyz[3] == 0.0 || xyz[3] == stripped_pipe_length))
getfacetset(grid, "pipe_endcaps_to_air")

n_cells = length(grid.cells)

struct Fluid <: AbstractPhysics end
struct Solid <: AbstractPhysics end
#struct CustomFlux <: AbstractPhysics end #this is for when we want to only apply thermal resistance to some interfaces
struct ThermocoupleAndJacket <: AbstractPhysics end
struct HeatingWire <: AbstractPhysics end

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "physics", "multiphase.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "physics", "energy.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "..", "physics", "momentum.jl"))

#we use a polynomial due to the fact that the endcaps allow for a lot more heat to escape the ends of the reactor 
#than the cells closer to the center of the reactor 
#we don't need p1

Revise.includet(joinpath(@__DIR__, "..", "..", "pseudo_2D_reactor", "packing_and_fluid_property_merging.jl"))

function update_copper_mesh_properties!(du, u, p, t, cell_id, vol, system)
    mw_avg!(u, cell_id)
    get_fluid_and_copper_mesh_packing_rho!(du, u, cell_id, vol)
    get_fluid_and_copper_mesh_packing_k!(du, u, cell_id, vol)
    get_fluid_and_copper_mesh_packing_cp!(du, u, cell_id, vol)
    molar_concentrations!(u, cell_id)
end

function update_silicon_carbide_bed_packing_properties!(du, u, p, t, cell_id, vol, system)
    mw_avg!(u, cell_id)
    get_fluid_and_silicon_carbide_bed_packing_rho!(du, u, cell_id, vol)
    get_fluid_and_silicon_carbide_bed_packing_k!(du, u, cell_id, vol)
    get_fluid_and_silicon_carbide_bed_packing_cp!(du, u, cell_id, vol)
    molar_concentrations!(u, cell_id)
end

function update_fluid_properties!(du, u, p, t, cell_id, vol, system)
    mw_avg!(u, cell_id)

    properties = ComponentVector(system.properties_vec, system.properties_axes)
    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.rho[cell_id]

    molar_concentrations!(u, cell_id)
end

function update_solid_properties!(du, u, p, t, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id]
    u.rho[cell_id] = properties.rho[cell_id]
    #we still have to come up with a way to automatically do this
    #we could probably make a generated function that automatically applies fixed values to cached values 
    #I checked and this is working for now
end

function update_steel_pipe_properties!(du, u, p, t, cell_id, vol, system)
    properties = ComponentVector(system.properties_vec, system.properties_axes)

    u.k[cell_id] = properties.k[cell_id]
    u.cp[cell_id] = properties.cp[cell_id] * u.steel_thermal_mass_multiplier[cell_id]
    u.rho[cell_id] = properties.rho[cell_id]
end

#you could also add any of these functions individually

function fluid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    sum_mass_flux_face_to_cell!(du, u, cell_id, vol) #this always has to go before cap_mass_flux_to_pressure_change!

    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    cap_mass_flux_to_pressure_change!(du, u, cell_id, vol)
    cap_species_mass_flux_to_mass_fraction_change!(du, u, cell_id, vol)
end

function solid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
end

u_proto = ComponentVector(
    mass_fractions = (
        #methanol = zeros(n_cells)u"kg/kg",
        water = zeros(n_cells)u"kg/kg",
        #carbon_monoxide = zeros(n_cells)u"kg/kg",
        #hydrogen = zeros(n_cells)u"kg/kg",
        #carbon_dioxide = zeros(n_cells)u"kg/kg",
        air = zeros(n_cells)u"kg/kg"
    ),
    pressure = zeros(n_cells)u"Pa",
    temp = zeros(n_cells)u"K",
    #liquid_holdup = zeros(n_cells)u"m^3/m^3",
    #gas_holdup = zeros(n_cells)u"m^3/m^3"
)

config = create_fvm_config(grid, u_proto);

cell_lengths_along_pipe = [config.geo.cell_centroids[i][3]u"m" for i in 1:length(config.geo.cell_centroids)]

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "5_layer_common_properties.jl"))
common_properties = get_5_layer_common_properties(pipe_length, cell_lengths_along_pipe, n_cells_axial)

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "packed_media_properties", "copper_mesh_reformer_properties.jl"))
pre_copper_mesh_reformer_properties = get_copper_mesh_reformer_properties()
copper_mesh_reformer_properties = merge_properties(pre_copper_mesh_reformer_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "packed_media_properties", "silicon_carbide_sand_properties.jl"))
pre_silicon_carbide_sand_properties = get_silicon_carbide_sand_properties()
silicon_carbide_sand_properties = merge_properties(pre_silicon_carbide_sand_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "steel_pipe_wall_properties", "steel_pipe_wall_properties.jl"))
pre_steel_pipe_wall_properties = get_steel_pipe_wall_properties()
steel_pipe_wall_properties = merge_properties(pre_steel_pipe_wall_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "thermocouple_and_heating_wire_properties", "divided_thermocouple_and_heating_wire_properties", "thermocouple_and_jacket_properties.jl"))
pre_thermocouple_and_jacket_properties = get_thermocouple_and_jacket_properties(grid, pipe_length, n_cells_axial, cell_lengths_along_pipe)
thermocouple_and_jacket_properties = merge_properties(pre_thermocouple_and_jacket_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "thermocouple_and_heating_wire_properties", "divided_thermocouple_and_heating_wire_properties", "heating_wire_properties.jl"))
pre_heating_wire_properties = get_heating_wire_properties(grid, pipe_length, common_properties)
heating_wire_properties = merge_properties(pre_heating_wire_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "properties", "insulation_properties", "insulation_properties.jl"))
pre_insulation_properties = get_insulation_properties()
insulation_properties = merge_properties(pre_insulation_properties, common_properties)

Revise.includet(joinpath(@__DIR__, "pseudo_2D_geometry_editing_5_layers.jl"))
edit_pseudo_2D_geometry_5_layers!(config.geo, grid, common_properties)

n_cells = length(config.geo.cell_volumes)
n_faces = length(config.geo.cell_neighbor_areas[1])
reaction_names = keys(copper_mesh_reformer_properties.reactions.reforming_reactions)
species_names = keys(copper_mesh_reformer_properties.molecular_weights)

add_setup_syms!(config;
    cache_syms_and_units = (
        heat = u"J",
        mw_avg = u"kg/mol",
        k = u"W/(m*K)",
        cp = u"J/(kg*K)",
        rho = u"kg/m^3",
        insulation_to_air_overall_heat_transfer_coefficient_to_environment = u"W/(m^2*K)",
        pipe_endcaps_to_air_thermal_conductance = u"W/K",
        heater_weight_1 = u"1",
        heater_weight_2 = u"1",
        heater_weight_3 = u"1",
        heater_weight_4 = u"1",
        fluid_to_steel_pipe_convective_heat_transfer_coefficient = u"W/(m^2*K)",
        steel_thermal_mass_multiplier = u"1",
        thermocouple_to_heating_wire_thermal_resistance = u"K/W",
        wattage_received_per_m = u"W/m",
        molar_concentrations = u"mol/m^3",
        species_masses = u"kg",
        net_rates = u"mol/s",
        mass = u"kg",
        mass_face = u"kg",
        mass_evaporated = u"kg",
        superficial_velocity = u"m/s",
        measured_heater_wattage_per_cell = u"W", 
        pipe_mass_flow = u"kg/s"
    ),
    special_caches = ComponentArray(
        mass_face = zeros(n_cells, n_faces)u"kg",
        net_rates = (
            reforming_reactions = NamedTuple{reaction_names}(
                Tuple(zeros(n_cells)u"mol/s" for _ in 1:length(reaction_names))
            ), #don't forget the comma!
        ), 
        molar_concentrations = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"mol/m^3" for _ in 1:length(species_names))
        ),
        species_masses = NamedTuple{species_names}(
            Tuple(zeros(n_cells)u"kg" for _ in 1:length(species_names))
        )
    ),
    optimized_parameters = ComponentVector(
        insulation_to_air_overall_heat_transfer_coefficient_to_environment = 0.0u"W/(m^2*K)",
        pipe_endcaps_to_air_thermal_conductance = 0.0u"W/K",
        heater_weight_1 = 0.0u"1",
        heater_weight_2 = 0.0u"1",
        heater_weight_3 = 0.0u"1",
        heater_weight_4 = 0.0u"1",
        fluid_to_steel_pipe_convective_heat_transfer_coefficient = 0.0u"W/(m^2*K)",
        steel_thermal_mass_multiplier = 0.0,
        TC1_thermal_resistance = 0.0u"K/W",
        TC2_thermal_resistance = 0.0u"K/W",
        TC3_thermal_resistance = 0.0u"K/W",
        TC4_thermal_resistance = 0.0u"K/W",
        TC5_thermal_resistance = 0.0u"K/W"
    )
)

function fluid_physics_functions!(du, u, cell_id, vol)
end

function solid_physics_functions!(du, u, cell_id, vol)
end

#TODO: remember to turn off the heater and to set 
#the inlet mass flow to what was observed in the experiment

placeholder_temperature = 99999.9u"K" #we want this to error hard if it doesn't get updated by a trial's measured room temperature
placeholder_mass_fractions = ComponentVector(
    #methanol = 1e-99u"kg/kg",
    water = 1e-99u"kg/kg",
    #carbon_monoxide = 1e-99u"kg/kg",
    #hydrogen = 1e-99u"kg/kg",
    #carbon_dioxide = 1e-99u"kg/kg",
    air = 1e-99u"kg/kg"
)

add_region!(
    config, "pipe_inlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = silicon_carbide_sand_properties,
    property_update_function = update_silicon_carbide_bed_packing_properties!,
    region_function =
    function inlet!(du, u, p, t, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 1] += u.pipe_mass_flow[cell_id]
        
        du.heat[cell_id] *= 0.0
        #du.heat[cell_id] += u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id]

        fluid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)

        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
        end
    end
)

add_region!(
    config, "silicon_carbide_preheater";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = silicon_carbide_sand_properties,
    property_update_function = update_silicon_carbide_bed_packing_properties!,
    region_function =
    function evaporator!(du, u, p, t, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        fluid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    end
)

add_region!(
    config, "copper_mesh_reformer";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = copper_mesh_reformer_properties,
    property_update_function = update_copper_mesh_properties!,
    region_function =
    function reactor!(du, u, p, t, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        #PAM_reforming_react_cell!(du, u, cell_id, vol)

        fluid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    end
)

add_region!(
    config, "pipe_outlet";
    type = Fluid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = copper_mesh_reformer_properties,
    property_update_function = update_copper_mesh_properties!,
    region_function =
    function outlet!(du, u, p, t, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)
        
        du.mass_face[cell_id, 6] -= u.pipe_mass_flow[cell_id]

        #this is to prevent the concentration of all species from building up at the outlet
        for_fields!(u.mass_fractions, du.species_masses) do species, u_mass_fractions, du_species_mass
            du_species_mass[species[cell_id]] -= u.pipe_mass_flow[cell_id] * u_mass_fractions[species[cell_id]]
        end

        #du.heat[cell_id] *= 0.0
        du.heat[cell_id] -= u.pipe_mass_flow[cell_id] * u.fluid_cp[cell_id] * u.temp[cell_id] 
        #this is to prevent the temperature from building up at the outlet

        fluid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    end
)

add_region!(
    config, "steel_pipe_wall";
    type = Solid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = steel_pipe_wall_properties,
    property_update_function = update_steel_pipe_properties!,
    region_function =
    function steel_pipe_wall!(du, u, p, t, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    end
)

add_region!(
    config, "thermocouple_and_jacket";
    type = ThermocoupleAndJacket(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = thermocouple_and_jacket_properties,
    property_update_function = update_solid_properties!,
    region_function =
    function thermocouple_and_jacket!(du, u, p, t, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    end
)

add_patch!(
    config, "thermocouple_to_heating_wire";
    properties = ComponentVector(), #no new properties here
    patch_function =
    function thermocouple_to_heating_wire!(
        du, u, p, t,
        idx_a, idx_b, face_idx, #idx_b is not applicable here because it connects to nothing
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        du.heat[idx_a] += (u.temp[idx_b] - u.temp[idx_a]) / u.thermocouple_to_heating_wire_thermal_resistance[idx_a]
    end
)

add_region!(
    config, "heating_wire";
    type = HeatingWire(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = heating_wire_properties,
    property_update_function = update_solid_properties!,
    region_function =
    function heating_wire!(du, u, p, t, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    end
)

add_region!(
    config, "insulation";
    type = Solid(),
    initial_conditions = ComponentVector(
        mass_fractions = placeholder_mass_fractions,
        pressure = 1.0u"atm",
        temp = placeholder_temperature,
    ),
    properties = insulation_properties,
    property_update_function = update_solid_properties!,
    region_function =
    function insulation!(du, u, p, t, cell_id, vol)
        solid_physics_functions!(du, u, cell_id, vol)

        solid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    end
)

add_patch!(
    config, "insulation_to_air";
    properties = ComponentVector(),
    patch_function =
    function insulation_to_air!(
        du, u, p, t, 
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        du.heat[idx_a] += (u.insulation_to_air_overall_heat_transfer_coefficient_to_environment[idx_a] * u.insulation_to_environment_cell_areas[idx_a]) * (u.room_temperature[idx_a] - u.temp[idx_a])
    end
)

add_patch!(
    config, "pipe_endcaps_to_air";
    properties = ComponentVector(),
    patch_function =
    function pipe_endcaps_to_air!(
        du, u, p, t, 
        idx_a, idx_b, face_idx,
        cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
        cell_volumes
    )
        du.heat[idx_a] += (u.pipe_endcaps_to_air_thermal_conductance[idx_a]) * (u.room_temperature[idx_a] - u.temp[idx_a])
    end
)

#Connection functions
function fluid_fluid_flux!(
    du, u, p, t, 
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    #=new_pressure_driven_mass_flux!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
    )=#

    all_species_advection!(
        du, u, 
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
    
    enthalpy_advection!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )

    heat_diffusion!(
        du, u,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end

Revise.includet(joinpath(@__DIR__, "..", "..", "..", "physics", "convective_heat_transfer.jl"))

function fluid_solid_flux!(
    du, u, p, t, 
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    fluid_to_steel_pipe_convective_heat_transfer!(
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

function no_flux!(
    du, u, p, t, 
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    #we do nothing here because we're going to apply thermal resistance based model somewhere else
end

function connection_map_function(phys_a, phys_b)
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Fluid && return fluid_fluid_flux!
    typeof(phys_a) <: Fluid && typeof(phys_b) <: Solid && return fluid_solid_flux!

    typeof(phys_a) <: Solid && typeof(phys_b) <: Fluid && return fluid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: Solid && return solid_solid_flux!

    typeof(phys_a) <: ThermocoupleAndJacket && typeof(phys_b) <: Solid && return solid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: ThermocoupleAndJacket && return solid_solid_flux!
    typeof(phys_a) <: ThermocoupleAndJacket && typeof(phys_b) <: ThermocoupleAndJacket && return no_flux! 
    #this is made almost completely of fiberglass tape, so axial heat transfer along the pipe between thermocouple cells can be ignored
    
    typeof(phys_a) <: HeatingWire && typeof(phys_b) <: Solid && return solid_solid_flux!
    typeof(phys_a) <: Solid && typeof(phys_b) <: HeatingWire && return solid_solid_flux!
    typeof(phys_a) <: HeatingWire && typeof(phys_b) <: HeatingWire && return no_flux!
    #this is made of heating wire covered with ceramic fiber, so axial heat transfer along the pipe between thermocouple cells can be ignored
    
    typeof(phys_a) <: HeatingWire && typeof(phys_b) <: ThermocoupleAndJacket && return no_flux! #we apply custom thermal resistances between these two, so no flux should occur
    typeof(phys_a) <: ThermocoupleAndJacket && typeof(phys_b) <: HeatingWire && return no_flux!
end

#heated cells
heater_start = 1.0u"inch"
heater_end = 11.0u"inch"

heating_wire_properties.cell_lengths_along_pipe

function get_heated_cells(heater_start, heater_end)
    heated_cells = Int64[]
    for cell_id in collect(grid.cellsets["heating_wire"])
        distance_along_pipe = heating_wire_properties.cell_lengths_along_pipe[cell_id]
        if upreferred(heater_start) <= distance_along_pipe <= upreferred(heater_end)
            push!(heated_cells, cell_id)
        end
    end
    return heated_cells
end

total_heater_length = 10u"inch"
heater_start = 1u"inch"
n_heater_regions = 5

heater_1_start = heater_start + (total_heater_length / n_heater_regions) * 0
heater_1_end = heater_start + (total_heater_length / n_heater_regions) * 1

heater_2_start = heater_start + (total_heater_length / n_heater_regions) * 1
heater_2_end = heater_start + (total_heater_length / n_heater_regions) * 2

heater_3_start = heater_start + (total_heater_length / n_heater_regions) * 2
heater_3_end = heater_start + (total_heater_length / n_heater_regions) * 3

heater_4_start = heater_start + (total_heater_length / n_heater_regions) * 3
heater_4_end = heater_start + (total_heater_length / n_heater_regions) * 4

heater_5_start = heater_start + (total_heater_length / n_heater_regions) * 4
heater_5_end = heater_start + (total_heater_length / n_heater_regions) * 5

heater_1_cells = get_heated_cells(heater_1_start, heater_1_end)
n_heater_1_cells = length(heater_1_cells)

heater_2_cells = get_heated_cells(heater_2_start, heater_2_end)
n_heater_2_cells = length(heater_2_cells)

heater_3_cells = get_heated_cells(heater_3_start, heater_3_end)
n_heater_3_cells = length(heater_3_cells)

heater_4_cells = get_heated_cells(heater_4_start, heater_4_end)
n_heater_4_cells = length(heater_4_cells)

heater_5_cells = get_heated_cells(heater_5_start, heater_5_end)
n_heater_5_cells = length(heater_5_cells)

n_total_heated_cells = n_heater_1_cells + n_heater_2_cells + n_heater_3_cells + n_heater_4_cells + n_heater_5_cells


#thermal reistance interp
TC1_TC2_midpoint = (thermocouple_and_jacket_properties.TC1_length_along_pipe + thermocouple_and_jacket_properties.TC2_length_along_pipe) / 2
TC2_TC3_midpoint = (thermocouple_and_jacket_properties.TC2_length_along_pipe + thermocouple_and_jacket_properties.TC3_length_along_pipe) / 2
TC3_TC4_midpoint = (thermocouple_and_jacket_properties.TC3_length_along_pipe + thermocouple_and_jacket_properties.TC4_length_along_pipe) / 2
TC4_TC5_midpoint = (thermocouple_and_jacket_properties.TC4_length_along_pipe + thermocouple_and_jacket_properties.TC5_length_along_pipe) / 2

TC1_area_start = 0.0u"inch"
TC1_area_end = TC1_TC2_midpoint

TC2_area_start = TC1_TC2_midpoint
TC2_area_end = TC2_TC3_midpoint

TC3_area_start = TC2_TC3_midpoint
TC3_area_end = TC3_TC4_midpoint

TC4_area_start = TC3_TC4_midpoint
TC4_area_end = TC4_TC5_midpoint

TC5_area_start = TC4_TC5_midpoint
TC5_area_end = common_properties.pipe_length

function get_cells_along_pipe_range(common_properties, region_names, start_point, end_point)
    cells_along_pipe_range = Int64[]

    for region_name in region_names
        region_cell_ids = grid.cellsets[region_name]
        for cell_id in region_cell_ids
            if start_point <= common_properties.cell_lengths_along_pipe[cell_id] < end_point
                push!(cells_along_pipe_range, cell_id)
            end
        end
    end

    return cells_along_pipe_range
end

TC1_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC1_area_start, TC1_area_end)
TC2_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC2_area_start, TC2_area_end)
TC3_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC3_area_start, TC3_area_end)
TC4_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC4_area_start, TC4_area_end)
TC5_cells = get_cells_along_pipe_range(common_properties, ["thermocouple_and_jacket", "heating_wire"], TC5_area_start, TC5_area_end)

Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_dry", "empty_trial_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_dry", "air_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "thermocouple_data_processing", "thermocouple_data_with_wattage.jl"))

function dry_run_trial()
    dry_run_config = deepcopy(config)
    dry_run_properties = deepcopy(common_properties)

    general_trial_properties = get_empty_trial_properties() 
    dry_run_properties = merge_properties(dry_run_properties, general_trial_properties)

    air_properties = get_air_properties()
    dry_run_properties = merge_properties(dry_run_properties, air_properties)

    region_names = [dry_run_config.regions[i].name for i in eachindex(dry_run_config.regions)]
    regions = Dict(region_names .=> dry_run_config.regions)

    function dry_run_inlet!(du, u, p, t, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        fluid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)
    end

    regions["pipe_inlet"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["pipe_inlet"].initial_conditions.mass_fractions = dry_run_properties.inlet_mass_fractions
    regions["pipe_inlet"].region_function = dry_run_inlet!
    regions["pipe_inlet"].properties = merge_properties(regions["pipe_inlet"].properties, air_properties)
    regions["pipe_inlet"].properties = merge_properties(regions["pipe_inlet"].properties, general_trial_properties)
    update_region!(dry_run_config, "pipe_inlet")

    regions["silicon_carbide_preheater"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["silicon_carbide_preheater"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions
    regions["silicon_carbide_preheater"].properties = merge_properties(regions["silicon_carbide_preheater"].properties, air_properties)
    regions["silicon_carbide_preheater"].properties = merge_properties(regions["silicon_carbide_preheater"].properties, general_trial_properties)
    update_region!(dry_run_config, "silicon_carbide_preheater")

    regions["copper_mesh_reformer"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["copper_mesh_reformer"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions
    regions["copper_mesh_reformer"].properties = merge_properties(regions["copper_mesh_reformer"].properties, air_properties)
    regions["copper_mesh_reformer"].properties = merge_properties(regions["copper_mesh_reformer"].properties, general_trial_properties)
    update_region!(dry_run_config, "copper_mesh_reformer")

    regions["pipe_outlet"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["pipe_outlet"].initial_conditions.mass_fractions = dry_run_properties.empty_mass_fractions
    regions["pipe_outlet"].properties = merge_properties(regions["pipe_outlet"].properties, air_properties)
    regions["pipe_outlet"].properties = merge_properties(regions["pipe_outlet"].properties, general_trial_properties)
    update_region!(dry_run_config, "pipe_outlet")

    regions["steel_pipe_wall"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["steel_pipe_wall"].properties = merge_properties(regions["steel_pipe_wall"].properties, general_trial_properties)
    update_region!(dry_run_config, "steel_pipe_wall")

    regions["thermocouple_and_jacket"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["thermocouple_and_jacket"].properties = merge_properties(regions["thermocouple_and_jacket"].properties, general_trial_properties)
    update_region!(dry_run_config, "thermocouple_and_jacket")

    regions["heating_wire"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["heating_wire"].properties = merge_properties(regions["heating_wire"].properties, general_trial_properties)
    update_region!(dry_run_config, "heating_wire")

    regions["insulation"].initial_conditions.temp = dry_run_properties.room_temperature
    regions["insulation"].properties = merge_properties(regions["insulation"].properties, general_trial_properties)
    update_region!(dry_run_config, "insulation")

    thermocouple_data = get_dry_run_thermocouple_data()

    function heater_1_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * p.heater_weight_1[1]) / n_heater_1_cells
    end

    function heater_2_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * p.heater_weight_2[1]) / n_heater_2_cells
    end

    function heater_3_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * p.heater_weight_3[1]) / n_heater_3_cells
    end

    function heater_4_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * (1 - p.heater_weight_3[1]) * p.heater_weight_4[1]) / n_heater_4_cells
    end

    function heater_5_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * (1 - p.heater_weight_3[1]) * (1 - p.heater_weight_4[1])) / n_heater_5_cells
    end

    return dry_run_config, dry_run_properties, thermocouple_data, heater_1_wattage_per_cell, heater_2_wattage_per_cell, heater_3_wattage_per_cell, heater_4_wattage_per_cell, heater_5_wattage_per_cell
end

Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_water_flow_trial", "hot_water_trial_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_water_flow_trial", "water_properties.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_water_flow_trial", "recorded_data_processing", "inlet_and_outlet_temperatures.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "packed_bed_reactor_water_flow_trial", "recorded_data_processing", "values_of_note.jl"))
Revise.includet(joinpath(@__DIR__, "..", "..", "thermocouple_data_processing", "thermocouple_data.jl"))

function hot_water_trial()
    hot_water_config = deepcopy(config)
    hot_water_trial_properties = deepcopy(common_properties)

    values_of_note = get_packed_bed_reactor_water_flow_trial_values_of_note()
    pump_shutoff_timestamp = ustrip(upreferred(values_of_note.pump_shut_off_time))

    general_trial_properties = get_hot_water_trial_properties(values_of_note)
    hot_water_trial_properties = merge_properties(hot_water_trial_properties, general_trial_properties)

    water_properties = get_water_properties()
    hot_water_trial_properties = merge_properties(hot_water_trial_properties, water_properties)

    inlet_and_outlet_temperatures = get_inlet_and_outlet_temperature_correlations()
    inlet_temp_interp = inlet_and_outlet_temperatures.inlet_temp_interp

    region_names = [hot_water_config.regions[i].name for i in eachindex(hot_water_config.regions)]
    regions = Dict(region_names .=> hot_water_config.regions)

    function hot_water_inlet!(du, u, p, t, cell_id, vol)
        fluid_physics_functions!(du, u, cell_id, vol)

        du.mass_face[cell_id, 1] += u.pipe_mass_flow[cell_id]
        
        du.heat[cell_id] *= 0.0
        #du.heat[cell_id] += u.pipe_mass_flow[cell_id] * u.cp[cell_id] * u.temp[cell_id]
        #the only reason we don't have to do this for the outlet is because the outlet actually does the du.heat[cell_id] += etc..

        fluid_sum_and_cap_fluxes!(du, u, p, t, cell_id, vol)

        for_fields!(du.mass_fractions) do species, du_mass_fractions
            du_mass_fractions[species[cell_id]] *= 0.0
        end
    end

    regions["pipe_inlet"].initial_conditions.temp = inlet_temp_interp(0.0)
    regions["pipe_inlet"].initial_conditions.mass_fractions = hot_water_trial_properties.inlet_mass_fractions
    regions["pipe_inlet"].region_function = hot_water_inlet!
    regions["pipe_inlet"].properties = merge_properties(regions["pipe_inlet"].properties, water_properties)
    regions["pipe_inlet"].properties = merge_properties(regions["pipe_inlet"].properties, general_trial_properties)
    update_region!(hot_water_config, "pipe_inlet")

    regions["silicon_carbide_preheater"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["silicon_carbide_preheater"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions
    regions["silicon_carbide_preheater"].properties = merge_properties(regions["silicon_carbide_preheater"].properties, water_properties)
    regions["silicon_carbide_preheater"].properties = merge_properties(regions["silicon_carbide_preheater"].properties, general_trial_properties)
    update_region!(hot_water_config, "silicon_carbide_preheater")

    regions["copper_mesh_reformer"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["copper_mesh_reformer"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions
    regions["copper_mesh_reformer"].properties = merge_properties(regions["copper_mesh_reformer"].properties, water_properties)
    regions["copper_mesh_reformer"].properties = merge_properties(regions["copper_mesh_reformer"].properties, general_trial_properties)
    update_region!(hot_water_config, "copper_mesh_reformer")

    regions["pipe_outlet"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["pipe_outlet"].initial_conditions.mass_fractions = hot_water_trial_properties.empty_mass_fractions
    regions["pipe_outlet"].properties = merge_properties(regions["pipe_outlet"].properties, water_properties)
    regions["pipe_outlet"].properties = merge_properties(regions["pipe_outlet"].properties, general_trial_properties)
    update_region!(hot_water_config, "pipe_outlet")

    regions["steel_pipe_wall"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["steel_pipe_wall"].properties = merge_properties(regions["steel_pipe_wall"].properties, general_trial_properties)
    update_region!(hot_water_config, "steel_pipe_wall")

    regions["thermocouple_and_jacket"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["thermocouple_and_jacket"].properties = merge_properties(regions["thermocouple_and_jacket"].properties, general_trial_properties)
    update_region!(hot_water_config, "thermocouple_and_jacket")

    regions["heating_wire"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["heating_wire"].properties = merge_properties(regions["heating_wire"].properties, general_trial_properties)
    update_region!(hot_water_config, "heating_wire")

    regions["insulation"].initial_conditions.temp = hot_water_trial_properties.room_temperature
    regions["insulation"].properties = merge_properties(regions["insulation"].properties, general_trial_properties)
    update_region!(hot_water_config, "insulation")

    inlet_cell_id = collect(grid.cellsets["pipe_inlet"])[1]

    function pump_shut_off(du, u, cell_id, t)
        if t <= pump_shutoff_timestamp #pump on
            #u.temp[inlet_cell_id] = inlet_temp_interp(ForwardDiff.value(t)) #this is a lot slower
            du.temp[inlet_cell_id] *= 0.0
            du.temp[inlet_cell_id] += DataInterpolations.derivative(inlet_temp_interp, ForwardDiff.value(t))

            u.pipe_mass_flow[cell_id] = ustrip(upreferred(hot_water_trial_properties.pipe_mass_flow))
        else #pump shut off
            u.pipe_mass_flow[cell_id] = 0.0
        end
    end

    thermocouple_data = get_hot_water_thermocouple_data()

    function heater_1_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * p.heater_weight_1[1]) / n_heater_1_cells
    end

    function heater_2_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * p.heater_weight_2[1]) / n_heater_2_cells
    end

    function heater_3_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * p.heater_weight_3[1]) / n_heater_3_cells
    end

    function heater_4_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * (1 - p.heater_weight_3[1]) * p.heater_weight_4[1]) / n_heater_4_cells
    end

    function heater_5_wattage_per_cell(p, t)
        return (thermocouple_data.heater_power_interp(ForwardDiff.value(t)) * (1 - p.heater_weight_1[1]) * (1 - p.heater_weight_2[1]) * (1 - p.heater_weight_3[1]) * (1 - p.heater_weight_4[1])) / n_heater_5_cells
    end

    return hot_water_config, hot_water_trial_properties, pump_shut_off, thermocouple_data, heater_1_wattage_per_cell, heater_2_wattage_per_cell, heater_3_wattage_per_cell, heater_4_wattage_per_cell, heater_5_wattage_per_cell
end

#TODO: I think we could have an integrated way that allows for the properties of many trials to be sequentially added by doing something like in the initial conditions:
#=
    initial_conditions = (
        mass_fractions = TrialSelector(
            dry_run = empty_mass_fractions,
            hot_water = common_properties.inlet_mass_fractions,
        ),
        temp = TrialSelector(
            dry_run = common_properties.room_temperature,
            hot_water = inlet_temp_interp(0.0u"s")
        ),
        pressure = 1.0u"bar", #this is automatically applied to both trials
    )
=#
#we'll do that next time if it's necessary

fluid_regions = ["pipe_inlet", "silicon_carbide_preheater", "copper_mesh_reformer", "pipe_outlet"]
advecting_fluid_cells = vcat(collect(grid.cellsets["pipe_inlet"]), collect(grid.cellsets["silicon_carbide_preheater"]), collect(grid.cellsets["copper_mesh_reformer"]), collect(grid.cellsets["pipe_outlet"]))

function trial_independent_solve_system!(du, u, p_vec, t, geo, system)
    p = ComponentVector(p_vec, system.p_axes)

    #TODO: Automate this process of applying the guessed values to the actual variables
    for cell_id in eachindex(geo.cell_volumes)
        u.insulation_to_air_overall_heat_transfer_coefficient_to_environment[cell_id] = p.insulation_to_air_overall_heat_transfer_coefficient_to_environment[1]
        u.pipe_endcaps_to_air_thermal_conductance[cell_id] = p.pipe_endcaps_to_air_thermal_conductance[1]
        u.fluid_to_steel_pipe_convective_heat_transfer_coefficient[cell_id] = p.fluid_to_steel_pipe_convective_heat_transfer_coefficient[1]
        u.steel_thermal_mass_multiplier[cell_id] = p.steel_thermal_mass_multiplier[1]
    end

    update_region_groups!(du, u, p, t, geo, system)

    for cell_id in TC1_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC1_thermal_resistance[1]
    end
    for cell_id in TC2_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC2_thermal_resistance[1]
    end
    for cell_id in TC3_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC3_thermal_resistance[1]
    end
    for cell_id in TC4_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC4_thermal_resistance[1]
    end
    for cell_id in TC5_cells
        u.thermocouple_to_heating_wire_thermal_resistance[cell_id] = p.TC5_thermal_resistance[1]
    end

    #VERY IMPORTANT: since most software uses 0-based indexing, you need to adjust the cell id by +1
    #for example, if you mouse over cell_id 5161 in ParaView, you need to use 5162 in the code because julia uses 1-based indexing 

    for i in 1:length(advecting_fluid_cells) - 1 #we don't take mass out of the outlet
        idx_a = advecting_fluid_cells[i]
        idx_b = advecting_fluid_cells[i + 1]
        
        du.mass_face[idx_a, 6] -= u.pipe_mass_flow[idx_a]
        du.mass_face[idx_b, 1] += u.pipe_mass_flow[idx_a]
    end

    solve_connection_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end

dry_run_config, dry_run_properties, dry_run_thermocouple_data, 
dry_run_heater_1_wattage_per_cell, dry_run_heater_2_wattage_per_cell, dry_run_heater_3_wattage_per_cell, 
dry_run_heater_4_wattage_per_cell, dry_run_heater_5_wattage_per_cell = dry_run_trial();
dry_run_du0_vec, dry_run_u0_vec, dry_run_geo, dry_run_system = finish_fvm_config(dry_run_config, connection_map_function, check_units = false);
#TODO: the time finish_fvm_config(...) takes to complete needs to be reduced

hot_water_config, hot_water_properties, hot_water_pump_shut_off, hot_water_thermocouple_data, 
hot_water_heater_1_wattage_per_cell, hot_water_heater_2_wattage_per_cell, hot_water_heater_3_wattage_per_cell, 
hot_water_heater_4_wattage_per_cell, hot_water_heater_5_wattage_per_cell = hot_water_trial();
hot_water_du0_vec, hot_water_u0_vec, hot_water_geo, hot_water_system = finish_fvm_config(hot_water_config, connection_map_function, check_units = false);

p_axes = hot_water_system.p_axes

function dry_run_solve_system!(du, u, p_vec, t, geo, system)
    p = ComponentVector(p_vec, p_axes)

    for cell_id in heater_1_cells
        du.heat[cell_id] += dry_run_heater_1_wattage_per_cell(p, t)
    end

    for cell_id in heater_2_cells
        du.heat[cell_id] += dry_run_heater_2_wattage_per_cell(p, t)
    end

    for cell_id in heater_3_cells
        du.heat[cell_id] += dry_run_heater_3_wattage_per_cell(p, t)
    end

    for cell_id in heater_4_cells
        du.heat[cell_id] += dry_run_heater_4_wattage_per_cell(p, t)
    end

    for cell_id in heater_5_cells
        du.heat[cell_id] += dry_run_heater_5_wattage_per_cell(p, t)
    end

    trial_independent_solve_system!(du, u, p_vec, t, geo, system)
end

function hot_water_solve_system!(du, u, p_vec, t, geo, system)
    p = ComponentVector(p_vec, p_axes)

    for cell_id in advecting_fluid_cells
        hot_water_pump_shut_off(du, u, cell_id, t)
    end

    for cell_id in heater_1_cells
        du.heat[cell_id] += hot_water_heater_1_wattage_per_cell(p, t)
    end

    for cell_id in heater_2_cells
        du.heat[cell_id] += hot_water_heater_2_wattage_per_cell(p, t)
    end

    for cell_id in heater_3_cells
        du.heat[cell_id] += hot_water_heater_3_wattage_per_cell(p, t)
    end

    for cell_id in heater_4_cells
        du.heat[cell_id] += hot_water_heater_4_wattage_per_cell(p, t)
    end

    for cell_id in heater_5_cells
        du.heat[cell_id] += hot_water_heater_5_wattage_per_cell(p, t)
    end
    
    trial_independent_solve_system!(du, u, p_vec, t, geo, system)
end

function build_trial_implicit_prob(f_closure, du0_vec, u0_vec, thermocouple_data, p_guess)
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
    )

    ode_func = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

    t0 = 0.0
    tMax = ustrip(upreferred(thermocouple_data.timestamps[end]))
    tspan = (t0, tMax)

    implicit_prob = ODEProblem(ode_func, u0_vec, tspan, p_guess)

    return implicit_prob
end

f_closure_dry_run = (du, u, p, t) -> fvm_operator!(du, u, p, t, dry_run_solve_system!, dry_run_geo, dry_run_system)
f_closure_hot_water = (du, u, p, t) -> fvm_operator!(du, u, p, t, hot_water_solve_system!, hot_water_geo, hot_water_system)

first_p_guess_init = ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 0.6660961855285039u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 0.2495764613764466u"W/K",
    heater_weight_1 = 0.07428260028432067,
    heater_weight_2 = 0.3499726682885098,
    heater_weight_3 = 0.26917982171867927,
    heater_weight_4 = 0.5781643096813001,
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 34.02609390801851u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 2.903109036597523,
    TC1_thermal_resistance = 68.52070429001392u"K/W",
    TC2_thermal_resistance = 195.16398786764938u"K/W",
    TC3_thermal_resistance = 129.78394261858463u"K/W",
    TC4_thermal_resistance = 157.8874124214126u"K/W",
    TC5_thermal_resistance = 21.669354950158237u"K/W",
)

dry_run_best_params = ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 1.3613733508973453u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 0.046224916624426536u"W/K",
    heater_weight_1 = 0.04233643564131222,
    heater_weight_2 = 0.6599848372698764,
    heater_weight_3 = 0.4480116180759421,
    heater_weight_4 = 0.007962092313654773,
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 1194.6515661535152u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 2.1703947052768693,
    TC1_thermal_resistance = 28.129313413301414u"K/W",
    TC2_thermal_resistance = 97.1810904949918u"K/W",
    TC3_thermal_resistance = 91.60608551471157u"K/W",
    TC4_thermal_resistance = 57.90696805772799u"K/W",
    TC5_thermal_resistance = 67.15056955503321u"K/W",
)

hot_water_best_params = ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 3.9215249544953634u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 0.21843778568123667u"W/K",
    heater_weight_1 = 0.3211331590638462,
    heater_weight_2 = 0.14535035339598748,
    heater_weight_3 = 0.6840217991331898,
    heater_weight_4 = 0.9984762651258827,
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 61.49737866475787u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 3.7753230723789706,
    TC1_thermal_resistance = 31.29201644638937u"K/W",
    TC2_thermal_resistance = 2.204398102203747u"K/W",
    TC3_thermal_resistance = 4.116298695255978u"K/W",
    TC4_thermal_resistance = 89.06507428995272u"K/W",
    TC5_thermal_resistance = 185.05463278797103u"K/W"
)

p_axes = getaxes(first_p_guess_init)
p_guess = ustrip.(upreferred.(Vector(first_p_guess_init)))

dry_run_p_best = ustrip.(upreferred.(Vector(dry_run_best_params)))
hot_water_p_best = ustrip.(upreferred.(Vector(hot_water_best_params)))

dry_run_implicit_prob = build_trial_implicit_prob(f_closure_dry_run, dry_run_du0_vec, dry_run_u0_vec, dry_run_thermocouple_data, dry_run_p_best);
hot_water_implicit_prob = build_trial_implicit_prob(f_closure_hot_water, hot_water_du0_vec, hot_water_u0_vec, hot_water_thermocouple_data, hot_water_p_best);

#NOTE!!: since we use an interpolation for heater wattage and nothing else happens in dry_run_saveat, 
#the solver will skip over times where the heater is on, resulting in no change
#this obviously isn't a problem for the hot water trial because the solver has to slow down because of all the physics that do not rely on interpolated values
desired_steps = 100
dry_run_saveat = ustrip(upreferred(dry_run_thermocouple_data.timestamps[end])) / desired_steps
hot_water_saveat = ustrip(upreferred(hot_water_thermocouple_data.timestamps[end])) / desired_steps

@time dry_run_test_sol = solve(dry_run_implicit_prob, FBDF(linsolve = SparspakFactorization()), dtmax = dry_run_saveat, callback = approximate_time_to_finish_cb)
@time hot_water_test_sol = solve(hot_water_implicit_prob, FBDF(linsolve = SparspakFactorization()), callback = approximate_time_to_finish_cb)

#=
dry_run_du_complete, dry_run_u_complete = regenerate_fvm_state(dry_run_test_sol, dry_run_system, dry_run_solve_system!, dry_run_geo, p_guess, u_additional_information = ComponentVector());
hot_water_du_complete, hot_water_u_complete = regenerate_fvm_state(hot_water_test_sol, hot_water_system, hot_water_solve_system!, hot_water_geo, p_guess, u_additional_information = ComponentVector());

root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files"

sol_to_vtk(dry_run_test_sol, dry_run_du_complete, dry_run_u_complete, grid, dry_run_geo, @__FILE__, root_dir; include_zeros_fields = false)
sol_to_vtk(hot_water_test_sol, hot_water_du_complete, hot_water_u_complete, grid, hot_water_geo, @__FILE__, root_dir; include_zeros_fields = false)
=#

#@time hot_water_krylov_sol = solve(hot_water_implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), callback = approximate_time_to_finish_cb)

dry_run_u_named = [ComponentVector(dry_run_test_sol.u[i], dry_run_system.state_axes) for i in 1:length(dry_run_test_sol.u)];
hot_water_u_named = [ComponentVector(hot_water_test_sol.u[i], hot_water_system.state_axes) for i in 1:length(hot_water_test_sol.u)];

sim_file = @__FILE__
dry_run_root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files\\dry_run_trial_output"
hot_water_root_dir = "C:\\Users\\wille\\Desktop\\Julia_cfd_output_files\\hot_water_flow_trial_output"

#sol_to_vtk(dry_run_test_sol, dry_run_u_named, grid, sim_file, dry_run_root_dir)
#sol_to_vtk(hot_water_test_sol, hot_water_u_named, grid, sim_file, hot_water_root_dir)

dry_run_timestamps = ustrip.(dry_run_thermocouple_data.timestamps)
hot_water_timestamps = ustrip.(hot_water_thermocouple_data.timestamps)

TC1_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC1_closest_cell_id))
TC2_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC2_closest_cell_id))
TC3_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC3_closest_cell_id))
TC4_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC4_closest_cell_id))
TC5_cell_id = Int(ustrip(thermocouple_and_jacket_properties.TC5_closest_cell_id))

n_saves = 100
dry_run_save_interval = ustrip(upreferred(dry_run_thermocouple_data.timestamps[end])) / n_saves
hot_water_save_interval = ustrip(upreferred(hot_water_thermocouple_data.timestamps[end])) / n_saves

dry_run_heater_powers = dry_run_thermocouple_data.heater_power_interp.(ustrip.(dry_run_thermocouple_data.timestamps))
dry_run_tstops = ustrip.(upreferred.(dry_run_thermocouple_data.timestamps[findall(!iszero, dry_run_heater_powers)]))

function get_trial_temperature_variance(thermocouple_data)
    n_thermocouples = 5
    all_recorded_temperatures = zeros(Float64, n_thermocouples, length(thermocouple_data.timestamps))

    for (time_index, t) in enumerate(ustrip.(thermocouple_data.timestamps))
        all_recorded_temperatures[1, time_index] = thermocouple_data.TC1_temps_interp(t)
        all_recorded_temperatures[2, time_index] = thermocouple_data.TC2_temps_interp(t)
        all_recorded_temperatures[3, time_index] = thermocouple_data.TC3_temps_interp(t)
        all_recorded_temperatures[4, time_index] = thermocouple_data.TC4_temps_interp(t)
        all_recorded_temperatures[5, time_index] = thermocouple_data.TC5_temps_interp(t)
    end

    mean_recorded_temperatures = sum(all_recorded_temperatures) / (n_thermocouples * length(thermocouple_data.timestamps))

    sum_of_squares = 0.0

    for i in eachindex(all_recorded_temperatures)
        sum_of_squares += (all_recorded_temperatures[i] - mean_recorded_temperatures)^2
    end

    return sum_of_squares / (n_thermocouples * length(thermocouple_data.timestamps))
end

dry_run_temperature_variance = get_trial_temperature_variance(dry_run_thermocouple_data)
hot_water_temperature_variance = get_trial_temperature_variance(hot_water_thermocouple_data)

function dry_run_loss(θ)
    prob = get!(task_local_storage(), :dry_run_implicit_prob) do
        # Deepcopy the mutable state and geometry objects
        system_copy = deepcopy(dry_run_system)
        geo_copy = deepcopy(dry_run_geo)
        
        # Build a new closure bound to the thread-isolated copies
        f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, dry_run_solve_system!, geo_copy, system_copy)
        
        # Re-generate the Jacobian sparsity and the implicit ODEProblem
        build_trial_implicit_prob(f_closure, dry_run_du0_vec, dry_run_u0_vec, dry_run_thermocouple_data, θ)
    end
    
    dry_run_loss_prob = remake(prob, p = θ)

    dry_run_sol = solve(
        dry_run_loss_prob, 
        FBDF(linsolve = SparspakFactorization(),), 
        sensealg = ForwardSensitivity(),
        saveat = dry_run_save_interval,
        tstops = dry_run_tstops,
    )

    if length(dry_run_sol.t) < n_saves + 1
        return nothing, nothing, 1e10
    end

    u_named = [ComponentVector(dry_run_sol.u[i], dry_run_system.state_axes) for i in eachindex(dry_run_sol.u)]

    mean_squared_error = 0.0

    for i in eachindex(dry_run_sol.t)
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC1_cell_id]) - ustrip(dry_run_thermocouple_data.TC1_temps_interp(dry_run_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC2_cell_id]) - ustrip(dry_run_thermocouple_data.TC2_temps_interp(dry_run_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC3_cell_id]) - ustrip(dry_run_thermocouple_data.TC3_temps_interp(dry_run_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC4_cell_id]) - ustrip(dry_run_thermocouple_data.TC4_temps_interp(dry_run_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC5_cell_id]) - ustrip(dry_run_thermocouple_data.TC5_temps_interp(dry_run_sol.t[i])))
        #mean_squared_error += abs2(ustrip(u_named[i].temp[end]) - ustrip(outlet_temp_interp(sol.t[i])))
        #I don't think using the measured outlet temperature is actually useful
    end

    loss = (mean_squared_error / dry_run_temperature_variance) / length(dry_run_sol.t)

    return dry_run_sol.t, u_named, loss
end

function viewable_dry_run_loss(θ)
    sol_t, u_named, loss = dry_run_loss(θ)

    simulated_thermocouple_data = (
        dry_run_times = sol_t,
        dry_run_TC1 = [u_named[i].temp[TC1_cell_id] for i in eachindex(sol_t)],
        dry_run_TC2 = [u_named[i].temp[TC2_cell_id] for i in eachindex(sol_t)],
        dry_run_TC3 = [u_named[i].temp[TC3_cell_id] for i in eachindex(sol_t)],
        dry_run_TC4 = [u_named[i].temp[TC4_cell_id] for i in eachindex(sol_t)],
        dry_run_TC5 = [u_named[i].temp[TC5_cell_id] for i in eachindex(sol_t)],
        loss = loss
    )

    return simulated_thermocouple_data
end

function pure_dry_run_loss(θ)
    sol_t, u_named, loss = dry_run_loss(θ)

    return loss 
end

function hot_water_loss(θ)
    prob = get!(task_local_storage(), :hot_water_implicit_prob) do
        # Deepcopy the mutable state and geometry objects
        system_copy = deepcopy(hot_water_system)
        geo_copy = deepcopy(hot_water_geo)
        
        # Build a new closure bound to the thread-isolated copies
        f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, hot_water_solve_system!, geo_copy, system_copy)
        
        # Re-generate the Jacobian sparsity and the implicit ODEProblem
        build_trial_implicit_prob(f_closure, hot_water_du0_vec, hot_water_u0_vec, hot_water_thermocouple_data, θ)
    end
    
    hot_water_loss_prob = remake(prob, p = θ)

    hot_water_sol = solve(
        hot_water_loss_prob, 
        FBDF(linsolve = SparspakFactorization(),), 
        sensealg = ForwardSensitivity(),
        saveat = hot_water_save_interval
    )

    # this is to prevent the optimizer from crashing the simulation to get a lower mean_squared_error
    if length(hot_water_sol.t) < n_saves + 1
        return nothing, nothing, 1e10
    end

    u_named = [ComponentVector(hot_water_sol.u[i], hot_water_system.state_axes) for i in eachindex(hot_water_sol.u)]

    mean_squared_error = 0.0

    for i in eachindex(hot_water_sol.t)
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC1_cell_id]) - ustrip(hot_water_thermocouple_data.TC1_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC2_cell_id]) - ustrip(hot_water_thermocouple_data.TC2_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC3_cell_id]) - ustrip(hot_water_thermocouple_data.TC3_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC4_cell_id]) - ustrip(hot_water_thermocouple_data.TC4_temps_interp(hot_water_sol.t[i])))
        mean_squared_error += abs2(ustrip(u_named[i].temp[TC5_cell_id]) - ustrip(hot_water_thermocouple_data.TC5_temps_interp(hot_water_sol.t[i])))
        #mean_squared_error += abs2(ustrip(u_named[i].temp[end]) - ustrip(outlet_temp_interp(sol.t[i])))
        #I don't think using the measured outlet temperature is actually useful
    end

    loss = (mean_squared_error / hot_water_temperature_variance) / length(hot_water_sol.t)

    return hot_water_sol.t, u_named, loss
end

function viewable_hot_water_loss(θ)
    sol_t, u_named, loss = hot_water_loss(θ)

    simulated_thermocouple_data = (
        hot_water_times = sol_t,
        hot_water_TC1 = [u_named[i].temp[TC1_cell_id] for i in eachindex(sol_t)],
        hot_water_TC2 = [u_named[i].temp[TC2_cell_id] for i in eachindex(sol_t)],
        hot_water_TC3 = [u_named[i].temp[TC3_cell_id] for i in eachindex(sol_t)],
        hot_water_TC4 = [u_named[i].temp[TC4_cell_id] for i in eachindex(sol_t)],
        hot_water_TC5 = [u_named[i].temp[TC5_cell_id] for i in eachindex(sol_t)],
        loss = loss
    )
    
    return simulated_thermocouple_data
end

function pure_hot_water_loss(θ)
    sol_t, u_named, loss = hot_water_loss(θ)

    return loss
end

function loss(θ)
    return (pure_dry_run_loss(θ) + pure_hot_water_loss(θ)) / 2.0
end

loss(p_guess)

hot_water_output_dir = joinpath(@__DIR__, "..", "graphs", "hot_water_graphs")
dry_run_output_dir = joinpath(@__DIR__, "..", "graphs", "dry_run_graphs")
combined_parameters_output_dir = joinpath(@__DIR__, "..", "graphs", "combined_parameters_loss_graphs")

function generate_dry_run_plots(dry_run_p_best, results_path)
    dry_run_simulated = viewable_dry_run_loss(dry_run_p_best)
    dry_run_times = dry_run_simulated.dry_run_times

    TC1_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC1, label="Sim TC1", linewidth=2)
    plot!(TC1_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC1_temps_interp.(dry_run_times)), label="Exp TC1", linewidth=2)
    savefig(TC1_plt, joinpath(results_path, "TC1_dry_run_loss.png"))

    TC2_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC2, label="Sim TC2", linewidth=2)
    plot!(TC2_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC2_temps_interp.(dry_run_times)), label="Exp TC2", linewidth=2)
    savefig(TC2_plt, joinpath(results_path, "TC2_dry_run_loss.png"))

    TC3_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC3, label="Sim TC3", linewidth=2)
    plot!(TC3_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC3_temps_interp.(dry_run_times)), label="Exp TC3", linewidth=2)
    savefig(TC3_plt, joinpath(results_path, "TC3_dry_run_loss.png"))

    TC4_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC4, label="Sim TC4", linewidth=2)
    plot!(TC4_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC4_temps_interp.(dry_run_times)), label="Exp TC4", linewidth=2)
    savefig(TC4_plt, joinpath(results_path, "TC4_dry_run_loss.png"))

    TC5_plt = plot(dry_run_times, dry_run_simulated.dry_run_TC5, label="Sim TC5", linewidth=2)
    plot!(TC5_plt, dry_run_times, ustrip(dry_run_thermocouple_data.TC5_temps_interp.(dry_run_times)), label="Exp TC5", linewidth=2)
    savefig(TC5_plt, joinpath(results_path, "TC5_dry_run_loss.png"))

    overall_plot = plot(TC1_plt, TC2_plt, TC3_plt, TC4_plt, TC5_plt, layout=(5, 1), size=(1000, 250*5), ylims=(250, 600))
    savefig(overall_plot, joinpath(results_path, "overall_dry_run_loss.png"))
end

function generate_hot_water_plots(hot_water_p_best, results_path)
    hot_water_simulated = viewable_hot_water_loss(hot_water_p_best)
    hot_water_times = hot_water_simulated.hot_water_times

    TC1_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC1, label="Sim TC1", linewidth=2)
    plot!(TC1_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC1_temps_interp.(hot_water_times)), label="Exp TC1", linewidth=2)
    savefig(TC1_hw_plt, joinpath(results_path, "TC1_hot_water_loss.png"))

    TC2_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC2, label="Sim TC2", linewidth=2)
    plot!(TC2_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC2_temps_interp.(hot_water_times)), label="Exp TC2", linewidth=2)
    savefig(TC2_hw_plt, joinpath(results_path, "TC2_hot_water_loss.png"))

    TC3_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC3, label="Sim TC3", linewidth=2)
    plot!(TC3_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC3_temps_interp.(hot_water_times)), label="Exp TC3", linewidth=2)
    savefig(TC3_hw_plt, joinpath(results_path, "TC3_hot_water_loss.png"))

    TC4_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC4, label="Sim TC4", linewidth=2)
    plot!(TC4_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC4_temps_interp.(hot_water_times)), label="Exp TC4", linewidth=2)
    savefig(TC4_hw_plt, joinpath(results_path, "TC4_hot_water_loss.png"))

    TC5_hw_plt = plot(hot_water_times, hot_water_simulated.hot_water_TC5, label="Sim TC5", linewidth=2)
    plot!(TC5_hw_plt, hot_water_times, ustrip(hot_water_thermocouple_data.TC5_temps_interp.(hot_water_times)), label="Exp TC5", linewidth=2)
    savefig(TC5_hw_plt, joinpath(results_path, "TC5_hot_water_loss.png"))

    overall_hw_plot = plot(TC1_hw_plt, TC2_hw_plt, TC3_hw_plt, TC4_hw_plt, TC5_hw_plt, layout=(5, 1), size=(1000, 250*5), ylims=(290, 340))
    savefig(overall_hw_plot, joinpath(results_path, "overall_hot_water_loss.png"))
end

function generate_combined_validation_plots(p_best, plot_output_dir)
    # Ensure plots directory exists
    mkpath(plot_output_dir)
    new_path = joinpath(plot_output_dir, Dates.format(now(), "yyyy-mm-dd_HH-MM-SS"))
    nested_dry_run_path = joinpath(new_path, "dry_run")
    nested_hot_water_path = joinpath(new_path, "hot_water")
    #parameters_path = joinpath(new_path, "parameters")
    mkpath(nested_dry_run_path)
    mkpath(nested_hot_water_path)
    #mkpath(parameters_path)

    CSV.write(joinpath(new_path, "parameters.csv"), DataFrame([NamedTuple(ComponentVector(p_best, p_axes))]))
    
    generate_dry_run_plots(p_best, nested_dry_run_path)
    generate_hot_water_plots(p_best, nested_hot_water_path)
end

generate_combined_validation_plots(p_guess, combined_parameters_output_dir)

#OPTIMIZATION
adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)

p_lower_bounds = ustrip.(upreferred.(Vector(ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 0.01u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 0.01u"W/K",
    heater_weight_1 = 0.05u"1",
    heater_weight_2 = 0.12u"1",
    heater_weight_3 = 0.20u"1",
    heater_weight_4 = 0.30u"1",
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 1.0u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 0.1,
    TC1_thermal_resistance = 0.1u"K/W",
    TC2_thermal_resistance = 0.1u"K/W",
    TC3_thermal_resistance = 0.1u"K/W",
    TC4_thermal_resistance = 0.1u"K/W",
    TC5_thermal_resistance = 0.1u"K/W",
))))

p_upper_bounds = ustrip.(upreferred.(Vector(ComponentVector(
    insulation_to_air_overall_heat_transfer_coefficient_to_environment = 100.0u"W/(m^2*K)",
    pipe_endcaps_to_air_thermal_conductance = 5.0u"W/K",
    heater_weight_1 = 0.3u"1",
    heater_weight_2 = 0.35u"1",
    heater_weight_3 = 0.50u"1",
    heater_weight_4 = 0.80u"1",
    fluid_to_steel_pipe_convective_heat_transfer_coefficient = 2000.0u"W/(m^2*K)",
    steel_thermal_mass_multiplier = 10.0,
    TC1_thermal_resistance = 200.0u"K/W",
    TC2_thermal_resistance = 200.0u"K/W",
    TC3_thermal_resistance = 200.0u"K/W",
    TC4_thermal_resistance = 200.0u"K/W",
    TC5_thermal_resistance = 200.0u"K/W",
))))

optprob = Optimization.OptimizationProblem(optf, p_guess, lb = p_lower_bounds, ub = p_upper_bounds)

function randomize(lower, upper)
    return lower + (upper - lower) * rand()
end

function randomize_list(lower, upper, length)
    return [lower + (upper - lower) * rand() for _ in 1:length]
end

p_ensemble = [[randomize(p_lower_bounds[i], p_upper_bounds[i]) for i in eachindex(p_lower_bounds)] for _ in 1:Sys.CPU_THREADS]

function prob_func(prob, i, repeat)
    return remake(prob, u0 = p_ensemble[i])
end

ensembleprob = EnsembleProblem(optprob; prob_func)

LOSS = Float64[]
PARS = []

# Ensure the directory exists and use a filename-safe date format (colons are invalid on Windows)
mkpath(joinpath(@__DIR__, "optimization_results"))
results_path = joinpath(@__DIR__, "optimization_results", "optimization_results_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).csv")

# Create the file and manually write the header string using propertynames
open(results_path, "w") do io
    header_str = "loss," * join(string.(propertynames(first_p_guess_init)), ",")
    println(io, header_str)
end

const cb_lock = ReentrantLock()

cb = function (state, l)
    display(l)
    display(state.u)
    
    lock(cb_lock) do
        push!(LOSS, l)
        push!(PARS, state.u)
        
        # Convert state.u to a named tuple using your p_axes so the CSV has nice column headers
        p_named = NamedTuple(ComponentVector(state.u, p_axes))
        row = merge((loss = l, ), p_named)
        
        # CSV.write with append=true automatically opens, appends, and closes (flushes) the file
        CSV.write(results_path, DataFrame([row]), append=true)
    end
    
    false
end

@time loss(p_guess)


#=@time res = Optimization.solve(
    ensembleprob,
    callback = cb,
    OptimizationOptimJL.LBFGS(),
    EnsembleThreads(),
    trajectories = 50,
    #Sys.CPU_THREADS,
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works kinda fine
    #f_abstol=1e-8,
    #g_abstol=1e-8,
)=#


@time res = Optimization.solve(
    optprob,
    callback = cb,
    OptimizationOptimJL.LBFGS(),
    #EnsembleThreads(),
    #trajectories = 50,
    #Sys.CPU_THREADS,
    #LBFGS, BFGS, and Fminbox don't work if the guess is very far away from the actual value
    #IPNewton works kinda fine
    #f_abstol=1e-8,
    #g_abstol=1e-8,
)

#=
@time res = solve(
    ensembleprob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(), 
    EnsembleThreads(),
    trajectories = Sys.CPU_THREADS,
    #PopulationSize = 100,
    #Method = :RandomSearcher,
    callback = cb,
)
=#
#=
@time res = solve(
    optprob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(),
    callback = cb,
    PoulationSize = 100,
    Method = :RandomSearcher
    #Method = :SepReal
)
    =#

optimization_results_file = joinpath(@__DIR__, "optimization_results", "optimization_results_2026-06-10_14-41-25.csv")
data = CSV.read(optimization_results_file, DataFrame)

n_best_solutions = 10

best_losses = sort(unique(data.loss))[1:n_best_solutions]

function find_best_parameters(best_losses, data, p_axes)
    corresponding_best_parameters = ComponentVector{Float64}[]

    for i in 1:n_best_solutions
        found_index = findlast(x -> x == best_losses[i], data.loss)
        best_parameters = collect(data[found_index, 2:end])
        @show best_parameters
        @show length(best_parameters)
        push!(corresponding_best_parameters, ComponentVector(best_parameters, p_axes))
    end
    
    return corresponding_best_parameters
end

corresponding_best_parameters = find_best_parameters(best_losses, data, p_axes)

corresponding_best_parameters[1]

#generate_combined_validation_plots(Vector(corresponding_best_parameters[1]), combined_parameters_output_dir)