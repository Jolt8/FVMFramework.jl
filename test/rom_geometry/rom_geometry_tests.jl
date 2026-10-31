#=
Tests of the ROM geometry loader and of cell/face sets on FVMGeometry.

Data (exported by the freecad-rom-geometry workbench):
  data/rom_fixture_export.json         its test fixture: cavity, liner, housing,
                                       end cap (split radially) and a nozzle
  data/rom_stepped_contact_export.json a plug seated in a cup: the two cells
                                       touch on a plane AND a cylinder

Run on its own with:
    julia --project=. test/rom_geometry/rom_geometry_tests.jl
=#

using Test
using FVMFramework
using ComponentArrays
using Unitful
using Ferrite

const ROM_TEST_DATA = joinpath(@__DIR__, "data")
const ROM_FIXTURE_PATH = joinpath(ROM_TEST_DATA, "rom_fixture_export.json")
const ROM_STEPPED_CONTACT_PATH = joinpath(ROM_TEST_DATA, "rom_stepped_contact_export.json")

struct ROMTestSolid <: AbstractPhysics end

no_property_update!(du, u, p, t, system, geo, cell_id) = nothing
no_region_function!(du, u, p, t, system, geo, cell_id) = nothing
no_patch_function!(du, u, p, t, system, geo, idx_a, face_a, idx_b, face_b) = nothing

function add_test_solid_region!(config, name, temperature, conductivity)
    add_region!(
        config, name;
        type = ROMTestSolid(),
        initial_conditions = ComponentVector(temp = temperature),
        properties = ComponentVector(k = conductivity, rho = 2700.0u"kg/m^3", cp = 900.0u"J/(kg*K)"),
        property_update_function = no_property_update!,
        region_function = no_region_function!,
    )
end

function total_face_count(facetset)
    face_count = 0
    for (idx_a, neighbor_list) in facetset
        face_count += length(neighbor_list)
    end
    return face_count
end

@testset "ROM geometry" begin

    @testset "fixture geometry matches the export" begin
        rom_data = read_rom_geometry_json(ROM_FIXTURE_PATH)
        geo = build_rom_fvm_geometry(rom_data)

        @test geo.n_cells == 11
        @test isempty(geo.node_coordinates)

        for rom_cell in rom_data["cells"]
            cell_id = rom_cell["id"]
            @test geo.cell_volumes[cell_id] ≈ rom_cell["volume"] rtol = 1e-14
            @test geo.cell_n_faces[cell_id] == length(rom_cell["faces"])
            for (face_idx, rom_face) in enumerate(rom_cell["faces"])
                @test geo.cell_face_areas[cell_id, face_idx] ≈ rom_face["area"] rtol = 1e-14
                @test geo.cell_face_distances[cell_id, face_idx] ≈ rom_face["half_cell_distance"] rtol = 1e-14
            end
        end

        # Every interface half-face points back at its partner.
        for (idx_a, neighbor_list) in geo.cell_neighbors
            for (idx_b, face_idx_a, face_idx_b) in neighbor_list
                @test (idx_a, face_idx_b, face_idx_a) in geo.cell_neighbors[idx_b][2]
                @test geo.cell_neighbor_distances[idx_a, face_idx_a] ≈ geo.cell_face_distances[idx_a, face_idx_a] + geo.cell_face_distances[idx_b, face_idx_b] rtol = 1e-12
            end
        end

        # Cell sets keep the model's region order.
        @test collect(keys(geo.cellsets)) == ["cavity", "liner", "housing", "end_cap_center", "end_cap_outer", "nozzle"]
        @test get_cellset(geo, "liner") == [4, 5, 6]
    end

    @testset "get_facetset by name and by region pair" begin
        geo = build_rom_fvm_geometry(read_rom_geometry_json(ROM_FIXTURE_PATH))

        cavity_to_liner = get_facetset(geo, "cavity", "liner")
        liner_to_cavity = get_facetset(geo, "liner", "cavity")
        @test [idx_a for (idx_a, neighbor_list) in cavity_to_liner] == [1, 2, 3]
        @test [idx_a for (idx_a, neighbor_list) in liner_to_cavity] == [4, 5, 6]

        # Same faces, seen from the other side.
        for (idx_a, neighbor_list) in cavity_to_liner
            for (idx_b, face_idx_a, face_idx_b) in neighbor_list
                liner_entry = only(neighbor_list_b for (idx_b_listed, neighbor_list_b) in liner_to_cavity if idx_b_listed == idx_b)
                @test (idx_a, face_idx_b, face_idx_a) in liner_entry
            end
        end

        # By name: a patch, and a stored region pair, give the same as above.
        @test get_facetset(geo, "cavity__liner") == cavity_to_liner
        @test all(idx_b == 0 for (idx_a, neighbor_list) in get_facetset(geo, "housing_ambient") for (idx_b, face_idx_a, face_idx_b) in neighbor_list)

        # Inside one region every axial face is listed once.
        @test total_face_count(get_facetset(geo, "liner", "liner")) == 2
        @test_throws ErrorException get_facetset(geo, "no_such_patch")
    end

    @testset "add_region! and add_patch! on a ROM config" begin
        rom_data = read_rom_geometry_json(ROM_FIXTURE_PATH)
        n_cells = rom_cell_count(rom_data)
        config = create_fvm_config_from_rom(rom_data, ComponentVector(temp = zeros(n_cells)u"K"))

        for (region_position, region_name) in enumerate(keys(config.geo.cellsets))
            add_test_solid_region!(config, region_name, (290.0 + region_position)u"K", 100.0u"W/(m*K)")
        end
        @test length(config.regions) == 6
        @test config.u_proto.temp[4] == 292.0u"K"

        housing_volume_in_export = sum(rom_cell["volume"] for rom_cell in rom_data["cells"] if rom_cell["region"] == "housing")
        @test get_cellset_volume(config, "housing") ≈ housing_volume_in_export rtol = 1e-14
        @test get_cellset_volume(config, "housing") ≈ π * (0.024^2 - 0.020^2) * 0.060 rtol = 1e-9

        # A boundary patch.
        add_patch!(config, "housing_ambient"; patch_function = no_patch_function!)
        @test get_facetset_area(config, "housing_ambient") ≈ 2π * 0.024 * 0.060 rtol = 1e-9

        # An interface set, asked for from the liner's side (stored the other way round).
        add_patch!(config, "liner__cavity"; owner_region = "liner", patch_function = no_patch_function!)
        liner_patch = config.patches[end]
        @test [idx_a for (idx_a, neighbor_list) in liner_patch.cell_neighbors] == [4, 5, 6]
        @test get_facetset_area(config, "liner__cavity") ≈ 2π * 0.018 * 0.060 rtol = 1e-9

        update_patch!(config, "liner__cavity")
        @test config.patches[end].cell_neighbors == liner_patch.cell_neighbors
    end

    @testset "add_virtual_cell!" begin
        rom_data = read_rom_geometry_json(ROM_FIXTURE_PATH)
        virtual_cell_id = add_virtual_cell!(
            rom_data;
            name = "recirculation_zone",
            volume = 2.0e-6,
            faces = [(neighbor_region = "cavity", neighbor_axial_index = :last, area = 1.0e-4, half_cell_distance = 1.0e-3, neighbor_half_cell_distance = 2.0e-3)],
        )
        geo = build_rom_fvm_geometry(rom_data)

        @test virtual_cell_id == 12
        @test geo.n_cells == 12
        @test geo.cell_n_faces[3] == 4  # the top cavity cell gained a face
        @test get_cellset(geo, "recirculation_zone") == [12]
        @test get_facetset(geo, "cavity__recirculation_zone") == [(3, [(12, 4, 1)])]
        @test get_facetset(geo, "recirculation_zone", "cavity") == [(12, [(3, 1, 4)])]
        @test geo.cell_neighbor_distances[12, 1] ≈ 3.0e-3
        @test geo.cell_volumes[12] == 2.0e-6
    end

    @testset "stepped contact: both faces conduct (finish_fvm_config + fvm_operator!)" begin
        rom_data = read_rom_geometry_json(ROM_STEPPED_CONTACT_PATH)
        config = create_fvm_config_from_rom(rom_data, ComponentVector(temp = zeros(2)u"K"))
        add_setup_syms!(
            config;
            cache_syms_and_units = (heat = u"J",),
            special_caches = ComponentVector(),
            optimized_parameters = ComponentVector(unused_parameter = 1.0),
        )

        cup_conductivity = 200.0
        plug_conductivity = 50.0
        add_test_solid_region!(config, "cup", 300.0u"K", cup_conductivity * u"W/(m*K)")
        add_test_solid_region!(config, "plug", 400.0u"K", plug_conductivity * u"W/(m*K)")

        solid_connection_map(type_a, type_b) = heat_diffusion_series_resistance!
        du0_vec, u0_vec, system, geo = finish_fvm_config(config, solid_connection_map, NamedTuple(); check_units = false)

        # The plug touches the cup on a plane AND a cylinder: two faces, two fluxes.
        @test total_face_count(only(system.connection_groups).cell_neighbors) == 2

        recorded_heat = zeros(2)
        function solve_and_record_heat!(du, u, p, t, system, geo)
            solve_connection_groups!(du, u, p, t, system, geo)
            recorded_heat[1] = du.heat[1]
            recorded_heat[2] = du.heat[2]
        end
        fvm_operator!(zero(du0_vec), u0_vec, system.p_vec, 0.0, system, geo, solve_and_record_heat!)

        # Hand calculation from the export: per face, two half-cell
        # resistances in series, R = d_cup / (k_cup A) + d_plug / (k_plug A).
        cup_cell, plug_cell = rom_data["cells"]
        expected_heat_into_cup = 0.0
        for cup_face in cup_cell["faces"]
            if cup_face["kind"] != "interface"
                continue
            end
            plug_face = plug_cell["faces"][cup_face["neighbor_local_idx"]]
            half_cell_resistance_cup = cup_face["half_cell_distance"] / (cup_conductivity * cup_face["area"])
            half_cell_resistance_plug = plug_face["half_cell_distance"] / (plug_conductivity * plug_face["area"])
            expected_heat_into_cup += (400.0 - 300.0) / (half_cell_resistance_cup + half_cell_resistance_plug)
        end

        @test recorded_heat[1] ≈ expected_heat_into_cup rtol = 1e-12
        @test recorded_heat[2] ≈ -expected_heat_into_cup rtol = 1e-12
    end

    @testset "series resistance equals heat_diffusion! for equal halves" begin
        geo = build_rom_fvm_geometry(read_rom_geometry_json(ROM_FIXTURE_PATH))
        # Cavity cells 1 and 2 meet on an axial plane, 10 mm from each centroid.
        idx_a, face_idx_a, face_idx_b = 1, 1, 1
        idx_b = geo.cell_neighbors[idx_a][2][face_idx_a][1]
        @test geo.cell_face_distances[idx_a, face_idx_a] ≈ geo.cell_face_distances[idx_b, face_idx_b]

        u = (k = [20.0, 80.0], temp = [300.0, 350.0])
        series_du = (heat = zeros(2),)
        harmonic_du = (heat = zeros(2),)
        heat_diffusion_series_resistance!(series_du, u, nothing, 0.0, nothing, geo, idx_a, face_idx_a, idx_b, face_idx_b)
        heat_diffusion!(harmonic_du, u, nothing, 0.0, nothing, geo, idx_a, face_idx_a, idx_b, face_idx_b)
        @test series_du.heat ≈ harmonic_du.heat rtol = 1e-12
        @test series_du.heat[1] ≈ -series_du.heat[2]
    end

    @testset "Ferrite mesh path: sets copied into geo, region pairs derived" begin
        grid = generate_grid(Hexahedron, (4, 1, 1))
        addcellset!(grid, "left_cells", x -> x[1] <= 1e-9)
        addcellset!(grid, "right_cells", x -> x[1] >= -1e-9)
        config = create_fvm_config(grid, ComponentVector(temp = zeros(4)u"K"))
        geo = config.geo

        @test get_cellset(geo, "left_cells") == [1, 2]
        @test only(get_facetset(geo, "left"))[1] == 1
        @test all(idx_b == 0 for (idx_a, neighbor_list) in get_facetset(geo, "left") for (idx_b, face_idx_a, face_idx_b) in neighbor_list)

        left_to_right = get_facetset(geo, "left_cells", "right_cells")
        right_to_left = get_facetset(geo, "right_cells", "left_cells")
        @test [(idx_a, [idx_b for (idx_b, face_idx_a, face_idx_b) in neighbor_list]) for (idx_a, neighbor_list) in left_to_right] == [(2, [3])]
        @test [(idx_a, [idx_b for (idx_b, face_idx_a, face_idx_b) in neighbor_list]) for (idx_a, neighbor_list) in right_to_left] == [(3, [2])]
        @test total_face_count(get_facetset(geo, "left_cells", "left_cells")) == 1

        add_patch!(config, "left"; patch_function = no_patch_function!)
        @test get_facetset_area(config, "left") ≈ 4.0  # the x = -1 end of the [-1, 1]^3 domain
        add_patch!(config, "left_cells__right_cells"; patch_function = no_patch_function!)
        @test total_face_count(config.patches[end].cell_neighbors) == 1
    end
end
