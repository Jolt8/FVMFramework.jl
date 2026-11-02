# Handoff: ROM Geometry Workbench for FreeCAD + FVMFramework loader

**Owner:** Will Martin
**Written:** 2026-10-30
**For:** Claude Code. Read all of this before you write any code.

---

## 1. Context

I'm building a reduced-order "digital twin" of a nitrous-oxide hybrid rocket engine (USST) with my Julia finite-volume framework, **FVMFramework.jl**. The model is a lumped thermal/flow network. Each node is a physical region, for example `injector_bulkhead_top`, `phenolic_liner_top` or `nozzle_top`. Some regions are split into 1..n axial cells along the fuel grain: `fluid_domain`, `fuel_grain`, `phenolic_liner` and `aluminum_housing`.

The geometry already exists as a parametric FreeCAD assembly. The hard part is turning that CAD into the volumes, interface areas and conduction distances the model needs. Splitting the engine into lumps is **my modelling decision, not a property of the CAD**. So I want a GUI tool where I define the regions myself, rather than a script that guesses them.

The key realization: **a reduced-order model is just a very coarse unstructured mesh whose cells have unusual shapes.** FVMFramework's `CellFaceLayout` already allows a different number of faces for each cell. If the workbench exports cells, half-faces, neighbours, areas, distances, cell sets and face sets in the same meaning FVMFramework already uses, then `add_region!`, `add_patch!`, `interface_geometry`, `boundary_geometry`, `heat_diffusion!` and the connection-group machinery all work **unchanged**. Do not redesign FVMFramework's internals.

### Files to read first

| Path | Why |
|---|---|
| `FVMFramework/examples/1D_hybrid_rocket_engine/CLAUDE.md` | My code-style rules. They apply to Python as well as Julia. |
| `FVMFramework/examples/1D_hybrid_rocket_engine/hybrid_digital_twin/hybrid_digital_twin.jl` | Draft model and region list. The `add_patch!` calls there are pseudocode; the real signature is in `src/setup/sim_config_patches.jl`. |
| `FVMFramework/src/geometry/geometry_building.jl` | `FVMGeometry` struct. The export must be able to fill every field. |
| `FVMFramework/src/geometry/cell_face_layout.jl` | Per-cell variable face counts (half-face storage). |
| `FVMFramework/src/setup/sim_config_base.jl`, `sim_config_regions.jl`, `sim_config_patches.jl`, `sim_config_finishing.jl` | How regions and patches look up `grid.cellsets`, `grid.facetsets` and `top.face_face_neighbor`. |
| `FVMFramework/src/geometry/get_x_geometry_functions/*.jl`, `src/physics/flux_physics/heat_transfer.jl` | How physics reads geometry. |
| `FVMFramework/examples/1D_hybrid_rocket_engine/fvm/meshes/gmsh_and_freecad_findings.md` | FreeCAD 1.1 API findings, including headless loading and object structure. |
| `FVMFramework/examples/1D_hybrid_rocket_engine/fvm/meshes/new_meshing_script.py` | A working way to re-run a script under FreeCAD's bundled Python. Reuse it. |
| `hybrid-rocket-CAD/src/hybrid_assembly_for_cfd.FCStd` | The real assembly. Read only; see the non-goals. |
| `hybrid-rocket-CAD/src/table_hybrid.FCStd` | Master parameter spreadsheet. Every key dimension has an alias. |

### Environment

- Windows. FreeCAD **1.1** at `C:\Program Files\FreeCAD 1.1\` (bundled Python **3.11**). Headless: `FreeCADCmd.exe`, or the re-exec trick in `new_meshing_script.py`.
- FreeCAD's internal length unit is mm. **Every exported value is SI** (m, m², m³, rad).
- Engine axis: global **Z** in the current assembly. The CFD slice spans roughly Z = −35 mm to +390 mm, and there is an `LCS_Origin` in the assembly. Make the axis configurable, with global Z as the default.
- Julia side: FVMFramework uses Ferrite, ComponentArrays, StaticArrays and OrderedCollections.

---

## 2. Goals and non-goals

### Goals

1. A FreeCAD workbench in a **new, separate git repository**. Suggested names: repo `freecad-rom-geometry`, workbench `ROMGeometry`; check with me before you finalize them. It lets me:
   - define **regions** as solids cut to a z-range (optional radial range), each split into N equal axial cells;
   - define **boundary patches** by picking CAD faces and naming them (e.g. `flat_plate_ambient_heat_exchange`);
   - preview regions and patches in colour in the 3D view;
   - review the interfaces between regions that were found automatically, and ignore any I don't want;
   - validate, then **export one JSON file** that FVMFramework can load.
2. All geometry logic lives in a **pure-Python core module with no GUI imports**, so the export also runs headless (`FreeCADCmd`) and can be tested.
3. A Julia loader in FVMFramework (phase 4) that turns the JSON into a `SimulationConfigInfo` so the existing `add_region!` and `add_patch!` work.

### Non-goals

- **Virtual regions** (e.g. `recirculation_zone`, which has no physical volume) are *not* defined in the workbench. I define them in the Julia script. The Julia loader must accept them (section 7).
- Do not define the hybrid engine's regions for me. I'll do that in the GUI. You may build a small fixture of your own.
- **Do not modify any file in `hybrid-rocket-CAD`.** Work on copies if you need to test against the real assembly.
- Don't touch the existing FVM/CFD example code or the meshing scripts.
- No regression or moving-geometry logic. The fuel port radius is a state in the Julia model, and patch functions compute the changing areas from `u`. The export describes the **initial** geometry only.

---

## 3. Repository layout (suggested)

```
freecad-rom-geometry/
  package.xml                    # Addon Manager metadata
  README.md                      # install = clone into %APPDATA%/FreeCAD/Mod (check the 1.1 path), usage
  CLAUDE.md                      # copy the style rules from my CLAUDE.md section 1
  Init.py
  InitGui.py                     # workbench class, toolbar, menu
  rom_geometry/
    core/                        # NO FreeCADGui imports anywhere in here
      regions.py                 # body ∩ slab → region solid → axial cell solids
      interfaces.py              # shared-face detection between cells, boundary faces
      distances.py               # half-cell distance rules (section 5.4)
      boundary_patches.py        # face fingerprints, re-matching, assignment
      parameters.py              # spreadsheet aliases → SI values
      validation.py              # all checks in section 6
      export_json.py             # writes the format in section 5
    objects/                     # FeaturePython proxies (data + recompute)
      rom_model.py
      rom_region.py
      rom_boundary_patch.py
    gui/                         # ViewProviders, commands, task panels
      commands.py
      region_task_panel.py
      boundary_patch_task_panel.py
      interface_review_panel.py
    resources/icons/
  scripts/
    export_headless.py           # FreeCADCmd export_headless.py model.FCStd out.json
  tests/
    build_fixture.py             # builds the test .FCStd from code (don't commit hand-made binaries)
    test_core.py                 # headless tests run under FreeCAD's Python
    expected/fixture_export.json
  docs/
    export_format.md             # the format spec; this repo owns it
```

---

## 4. Workbench objects and commands

### 4.1 `ROMModel` (one per document; a group container)

Properties:
- `Axis`: link to an LCS (optional; default global Z through the origin).
- `ParameterSheets`: link list to Spreadsheet objects (e.g. the sheet inside `table_hybrid.FCStd`). All their aliases are exported.
- `ExportPath`: file path.
- `FuzzyTolerance`: length, default about 0.01 mm. Used for the boolean operations.
- `NearContactTolerance`: length, default about 0.5 mm. Used for gap warnings.
- `IgnoredInterfaces`: string list of `"region_a|region_b"` pairs I've chosen to ignore.
- The **order of regions inside the group** is meaningful: it sets cell numbering and the owner direction of interface face sets (section 5.3).

### 4.2 `ROMRegion` (Part::FeaturePython so it renders)

Properties:
- `Name` (string, the name the model uses; must be unique), `Material` (string, free text, e.g. `aluminum_6061`).
- `SourceShapes`: link list of bodies, links or solids. Fluid volumes are ordinary solids too (section 4.4).
- `ZMin`, `ZMax`: Length, measured along `ROMModel.Axis`. **These must accept FreeCAD expressions** so I can bind them to spreadsheet aliases, e.g. `=table_hybrid#Spreadsheet.pre_combustion_chamber_length`. Check that cross-document expressions work in 1.1 and tell me if they don't.
- `RMin`, `RMax`: optional Length (0 / empty = unbounded).
- `NAxialCells`: integer ≥ 1. Equal-length z-slabs. Cells are numbered in increasing z.
- `Color`.
- `Shape`: the computed region solid (the union of its cells), shown semi-transparent.

On recompute: fuse `SourceShapes`, intersect with the slab (and the radial annulus if one is set), then split into `NAxialCells` cell solids. Warn if the result is empty or has more than one disconnected solid.

### 4.3 `ROMBoundaryPatch`

Properties:
- `Name` (e.g. `flat_plate_ambient_heat_exchange`), `Color`.
- `Faces`: LinkSubList to the CAD faces I picked.
- `Fingerprints`: for each picked face, store its centre of mass, area, surface type, and normal at the centre (plus radius and axis for cylinders). Stored as a JSON string property.

**Topological naming:** the `FaceN` indices change when the CAD is edited upstream. On every recompute or export, re-match each picked face by its fingerprint, not by index. If a match is missing or ambiguous, **stop with a clear error** that names the patch and the face. Do not silently use a different face. Provide a "re-pick" action in the task panel.

### 4.4 Commands (toolbar + menu)

1. **Create ROM Model**
2. **Add Region**: from the selected solid(s); opens the task panel (name, material, z/r limits, N cells).
3. **Make Cavity Solid**: helper for fluid regions such as `fluid_domain`, `pre_chamber` and `mid_section`. Select an enclosing solid (or create a cylinder on the axis), subtract the selected solids, and produce a normal Part feature I can then use as a region source. Fluid volumes are the empty space inside the engine, not CAD bodies, so I need this.
4. **Add Boundary Patch**: from the selected faces.
5. **Review Interfaces**: a panel listing every detected region-pair interface with its total area. Each has an ignore checkbox; clicking a row highlights it. It also lists near-contacts (gaps below `NearContactTolerance`) with their gap size, and unassigned outer surfaces.
6. **Validate**: runs section 6 and shows the report.
7. **Export**: validate, then write the JSON. Refuse to write if there are errors; warnings are allowed.

---

## 5. Geometry algorithms and export format

### 5.1 Cells and interfaces

- Build every cell solid of every region.
- Run **one `generalFuse`** over all cell solids together (`Part.Shape.generalFuse(shapes, fuzzy)` or BOPTools), using `FuzzyTolerance`. In the fused result:
  - a sub-face shared by two cells is an **interface half-face pair**;
  - a sub-face belonging to only one cell is a **boundary half-face** of that cell.
- Merge sub-faces between the same two cells into **one face per neighbouring cell** where they are the same kind (all planar, or all coaxial cylinders of the same radius). If merging would mix a flat face and a cylindrical one, keep them as separate faces and explain why in the export (`notes`). The aim is for each cell to have about one face per neighbour, because physics functions are written per face.
- Interfaces listed in `IgnoredInterfaces` are exported as **boundary faces with no patch** and `kind = "ignored_interface"`, so their area is still accounted for.
- **Pairwise overlap check:** the `common()` volume between any two cells of *different* regions must be about 0. If it isn't, it's an error (interfering CAD, e.g. threads or press fits). Report the overlap volume.
- **Near-contacts:** pairs of cells whose minimum `distToShape` is below `NearContactTolerance` but which share no face. Report them as warnings with the gap size. A gap is a physical thermal resistance, so I decide what to do about it.

### 5.2 Boundary patch assignment

For each boundary half-face, it belongs to patch P if it lies on one of P's re-matched faces: its centroid is on that surface within tolerance, and the `common` area is at least 99% of the sub-face area. A half-face may belong to **at most one** patch (error otherwise). Boundary half-faces with no patch go in `unassigned_boundary_faces`, with area and centroid, so adiabatic surfaces are a choice I can see.

### 5.3 Face-set naming and direction (this matters)

`add_patch!` loops over every `(cell, face)` in a face set and looks up the neighbour. **If a face set contains both halves of an interface, the flux is applied twice.** So:

- Region-pair interface set: name `"<region_a>__<region_b>"`. It contains **only region_a's half-faces**, where region_a is the region that comes *earlier* in the `ROMModel` group order.
- Internal axial faces within one region: name `"<region>__<region>"`, containing each internal face once (the lower-z cell's half-face).
- Boundary patch set: the patch name, containing that patch's boundary half-faces.
- Cell sets: one per region name, with cell ids in increasing z.

### 5.4 Half-cell distances (`cell_face_distances`)

These set the conduction resistance, so get them right and document them in `docs/export_format.md`:

- **Planar faces:** `|dot(face_centroid − cell_centroid, face_normal)|`. This is the same definition FVMFramework uses.
- **Coaxial cylindrical faces** (radial conduction through shells): use the distance that makes a flat-wall conduction formula exact for a cylindrical shell:
  `d_eq = r_f * |ln(r_f / r_c)|`
  where `r_f` is the face radius and `r_c` is the cell's radial node position. Default `r_c = sqrt(r_inner * r_outer)` (geometric mean, which makes the inner and outer half-resistances equal). Use the cell's own inner and outer radii. If the cell is a solid cylinder (`r_inner = 0`), fall back to `r_c = r_outer / 2` and say so in `notes`. Make the choice of `r_c` a ROMModel property so I can change it.
- **Any other surface** (cones in the nozzle, fillets): use the normal distance at the face's centroid, and set `distance_method = "approximate"` so I can find those faces.
- `cell_neighbor_distances[a, face] = d_a + d_b` (the sum of the two half-cell distances). **Also** export the straight centroid-to-centroid distance as `centroid_distance` for reference.
- **Normals:** planar faces get their unit normal pointing out of the cell. A full 360° coaxial cylindrical face has no single normal (its centroid is on the axis), so export `normal: null`, `surface: "cylinder"`, `radius`, and `outward: "+r"` or `"-r"`.

### 5.5 Export JSON (format_version 1)

All ids are **1-based** (Julia). All values are SI.

```json
{
  "format_version": 1,
  "generator": {"name": "ROMGeometry", "version": "0.1.0", "freecad_version": "1.1.x"},
  "source_document": "C:/.../hybrid_assembly_for_cfd.FCStd",
  "exported_at": "2026-10-30T22:40:00-06:00",
  "axis": {"origin": [0.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]},
  "radial_node_rule": "geometric_mean",

  "parameters": {
    "phenolic_liner_inner_diameter": {"value": 0.03332, "unit": "m"},
    "converging_half_angle": {"value": 0.5236, "unit": "rad"},
    "flange_n_holes": {"value": 8, "unit": ""}
  },

  "regions": [
    {"name": "aluminum_housing_top", "material": "aluminum_6061",
     "source_shapes": ["aluminum_housing"],
     "z_min": 0.0, "z_max": 0.015, "r_min": null, "r_max": null,
     "n_axial_cells": 1, "cell_ids": [1],
     "z_min_expression": "0 mm", "z_max_expression": "table_hybrid#Spreadsheet.pre_combustion_chamber_length"}
  ],

  "cells": [
    {"id": 1, "region": "aluminum_housing_top", "axial_index": 1,
     "volume": 1.23e-6, "centroid": [0.0, 0.0, 0.0075],
     "z_min": 0.0, "z_max": 0.015, "r_inner": 0.017465, "r_outer": 0.020965,
     "faces": [
       {"local_idx": 1, "kind": "interface",
        "neighbor_cell": 7, "neighbor_local_idx": 3,
        "area": 1.6e-3, "surface": "cylinder", "radius": 0.017465, "outward": "-r",
        "normal": null, "centroid": [0.0, 0.0, 0.0075],
        "half_cell_distance": 0.00088, "neighbor_distance": 0.0017,
        "centroid_distance": 0.0, "distance_method": "cylindrical_log",
        "notes": ""},
       {"local_idx": 2, "kind": "boundary",
        "neighbor_cell": 0, "neighbor_local_idx": 0,
        "patch": "flat_plate_ambient_heat_exchange",
        "area": 1.97e-3, "surface": "cylinder", "radius": 0.020965, "outward": "+r",
        "normal": null, "centroid": [0.0, 0.0, 0.0075],
        "half_cell_distance": 0.00088, "distance_method": "cylindrical_log",
        "notes": ""}
     ]}
  ],

  "cellsets": {"aluminum_housing_top": [1]},
  "facetsets": {
    "aluminum_housing_top__phenolic_liner_top": [[1, 1]],
    "flat_plate_ambient_heat_exchange": [[1, 2]]
  },

  "unassigned_boundary_faces": [{"cell": 1, "local_idx": 3, "area": 4.1e-4, "centroid": [0, 0, 0]}],
  "near_contacts": [{"region_a": "...", "region_b": "...", "gap": 0.0003}],
  "validation": {"errors": [], "warnings": ["..."]}
}
```

`parameters`: every alias from every linked spreadsheet, evaluated by FreeCAD after a recompute. Lengths → m, areas → m², volumes → m³, angles → rad; dimensionless values stay as they are. Record the SI unit string with each value.

---

## 6. Validation (run before every export; also a separate command)

Errors (export refused):
1. Each region's cell volumes sum to its region solid volume (relative tolerance 1e-6).
2. Each interface pair is reciprocal: `neighbor_cell` / `neighbor_local_idx` point back, and `area_a ≈ area_b` (relative 1e-6).
3. No overlap volume between cells of different regions (above a tolerance relative to the smaller cell).
4. Region names are unique, and no region is empty.
5. Each boundary half-face belongs to at most one patch. Every picked patch face is re-matched uniquely.
6. In every face set, no interface appears from both sides.
7. Closure: for each cell, the sum of its face area vectors ≈ 0. Planar faces only; skip cylinders, or use the exact cylinder contribution.

Warnings:
- For each source body: the total volume of the regions that use it vs the body's volume. Shows parts of a body I left out (fine if intentional; report the percentage).
- Near-contacts, unassigned boundary area (per region, in m²), faces marked `approximate`, and any face merges that were refused.

---

## 7. Julia side: FVMFramework loader (phase 4)

Put this in FVMFramework (e.g. `src/geometry/rom_geometry_loading.jl`), not in the workbench repo.

```julia
rom_data = read_rom_geometry_json(path)            # checks format_version

add_virtual_cell!(rom_data;
    name = "recirculation_zone",
    volume = recirculation_zone_volume,
    faces = [(neighbor_region = "fluid_domain", neighbor_axial_index = :last,
              area = ..., half_cell_distance = ..., neighbor_half_cell_distance = ...)]
)
# Virtual cells must be added BEFORE geo is built: CellFaceLayout offsets are fixed.
# Adding a virtual cell also adds the matching half-face to the real neighbour cell,
# and creates a cell set (its name) and a face set "<real_region>__<virtual_name>".

config = create_fvm_config_from_rom(rom_data, u_proto)
```

- Fill `FVMGeometry` directly from the JSON. Node-based fields (`cell_node_ids`, `cell_faces_node_ids`, `node_coordinates`, `cell_neighbors_node_ids`) can be empty or placeholders; check that nothing in the ROM path reads them. `null` normals become zero vectors. **Don't** call `update_fvm_geometry!` on ROM geometry, because it would recompute everything from nodes.
- **`SimulationConfigInfo` currently declares `grid::Ferrite.Grid` and `top::ExclusiveTopology`.** A Ferrite cell type has a fixed number of faces, so a lumped region can't be one. Proposed minimal change: add `ROMGrid` (fields `cells` or `n_cells`, `cellsets::Dict{String, OrderedSet{Int}}`, `facetsets::Dict{String, OrderedSet{Ferrite.FacetIndex}}`) and `ROMTopology` (`face_face_neighbor`, indexable as `[cell_id, face_idx]` and returning a collection whose `first(...).idx == (neighbor_cell, neighbor_face)`, empty for boundary faces). Widen the two field types (e.g. to a Union, or drop the annotations), and add a `getcellset` method for `ROMGrid`. Look at every use of `config.grid`, `config.top` and `getcellset` in `src/` first, then **show me the exact diff and wait for my OK** before changing `SimulationConfigInfo`. If there's a less invasive way, propose it instead.
- Round-trip test: load `tests/expected/fixture_export.json` (copy it into FVMFramework's tests), build a config, run `add_region!` on every region and `add_patch!` on one interface set and one boundary set, and check total volumes and face-set areas against the JSON.
- **Optional; ask me first.** `heat_diffusion!` uses the harmonic mean of k with the total centroid distance. That is exact only when both half-cell distances are equal, which is badly wrong for, say, a thin liner against a thick wall. A per-side series-resistance variant (a *new* function; don't change the existing one):

  ```julia
  half_cell_resistance_a = face_distance_a / (u.k[idx_a] * face_area_a)

  half_cell_resistance_b = face_distance_b / (u.k[idx_b] * face_area_b)

  conduction_heat_flow = (u.temp[idx_b] - u.temp[idx_a]) / (half_cell_resistance_a + half_cell_resistance_b)

  du.heat[idx_a] += conduction_heat_flow
  du.heat[idx_b] -= conduction_heat_flow
  ```

---

## 8. Phases and acceptance criteria

Finish each phase, run its tests, and show me the result before starting the next.

**Phase 1: headless core + fixture**
- `tests/build_fixture.py` builds a small `.FCStd` from code: a stack of coaxial cylinders and annuli imitating a liner inside a housing, with a solid end cap, a cavity, a deliberate 0.2 mm gap somewhere, and one conical face. It also includes a tiny spreadsheet whose aliases drive at least one region's z-limits.
- The core produces cells, interfaces, distances and the JSON from code-defined region specs (no GUI objects yet).
- Tests cover: analytic volumes and areas of the cylinders to 1e-9 relative; `d_eq` against the hand formula; every check in section 6 passing; the gap reported as a near-contact; the cone marked `approximate`; and a deliberately overlapping case raising the overlap error.

**Phase 2: FeaturePython objects + GUI preview**
- `ROMModel`, `ROMRegion`, task panel, coloured preview, save/reopen survives, expression binding to spreadsheet aliases works.
- `scripts/export_headless.py` exports a document containing these objects without the GUI.

**Phase 3: boundary patches, interface review, validation UI**
- Face picking and fingerprints. Test: edit the fixture so face indices shift, and the patch still re-matches. Make a face ambiguous and it fails loudly.
- Interface review panel with ignore toggles; validate and export commands.

**Phase 4: Julia loader in FVMFramework** (section 7), including the round-trip test.

---

## 9. Code style (from my CLAUDE.md; applies to Python and Julia)

- Engineering first and physically explicit. Keep the main geometric or physical calculation visible in the function that does it. Avoid deep helper chains.
- Descriptive, domain-specific names (`cell_face_half_cell_distance`, not `d`). Long names are fine.
- Explicit intermediate variables and plain loops rather than compact idioms.
- Comments explain geometric assumptions and why (e.g. why `r_c` is the geometric mean).
- Keep assignments on one line. Moderate indentation, and blank lines between separate calculations. Don't wrap `x =` onto a new line.
- Readability beats clever architecture. A few plain modules beat a deep class hierarchy.

---

## 10. Ask me (don't guess) about

- Final repo and workbench names.
- The `SimulationConfigInfo` change (section 7). Show the diff first.
- Whether to add the series-resistance heat flux function.
- Anything in my real assembly that doesn't fit this design, e.g. regions that need a cut other than z/r, or bodies that interfere.
