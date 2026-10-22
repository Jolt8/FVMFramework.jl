# FreeCAD API and Gmsh Integration Findings for Hybrid Rocket Engine CFD

## 1. Executive Summary

This report documents the architectural, geometric, and programmatic findings from investigating the **FreeCAD 1.1 Python API** and its integration with **Gmsh** for generating finite-volume computational fluid dynamics (CFD) meshes.

The automated workflow developed in [`meshing_script.py`](file:///c:/Users/wille/Desktop/FVMFramework/examples/1D_hybrid_rocket_engine/meshes/meshing_script.py):
1. Solves the Python interpreter mismatch between system Python (3.14/3.12) and FreeCAD 1.1's embedded Python (3.11) via transparent process re-execution.
2. Headlessly loads the full rocket CAD assembly [`hybrid_assembly_for_cfd.FCStd`](file:///C:/Users/wille/Desktop/hybrid-rocket-CAD/src/hybrid_assembly_for_cfd.FCStd).
3. Inspects and extracts geometric and topological metadata (names, labels, placements, shapes, bounding boxes, areas, and volumes).
4. Matches the 2D cross-section slice faces (`Exploded Slice`) to their corresponding 3D assembly solid bodies via boolean intersection (`face.common(body.Shape).Area`), automatically segregating structural bodies from the CFD fluid chamber.
5. Exports the sliced faces into a standardized STEP file (`exploded_slice.step`) while **preserving entity labels** via FreeCAD's `Import.export` (Open CASCADE XDE).
6. Imports the STEP model directly into Gmsh, maps the embedded labels to **Gmsh 2D Physical Groups**, recombines triangular elements into high-quality quadrilaterals via Blossom recombination, scales dimensions to SI meters, and generates the final mesh (`hybrid_rocket_cfd.msh`) and a companion `.geo` viewing script.

---

## 2. FreeCAD API Findings

### 2.1 Headless Document Loading
FreeCAD assemblies are opened headlessly using:
```python
import FreeCAD as App

doc = App.openDocument(path_to_cfd_assembly)
```
- `App.openDocument(filepath)` deserializes `.FCStd` files into memory without requiring any GUI or display server.
- The document object provides `doc.Objects` (flat list of all features), `doc.RootObjects` (top-level tree nodes), and direct dictionary-style indexing `doc.getObject("Name")` or `doc.Name`.

### 2.2 Object Hierarchy & Naming Conventions
Every object in FreeCAD has two fundamental identifier properties:
| Property | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `obj.Name` | String | Immutable, unique internal document ID | `'Body'`, `'Body004'`, `'Slice_child0'` |
| `obj.Label` | String | User-defined, mutable display label shown in the tree | `'aluminum_housing'`, `'fuel_grain'` |
| `obj.TypeId` | String | C++ type identifier | `'App::Link'`, `'Part::Compound'`, `'App::Part'` |

### 2.3 Assembly Structure in `hybrid_assembly_for_cfd.FCStd`
The CFD assembly document contains **37 document objects**:
- **Assembly Component Links (`App::Link`)**:
  - `Body` $\rightarrow$ `aluminum_housing` (linked from `aluminum_housing.FCStd`)
  - `Body001` $\rightarrow$ `injector_bulkhead` (linked from `injector_assembly.FCStd`)
  - `Body002` $\rightarrow$ `injector_disk` (linked from `injector_disk.FCStd`)
  - `ASM0002_ASM` $\rightarrow$ `ASM0002_ASM` (linked from `fitting_with_lcs.FCStd`)
  - `Body003` $\rightarrow$ `phenolic_liner` (linked from `phenolic_liner.FCStd`)
  - `Body004` $\rightarrow$ `fuel_grain` (linked from `fuel_grain.FCStd`)
  - `Body005` $\rightarrow$ `inner_graphite_nozzle` (linked from `inner_graphite_nozzle.FCStd`)
  - `Body006` $\rightarrow$ `outer_aluminum_nozzle_casing` (linked from `outer_aluminum_nozzle_casing.FCStd`)
- **Containers**: `hybrid_assembly` (`App::Part`), `Parts` (`App::DocumentObjectGroup`), `Constraints`, `Configurations`.
- **CFD Cross-Section Slicing Pipeline**:
  - `Compound` (`Part::Compound`): Groups all solid bodies.
  - `Sketch001` & `Face` (`Part::Face`): A large bounding planar rectangle ($9,454.33\text{ mm}^2$) covering the longitudinal half of the engine from $Z = -35\text{ mm}$ to $Z = +390\text{ mm}$.
  - `Slice` (`Part::FeaturePython`, BOPTools): Slices the planar cutting face against the `Compound` solids.
  - `GrExplode_Slice` (`App::DocumentObjectGroup`): Explodes `Slice` into individual topological faces (`Slice.0` through `Slice.8`).

### 2.4 Queryable Geometric & Physical Properties
For any object `obj` with a `.Shape` attribute:
- **Topology**:
  - `obj.Shape.ShapeType`: `'Solid'`, `'Compound'`, `'Face'`, `'Edge'`, `'Vertex'`
  - `obj.Shape.Solids`, `Faces`, `Edges`, `Vertexes`: Iterators over topological sub-elements
  - `obj.Shape.isValid()`: Boolean verifying topological integrity
- **Measurements**:
  - `obj.Shape.Volume`: Volume in $\text{mm}^3$ (applicable to solids and closed shells)
  - `obj.Shape.Area`: Surface area in $\text{mm}^2$
  - `obj.Shape.BoundBox`: Min/max extents (`XMin, XMax, YMin, YMax, ZMin, ZMax`)
  - `obj.Shape.CenterOfMass`: Center of mass vector $(x, y, z)$
- **Coordinate Transformations**:
  - `obj.Placement`: Local position vector `Base` and rotation quaternion / Euler angles
  - `obj.getGlobalPlacement()`: Absolute transformation matrix through the assembly tree
- **API Gotcha**:
  - `obj.Shape.CenterOfMass` is defined on `Part.Solid` and `Part.Face`, but raises an `AttributeError` on raw `Part.Compound` objects. Use `getattr(shape, "CenterOfMass", None)` for safe defensive querying.

---

## 3. Exploded Slice Mapping & Analysis

The 2D longitudinal cutting face has an initial area of **$9,454.33\text{ mm}^2$**. The exploded slice partitions this area into **9 disjoint faces**:

| Slice Child | Label | Area ($\text{mm}^2$) | Extent $Z$ ($\text{mm}$) | Extent $X$ ($\text{mm}$) | Matched Physical Region | Notes |
| :--- | :--- | :---: | :---: | :---: | :--- | :--- |
| `Slice_child0` | `Slice.0` | 5,297.60 | $[-35.0, 390.0]$ | $[0.0, 36.9]$ | `fluid_domain` | Core chamber, grain port, nozzle throat, plume exit |
| `Slice_child1` | `Slice.1` | 67.18 | $[-5.0, 10.0]$ | $[14.8, 23.8]$ | `outer_aluminum_nozzle_casing` | Aft nozzle retaining casing |
| `Slice_child2` | `Slice.2` | 1,394.10 | $[0.0, 390.0]$ | $[17.2, 23.8]$ | `aluminum_housing` | Outer motor casing pressure vessel |
| `Slice_child3` | `Slice.3` | 234.95 | $[0.0, 41.5]$ | $[6.9, 17.2]$ | `inner_graphite_nozzle` | Converging-diverging graphite throat insert |
| `Slice_child4` | `Slice.4` | 13.80 | $[32.9, 50.0]$ | $[16.4, 17.2]$ | `phenolic_liner` | Aft liner section before fuel grain |
| `Slice_child5` | `Slice.5` | 249.55 | $[50.0, 360.0]$ | $[16.4, 17.2]$ | `phenolic_liner` | Liner section surrounding fuel grain |
| `Slice_child6` | `Slice.6` | 24.15 | $[360.0, 390.0]$ | $[16.4, 17.2]$ | `phenolic_liner` | Forward liner section before injector bulkhead |
| `Slice_child7` | `Slice.7` | 2,170.00 | $[50.0, 360.0]$ | $[9.5, 16.4]$ | `fuel_grain` | Solid fuel port annulus ($r_\text{in} \approx 9.5\text{ mm}$, $r_\text{out} \approx 16.4\text{ mm}$) |
| `Slice_child8` | `Slice.8` | 3.00 | $[388.5, 390.0]$ | $[19.2, 21.1]$ | `fluid_domain` | Bulkhead seal groove / clearance void |
| **Sum** | | **9,454.33** | | | | **Exact conservation of sliced area** |

### Key Observations:
1. **Fluid Domain Identification**: `Slice.0` and `Slice.8` have **0.0 mm² overlap** with any solid body. They represent the interior gas passage and voids. In particular, `Slice.0` spans from $X = 0$ (the axis of symmetry) across the pre-combustion chamber, through the fuel port, and out the supersonic nozzle expansion.
2. **Phenolic Liner Partitioning**: The phenolic liner is partitioned into 3 separate faces (`Slice.4`, `Slice.5`, `Slice.6`) because the solid bodies and step cuts divide the cross-section longitudinally. In Gmsh, these are unified under a single physical group `phenolic_liner`.

---

## 4. Bridging FreeCAD Labels to Gmsh Physical Groups

### 4.1 The Metadata Loss Problem with Standard STEP Export
- Calling `Part.export(objects, filename)` exports raw Open CASCADE geometric surfaces without XDE metadata.
- When imported into Gmsh (`gmsh.open()`), `gmsh.model.getEntityName(dim, tag)` returns empty strings `""`, forcing error-prone manual tag mapping.

### 4.2 Solution: FreeCAD `Import.export` (Open CASCADE XDE)
- FreeCAD includes a high-level STEP translator module: `import Import`.
- Calling `Import.export(objects, filename)` builds an Open CASCADE XDE document structure and embeds each object's `Label` directly into the STEP AP214/AP242 file.
- When Gmsh opens this file, Gmsh reads the embedded labels:
  ```
  Shapes/hybrid_assembly_for_cfd/fluid_domain/fluid_domain
  Shapes/hybrid_assembly_for_cfd/fuel_grain001/fuel_grain001
  Shapes/hybrid_assembly_for_cfd/inner_graphite_nozzle001/inner_graphite_nozzle001
  Shapes/hybrid_assembly_for_cfd/aluminum_housing001/aluminum_housing001
  Shapes/hybrid_assembly_for_cfd/outer_aluminum_nozzle_casing001/outer_aluminum_nozzle_casing001
  Shapes/hybrid_assembly_for_cfd/phenolic_liner001/phenolic_liner001
  ```

### 4.3 Automated Physical Group Assignment in Python
In [`meshing_script.py`](file:///c:/Users/wille/Desktop/FVMFramework/examples/1D_hybrid_rocket_engine/meshes/meshing_script.py), slice faces are matched against bodies using boolean overlap:
```python
slice_face_body_names = {}
for child in exploded_slice.Group:
    face = child.Shape
    matched_name = "fluid_domain"
    for body in solid_bodies:
        if face.common(body.Shape).Area > 0.5 * face.Area:
            matched_name = body.Label
            break
    slice_face_body_names[child.Name] = matched_name
```
After exporting to STEP and importing into Gmsh:
```python
physical_groups = {}
for dim, tag in gmsh.model.getEntities(2):
    raw_name = gmsh.model.getEntityName(dim, tag)
    clean_name = raw_name.split("/")[-1].rstrip("0123456789")
    physical_groups.setdefault(clean_name, []).append(tag)

for group_name, tags in physical_groups.items():
    gmsh.model.addPhysicalGroup(2, tags, name=group_name)
```
This produces the physical group structure:
- `fluid_domain` $\rightarrow$ Surface tags `[1, 9]`
- `outer_aluminum_nozzle_casing` $\rightarrow$ Surface tag `[2]`
- `aluminum_housing` $\rightarrow$ Surface tag `[3]`
- `inner_graphite_nozzle` $\rightarrow$ Surface tag `[4]`
- `phenolic_liner` $\rightarrow$ Surface tags `[5, 6, 7]`
- `fuel_grain` $\rightarrow$ Surface tag `[8]`

---

## 5. Gmsh Meshing Configuration

### 5.1 Unit Scaling (SI Meters)
FreeCAD operates in **millimeters**, while fluid solvers (FVM, OpenFOAM, SU2) expect **meters**:
```python
gmsh.option.setNumber("Mesh.ScalingFactor", 0.001)  # Convert mm to meters
```

### 5.2 Quadrilateral Blossom Recombination
To generate quad-dominant or pure quad meshes suitable for boundary-layer CFD:
```python
gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.05)
gmsh.option.setNumber("Mesh.RecombineAll", 1)       # Blossom full quad recombination
gmsh.model.mesh.generate(2)
```
**Mesh Generation Results**:
- **Nodes**: 3,668
- **Elements**: 4,519 (pure quadrilaterals across surfaces)
- **Minimum Quad Quality ($Q_\text{min}$)**: $> 0.33$
- **Average Quad Quality ($Q_\text{avg}$)**: $> 0.82$

---

## 6. Gmsh `.geo` Script for Interactive Viewing

To allow one-click viewing of both geometry and mesh in the Gmsh GUI, [`meshing_script.py`](file:///c:/Users/wille/Desktop/FVMFramework/examples/1D_hybrid_rocket_engine/meshes/meshing_script.py) generates `view_hybrid_rocket.geo`:

```gmsh
// Auto-generated Gmsh script for hybrid rocket CFD visualization
SetFactory("OpenCASCADE");

// Merge CAD geometry and generated mesh
Merge "exploded_slice.step";
Merge "hybrid_rocket_cfd.msh";

// Visual display options
General.Trackball = 0;
Mesh.SurfaceFaces = 1;
Mesh.Points = 0;
Mesh.Lines = 1;
Mesh.ColorCarousel = 1;
```
Opening `view_hybrid_rocket.geo` in Gmsh immediately displays all labeled physical surfaces with color differentiation and quadrilateral element outlines.

---

## 7. Recommended Next Steps for 1D/2D FVM CFD

1. **Boundary Patch Extraction (1D Physical Curves)**:
   Extract and label key boundary curves for the fluid solver:
   - `inlet_patch`: Upstream fluid boundary at $Z = 390\text{ mm}$
   - `symmetry_patch`: Axis of symmetry along $X = 0$
   - `fuel_regression_wall`: Interface curve between `fluid_domain` and `fuel_grain`
   - `nozzle_wall`: Interface curve between `fluid_domain` and `inner_graphite_nozzle`
   - `outlet_patch`: Downstream supersonic exit boundary at $Z = -35\text{ mm}$
2. **Axisymmetric Wedge Extrusion**:
   If running an axisymmetric FVM solver (e.g. OpenFOAM 2D wedge of $5^\circ$), the 2D planar surfaces can be extruded around the Z-axis via Gmsh's `Extrude` or OpenFOAM's `makeAxialMesh` / `extrudeMesh`.
