import os
import sys
import subprocess

# FreeCAD 1.1 requires Python 3.11 (python311.dll).
# If launched with a different Python interpreter (e.g. Python 3.14),
# automatically re-launch this script using FreeCAD's bundled Python executable.
freecad_bin = r"C:\Program Files\FreeCAD 1.1\bin"
freecad_lib = r"C:\Program Files\FreeCAD 1.1\lib"
freecad_python = os.path.join(freecad_bin, "python.exe")

if os.path.normcase(sys.executable) != os.path.normcase(freecad_python):
    if os.path.exists(freecad_python):
        result = subprocess.run([freecad_python, __file__] + sys.argv[1:])
        sys.exit(result.returncode)
    else:
        raise RuntimeError(
            f"FreeCAD's Python executable was not found at '{freecad_python}'. "
            "Please ensure FreeCAD 1.1 is installed."
        )

import gmsh

if hasattr(os, "add_dll_directory") and os.path.isdir(freecad_bin):
    os.add_dll_directory(freecad_bin)

sys.path.append(freecad_bin)
sys.path.append(freecad_lib)
import FreeCAD as App
import Part
import Import

path_to_hybrid_cad_files = r"C:\Users\wille\Desktop\hybrid-rocket-CAD\src"

# ---------------------------------------------------------------------------
# 1. Path to assembly meant for CFD
# ---------------------------------------------------------------------------
path_to_cfd_assembly = os.path.join(path_to_hybrid_cad_files, "hybrid_assembly_for_cfd.FCStd")

# ---------------------------------------------------------------------------
# 2. Open the assembly
# ---------------------------------------------------------------------------
# In FreeCAD's Python API, App.openDocument(path) loads an .FCStd document headlessly.
if not os.path.exists(path_to_cfd_assembly):
    raise FileNotFoundError(f"CFD assembly file not found at: {path_to_cfd_assembly}")

doc = App.openDocument(path_to_cfd_assembly)
print("=" * 80)
print(f"Opened FreeCAD Document: '{doc.Label}' (Internal Name: '{doc.Name}')")
print("=" * 80)

# ---------------------------------------------------------------------------
# 3. Report what you can see inside the assembly
# ---------------------------------------------------------------------------
# In FreeCAD:
# - doc.Objects contains all objects in the document hierarchy.
# - doc.RootObjects returns the top-level objects in the tree view.
print("\n[Assembly Structure Overview]")
print(f"Total objects in document: {len(doc.Objects)}")

links = [o for o in doc.Objects if o.isDerivedFrom("App::Link")]
groups = [o for o in doc.Objects if o.isDerivedFrom("App::DocumentObjectGroup")]
sketches = [o for o in doc.Objects if o.isDerivedFrom("Sketcher::SketchObject")]

print(f"- Component Links (App::Link): {len(links)}")
print(f"- Object Groups (App::DocumentObjectGroup): {len(groups)}")
print(f"- Sketches: {len(sketches)}")

print("\nKey Assembly Components & Slicing Objects:")
for obj in doc.Objects:
    if obj.isDerivedFrom("App::Link") or obj.Name in ["Compound", "Slice", "GrExplode_Slice"]:
        print(f"  * Name: {obj.Name:<18} | Label: {obj.Label:<28} | Type: {obj.TypeId}")

# ---------------------------------------------------------------------------
# 4. Can you see the name of each body?
#    What other information is available about each body?
# ---------------------------------------------------------------------------
# Answers:
# 1) Can you see the name of each body?
#    YES. FreeCAD provides:
#    - obj.Name: The unique, immutable internal identifier (e.g. 'Body', 'Body004')
#    - obj.Label: The editable user display name shown in the tree (e.g. 'fuel_grain')
#    - For App::Link objects, obj.LinkedObject points to the actual underlying body
#      and obj.LinkedObject.Document gives its source CAD file.
#
# 2) What other information is available about each body?
#    Via the FreeCAD Python API:
#    - Geometry (obj.Shape):
#        * obj.Shape.ShapeType: Solid, Compound, Face, Shell, etc.
#        * obj.Shape.Volume: Volume (in mm^3)
#        * obj.Shape.Area: Total surface area (in mm^2)
#        * obj.Shape.CenterOfMass: 3D vector (x, y, z)
#        * obj.Shape.BoundBox: Bounding box extents (XMin, XMax, YMin, YMax, ZMin, ZMax)
#        * obj.Shape.Solids, obj.Shape.Faces, obj.Shape.Edges, obj.Shape.Vertexes
#        * obj.Shape.MatrixOfInertia, moments of inertia, principal axes
#    - Spatial Placement:
#        * obj.Placement: Local position vector Base and Rotation quaternion/Euler angles
#        * obj.getGlobalPlacement(): World coordinates transformed through assembly
#    - Properties & Materials:
#        * obj.ShapeMaterial: Assigned physical material (density, elasticity, etc.)
#        * obj.PropertiesList: All readable/writable properties on the object
#        * obj.Visibility: Whether the object is visible or hidden

print("\n" + "=" * 80)
print("BODY INFORMATION & GEOMETRIC PROPERTIES")
print("=" * 80)

for link in links:
    shape = getattr(link, "Shape", None)
    target = getattr(link, "LinkedObject", None)
    target_doc = target.Document.Name if (target and target.Document) else "N/A"
    
    print(f"\nBody Label: '{link.Label}' (Internal Name: '{link.Name}')")
    print(f"  Source Document: {target_doc} -> Target Label: {target.Label if target else 'N/A'}")
    print(f"  Placement Base:  {link.Placement.Base}")
    
    if shape and shape.isValid():
        bb = shape.BoundBox
        print(f"  Shape Type:      {shape.ShapeType}")
        if hasattr(shape, "Volume") and abs(shape.Volume) < 1e12:
            print(f"  Volume:          {shape.Volume:,.2f} mm^3")
        print(f"  Surface Area:    {shape.Area:,.2f} mm^2")
        com = getattr(shape, "CenterOfMass", None)
        if com is not None:
            print(f"  Center of Mass:  ({com.x:.2f}, {com.y:.2f}, {com.z:.2f}) mm")
        if abs(bb.XMin) < 1e6 and abs(bb.XMax) < 1e6:
            print(f"  Bounding Box:    X[{bb.XMin:.1f}, {bb.XMax:.1f}], Y[{bb.YMin:.1f}, {bb.YMax:.1f}], Z[{bb.ZMin:.1f}, {bb.ZMax:.1f}] mm")

# ---------------------------------------------------------------------------
# 5. Is it possible to make these body names available to gmsh to define
#    physical groups like volumes and surfaces?
# ---------------------------------------------------------------------------
# Answer:
# YES, in two complementary ways:
#
# Method 1 (Direct STEP Label Metadata via FreeCAD `Import.export`):
# - Unlike `Part.export()`, FreeCAD's `Import.export()` writes Open CASCADE XDE
#   shape metadata and embeds each object's `Label` directly into the STEP AP214/AP242 file.
# - When Gmsh opens this STEP file (`gmsh.open()`), Gmsh automatically reads these
#   names into its geometry model! You can query them with `gmsh.model.getEntityName(dim, tag)`.
#
# Method 2 (Geometric Boolean Matching in Python):
# - Because FreeCAD and Gmsh both run in this script, we can match each 2D slice face
#   to its corresponding assembly body by checking geometric overlap (`face.common(body.Shape).Area`).
# - Any sliced face with zero solid overlap is identified as the CFD fluid domain!

exploded_slice = doc.getObject("GrExplode_Slice")
compound = doc.getObject("Compound")
solid_bodies = list(compound.Links) if compound else links

print("\n" + "=" * 80)
print("MAPPING EXPLODED SLICE FACES TO PHYSICAL CFD / SOLID REGIONS")
print("=" * 80)

slice_face_body_names = {}
for child in exploded_slice.Group:
    face = child.Shape
    matched_name = "fluid_domain"
    
    for body in solid_bodies:
        try:
            overlap = face.common(body.Shape).Area
            # If the slice face overlaps with this body
            if overlap > 0.5 * face.Area:
                matched_name = body.Label
                break
        except Exception:
            pass
            
    slice_face_body_names[child.Name] = matched_name
    print(f"  {child.Name} ({child.Label:<8}) Area = {face.Area:>8.2f} mm^2 -> Physical Region: '{matched_name}'")

# ---------------------------------------------------------------------------
# 6. Export the Exploded Slice as a .step file (preserving body names)
# ---------------------------------------------------------------------------
mesh_output_dir = os.path.dirname(os.path.abspath(__file__))
step_export_path = os.path.join(mesh_output_dir, "exploded_slice.step")

# Temporarily assign descriptive labels to slice children so FreeCAD embeds
# these names directly into the STEP file
original_labels = {child.Name: child.Label for child in exploded_slice.Group}
for child in exploded_slice.Group:
    child.Label = slice_face_body_names[child.Name]

print(f"\nExporting {len(exploded_slice.Group)} slice faces to STEP: '{step_export_path}'...")
Import.export(exploded_slice.Group, step_export_path)
print("STEP export completed successfully.")

# Restore original FreeCAD labels
for child in exploded_slice.Group:
    child.Label = original_labels[child.Name]

# ---------------------------------------------------------------------------
# 7. Import the .step file in gmsh and define Physical Groups
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("IMPORTING STEP INTO GMSH & DEFINING PHYSICAL GROUPS")
print("=" * 80)

gmsh.initialize()
gmsh.model.add("HybridRocketEngine2D")

# Import the STEP geometry
gmsh.open(step_export_path)

# Query 2D surface entities
surfaces = gmsh.model.getEntities(2)
print(f"Total 2D surface entities imported into Gmsh: {len(surfaces)}")

physical_groups = {}
for dim, tag in surfaces:
    # Gmsh retrieves the label embedded in the STEP file by FreeCAD
    entity_name = gmsh.model.getEntityName(dim, tag)
    # Extract base name (e.g. 'Shapes/hybrid_assembly_for_cfd/fuel_grain001/fuel_grain001' -> 'fuel_grain')
    clean_name = entity_name.split("/")[-1].rstrip("0123456789")
    if not clean_name:
        clean_name = "unassigned"
    physical_groups.setdefault(clean_name, []).append(tag)
    print(f"  Surface Tag {tag}: STEP Label = '{entity_name}' -> Group '{clean_name}'")

print("\nCreated Gmsh 2D Physical Groups:")
for group_name, tags in physical_groups.items():
    phys_tag = gmsh.model.addPhysicalGroup(2, tags, name=group_name)
    print(f"  * Physical Surface '{group_name}' [Physical ID: {phys_tag}]: Surface Tags {tags}")

# Generate a 2D mesh to verify geometry readiness
gmsh.option.setNumber("Mesh.ScalingFactor", 0.001) #Convert mm to meters
gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.05)
gmsh.option.setNumber("Mesh.RecombineAll", 1) #combined all wedge volumes into hexahedra

gmsh.model.mesh.generate(2)

nodes = gmsh.model.mesh.getNodes()
elements = gmsh.model.mesh.getElements()
num_elements = sum(len(e) for e in elements[1])
print(f"\nMesh Generation Summary:")
print(f"  - Total Mesh Nodes:    {len(nodes[0])}")
print(f"  - Total Mesh Elements: {num_elements}")

# Save the generated mesh file
msh_output_path = os.path.join(mesh_output_dir, "hybrid_rocket_cfd.msh")
gmsh.write(msh_output_path)
print(f"Saved Gmsh mesh to: '{msh_output_path}'")

# Create a .geo file that includes both the STEP file and the mesh for easy viewing in Gmsh
geo_output_path = os.path.join(mesh_output_dir, "view_hybrid_rocket.geo")
step_filename = os.path.basename(step_export_path)
msh_filename = os.path.basename(msh_output_path)

with open(geo_output_path, "w") as f:
    f.write('// Auto-generated Gmsh script for hybrid rocket CFD\n')
    f.write('SetFactory("OpenCASCADE");\n')
    f.write(f'Merge "{step_filename}";\n')
    f.write('// Scale the CAD geometry from mm to meters so it matches the scaled mesh\n')
    f.write('Dilate {{0, 0, 0}, 0.001} { Surface{:}; }\n')
    f.write(f'Merge "{msh_filename}";\n')
    f.write('Mesh.SurfaceFaces = 1;\n')
    f.write('Mesh.Lines = 1;\n')
    f.write('Mesh.ColorCarousel = 1;\n')

print(f"Created Gmsh .geo viewer script: '{geo_output_path}'")


gmsh.finalize()
print("\nWorkflow completed successfully!")
