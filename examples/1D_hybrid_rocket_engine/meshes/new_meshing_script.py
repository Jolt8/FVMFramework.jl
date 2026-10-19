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

# ---------------------------------------------------------------------------
# #import gmsh and freecad
# ---------------------------------------------------------------------------
if hasattr(os, "add_dll_directory") and os.path.isdir(freecad_bin):
    os.add_dll_directory(freecad_bin)

sys.path.append(freecad_bin)
sys.path.append(freecad_lib)
import FreeCAD as App
import Part
import Import
import BOPTools.SplitFeatures
import CompoundTools.Explode
import gmsh

# ---------------------------------------------------------------------------
# #import the assembly:
# ---------------------------------------------------------------------------
path_to_hybrid_cad_files = r"C:\Users\wille\Desktop\hybrid-rocket-CAD\src"
path_to_cfd_assembly = os.path.join(path_to_hybrid_cad_files, "hybrid_assembly_for_cfd.FCStd")

if not os.path.exists(path_to_cfd_assembly):
    raise FileNotFoundError(f"Assembly file not found at: {path_to_cfd_assembly}")

print("=" * 80)
print(f"Opening assembly document: '{path_to_cfd_assembly}'")
doc = App.openDocument(path_to_cfd_assembly)
print(f"Document opened: {doc.Label} ({doc.Name})")
print("=" * 80)

# ---------------------------------------------------------------------------
# #use the sketch named '1_sided_axisymmetric_clipping_plane_sketch' to make a plane through the part workbench
# ---------------------------------------------------------------------------
sk_list = doc.getObjectsByLabel("1_sided_axisymmetric_clipping_plane_sketch")
if not sk_list:
    raise RuntimeError("Sketch '1_sided_axisymmetric_clipping_plane_sketch' not found in document.")

clipping_sketch = sk_list[0]
print(f"Found clipping sketch: '{clipping_sketch.Label}' ({clipping_sketch.Name})")

clipping_face = doc.addObject("Part::Face", "CFD_ClippingFace")
clipping_face.Sources = [clipping_sketch]
doc.recompute()
print(f"Created clipping plane face: Area = {clipping_face.Shape.Area:.2f} mm^2")

# ---------------------------------------------------------------------------
# #make a compound using every body except for the part 'ASM0002_ASM'
# ---------------------------------------------------------------------------
solid_bodies = [
    obj for obj in doc.Objects
    if obj.isDerivedFrom("App::Link") and "ASM0002_ASM" not in (obj.Name, obj.Label)
]

print("\nSolid bodies included in CFD compound:")
for b in solid_bodies:
    print(f"  * {b.Label} ({b.Name})")

compound_obj = doc.addObject("Part::Compound", "CFD_Compound")
compound_obj.Links = solid_bodies
doc.recompute()

# ---------------------------------------------------------------------------
# #slice the compound using the plane created from the sketch and this new compound
# ---------------------------------------------------------------------------
slice_obj = BOPTools.SplitFeatures.makeSlice(name="CFD_Slice")
slice_obj.Base = clipping_face
slice_obj.Tools = [compound_obj]
slice_obj.Mode = "Split"
doc.recompute()
print(f"\nSliced compound with clipping plane: {len(slice_obj.Shape.Faces)} 2D faces produced.")

# ---------------------------------------------------------------------------
# #Import.export the resulting 'Exploded Slice'
# ---------------------------------------------------------------------------
# Explode the sliced compound into individual topological faces
exploded_res = CompoundTools.Explode.explodeCompound(slice_obj)
exploded_group = exploded_res[0]
exploded_children = exploded_res[1]
doc.recompute()

# Match each sliced child face to its corresponding solid body (or fluid domain)
slice_face_body_names = {}
for child in exploded_children:
    face = child.Shape
    matched_name = "fluid_domain"
    
    for body in solid_bodies:
        try:
            overlap = face.common(body.Shape).Area
            if overlap > 0.5 * face.Area:
                matched_name = body.Label
                break
        except Exception:
            pass
            
    slice_face_body_names[child.Name] = matched_name
    print(f"  Child {child.Name:<14} (Area = {face.Area:>8.2f} mm^2) -> '{matched_name}'")

# Temporarily assign descriptive labels so FreeCAD's Import.export writes them into STEP AP214/AP242 XDE
orig_labels = {c.Name: c.Label for c in exploded_children}
for c in exploded_children:
    c.Label = slice_face_body_names[c.Name]

mesh_output_dir = os.path.dirname(os.path.abspath(__file__))
step_export_path = os.path.join(mesh_output_dir, "exploded_slice.step")
Import.export(exploded_children, step_export_path)
print(f"\nExported Exploded Slice to STEP with labels: '{step_export_path}'")

# Restore original FreeCAD labels
for c in exploded_children:
    c.Label = orig_labels[c.Name]

# ===========================================================================
# GMSH 3D EXTRUSION & STRUCTURED TRANSFINITE MESHING PIPELINE
# ===========================================================================
gmsh.initialize()
gmsh.model.add("HybridRocketCFD_3D")
gmsh.open(step_export_path)

# ---------------------------------------------------------------------------
# #create a dedicated and very readable section where I can input all the gmsh options
# #as an example of one such option: 
# #gmsh.option.setNumber("Mesh.ScalingFactor", 0.001)
# ---------------------------------------------------------------------------
GMSH_OPTIONS = {
    # Unit scaling: convert mm to meters on output mesh file
    "Mesh.ScalingFactor": 0.001,
    
    # Global mesh sizing control:
    # Set to 1.0 to match the baseline element sizing from 'hybrid_assembly_for_cfd-Exploded Slice.geo'.
    # To globally refine the unstructured fluid and nozzle domains, lower this value (e.g. 0.5, 0.2, etc.).
    "Mesh.MeshSizeFactor": 0.05,
    #"Mesh.MeshSizeMax": 100.0,
    
    # 2D surface recombination:
    # Set to 0 to keep triangular 2D surface elements (which extrude into clean triangular prisms / wedges).
    # Set to 1 to recombine triangles into quadrilaterals (which extrude into hexahedra).
    "Mesh.RecombineAll": 1,
    
    # 2D meshing algorithm (1: MeshAdapt, 5: Delaunay, 6: Frontal-Delaunay, 8: Frontal-Delaunay for Quads)
    "Mesh.Algorithm": 6,
    #"Mesh.RecombinationAlgorithm": 3,
    
    # 3D meshing algorithm
    "Mesh.Algorithm3D": 1,
}

print("\n" + "=" * 80)
print("CONFIGURING GMSH OPTIONS")
print("=" * 80)
for opt_name, opt_val in GMSH_OPTIONS.items():
    gmsh.option.setNumber(opt_name, opt_val)
    print(f"  * {opt_name:<30} = {opt_val}")

# ---------------------------------------------------------------------------
# #Extrude the 2D model into a flat 3D model
# #Before doing that, call the API equivalent of Coherance;
# #Basically just do the python API equivalent of this:
# ''' 
# # This extrudes by 1 mm:
# Extrude {0, -1.0, 0} {
#   Surface{:}; Layers {1}; Recombine;
# }
# '''
# #then call the API equivalent of Coherance; again
# #Note, I know that using Cartesian coordinates for an axisymmetric case is inaccurate
# #I'm going to scale all of the cell volumes, axial face areas, and radial face areas by the approximate geometry conversion later:
# # Cell Volume: $V_i = 2\pi x_i \cdot \Delta x \cdot \Delta z$
# # Axial Face Area: $A_z = 2\pi x_{\text{cell}} \cdot \Delta x$
# # Radial Face Area: $A_x = 2\pi x_{\text{face}} \cdot \Delta z$
# ---------------------------------------------------------------------------

# Before doing extrusion: call API equivalent of Coherence;
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

initial_surfs = gmsh.model.getEntities(2)

# Build map from initial 2D surface tag to part label
surf_to_part = {}
for dim, tag in initial_surfs:
    raw_name = gmsh.model.getEntityName(dim, tag)
    clean_name = raw_name.split("/")[-1].rstrip("0123456789")
    if not clean_name:
        clean_name = "fluid_domain"
    surf_to_part[tag] = clean_name

# (Note: Transfinite structured meshing for fuel grain and liner is applied
#  after extrusion and coherence below, ensuring OpenCASCADE does not discard tags)

# Extrude by 1 mm along (0, -1.0, 0):
out_of_plane_surfaces = []
vol_to_part = {}

recombine_extrusion = bool(GMSH_OPTIONS.get("Mesh.RecombineAll", 0))
for dim, s_tag in initial_surfs:
    # Extrude by 1 layer (recombine into prisms/hexes only if RecombineAll is enabled)
    ext = gmsh.model.occ.extrude([(dim, s_tag)], 0, -1.0, 0, numElements=[1], recombine=recombine_extrusion)
    
    # In Gmsh OCC, extruding a 2D surface s_tag returns:
    # ext[0]: top surface at Y = -1.0
    # ext[1]: extruded 3D volume
    # ext[2:]: lateral boundary surfaces
    out_of_plane_surfaces.append(s_tag)       # Front face (Y = 0)
    out_of_plane_surfaces.append(ext[0][1])    # Back face (Y = -1.0)
    
    vol_tag = ext[1][1]
    vol_to_part[vol_tag] = surf_to_part[s_tag]

# Then call the API equivalent of Coherence; again
gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

print("\nExtrusion and Coherence completed successfully.")

# ---------------------------------------------------------------------------
# PHYSICAL GROUPS DEFINITIONS
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("DEFINING PHYSICAL VOLUMES & SURFACES")
print("=" * 80)

# General Rules:
#     - For all volumes: Call them '{part_name}_volume'
#     - Special Rule: Call the fluid domain 'fluid_domain'
volumes_by_name = {}
for v_dim, v_tag in gmsh.model.getEntities(3):
    pname = vol_to_part.get(v_tag, "fluid_domain")
    if pname == "fluid_domain":
        v_group_name = "fluid_domain"
    else:
        v_group_name = f"{pname}_volume"
    volumes_by_name.setdefault(v_group_name, []).append(v_tag)

for vname, vtags in sorted(volumes_by_name.items()):
    phys_vtag = gmsh.model.addPhysicalGroup(3, vtags, name=vname)
    print(f"  [Physical Volume] '{vname}' [Tag {phys_vtag}]: Volumes {vtags}")

# Classify 2D Surfaces
all_surfs_3d = gmsh.model.getEntities(2)
surface_groups = {}

for dim, s_tag in all_surfs_3d:
    bounds = gmsh.model.getParametrizationBounds(dim, s_tag)
    u_mid = 0.5 * (bounds[0][0] + bounds[0][1])
    v_mid = 0.5 * (bounds[1][0] + bounds[1][1])
    norm = gmsh.model.getNormal(s_tag, [u_mid, v_mid])
    bb = gmsh.model.getBoundingBox(dim, s_tag)
    up, down = gmsh.model.getAdjacencies(dim, s_tag)
    
    # -----------------------------------------------------------------------
    # Special Rule:
    #     - For Planes that:
    #         - Are out of plane patches whose normal vector is aligned with the extrusion direction
    #         - Call them 'out_of_plane_patch'
    # -----------------------------------------------------------------------
    if abs(norm[1]) > 0.9 and abs(norm[2]) < 1e-3:
        surface_groups.setdefault("out_of_plane_patch", []).append(s_tag)
        
    # -----------------------------------------------------------------------
    # Internal Interface between two different parts/volumes:
    #     - Call them '{part_name_1}_to_{part_name_2}_interface'
    #     - Special Rule: Interface between fluid_domain and fuel grain -> 'fuel_wall_patch'
    # -----------------------------------------------------------------------
    elif len(up) == 2:
        v1, v2 = up[0], up[1]
        name1, name2 = vol_to_part[v1], vol_to_part[v2]
        
        # Special Rule: fluid_domain and fuel grain interface
        if (name1 == "fluid_domain" and name2 == "fuel_grain") or (name2 == "fluid_domain" and name1 == "fuel_grain"):
            surface_groups.setdefault("fuel_wall_patch", []).append(s_tag)
        elif name1 != name2:
            sorted_names = sorted([name1, name2])
            interface_name = f"{sorted_names[0]}_to_{sorted_names[1]}_interface"
            surface_groups.setdefault(interface_name, []).append(s_tag)
            
    # -----------------------------------------------------------------------
    # External Boundary Surfaces (interface with nothing):
    # -----------------------------------------------------------------------
    elif len(up) == 1:
        v = up[0]
        vname = vol_to_part[v]
        
        if vname == "fluid_domain":
            # Special Rule: Essentially defines the z axis -> symmetry_patch
            if abs(bb[0]) < 1e-3 and abs(bb[3]) < 1e-3:
                surface_groups.setdefault("symmetry_patch", []).append(s_tag)
                
            # Special Rule: Facing upwards at highest z coordinate (Z ~ 390 mm), closest to center z axis
            elif abs(bb[5] - 390.0) < 1.0 and abs(bb[2] - 390.0) < 1.0 and bb[0] < 2.0:
                surface_groups.setdefault("injector_inlet_patch", []).append(s_tag)
                
            # Special Rule: Lower than the face that defines the bottom of the fuel grain (Z < 50 mm)
            elif bb[5] < 50.0:
                # Face normal parallel to z axis (flat horizontal bottom at Z = -35 mm) -> outlet_patch
                if abs(bb[2] - bb[5]) < 1e-3:
                    surface_groups.setdefault("outlet_patch", []).append(s_tag)
                # Face normal perpendicular to z axis (vertical outer boundary at X ~ 36.9 mm) -> far_field_patch
                else:
                    surface_groups.setdefault("far_field_patch", []).append(s_tag)
                    
            # All other fluid walls with nothing (e.g. top bulkhead contact or casing step)
            else:
                surface_groups.setdefault("fluid_domain_wall_patch", []).append(s_tag)
                
        # General Rule: For all surfaces that are not out of plane patches but still form a wall with nothing
        #     - Call them '{part_name}_wall_patch'
        else:
            surface_groups.setdefault(f"{vname}_wall_patch", []).append(s_tag)

# Register physical surface groups in Gmsh
for sname, stags in sorted(surface_groups.items()):
    phys_stag = gmsh.model.addPhysicalGroup(2, stags, name=sname)
    print(f"  [Physical Surface] '{sname}' [Tag {phys_stag}]: Surfaces {stags}")

# ---------------------------------------------------------------------------
# Special meshing instructions:
#     - For the fuel grain and phenolic liner, do some resarch into how to use Transfinite Curves, Surfaces, and Volumes to 
#     create a fully structured mesh for these rectancular primsm sections
# ---------------------------------------------------------------------------
# Special meshing instructions:
#     - For the fuel grain, use Transfinite Curves, Surfaces, and Volume to 
#       create a fully structured hexahedral mesh.
#     - The axial and radial subdivisions are determined dynamically from Gmsh's
#       global sizing pass (Mesh.MeshSizeFactor / Characteristic Length), so the
#       fuel grain cells automatically scale consistently with the surrounding mesh.
#     - The phenolic liner is meshed regularly (no transfinite volume).
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("CONFIGURING STRUCTURED TRANSFINITE MESHING (FUEL GRAIN)")
print("=" * 80)

# Dynamically locate fuel grain volume and determine its dimensions from CAD bounding box
v_grain = [v_tag for v_dim, v_tag in gmsh.model.getEntities(3) if vol_to_part.get(v_tag) == "fuel_grain"][0]
bb_grain = gmsh.model.getBoundingBox(3, v_grain)

L_axial = bb_grain[5] - bb_grain[2]   # Z length from bounding box
L_radial = bb_grain[3] - bb_grain[0]  # X radial thickness from bounding box

print(f"  * Dynamically queried Fuel Grain dimensions:")
print(f"      - Axial Length:     {L_axial:.2f} mm")
print(f"      - Radial Thickness: {L_radial:.2f} mm")

# Query bounding surfaces and curves of the fuel grain
surfs_grain = gmsh.model.getAdjacencies(3, v_grain)[1]
grain_curves = set()
for s_tag in surfs_grain:
    grain_curves.update(gmsh.model.getAdjacencies(2, s_tag)[1])

# Classify curves of fuel grain dynamically based on orientation
axial_curves = []
radial_curves = []
thickness_curves = []

for c_tag in grain_curves:
    bb = gmsh.model.getBoundingBox(1, c_tag)
    dx = bb[3] - bb[0]
    dy = bb[4] - bb[1]
    dz = bb[5] - bb[2]
    length = (dx**2 + dy**2 + dz**2)**0.5
    
    # Axial curve aligned with rocket axis (Z):
    if abs(dz) > 0.8 * L_axial:
        axial_curves.append(c_tag)
    # Extrusion thickness curve along Y:
    elif abs(dy - 1.0) < 1e-3 and abs(dz) < 1e-3 and abs(dx) < 1e-3:
        thickness_curves.append(c_tag)
    # Cross-sectional radial curve:
    else:
        radial_curves.append(c_tag)

# Perform 1D discretization pass to sample Gmsh's nominal element size based on GMSH_OPTIONS
gmsh.model.mesh.generate(1)
sampled_nodes = gmsh.model.mesh.getNodes(1, axial_curves[0])[0]
num_elems_axial = max(1, len(sampled_nodes) - 1)
h_mesh = L_axial / num_elems_axial
num_elems_radial = max(1, int(round(L_radial / h_mesh)))

N_THICKNESS = 2                  # 1 cell layer across 1.0 mm extrusion thickness (Y-axis)
N_AXIAL_FUEL_GRAIN = num_elems_axial + 1
N_RADIAL_FUEL_GRAIN = num_elems_radial + 1

print(f"  * Sampled Gmsh nominal element size h: {h_mesh:.2f} mm")
print(f"  * Dynamically scaled N_AXIAL_FUEL_GRAIN:  {N_AXIAL_FUEL_GRAIN} nodes ({num_elems_axial} elements)")
print(f"  * Dynamically scaled N_RADIAL_FUEL_GRAIN: {N_RADIAL_FUEL_GRAIN} nodes ({num_elems_radial} elements)")

# 1. Apply Transfinite Curve constraints ONLY to fuel grain edges
for c_tag in thickness_curves:
    gmsh.model.mesh.setTransfiniteCurve(c_tag, N_THICKNESS)
for c_tag in axial_curves:
    gmsh.model.mesh.setTransfiniteCurve(c_tag, N_AXIAL_FUEL_GRAIN)
for c_tag in radial_curves:
    gmsh.model.mesh.setTransfiniteCurve(c_tag, N_RADIAL_FUEL_GRAIN)

# 2. Configure 2D surfaces and 3D volume of fuel grain as Transfinite
# If Mesh.RecombineAll is enabled, recombine into pure Hexahedra.
# If Mesh.RecombineAll is disabled (0), generate pure structured Tetrahedra to prevent mixed-element issues in FVMFramework.
recombine_mesh = bool(GMSH_OPTIONS.get("Mesh.RecombineAll", 0))

for s_tag in surfs_grain:
    gmsh.model.mesh.setTransfiniteSurface(s_tag)
    if recombine_mesh:
        gmsh.model.mesh.setRecombine(2, s_tag)

gmsh.model.mesh.setTransfiniteVolume(v_grain)
if recombine_mesh:
    gmsh.model.mesh.setRecombine(3, v_grain)
    print(f"  * Configured Transfinite Volume & Surfaces on 'fuel_grain' (Hexahedra) (Volume {v_grain})")
else:
    print(f"  * Configured Transfinite Volume & Surfaces on 'fuel_grain' (Tetrahedra) (Volume {v_grain})")

# ---------------------------------------------------------------------------
# GENERATE MESH & COHERENCE MESH
# ---------------------------------------------------------------------------
print("\nGenerating 3D mesh...")
gmsh.model.mesh.generate(3)

# Other notes:
#     - After doing the 3D mesh, make sure to call the API equivalent of Coherence Mesh;
gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.mesh.removeDuplicateElements()

nodes = gmsh.model.mesh.getNodes()
elements = gmsh.model.mesh.getElements()
num_elements = sum(len(e) for e in elements[1])
print(f"\n3D Mesh Generation Summary:")
print(f"  - Total Mesh Nodes:    {len(nodes[0])}")
print(f"  - Total Mesh Elements: {num_elements}")

n_hexahedrons = gmsh.model.mesh.getElementsByType(5)
n_tetrahedrons = gmsh.model.mesh.getElementsByType(4)
n_wedges = gmsh.model.mesh.getElementsByType(6)

print(f"  - Total Mesh Hexahedrons: {len(n_hexahedrons[0])}")
print(f"  - Total Mesh Tetrahedrons: {len(n_tetrahedrons[0])}")
print(f"  - Total Mesh Wedges: {len(n_wedges[0])}")

# Save the generated 3D mesh file
msh_output_path = os.path.join(mesh_output_dir, "hybrid_rocket_cfd_3d.msh")
gmsh.write(msh_output_path)
print(f"Saved Gmsh 3D mesh to: '{msh_output_path}'")

# Create a companion .geo file for easy viewing in the Gmsh GUI
geo_viewer_path = os.path.join(mesh_output_dir, "view_hybrid_rocket_3d.geo")
step_filename = os.path.basename(step_export_path)
with open(geo_viewer_path, "w") as f:
    f.write('// Auto-generated Gmsh viewer script for 3D Hybrid Rocket CFD\n')
    f.write('SetFactory("OpenCASCADE");\n')
    f.write(f'Merge "{step_filename}";\n')
    f.write('// Scale the CAD geometry from mm to meters so it matches the scaled mesh\n')
    f.write('Dilate {{0, 0, 0}, 0.001} { Surface{:}; }\n')
    f.write(f'Merge "{msh_output_path}";\n')
    f.write('Mesh.SurfaceFaces = 1;\n')
    f.write('Mesh.Lines = 1;\n')
    f.write('Mesh.ColorCarousel = 1;\n')

print(f"Created Gmsh .geo viewer script: '{geo_viewer_path}'")

gmsh.finalize()
print("\n" + "=" * 80)
print("Meshing script workflow completed successfully!")
print("=" * 80)
