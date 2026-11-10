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
candidate_paths = [
    r"C:\Users\wille\Desktop\hybrid-rocket-CAD\src",
    r"C:\Users\wille\OneDrive\Desktop\hybrid-rocket-CAD\src",
]
path_to_hybrid_cad_files = next((p for p in candidate_paths if os.path.exists(p)), candidate_paths[0])
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

# Detect fuel grain bounds on the clipping plane to partition the fluid port box [0, R_in] x [Z_min, Z_max]
fg_bodies = [b for b in solid_bodies if "fuel_grain" in b.Label]
if fg_bodies:
    fg_body = fg_bodies[0]
    fg_common = clipping_face.Shape.common(fg_body.Shape)
    fg_r_in = fg_common.BoundBox.XMin
    fg_z_min = fg_common.BoundBox.ZMin
    fg_z_max = fg_common.BoundBox.ZMax
    print(f"\nFuel grain bore region detected: R_inner={fg_r_in:.2f} mm, Z=[{fg_z_min:.2f}, {fg_z_max:.2f}] mm")
    
    # Create 3D solid tool to partition the fluid port region [0, fg_r_in] x [-10, 10] x [fg_z_min, fg_z_max]
    V = App.Vector
    port_poly = Part.makePolygon([
        V(0, 0, fg_z_min),
        V(fg_r_in, 0, fg_z_min),
        V(fg_r_in, 0, fg_z_max),
        V(0, 0, fg_z_max),
        V(0, 0, fg_z_min)
    ])
    port_face = Part.Face(port_poly)
    port_solid = port_face.extrude(V(0, 20.0, 0))
    port_solid.translate(V(0, -10.0, 0))
    port_tool = doc.addObject("Part::Feature", "CFD_PortSolidTool")
    port_tool.Shape = port_solid
    
    slice_compound = doc.addObject("Part::Compound", "CFD_SliceToolsCompound")
    slice_compound.Links = [compound_obj, port_tool]
    doc.recompute()
    tools_for_slice = [slice_compound]
else:
    tools_for_slice = [compound_obj]

# ---------------------------------------------------------------------------
# #slice the compound using the plane created from the sketch and this new compound
# ---------------------------------------------------------------------------
slice_obj = BOPTools.SplitFeatures.makeSlice(name="CFD_Slice")
slice_obj.Base = clipping_face
slice_obj.Tools = tools_for_slice
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
gmsh.logger.start()
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

# Pre-identify fuel grain and fluid port box
v_grain = [v_tag for v_dim, v_tag in gmsh.model.getEntities(3) if vol_to_part.get(v_tag) == "fuel_grain"][0]
bb_grain = gmsh.model.getBoundingBox(3, v_grain)
L_axial = bb_grain[5] - bb_grain[2]
L_radial_grain = bb_grain[3] - bb_grain[0]
r_inner_grain = bb_grain[0]

v_port = None
for v_dim, v_tag in gmsh.model.getEntities(3):
    if vol_to_part.get(v_tag) == "fluid_domain":
        bb = gmsh.model.getBoundingBox(3, v_tag)
        if abs(bb[0]) < 1e-3 and abs(bb[3] - r_inner_grain) < 1.0 and abs(bb[2] - bb_grain[2]) < 1.0 and abs(bb[5] - bb_grain[5]) < 1.0:
            v_port = v_tag
            break

bb_port = gmsh.model.getBoundingBox(3, v_port) if v_port else None
L_radial_port = (bb_port[3] - bb_port[0]) if bb_port else r_inner_grain

# General Rules:
#     - For all volumes: Call them '{part_name}_volume'
#     - Special Rule: Call the fluid domain 'fluid_domain'
volumes_by_name = {}
special_volumes = set()
non_special_volumes = set()

for v_dim, v_tag in gmsh.model.getEntities(3):
    pname = vol_to_part.get(v_tag, "fluid_domain")
    if pname == "fluid_domain":
        v_group_name = "fluid_domain"
        special_volumes.add(v_group_name)
    else:
        v_group_name = f"{pname}_volume"
        non_special_volumes.add(v_group_name)
    volumes_by_name.setdefault(v_group_name, []).append(v_tag)

# Special Rule: Add dedicated 'fluid_port' physical volume containing only the port box
if v_port is not None:
    volumes_by_name["fluid_port"] = [v_port]
    special_volumes.add("fluid_port")

for vname, vtags in sorted(volumes_by_name.items()):
    phys_vtag = gmsh.model.addPhysicalGroup(3, vtags, name=vname)
    print(f"  [Physical Volume] '{vname}' [Tag {phys_vtag}]: Volumes {vtags}")

# Classify 2D Surfaces
all_surfs_3d = gmsh.model.getEntities(2)
surface_groups = {}
special_interfaces = set()
non_special_interfaces = set()
special_wall_patches = set()
non_special_wall_patches = set()

# Pre-identify injector_inlet_patch:
# - Has a face normal that's parallel to the z axis (abs(norm[2]) > 0.8)
# - Closest to the center z axis (minimum X distance)
# - Has a node that's on the center z axis (abs(bb[0]) < 1e-3)
# - Higher than the bottom of the fuel grain (so it is not the outlet_patch at Z = -35 mm)
inlet_candidates = []
for dim, s_tag in all_surfs_3d:
    up, _ = gmsh.model.getAdjacencies(dim, s_tag)
    if len(up) == 1:
        bounds = gmsh.model.getParametrizationBounds(dim, s_tag)
        u_mid = 0.5 * (bounds[0][0] + bounds[0][1])
        v_mid = 0.5 * (bounds[1][0] + bounds[1][1])
        norm = gmsh.model.getNormal(s_tag, [u_mid, v_mid])
        if abs(norm[1]) < 0.5 and abs(norm[2]) > 0.8:
            bb = gmsh.model.getBoundingBox(dim, s_tag)
            # Must have a node on the center Z axis (X_min ~ 0) and be above fuel grain bottom (not the nozzle outlet)
            if abs(bb[0]) < 1e-3 and bb[2] >= 50.0:
                inlet_candidates.append((s_tag, bb))

if inlet_candidates:
    # Closest to center Z axis: pick candidate with smallest radial extent (min X_max)
    injector_inlet_tag = min(inlet_candidates, key=lambda cand: cand[1][3])[0]
else:
    injector_inlet_tag = None

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
        special_wall_patches.add("out_of_plane_patch")
        
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
            special_interfaces.add("fuel_wall_patch")
        elif name1 != name2:
            sorted_names = sorted([name1, name2])
            interface_name = f"{sorted_names[0]}_to_{sorted_names[1]}_interface"
            surface_groups.setdefault(interface_name, []).append(s_tag)
            non_special_interfaces.add(interface_name)
            
    # -----------------------------------------------------------------------
    # External Boundary Surfaces (interface with nothing):
    # -----------------------------------------------------------------------
    elif len(up) == 1:
        v = up[0]
        vname = vol_to_part[v]
        
        # Special Rule: Facing upwards at highest z coordinate of model, closest to center z axis
        if s_tag == injector_inlet_tag:
            surface_groups.setdefault("injector_inlet_patch", []).append(s_tag)
            special_wall_patches.add("injector_inlet_patch")
            
        elif vname == "fluid_domain":
            # Special Rule: Essentially defines the z axis -> symmetry_patch
            if abs(bb[0]) < 1e-3 and abs(bb[3]) < 1e-3:
                surface_groups.setdefault("symmetry_patch", []).append(s_tag)
                special_wall_patches.add("symmetry_patch")
                
            # Special Rule: Lower than the face that defines the bottom of the fuel grain (Z < 50 mm)
            elif bb[5] < 50.0:
                # Face normal parallel to z axis (flat horizontal bottom at Z = -35 mm) -> outlet_patch
                if abs(bb[2] - bb[5]) < 1e-3:
                    surface_groups.setdefault("outlet_patch", []).append(s_tag)
                    special_wall_patches.add("outlet_patch")
                # Face normal perpendicular to z axis (vertical outer boundary at X ~ 36.9 mm) -> far_field_patch
                else:
                    surface_groups.setdefault("far_field_patch", []).append(s_tag)
                    special_wall_patches.add("far_field_patch")
                    
            # All other fluid walls with nothing (e.g. top bulkhead contact or casing step)
            else:
                surface_groups.setdefault("fluid_domain_wall_patch", []).append(s_tag)
                non_special_wall_patches.add("fluid_domain_wall_patch")
                
        # General Rule: For all surfaces that are not out of plane patches but still form a wall with nothing
        #     - Call them '{part_name}_wall_patch'
        else:
            patch_name = f"{vname}_wall_patch"
            surface_groups.setdefault(patch_name, []).append(s_tag)
            non_special_wall_patches.add(patch_name)

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
print("CONFIGURING STRUCTURED TRANSFINITE MESHING (FUEL GRAIN & FLUID PORT)")
print("=" * 80)

print(f"  * Dynamically queried Fuel Grain dimensions:")
print(f"      - Axial Length:     {L_axial:.2f} mm")
print(f"      - Radial Thickness: {L_radial_grain:.2f} mm (Fuel Grain) [Volume {v_grain}]")
if v_port:
    print(f"      - Radial Thickness: {L_radial_port:.2f} mm (Fluid Port) [Volume {v_port}]")

# Helper function to classify curves of a structured box volume
def classify_box_curves(v_tag, L_ax):
    surfs = gmsh.model.getAdjacencies(3, v_tag)[1]
    curves = set()
    for s_tag in surfs:
        curves.update(gmsh.model.getAdjacencies(2, s_tag)[1])
    ax_c, rad_c, th_c = [], [], []
    for c_tag in curves:
        bb = gmsh.model.getBoundingBox(1, c_tag)
        dx = bb[3] - bb[0]
        dy = bb[4] - bb[1]
        dz = bb[5] - bb[2]
        if abs(dz) > 0.8 * L_ax:
            ax_c.append(c_tag)
        elif abs(dy - 1.0) < 1e-3 and abs(dz) < 1e-3 and abs(dx) < 1e-3:
            th_c.append(c_tag)
        else:
            rad_c.append(c_tag)
    return surfs, ax_c, rad_c, th_c

surfs_grain, ax_curves_grain, rad_curves_grain, th_curves_grain = classify_box_curves(v_grain, L_axial)
if v_port:
    surfs_port, ax_curves_port, rad_curves_port, th_curves_port = classify_box_curves(v_port, L_axial)
else:
    surfs_port, ax_curves_port, rad_curves_port, th_curves_port = [], [], [], []

# Perform 1D discretization pass to sample Gmsh's nominal element size based on GMSH_OPTIONS
gmsh.model.mesh.generate(1)
sampled_nodes = gmsh.model.mesh.getNodes(1, ax_curves_grain[0])[0]
num_elems_axial = max(1, len(sampled_nodes) - 1)
h_mesh = L_axial / num_elems_axial
num_elems_radial_grain = max(1, int(round(L_radial_grain / h_mesh)))
num_elems_radial_port = max(1, int(round(L_radial_port / h_mesh)))
# Clear temporary 1D sampling mesh so Gmsh rebuilds 1D with full transfinite constraints
gmsh.model.mesh.clear()

N_THICKNESS = 2                  # 1 cell layer across 1.0 mm extrusion thickness (Y-axis)
N_AXIAL = num_elems_axial + 1
N_RADIAL_GRAIN = num_elems_radial_grain + 1
N_RADIAL_PORT = num_elems_radial_port + 1

print(f"  * Sampled Gmsh nominal element size h: {h_mesh:.2f} mm")
print(f"  * Dynamically scaled N_AXIAL:        {N_AXIAL} nodes ({num_elems_axial} elements)")
print(f"  * Dynamically scaled N_RADIAL_GRAIN: {N_RADIAL_GRAIN} nodes ({num_elems_radial_grain} elements)")
if v_port:
    print(f"  * Dynamically scaled N_RADIAL_PORT:  {N_RADIAL_PORT} nodes ({num_elems_radial_port} elements)")

recombine_mesh = bool(GMSH_OPTIONS.get("Mesh.RecombineAll", 0))

# 1. Apply Transfinite Curve constraints
for c_tag in set(th_curves_grain + th_curves_port):
    gmsh.model.mesh.setTransfiniteCurve(c_tag, N_THICKNESS)
for c_tag in set(ax_curves_grain + ax_curves_port):
    gmsh.model.mesh.setTransfiniteCurve(c_tag, N_AXIAL)

def apply_radial_curve_transfinite(v_tag, rad_curves, num_elems_rad, L_rad):
    bb = gmsh.model.getBoundingBox(3, v_tag)
    bot_rad = [c for c in rad_curves if abs(gmsh.model.getBoundingBox(1, c)[2] - bb[2]) < 1.0]
    top_rad = [c for c in rad_curves if abs(gmsh.model.getBoundingBox(1, c)[2] - bb[5]) < 1.0]
    
    # Bottom full-width radial curves
    for c in bot_rad:
        gmsh.model.mesh.setTransfiniteCurve(c, num_elems_rad + 1)
        
    # Top curves (handles split curves at interfaces if present)
    for y_val in [0.0, -1.0]:
        segs = [c for c in top_rad if abs(gmsh.model.getBoundingBox(1, c)[1] - y_val) < 1e-3]
        if len(segs) == 1:
            gmsh.model.mesh.setTransfiniteCurve(segs[0], num_elems_rad + 1)
        elif len(segs) > 1:
            segs.sort(key=lambda c: gmsh.model.getBoundingBox(1, c)[0])
            lengths = [(gmsh.model.getBoundingBox(1, c)[3] - gmsh.model.getBoundingBox(1, c)[0]) for c in segs]
            allocated = []
            rem = num_elems_rad
            for l in lengths[:-1]:
                n = max(1, int(round(num_elems_rad * l / L_rad)))
                allocated.append(n)
                rem -= n
            allocated.append(max(1, rem))
            diff = sum(allocated) - num_elems_rad
            if diff != 0:
                allocated[-1] -= diff
            for c, n in zip(segs, allocated):
                gmsh.model.mesh.setTransfiniteCurve(c, n + 1)

apply_radial_curve_transfinite(v_grain, rad_curves_grain, num_elems_radial_grain, L_radial_grain)
if v_port:
    apply_radial_curve_transfinite(v_port, rad_curves_port, num_elems_radial_port, L_radial_port)

# 2. Configure 2D surfaces and 3D volume of fuel grain and fluid port as Transfinite
for s_tag in set(list(surfs_grain) + list(surfs_port)):
    bnd_pts = gmsh.model.getBoundary([(2, s_tag)], combined=False, oriented=False, recursive=True)
    pts = [p[1] for p in bnd_pts]
    
    if len(pts) == 5:
        # Option A: Deduce the 4 geometric outer corners (bottom-left, bottom-right, top-right, top-left)
        # by checking point coordinates relative to each other (using 2D bounding box in X-Z plane).
        # This safely ignores intermediate split nodes on composite boundary edges (e.g. at the injector interfaces).
        bb = gmsh.model.getBoundingBox(2, s_tag)
        target_corners = [
            (bb[0], bb[2]),  # Bottom-Left:  (X_min, Z_min)
            (bb[3], bb[2]),  # Bottom-Right: (X_max, Z_min)
            (bb[3], bb[5]),  # Top-Right:    (X_max, Z_max)
            (bb[0], bb[5]),  # Top-Left:     (X_min, Z_max)
        ]
        corners = []
        for tx, tz in target_corners:
            best_p = min(pts, key=lambda p: (gmsh.model.getValue(0, p, [])[0] - tx)**2 + (gmsh.model.getValue(0, p, [])[2] - tz)**2)
            corners.append(best_p)
            
        gmsh.model.mesh.setTransfiniteSurface(s_tag, cornerTags=corners)
    else:
        # Standard 4-vertex (or 3-vertex) surface: use the old method
        gmsh.model.mesh.setTransfiniteSurface(s_tag)
        
    if recombine_mesh:
        gmsh.model.mesh.setRecombine(2, s_tag)

for v_target, v_label in [(v_grain, "fuel_grain"), (v_port, "fluid_port")]:
    if v_target is None:
        continue
    gmsh.model.mesh.setTransfiniteVolume(v_target)
    if recombine_mesh:
        gmsh.model.mesh.setRecombine(3, v_target)
        print(f"  * Configured Transfinite Volume & Surfaces on '{v_label}' (Hexahedra) (Volume {v_target})")
    else:
        print(f"  * Configured Transfinite Volume & Surfaces on '{v_label}' (Tetrahedra) (Volume {v_target})")

# ---------------------------------------------------------------------------
# GENERATE MESH & COHERENCE MESH
# ---------------------------------------------------------------------------
print("\nGenerating 3D mesh...")
gmsh.model.mesh.generate(3)

gmsh_logs = gmsh.logger.get()
for l in gmsh_logs:
    if "warning" in l.lower() or "error" in l.lower() or "transfinite" in l.lower():
        print(f"  [Gmsh Log] {l}")

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
    f.write(f'Merge "exploded_slice.step";\n')
    f.write('// Scale the CAD geometry from mm to meters so it matches the scaled mesh\n')
    f.write('Dilate {{0, 0, 0}, 0.001} { Surface{:}; }\n')
    f.write('Extrude {0, -0.001, 0.0} { Surface{:}; }\n')
    f.write(f'Merge "{msh_output_path}";\n')
    f.write('Mesh.SurfaceFaces = 1;\n')
    f.write('Mesh.Lines = 1;\n')
    f.write('Mesh.ColorCarousel = 1;\n')

print(f"Created Gmsh .geo viewer script: '{geo_viewer_path}'")

# Create markdown summary file for special and non-special groups
md_summary_path = os.path.join(mesh_output_dir, "mesh_groups_summary.md")
with open(md_summary_path, "w", encoding="utf-8") as f:
    f.write("# Mesh Physical Groups Summary\n\n")
    
    f.write("## Special Volumes\n")
    for v in sorted(special_volumes):
        f.write(f"- {v}\n")
    f.write("\n")
    
    f.write("## Non-Special Volumes\n")
    for v in sorted(non_special_volumes):
        f.write(f"- {v}\n")
    f.write("\n")
    
    f.write("## Special Interfaces\n")
    for s in sorted(special_interfaces):
        f.write(f"- {s}\n")
    f.write("\n")
    
    f.write("## Non-Special Interfaces\n")
    for s in sorted(non_special_interfaces):
        f.write(f"- {s}\n")
    f.write("\n")
    
    f.write("## Special Wall Patches\n")
    for p in sorted(special_wall_patches):
        f.write(f"- {p}\n")
    f.write("\n")
    
    f.write("## Non-Special Wall Patches\n")
    for p in sorted(non_special_wall_patches):
        f.write(f"- {p}\n")
    f.write("\n")

print(f"Created mesh groups markdown summary: '{md_summary_path}'")

gmsh.finalize()
print("\n" + "=" * 80)
print("Meshing script workflow completed successfully!")
print("=" * 80)
