SetFactory("OpenCASCADE");
Merge "airfoil_mesh_flat-Compound.step";
Merge "airfoil_mesh_flat-Compound.geo.opt";
Mesh.ScalingFactor = 0.001; //Convert mm to meters
Coherence;
//+
Extrude {0, 0, -1.0} {
  Surface{1}; Layers {1}; Recombine;
}
//+
Physical Volume("surrounding_air_volume", 17) = {1};
//+
Physical Surface("inlet_patch", 18) = {2};
//+
Physical Surface("outlet_patch", 19) = {4};
//+
Physical Surface("top_patch", 20) = {5};
//+
Physical Surface("bottom_patch", 21) = {3};
//+
Physical Surface("side_patches", 22) = {1, 10};
//+
Physical Surface("airfoil_surface_patches", 25) = {7, 6, 9, 8};
//+
Coherence;
//+
Mesh 3;
Coherence Mesh;
Save "airfoil_mesh_flat-Compound.msh";//+
