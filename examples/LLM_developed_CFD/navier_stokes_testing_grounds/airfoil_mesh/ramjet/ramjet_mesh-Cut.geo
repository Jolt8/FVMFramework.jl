SetFactory("OpenCASCADE");
Merge "ramjet_mesh-Cut.step";
Merge "ramjet_mesh-Cut.geo.opt";
Mesh.ScalingFactor = 0.001; //Convert mm to meters
Coherence;
//+
Extrude {0, 0, -1.0} {
  Surface{1}; Layers {1}; Recombine;
}
//+
Physical Volume("surrounding_air_volume", 17) = {1};
//+
Physical Surface("inlet_patch", 18) = {5};
//+
Physical Surface("outlet_patch", 19) = {3};
//+
Physical Surface("top_patch", 20) = {4};
//+
Physical Surface("bottom_patch", 21) = {2};
//+
Physical Surface("side_patches", 22) = {1, 22};
//+
Physical Surface("airfoil_surface_patches", 25) = {6 : 21};
//+
Coherence;
//+
Mesh 3;
Coherence Mesh;
Save "ramjet_mesh-Cut.msh";
//+
