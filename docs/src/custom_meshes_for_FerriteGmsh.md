# Settings recommended for FerriteGmsh to work properly with this framework

## Example .geo file
```
SetFactory("OpenCASCADE");
Merge "dialysis_tubing_model_cone.step";
Mesh.ScalingFactor = 0.001;
Coherence;

Physical Volume("dialysis_tubing_interior", 15) = {3, 2};
//+
Physical Volume("surrounding_fluid", 16) = {1};
//+
Physical Surface("dialysis_tubing_surface", 17) = {4, 5};
//+
Coherence;
Mesh 3;
Coherence Mesh; // !!! IMPORTANT !!!
Save "dialysis_tubing_cone_output.msh";
//+
```

Some things to note:
- `Mesh.ScalingFactor = 0.001` is included to convert the native units of mm from Gmsh to the standard units of m used in this solver.
- `Coherance` and `Coherenece Mesh` are included to remove duplicate nodes and elements from the mesh. 
    - I'm not entirely sure if where `Coherance` is placed matters, but I would just put it before and after to be safe.
    - `Coherenece Mesh` MUST be placed after `Mesh 3`
    - if you want to check for this in any future gmsh file, just turn on node labels once meshed and check for duplicate node labels
- Make sure that the .step file is in the same directory