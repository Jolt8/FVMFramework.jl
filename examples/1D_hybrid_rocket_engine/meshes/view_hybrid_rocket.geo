// Auto-generated Gmsh script for hybrid rocket CFD
SetFactory("OpenCASCADE");
Merge "exploded_slice.step";
// Scale the CAD geometry from mm to meters so it matches the scaled mesh
Dilate {{0, 0, 0}, 0.001} { Surface{:}; }
Merge "hybrid_rocket_cfd.msh";
Mesh.SurfaceFaces = 1;
Mesh.Lines = 1;
Mesh.ColorCarousel = 1;
