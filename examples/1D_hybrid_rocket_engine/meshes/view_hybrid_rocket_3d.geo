// Auto-generated Gmsh viewer script for 3D Hybrid Rocket CFD
SetFactory("OpenCASCADE");
Merge "exploded_slice.step";
// Scale the CAD geometry from mm to meters so it matches the scaled mesh
Dilate {{0, 0, 0}, 0.001} { Surface{:}; }
Merge "c:\Users\wille\Desktop\FVMFramework\examples\1D_hybrid_rocket_engine\meshes\hybrid_rocket_cfd_3d.msh";
Mesh.SurfaceFaces = 1;
Mesh.Lines = 1;
Mesh.ColorCarousel = 1;
