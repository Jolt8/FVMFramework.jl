import gmsh

import os
import sys

freecad_bin = r"C:\Program Files\FreeCAD 1.1\bin"
freecad_lib = r"C:\Program Files\FreeCAD 1.1\lib"

if hasattr(os, "add_dll_directory") and os.path.isdir(freecad_bin):
    os.add_dll_directory(freecad_bin)

sys.path.append(freecad_bin)
sys.path.append(freecad_lib)
import FreeCAD as App
import Part

path_to_hybrid_cad_files = r"C:\Users\wille\OneDrive\Desktop\FVMFramework\examples\1D_hybrid_rocket_engine\parts"

#import path to assembly meant for cfd 

#export the Exploded Slice as a .step file

#import the .step file in gmsh

