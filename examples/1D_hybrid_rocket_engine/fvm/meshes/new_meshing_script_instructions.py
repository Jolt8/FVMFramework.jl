#import gmsh and freecad


#import the assembly:
path_to_hybrid_cad_files = r"C:\Users\wille\Desktop\hybrid-rocket-CAD\src"
path_to_cfd_assembly = os.path.join(path_to_hybrid_cad_files, "hybrid_assembly_for_cfd.FCStd")

#use the sketch named '1_sided_axisymmetric_clipping_plane_sketch' to make a plane through the part workbench

#make a compound using every body except for the part 'ASM0002_ASM'

#slice the compound using the plane created from the sketch and this new compound

#Import.export the resulting 'Exploded Slice'


#create a dedicated and very readable section where I can input all the gmsh options
#as an example of one such option: 
#gmsh.option.setNumber("Mesh.ScalingFactor", 0.001)

#Extrude the 2D model into a flat 3D model
#Before doing that, call the API equivalent of Coherance;
#Basically just do the python API equivalent of this:
''' 
# This extrudes by 1 mm:
Extrude {0, -1.0, 0} {
  Surface{:}; Layers {1}; Recombine;
}
'''
#then call the API equivalent of Coherance; again
#Note, I know that using Cartesian coordinates for an axisymmetric case is inaccurate
#I'm going to scale all of the cell volumes, axial face areas, and radial face areas by the approximate geometry conversion later:
# Cell Volume: $V_i = 2\pi x_i \cdot \Delta x \cdot \Delta z$
# Axial Face Area: $A_z = 2\pi x_{\text{cell}} \cdot \Delta x$
# Radial Face Area: $A_x = 2\pi x_{\text{face}} \cdot \Delta z$

"""
Rules for meshing:
    - Call the main fluid domain 'fluid_domain'

Special meshing instructions
    - For the fuel grain and fluid_domain adjacent to the fuel_grain, do some resarch into how to use Transfinite Curves, Surfaces, and Volumes to 
    create a fully structured mesh for these rectancular primsm sections

    - For Planes that:
        - Are out of plane patches whose normal vector is aligned with the extrusion direction

        - Call them 'out_of_plane_patch'

    - For the plane that's
        - Has a face normal that's parallel to the z axis 
        - Closest to the center z axis
        - Has a node that's on the center z axis


        - Call it the injector_inlet_patch

    - For the plane that:
        - Essentially defines the z axis:

        - Call it the symmetry_patch

    - For the plane that:
        - Provides an interface between the fluid_domain and nothing 
        - Is not the symmetry patch
        - Is lower than the face that defines the bottom of the fuel grain
        - Has a face normal that is parallel to the z axis
        
        - Call it the outlet_patch
    
    - For the plane that:
        - Provides an interface between the fluid_domain and nothing 
        - Is not the symmetry patch
        - Is lower than the face that defines the bottom of the fuel grain
        - Has a face normal that is perpendicular to the z axis

        - Call it the far_field_patch

    - For the planes that:
        - Define the interface between the fluid_domain and the fuel grain

        - Group them under a 'fuel_wall_patch'

        - Note to future self: treat far_field_patch as either a slip wall / symmetry wall or an extrapolation outflow to keep the solver stable.

General Rules (does not take priority over the special meshing instructions above):
    - For all volumes 
        - Call them '{part_name}_volume'

    - For all surfaces that create an interface between two different parts/volumes:
        - Call them '{part_name_1}_to_{part_name_2}_interface'

    - For all surfaces that are not out of plane patches but still form a wall with nothing
        - Call them '{part_name}_wall_patch'

Other notes:
    - Changing the Mesh.ScalingFactor only changes the dimensions of the outputted mesh file.
    Thus, if you want to extrude by 1mm, you must set the extrusion vector to (0, -1.0, 0) because the step and geo file still use mm.

    - After doing the 3D mesh, make sure to call the API equivalent of Coherence Mesh;
"""