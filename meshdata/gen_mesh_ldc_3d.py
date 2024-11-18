# Generates the 3d lid-driven cavity problem mesh (which is just a the
#   unit cube, possibly refined near the corner) in Gmsh for use 
#   in Fenics. Actually accepts a meshsize parameter h instead of 
#   whatever meshsize/mesh density parameter Fenics uses lol

import gmsh
import sys
import meshio
from fenics import * 
import matplotlib.pyplot as plt

def main(h):
    gmsh.initialize() # must be done first

    # Create new gmsh model
    filename = "lid_driven_cavity_3d_h_%.4e"%h
    filepath = "meshdata/" + filename

    gmsh.model.add(filename)

    # shortcut. Didn't know Python could do this lol
    factory = gmsh.model.geo
    
    # add 8 key points
    pt_000 = factory.addPoint(0.0, 0.0, 0.0, h, 1)
    pt_001 = factory.addPoint(0.0, 0.0, 1.0, h, 2)
    pt_010 = factory.addPoint(0.0, 1.0, 0.0, h, 3)
    pt_011 = factory.addPoint(0.0, 1.0, 1.0, h, 4)
    pt_100 = factory.addPoint(1.0, 0.0, 0.0, h, 5)
    pt_101 = factory.addPoint(1.0, 0.0, 1.0, h, 6)
    pt_110 = factory.addPoint(1.0, 1.0, 0.0, h, 7)
    pt_111 = factory.addPoint(1.0, 1.0, 1.0, h, 8)

    # 12 edges. Theres clever things you can do with gmsh like factory.copy and factory.translate, but the cube is simple enough
    factory.addLine(pt_000, pt_100, 1)
    factory.addLine(pt_100, pt_110, 2)
    factory.addLine(pt_110, pt_010, 3)
    factory.addLine(pt_010, pt_000, 4)

    factory.addLine(pt_001, pt_101, 5)
    factory.addLine(pt_101, pt_111, 6)
    factory.addLine(pt_111, pt_011, 7)
    factory.addLine(pt_011, pt_001, 8)

    factory.addLine(pt_000, pt_001, 9)
    factory.addLine(pt_100, pt_101, 10)
    factory.addLine(pt_110, pt_111, 11)
    factory.addLine(pt_010, pt_011, 12)
    
    # 6 faces
    # add closed loops describing/enclosing each of the 6 faces
    loop_z_0 = factory.addCurveLoop([1,2,3,4], 1)
    loop_z_1 = factory.addCurveLoop([5,6,7,8], 2)
    loop_y_0 = factory.addCurveLoop([1, 10, -5, -9], 3)
    loop_y_1 = factory.addCurveLoop([-3, 11, 7, -12], 4)
    loop_x_0 = factory.addCurveLoop([-4, 12, 8, -9], 5)
    loop_x_1 = factory.addCurveLoop([2, 11, -6, -10], 6)

    # add the 6 actual faces as PlaneSurfaces
    faces = [0,0,0,0,0,0]
    for i in range(6):
        faces[i] = factory.addPlaneSurface([i+1], i+1)
    
    # 1 volume
    # add the 6 PlaneSurfaces as a closed "loop" enclosing the volume
    boundary_surface_loop = factory.addSurfaceLoop(faces, 1)
    # add the actual volume as a volume
    factory.addVolume([boundary_surface_loop], 1)
    

    # the following works in place of lines 40-62, but I think adding physical labels more difficult
    """bottom_loop = factory.addCurveLoop([1,2,3,4], 1)
    bottom_face = factory.addPlaneSurface([1], 1)

    top_face_dimtag = factory.extrude([(2,1)], 0.0, 0.0, 1.0)"""

    # addPhysicalGroup takes in dimension and array of SurfaceLoops or CurveLoop
    fixed_boundary_grp = gmsh.model.addPhysicalGroup(2, [1, 3, 4, 5, 6], 1) # all but top lid, z=1
    gmsh.model.setPhysicalName(2, fixed_boundary_grp, "Fixed Boundary")

    top_lid_grp = gmsh.model.addPhysicalGroup(2, [2], 2)
    gmsh.model.setPhysicalName(2, top_lid_grp, "Top Boundary")

    domain_grp = gmsh.model.addPhysicalGroup(3, [boundary_surface_loop], 1)
    gmsh.model.setPhysicalName(3, domain_grp, "Domain")



    factory.synchronize()
    gmsh.model.mesh.generate(3)

    # ... and save it to disk
    gmsh.write(filepath + ".msh")

    # Visualize mesh
    """
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    """

    # Always run this at the end
    gmsh.finalize()
    
    # This is the end of the GMSH construction. Now we need to make it play nicely with Fenics

    # This function was recommended by Dokken in:
    # https://jsdokken.com/src/pygmsh_tutorial.html, "Mesh generation and conversion with GMSH and PYGMSH"
    def create_mesh(mesh, cell_type, prune_z=False):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        
        # Prune to 2D mesh if requested, ie ignore z component. mesh.prune_z_0() doesn't want to work
        points = mesh.points[:, :2] if prune_z else mesh.points
        
        out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={"name_to_read": [cell_data]})
        
        return out_mesh
    
    # Read the Gmsh mesh
    msh = meshio.read(filepath + ".msh")

    # Create 2D mesh. "False" flag since it's a 2D mesh
    triangle_mesh = create_mesh(msh, "tetra", prune_z=False)
    meshio.write(filepath + "_tetra.xdmf", triangle_mesh)

    # Create 1D mesh
    line_mesh = create_mesh(msh, "triangle", prune_z=False)
    meshio.write(filepath + "_triangle.xdmf", line_mesh)

    print(".xdmf files written successfully")
    
    # If you had had a 3D mesh, you would need a 3D mesh and a 2D mesh 
    #   - No 1D mesh needed for a 3D mesh, I don't think
    # Replace 'triangle' with 'tetra' and 'line' with 'triangle'. Do not prune

    # Bring it back into FEniCS
    mymesh = Mesh()

    # 3D triangles 
    with XDMFFile(filepath + "_tetra.xdmf") as infile:
        infile.read(mymesh)
    mvc_3d = MeshValueCollection("size_t", mymesh, 3) 

    with XDMFFile(filepath + "_tetra.xdmf") as infile:
        infile.read(mvc_3d, "name_to_read")
    mf_3d = cpp.mesh.MeshFunctionSizet(mymesh, mvc_3d)

    # 2D surfaces
    mvc_2d = MeshValueCollection("size_t", mymesh, 2)

    with XDMFFile(filepath + "_triangle.xdmf") as infile:
        infile.read(mvc_2d, "name_to_read")
    mf_2d = cpp.mesh.MeshFunctionSizet(mymesh, mvc_2d)

    # Save mesh as .h5 file for easy FEniCS access, filepath.h5/mesh
    outfile = HDF5File(MPI.comm_world, filepath + ".h5", 'w')
    outfile.write(mymesh, '/mesh')
    outfile.close()




# In case it's called from command line
if __name__ == '__main__':
    h = float(sys.argv[1])
    main(h)

