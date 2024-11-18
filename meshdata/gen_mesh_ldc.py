# Generates the lid-driven cavity problem mesh (which is just a the
#   unit square, possibly refined near the corner) in Gmsh for use 
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
    filename = "lid_driven_cavity_h_%.4e"%h
    filepath = "meshdata/" + filename
    
    gmsh.model.add(filename)

    # shortcut. Didn't know Python could do this lol
    factory = gmsh.model.geo
    
    # add key points
    nw_pt = factory.addPoint(0.0, 1.0, 0.0, h)
    ne_pt = factory.addPoint(1.0, 1.0, 0.0, h)
    se_pt = factory.addPoint(1.0, 0.0, 0.0, h)
    sw_pt = factory.addPoint(0.0, 0.0, 0.0, h)
    corners = [nw_pt, ne_pt, se_pt, sw_pt]

    edges = [0,0,0,0]
    # describe 4 edges or Lines by adjecent corners
    for i in range(4):
        edges[i] = factory.addLine(corners[i], corners[(i+1)%4])

    # gather the edges into one boundary or closed loop
    boundary_loop = factory.addCurveLoop(edges)
    
    # define the domain as a surface
    domain_surface = factory.addPlaneSurface([boundary_loop])

    # denote physical groups
    boundary_grp = gmsh.model.addPhysicalGroup(1, [boundary_loop]) # "1" for 1 dimensional surface/curve, the boundary
    gmsh.model.setPhysicalName(1, boundary_grp, "boundary")

    domain_grp = gmsh.model.addPhysicalGroup(2, [domain_surface])
    gmsh.model.setPhysicalName(2, domain_grp, "domain")

    # Synchronize the CAD (.geo) entities with the model
    gmsh.model.geo.synchronize()

    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)

    # ... and save it to disk
    gmsh.write(filepath + ".msh")

    # Visualize mesh
    """if '-nopopup' not in sys.argv:
        gmsh.fltk.run()"""

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

    # Create 2D mesh. "True" flag since it's a 2D mesh
    triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
    meshio.write(filepath + "_triangle.xdmf", triangle_mesh)

    # Create 1D mesh
    line_mesh = create_mesh(msh, "line", prune_z=True)
    meshio.write(filepath + "_line.xdmf", line_mesh)

    #print(".xdmf files written successfully")
    
    # If you had had a 3D mesh, you would need a 3D mesh and a 2D mesh 
    #   - No 1D mesh needed for a 3D mesh, I don't think
    # Replace 'triangle' with 'tetra' and 'line' with 'triangle'. Do not prune

    # Bring it back into FEniCS
    mymesh = Mesh()

    # 2D triangles 
    with XDMFFile(filepath + "_triangle.xdmf") as infile:
        infile.read(mymesh)
    mvc_2d = MeshValueCollection("size_t", mymesh, 2) 

    with XDMFFile(filepath + "_triangle.xdmf") as infile:
        infile.read(mvc_2d, "name_to_read")
    mf_2d = cpp.mesh.MeshFunctionSizet(mymesh, mvc_2d)

    # 1D lines
    mvc_1d = MeshValueCollection("size_t", mymesh, 1)

    with XDMFFile(filepath + "_line.xdmf") as infile:
        infile.read(mvc_1d, "name_to_read")
    mf_1d = cpp.mesh.MeshFunctionSizet(mymesh, mvc_1d)

    # Save mesh as .h5 file for easy FEniCS access, filepath.h5/mesh
    outfile = HDF5File(MPI.comm_world, filepath + ".h5", 'w')
    outfile.write(mymesh, '/mesh')
    outfile.close()





# In case it's called from command line
if __name__ == '__main__':
    main(float(sys.argv[1]))
    