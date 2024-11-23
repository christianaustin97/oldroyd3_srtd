from fenics import *
import numpy as np
import csv # Saving Results
import matplotlib.pyplot as plt

from meshdata import gen_mesh_jb
import os


s=1.0

# LDC 2D
table_file = open('table_dof/LDC2D.csv', 'w') 
writer = csv.writer(table_file)
writer.writerow(['h','elements', 'SRTD NSE DoF', 'SRTD NSE Total DoF', 'SRTD Press DoF', 'SRTD Stress DoF', 'EVSS DoF', 'EVSS Total DoF'])

h_array = -1.0*np.array([0, 1, 2, 3, 4]) 
h_array = 0.1 * (2 ** h_array)

for h in h_array:
    nx = round(1/h)
    mesh = UnitSquareMesh(nx, nx)

    # Element spaces for SRTD
    P_elem = FiniteElement("CG", triangle, 1) #pressure and auxiliary pressure, degree 1 elements
    V_elem = VectorElement("CG", triangle, 2) #velocity, degree 2 elements
    T_elem = TensorElement("CG", triangle, 2, symmetry=True) #stress tensor, degree 2 elements 
    
    W_elem = MixedElement([V_elem, P_elem]) # Mixed/Taylor Hood element space for Navier-Stokes type equations

    W = FunctionSpace(mesh, W_elem) # Taylor-Hood/mixed space
    P = FunctionSpace(mesh, P_elem) # true pressure space
    V = FunctionSpace(mesh, V_elem) # velocity space (not used)
    T = FunctionSpace(mesh, T_elem) # tensor space

    g = Constant((1.0, 1.0))
    bc = DirichletBC(W.sub(0), g, 'on_boundary')

    u = Function(W)
    bc.apply(u.vector())

    nse_int_dof = len(np.array(np.where(u.vector().get_local() == 0.0)[0]))
    

    #EVSS
    D_elem = VectorElement("CG", triangle, 1) # Deformation tensor, defined as VectorElement to exploit symmetry

    W_evss_elem = MixedElement([V_elem, P_elem, T_elem, D_elem]) # Mixed element (u, p, Sigma, D_Vec)
    
    # Function Spaces
    W_evss = FunctionSpace(mesh, W_evss_elem) # Only need this one for pure EVSS, the mixed Function Space

    bc_evss = DirichletBC(W_evss.sub(0), g, 'on_boundary')

    u_evss = Function(W_evss)
    bc.apply(u_evss.vector())

    #bndry_dof = len(np.array(np.where(u.vector().get_local() == 1.0)[0]))
    evss_int_dof = len(np.array(np.where(u_evss.vector().get_local() == 0.0)[0]))


    outrow = ['%.4e'%h, '%d'%mesh.num_cells(), '%d'%nse_int_dof, '%d'%W.dim(), '%d'%P.dim(), '%d'%T.dim(), '%d'%evss_int_dof, '%d'%W_evss.dim()]
    writer.writerow(outrow)

table_file.close()



######################################################################################3
# Journal bearing problem


# LDC 2D
table_file = open('table_dof/jb.csv', 'w') 
writer = csv.writer(table_file)
writer.writerow(['h','elements', 'SRTD NSE DoF', 'SRTD NSE Total DoF', 'SRTD Press DoF', 'SRTD Stress DoF', 'EVSS DoF', 'EVSS Total DoF'])

h_array = -1.0*np.array([0, 1, 2, 3, 4, 5]) 
h_array = 0.2 * (2 ** h_array)

rad = 0.5
ecc = 0.25

for h in h_array:
    nx = round(1/h)
    mesh = UnitSquareMesh(nx, nx)

    meshfile = "meshdata/journal_bearing_h_%.4e.h5"%h
    
    if not os.path.exists(meshfile):
        print("Creating mesh...")
        gen_mesh_jb.main(h, rad, ecc)
    
    #then, simply read the mesh in 
    mesh = Mesh() #empty mesh
    infile = HDF5File(MPI.comm_world, meshfile, 'r')
    infile.read(mesh, '/mesh', True) #for some reason, need this flag to import a mesh?
    infile.close()
    print("Mesh loaded into FEniCS")

    # Element spaces for SRTD
    P_elem = FiniteElement("CG", triangle, 1) #pressure and auxiliary pressure, degree 1 elements
    V_elem = VectorElement("CG", triangle, 2) #velocity, degree 2 elements
    T_elem = TensorElement("CG", triangle, 2, symmetry=True) #stress tensor, degree 2 elements 
    
    W_elem = MixedElement([V_elem, P_elem]) # Mixed/Taylor Hood element space for Navier-Stokes type equations

    W = FunctionSpace(mesh, W_elem) # Taylor-Hood/mixed space
    P = FunctionSpace(mesh, P_elem) # true pressure space
    V = FunctionSpace(mesh, V_elem) # velocity space (not used)
    T = FunctionSpace(mesh, T_elem) # tensor space

    g = Constant((1.0, 1.0))
    bc = DirichletBC(W.sub(0), g, 'on_boundary')

    u = Function(W)
    bc.apply(u.vector())

    nse_int_dof = len(np.array(np.where(u.vector().get_local() == 0.0)[0]))
    

    #EVSS
    D_elem = VectorElement("CG", triangle, 1) # Deformation tensor, defined as VectorElement to exploit symmetry

    W_evss_elem = MixedElement([V_elem, P_elem, T_elem, D_elem]) # Mixed element (u, p, Sigma, D_Vec)
    
    # Function Spaces
    W_evss = FunctionSpace(mesh, W_evss_elem) # Only need this one for pure EVSS, the mixed Function Space
    
    bc_evss = DirichletBC(W_evss.sub(0), g, 'on_boundary')

    u_evss = Function(W_evss)
    bc.apply(u_evss.vector())

    #bndry_dof = len(np.array(np.where(u.vector().get_local() == 1.0)[0]))
    evss_int_dof = len(np.array(np.where(u_evss.vector().get_local() == 0.0)[0]))


    outrow = ['%.4e'%h, '%d'%mesh.num_cells(), '%d'%nse_int_dof, '%d'%W.dim(), '%d'%P.dim(), '%d'%T.dim(), '%d'%evss_int_dof, '%d'%W_evss.dim()]
    writer.writerow(outrow)

table_file.close()



######################################################################################3
# LDC 3D

table_file = open('table_dof/LDC3D.csv', 'w') 
writer = csv.writer(table_file)
writer.writerow(['h','elements', 'SRTD NSE DoF', 'SRTD NSE Total DoF', 'SRTD Press DoF', 'SRTD Stress DoF', 'EVSS DoF', 'EVSS Total DoF'])

h_array = -1.0*np.array([0, 1, 2]) 
h_array = 0.25 * (2 ** h_array)

for h in h_array:
    nx = round(1/h)
    mesh = UnitCubeMesh(nx, nx, nx)

    # Element spaces for SRTD
    P_elem = FiniteElement("CG", tetrahedron, 1) #pressure and auxiliary pressure, degree 1 elements
    V_elem = VectorElement("CG", tetrahedron, 2) #velocity, degree 2 elements
    T_elem = TensorElement("CG", tetrahedron, 2, symmetry=True) #stress tensor, degree 2 elements 
    
    W_elem = MixedElement([V_elem, P_elem]) # Mixed/Taylor Hood element space for Navier-Stokes type equations

    W = FunctionSpace(mesh, W_elem) # Taylor-Hood/mixed space
    P = FunctionSpace(mesh, P_elem) # true pressure space
    V = FunctionSpace(mesh, V_elem) # velocity space (not used)
    T = FunctionSpace(mesh, T_elem) # tensor space

    g = Constant((1.0, 1.0, 1.0))
    bc = DirichletBC(W.sub(0), g, 'on_boundary')

    u = Function(W)
    bc.apply(u.vector())

    nse_int_dof = len(np.array(np.where(u.vector().get_local() == 0.0)[0]))
    

    #EVSS
    D_elem = VectorElement("CG", tetrahedron, 1, dim=5) # Deformation tensor, defined as VectorElement to exploit symmetry
    # D is Symmetric (9 -> 6 indep vars) and traceless (6 -> 5 vars)

    W_evss_elem = MixedElement([V_elem, P_elem, T_elem, D_elem]) # Mixed element (u, p, Sigma, D_Vec)
    
    # Function Spaces
    W_evss = FunctionSpace(mesh, W_evss_elem) # Only need this one for pure EVSS, the mixed Function Space
    
    bc_evss = DirichletBC(W_evss.sub(0), g, 'on_boundary')

    u_evss = Function(W_evss)
    bc.apply(u_evss.vector())

    #bndry_dof = len(np.array(np.where(u.vector().get_local() == 1.0)[0]))
    evss_int_dof = len(np.array(np.where(u_evss.vector().get_local() == 0.0)[0]))


    outrow = ['%.4e'%h, '%d'%mesh.num_cells(), '%d'%nse_int_dof, '%d'%W.dim(), '%d'%P.dim(), '%d'%T.dim(), '%d'%evss_int_dof, '%d'%W_evss.dim()]
    writer.writerow(outrow)

table_file.close()

