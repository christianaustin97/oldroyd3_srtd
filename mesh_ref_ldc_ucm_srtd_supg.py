"""
This performs mesh refinement experiment on the SRTD algorithm (with SUPG). We
assume no knowledge of the exact solution, and the convergence rate is approximated 
using successive mesh refinements

To change the solver or geometry, make sure you change the the following variables:
    table_file
    meshfile (and where you call gen_mesh_*.main())
    solution
"""

from fenics import *
from oldroyd_3_EVSS import *
from oldroyd_3_SRTD import *
from oldroyd_3_SRTD_SUPG import *

import os
import matplotlib.pyplot as plt
import numpy as np

from math import log2 as log2 # For computing the rate
import time # Timing the computations
import csv # Saving Results

# Fluid parameters
eta = 1.0
lambda1 = 1e-2
mu1 = lambda1

# SRTD algorithm parameters
max_srtd_iters = 20
srtd_tol = 1e-9

# meshsize for mesh refinement experiment
h_array = -1.0*np.array([0, 1, 2, 3, 4]) 
h_array = 0.1 * (2 ** h_array)

# For saving the info
############# CHANGE THIS DIFFERENT METHODS ##############
table_file = open('mesh_ref_ldc_ucm_l1=%.3e_srtd_supg.csv'%lambda1, 'w') 
writer = csv.writer(table_file)
writer.writerow(['h','elements','$L^{2}$ error $\\vu$','$L^{2}$ rate $\\vu$', '$H^{1}$ error $\\vu$', '$H^{1}$ rate $\\vu$', '$L^{2}$ error $p$', '$L^{2}$ rate $p$','time(s)'])

# Placeholder values for previous errors and functions
l2_diff_prev = 1.0 
h1_diff_prev = 1.0
l2_p_diff_prev = 1.0

u_prev = Constant((1.0, 1.0))
p_prev = Constant(1.0)

# mesh refinement experiment
for h in h_array:
    meshfile = "meshdata/lid_driven_cavity_h_%.4e.h5"%h
    
    if not os.path.exists(meshfile):
        print("Creating mesh...")
        gen_mesh_ldc.main(h)
    
    # Read the mesh in 
    mesh = Mesh() #empty mesh
    infile = HDF5File(MPI.comm_world, meshfile, 'r')
    infile.read(mesh, '/mesh', True) #for some reason, need this True flag to import a mesh?
    infile.close()
    
    num_els = mesh.num_cells()
    
    # Finite element spaces for comparisons
    V_elem = VectorElement("CG", triangle, 2) 
    P_elem = FiniteElement("CG", triangle, 1)
    V = FunctionSpace(mesh, V_elem)
    P = FunctionSpace(mesh, P_elem)

    # Interpolate previous solutions onto current mesh
    u_prev = interpolate(u_prev, V) 
    p_prev = interpolate(p_prev, P)

    # Start solve
    ############# CHANGE THIS DIFFERENT METHODS ##############
    start_solve = time.time()
    solution = oldroyd_3_LDC_SRTD_SUPG(h, eta, lambda1, mu1, max_srtd_iters, srtd_tol)
    end_solve = time.time()
    solve_time = end_solve - start_solve

    u = solution.velocity
    p = solution.pressure

    l2_diff = errornorm(u, u_prev, norm_type = "l2")
    h1_diff = errornorm(u, u_prev, norm_type = "h1")
    l2_p_diff = errornorm(p, p_prev, norm_type = "l2")

    l2_rate = log2(l2_diff_prev / l2_diff)
    h1_rate = log2(h1_diff_prev / h1_diff)
    l2_p_rate = log2(l2_p_diff_prev / l2_p_diff)

    # save data
    csv_str = ['%.3e'%h, num_els, '%.3e'%l2_diff, "%.3f"%l2_rate,'%.3e'%h1_diff, "%.3f"%h1_rate, '%.3e'%l2_p_diff, "%.3f"%l2_p_rate, "%.3f"%solve_time]
    writer.writerow(csv_str)

    # update for next refinement/iteration
    u_prev = u
    u_prev.set_allow_extrapolation(True)
    p_prev = p
    p_prev.set_allow_extrapolation(True)
    l2_diff_prev = l2_diff
    h1_diff_prev = h1_diff
    l2_p_diff_prev = l2_p_diff
    
#post outer loop code here
table_file.close()


