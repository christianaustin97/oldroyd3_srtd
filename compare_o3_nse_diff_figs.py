from fenics import *
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

from meshdata import gen_mesh_jb
from meshdata import gen_mesh_ldc

# all of these plots were generated with the UCM, and l1=0.01


# JB geometry parameters
rad = 0.5
ecc = 0.25

# Mesh size
h = 0.025


meshfile = "meshdata/lid_driven_cavity_h_%.4e.h5"%h
    
if not os.path.exists(meshfile):
    print("Creating mesh...")
    gen_mesh_ldc.main(h)

#then, simply read the mesh in 
mesh = Mesh() #empty mesh
infile = HDF5File(MPI.comm_world, meshfile, 'r')
infile.read(mesh, '/mesh', True) #for some reason, need this flag to import a mesh?
infile.close()
print("Mesh loaded into FEniCS")

P_elem = FiniteElement("CG", triangle, 1) #pressure and auxiliary pressure, degree 1 elements
V_elem = VectorElement("CG", triangle, 2) #velocity, degree 2 elements

P = FunctionSpace(mesh, P_elem) # true pressure space
V = FunctionSpace(mesh, V_elem) # velocity space (not used)


ldc_infile = HDF5File(MPI.comm_world, 'results_compare_ldc/ldc_functions.h5', 'r') # open file for reading

u_nse_ldc = Function(V)
u_o3_ldc = Function(V)
p_nse_ldc = Function(P)
p_o3_ldc = Function(P)

ldc_infile.read(u_nse_ldc, "/u_nse")
ldc_infile.read(u_o3_ldc, "/u_o3")
ldc_infile.read(p_nse_ldc, "/p_nse")
ldc_infile.read(p_o3_ldc, "/p_o3")

ldc_diff = project(u_o3_ldc - u_nse_ldc, V) # difference between two solutions
magnitude_diff = dot(ldc_diff, ldc_diff)
x_diff, y_diff = ldc_diff.split(deepcopy= True)


p_x_diff = plot(x_diff)
plt.title("Horizontal component of the Difference")
plt.colorbar(p_x_diff)
plt.savefig("results_compare_ldc/x_comp_ucm_nse_ldc.svg")
plt.close()

p_mag_diff = plot(magnitude_diff)
plt.title("Magnitude of the Difference")
plt.colorbar(p_mag_diff)
plt.savefig("results_compare_ldc/magnitude_ucm_nse_ldc.svg")
plt.close()



# JB Stuff 
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

P_elem = FiniteElement("CG", triangle, 1) #pressure and auxiliary pressure, degree 1 elements
V_elem = VectorElement("CG", triangle, 2) #velocity, degree 2 elements

P = FunctionSpace(mesh, P_elem) # true pressure space
V = FunctionSpace(mesh, V_elem) # velocity space (not used)

jb_infile = HDF5File(MPI.comm_world, 'results_compare_jb/jb_functions.h5', 'r')

u_nse_jb = Function(V)
u_o3_jb = Function(V)
p_nse_jb = Function(P)
p_o3_jb = Function(P)

jb_infile.read(u_nse_jb, "/u_nse")
jb_infile.read(u_o3_jb, "/u_o3")
jb_infile.read(p_nse_jb, "/p_nse")
jb_infile.read(p_o3_jb, "/p_o3")


jb_diff = project(u_o3_jb - u_nse_jb, V) # difference between two solutions
jb_magnitude_diff = dot(jb_diff, jb_diff)
jb_x_diff, jb_y_diff = jb_diff.split(deepcopy= True)

p_mag_diff = plot(jb_magnitude_diff)
plt.title("Magnitude of the Difference")
plt.colorbar(p_mag_diff)
plt.savefig("results_compare_jb/magnitude_ucm_nse_jb.svg")
plt.close()


# azimuthal component of difference
sctheta = Expression(("(x[1]+0.25)/sqrt(x[0]*x[0]+(x[1]+0.25)*(x[1]+0.25))" , "x[0]/sqrt(x[0]*x[0]+(x[1]+0.25)*(x[1]+0.25))"), degree=4)
sctheta = interpolate(sctheta, V)

sintheta, costheta = sctheta.split(deepcopy=True)
sintheta_vec = sintheta.vector().get_local()
costheta_vec = costheta.vector().get_local()


jb_x_diff_vec = jb_x_diff.vector().get_local() 
jb_y_diff_vec = jb_y_diff.vector().get_local()

u_theta_vector = (np.multiply(sintheta_vec, jb_x_diff_vec)) - (np.multiply(costheta_vec, jb_y_diff_vec))
u_theta = Function(jb_x_diff.function_space())
u_theta.vector().set_local(u_theta_vector)


fig = plot(u_theta)
plt.title("Azimuthal component of the difference")
plt.colorbar(fig, label = "magnitude")
plt.savefig("results_compare_jb/azimuthal_ucm_nse_jb.svg")
plt.close()



