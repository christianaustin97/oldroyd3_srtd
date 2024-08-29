from fenics import *
import ufl
import steady_nse_solver
import oldroyd_3_SRTD 

import os
import sys
import matplotlib.pyplot as plt
import numpy as np

from math import log2 as log2 # For computing the rate
import time # Timing the computations

h=0.025
eta=1.0
lambda1 = 1e-2
a=1.0
mu1 = lambda1*a
speed=1.0

# geometry parameters
rad = 0.5
ecc = 0.25

max_srtd_iters = 20
srtd_tol = 1e-8

start_solve_se = time.time()
solution_se = steady_nse_solver.navier_stokes_JB(h, rad, ecc, speed, eta)
end_solve_se = time.time()
se_time = end_solve_se - start_solve_se
print("solve time:")
print(se_time)

velocity = solution_se.velocity
u_x, u_y = velocity.split(deepcopy = True) # u.e1 and u.e2

V = u_x.function_space()
V_vectorspace = velocity.function_space()
mesh = V.mesh()
plot(mesh)
plt.show()

# what I eventually want is u.e_theta and u.e_r, the polar basis vectors
"""
bearing_radius = Expression("sqrt(x[0]*x[0]+(x[1]+0.25)*(x[1]+0.25))", degree=4)
bearing_radius = interpolate(bearing_radius, V)
fig = plot(bearing_radius)
plt.title("radius value relative to the journal center")
plt.colorbar(fig, label = "mag")
plt.show()

bearing_theta = Expression("atan2(x[1]+0.25, x[0])", degree=4)
bearing_theta = interpolate(bearing_theta, V)
fig = plot(bearing_theta)
plt.title("theta value relative to the journal center")
plt.colorbar(fig, label = "mag")
plt.show()
"""
sintheta = Expression("(x[1]+0.25)/sqrt(x[0]*x[0]+(x[1]+0.25)*(x[1]+0.25))", degree=4)
sintheta = interpolate(sintheta, V)
"""fig = plot(sintheta)
plt.title("sin value relative to the journal center")
plt.colorbar(fig, label = "mag")
plt.show()"""

costheta = Expression("x[0]/sqrt(x[0]*x[0]+(x[1]+0.25)*(x[1]+0.25))", degree=4)
costheta = interpolate(costheta, V)
"""fig = plot(costheta)
plt.title("cos value relative to the journal center")
plt.colorbar(fig, label = "mag")
plt.show()"""

sintheta_vec = sintheta.vector().get_local()
costheta_vec = costheta.vector().get_local()

e_theta = Function(V_vectorspace)


u_x_vec = u_x.vector().get_local() 
u_y_vec = u_y.vector().get_local()

u_theta_vector = (np.multiply(sintheta_vec, u_x_vec)) - (np.multiply(costheta_vec, u_y_vec))
u_theta = Function(V)
u_theta.vector().set_local(u_theta_vector)


fig = plot(u_x)
plt.title("x component of magnitude")
plt.colorbar(fig, label = "magnitude")
plt.show()

fig = plot(u_theta)
plt.title("Angular/azimuthal component of velocity")
plt.colorbar(fig, label = "magnitude")
plt.show()
