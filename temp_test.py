from fenics import *
import steady_nse_solver
import oldroyd_3_SRTD 

import os
import sys
import matplotlib.pyplot as plt
import numpy as np

from math import log2 as log2 # For computing the rate
import time # Timing the computations

h=0.01
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
print("Num iters:")
print(solution_se.iters)
plot(solution_se.velocity)
plt.title("NSE solution")
plt.show()
print(se_time)

start_solve_original = time.time()
solution_original = oldroyd_3_SRTD.oldroyd_3_JB_SRTD(h, rad, ecc, speed, eta, lambda1, mu1, max_srtd_iters, srtd_tol)
end_solve_original = time.time()
original_time = end_solve_original - start_solve_original

plot(solution_se.velocity)
plt.title("SE solution")
plt.show()
print(original_time)

print("norm diff:")
print(errornorm(solution_original.velocity, solution_se.velocity))

