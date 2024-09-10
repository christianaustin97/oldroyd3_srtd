from fenics import *
import ufl
import oldroyd_3_SRTD 
import pointer_test_oldroyd3

import os
import sys
import time
import matplotlib.pyplot as plt
import numpy as np


#h=0.0125
h=0.025
eta=1.0
lambda1 = 4e-2
a=1.0
mu1 = lambda1*a
speed=1.0

max_srtd_iters = 20
srtd_tol = 1e-8

# JB geometry parameters
rad = 0.5
ecc = 0.25


start_solve_original = time.time()
original_soln = oldroyd_3_SRTD.oldroyd_3_JB_SRTD(h, rad, ecc, speed, eta, lambda1, mu1, max_srtd_iters, srtd_tol)
end_solve_original = time.time()
original_time = end_solve_original - start_solve_original

start_solve_pointer = time.time()
pointer_soln = pointer_test_oldroyd3.oldroyd_3_JB_SRTD(h, rad, ecc, speed, eta, lambda1, mu1, max_srtd_iters, srtd_tol)
end_solve_pointer = time.time()
pointer_time = end_solve_pointer - start_solve_pointer


print("diff between 2 solves:")
print("%0.5e"%errornorm(pointer_soln.velocity, original_soln.velocity, 'l2', degree_rise=0))

plt.semilogy(list(original_soln.residuals.keys())[1:], list(original_soln.residuals.values())[1:])
plt.show()
plt.semilogy(list(pointer_soln.residuals.keys())[1:], list(pointer_soln.residuals.values())[1:])
plt.show()




print("original version solve time:")
print(original_time)
print("pointer version solve time:")
print(pointer_time)




plot(pointer_soln.velocity)
plt.show()