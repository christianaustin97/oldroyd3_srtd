"""
Does some comparisons of Oldroyd 3 solution to the NSE/Newtonian solution
"""

from fenics import *
from steady_nse_solver import *
from oldroyd_3_EVSS import *
from oldroyd_3_SRTD import *

import os
import sys
import csv
import matplotlib.pyplot as plt
import numpy as np

if not (os.path.isfile('results_compare_ldc/ucm_nse_ldc.csv')):
    # create file and write headers
    print("Creating ucm_nse_ldc.csv file ...")
    outfile = open('results_compare_ldc/ucm_nse_ldc.csv', 'a')
    writer = csv.writer(outfile)
    writer.writerow(['meshsize', 'eta', 'lambda1', 'u l2_diff', 'u h1_diff', 'p l2_diff'])
    outfile.close()

if not (os.path.isfile('results_compare_jb/ucm_nse_jb.csv')):
    # create file and write headers
    print("Creating ucm_nse_jb.csv file ...")
    outfile = open('results_compare_jb/ucm_nse_jb.csv', 'a')
    writer = csv.writer(outfile)
    writer.writerow(['meshsize', 'radius', 'eccentricity', 'eta', 'lambda1', 'u l2_diff', 'u h1_diff', 'p l2_diff'])
    outfile.close()

ucm_ldc_outfile = open('results_compare_ldc/ucm_nse_ldc.csv', 'a')
ucm_ldc_writer = csv.writer(ucm_ldc_outfile)

ucm_jb_outfile = open('results_compare_jb/ucm_nse_jb.csv', 'a')
ucm_jb_writer = csv.writer(ucm_jb_outfile)

# SRTD algorithm parameters
max_srtd_iters = 20
srtd_tol = 1e-9

# JB geometry parameters
rad = 0.5
ecc = 0.25

# Mesh size
h = 0.025

# characteristic speeds
s = 1.0
eta = 1.0

# Differing l1 values to compare to NSE
coeffs = np.arange(1,10,2) # 1,3,5,7,9
pows = np.power(10.0, np.arange(-6,-1)) # 10^-6, 10^-5, ..., 10^-2
l1vals = np.concatenate((coeffs*pows[0], coeffs*pows[1], coeffs*pows[2], coeffs*pows[3], coeffs[0:1]*pows[4])) #1e-6, 3e-6, ..., 9e-6, ..., 1e-3,...9e-3, 1e-2

print(l1vals)
#input('Press <ENTER> to continue')

# L2 and H1 differences in UCM and 
ldc_l2_differences = np.zeros(len(l1vals))
ldc_h1_differences = np.zeros(len(l1vals))
ldc_l2_p_differences = np.zeros(len(l1vals))
jb_l2_differences = np.zeros(len(l1vals))
jb_h1_differences = np.zeros(len(l1vals))
jb_l2_p_differences = np.zeros(len(l1vals))

for i in range(len(l1vals)):
    l1 = float(l1vals[i])
    mu1 = l1 # UCM

    # LDC solutions
    o3_ldc_soln = oldroyd_3_LDC_SRTD(h, s, eta, l1, mu1, max_srtd_iters, srtd_tol)
    nse_ldc_soln = navier_stokes_LDC(h, s, eta)

    ldc_l2_diff = errornorm(o3_ldc_soln.velocity, nse_ldc_soln.velocity, 'l2')
    ldc_h1_diff = errornorm(o3_ldc_soln.velocity, nse_ldc_soln.velocity, 'h1')
    ldc_l2_p_diff = errornorm(o3_ldc_soln.pressure, nse_ldc_soln.pressure, 'l2')
    ldc_l2_differences[i] = ldc_l2_diff
    ldc_h1_differences[i] = ldc_h1_diff
    ldc_l2_p_differences[i] = ldc_l2_p_diff

    ucm_ldc_writer.writerow(['%.3e'%h, '%.3e'%eta, '%.3e'%l1, '%.3e'%ldc_l2_diff, '%.3e'%ldc_h1_diff, '%.3e'%ldc_l2_p_diff])

    # JB solutions
    o3_jb_soln = oldroyd_3_JB_SRTD(h, rad, ecc, s, eta, l1, mu1, max_srtd_iters, srtd_tol)
    nse_jb_soln = navier_stokes_JB(h, rad, ecc, s, eta)

    jb_l2_diff = errornorm(o3_jb_soln.velocity, nse_jb_soln.velocity, 'l2')
    jb_h1_diff = errornorm(o3_jb_soln.velocity, nse_jb_soln.velocity, 'h1')
    jb_l2_p_diff = errornorm(o3_jb_soln.pressure, nse_jb_soln.pressure, 'l2')
    jb_l2_differences[i] = jb_l2_diff
    jb_h1_differences[i] = jb_h1_diff
    jb_l2_p_differences[i] = jb_l2_p_diff

    ucm_jb_writer.writerow(['%.3e'%h, '%.3e'%rad, '%.3e'%ecc, '%.3e'%eta, '%.3e'%l1, '%.3e'%ldc_l2_diff, '%.3e'%ldc_h1_diff, '%.3e'%ldc_l2_p_diff])

ucm_ldc_outfile.close()
ucm_jb_outfile.close()


# save highest l1 solution results
ldc_function_outfile = HDF5File(MPI.comm_world, 'results_compare_ldc/ldc_functions.h5', 'w') # write to file called 'u_function.h5'
ldc_function_outfile.write(o3_ldc_soln.velocity, "/u_o3") 
ldc_function_outfile.close()

ldc_function_outfile = HDF5File(MPI.comm_world, 'results_compare_ldc/ldc_functions.h5', 'a') # 'a' for 'append' to file
ldc_function_outfile.write(nse_ldc_soln.velocity, "/u_nse") 
ldc_function_outfile.write(o3_ldc_soln.pressure, "/p_o3")
ldc_function_outfile.write(nse_ldc_soln.pressure, "/p_nse")
ldc_function_outfile.close()



jb_function_outfile = HDF5File(MPI.comm_world, 'results_compare_jb/jb_functions.h5', 'w') # write to file called 'u_function.h5'
jb_function_outfile.write(o3_jb_soln.velocity, "/u_o3") 
jb_function_outfile.close()

jb_function_outfile = HDF5File(MPI.comm_world, 'results_compare_jb/jb_functions.h5', 'a') # 'a' for 'append' to file
jb_function_outfile.write(nse_jb_soln.velocity, "/u_nse") 
jb_function_outfile.write(o3_jb_soln.pressure, "/p_o3")
jb_function_outfile.write(nse_jb_soln.pressure, "/p_nse")
jb_function_outfile.close()


