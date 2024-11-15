"""
Generates some plots and tables for number of SRTD iters it takes to get to 
different tolerance values/residual values/difference between consecutive iterates
for the UCM model
"""

from fenics import *
from oldroyd_3_SRTD import *

import os
import sys
import matplotlib.pyplot as plt
import numpy as np

from math import log2 as log2 # For computing the rate
import time # Timing the computations
import csv # Saving Results


jb_file = open('results_num_iters/jb_num_iters.csv', 'w') 
jb_writer = csv.writer(jb_file)
jb_writer.writerow(["l1", "1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8", "1e-9"])

jb_residuals_file = open('results_num_iters/jb_residuals.csv', 'w') 
jb_residuals_writer = csv.writer(jb_residuals_file)
jb_residuals_writer.writerow(["l1", "tol_at_iter"])

# meshsize
h_ldc = 1.25e-2
h_jb = 2.5e-2

# characteristic speed
speed=1.0 

eta = 1.0

coeffs = np.array([0.25, 0.5, 0.75, 1.0])
bases = 10.0**np.array([-5, -4, -3, -2, -1, 0])
l1vals = np.array([coeffs*b for b in bases]).ravel() # opposite of unravel i guess lol. combines into 1d array

# geometry parameters for jb
rad = 0.5
ecc = 0.25

# SRTD algorithm parameters
max_srtd_iters = 25
srtd_tol = 1e-9 
#srtd_tol = -1.0 # error is capped

# tolerances to check, 1e-1, 1e-3, ..., 1e-9
tols = 10.0**np.array([-1, -2, -3, -4, -5, -6, -7, -8, -9])

# Journal-bearing problem first I guess. Record number of iters required to reach each tol
for l1 in l1vals:

    lambda1 = float(l1)
    out_row = ['%.3e'%lambda1] + ['']*9
    mu1 = lambda1

    solution = oldroyd_3_JB_SRTD(h_jb, rad, ecc, speed, eta, lambda1, mu1, max_srtd_iters, srtd_tol)

    res_dict = solution.residuals
    iterations = np.array(list(res_dict.keys()))
    residuals = np.array(list(res_dict.values()))

    print(residuals) #maybe do plots here

    for i in range(9):
        if(np.count_nonzero(residuals < tols[i])>0): # if there is a residual smaller than this tol,
            out_row[i+1] = str(iterations[residuals < tols[i]][0]) # returns first iteration for which residuals < that tol
        else:
            out_row[i+1] = " - " # otherwise, just put a dash
    
    jb_writer.writerow(out_row)
    residuals_out_row = np.concatenate([[l1], residuals])
    formatted_residuals_out_row = ["%.4e"%num for num in residuals_out_row]
    #print("Test:")
    #print(formatted_residuals_out_row)
    jb_residuals_writer.writerow(formatted_residuals_out_row)

    # plot residuals
    plt.semilogy(iterations, residuals)
    plt.xticks(range(1, max_srtd_iters+1))
    plt.xlabel("SRTD Iteration")
    plt.ylabel("Residual")
    plt.savefig("results_num_iters/jb_l1_%.3e_all_20.pdf"%lambda1)
    plt.close()

    # plot residuals
    plt.semilogy(iterations, residuals)
    plt.xticks(iterations)
    plt.xlabel("SRTD Iteration")
    plt.ylabel("Residual")
    plt.savefig("results_num_iters/jb_l1_%.3e.pdf"%lambda1)
    plt.close()

jb_file.close()
jb_residuals_file.close()





# LDC problem next
ldc_file = open('results_num_iters/ldc_num_iters.csv', 'w') 
ldc_writer = csv.writer(ldc_file)
ldc_writer.writerow(["l1", "1-e1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8", "1e-9"])

ldc_residuals_file = open('results_num_iters/ldc_residuals.csv', 'w') 
ldc_residuals_writer = csv.writer(ldc_residuals_file)
ldc_residuals_writer.writerow(["l1", "tol_at_iter"])

# ldc problem next. Record number of iters required to reach each tol
for l1 in l1vals:
    
    lambda1 = float(l1)
    out_row = ['%.3e'%lambda1] + ['']*9
    mu1 = lambda1

    solution_ldc = oldroyd_3_LDC_SRTD(h_ldc, speed, eta, lambda1, mu1, max_srtd_iters, srtd_tol)

    res_dict = solution_ldc.residuals
    iterations = np.array(list(res_dict.keys()))
    residuals = np.array(list(res_dict.values()))

    print(residuals)

    for i in range(9):
        if(np.count_nonzero(residuals < tols[i])>0): # if there is a residual smaller than this tol,
            out_row[i+1] = str(iterations[residuals < tols[i]][0]) # returns first iteration for which residuals < that tol
        else:
            out_row[i+1] = " - " # otherwise, just put a dash
    

    ldc_writer.writerow(out_row)
    residuals_out_row = np.concatenate([[l1], residuals])
    formatted_residuals_out_row = ["%.4e"%num for num in residuals_out_row]
    #print("Test:")
    #print(formatted_residuals_out_row)
    ldc_residuals_writer.writerow(formatted_residuals_out_row)

    # plot residuals
    plt.semilogy(iterations, residuals)
    plt.xticks(range(1, max_srtd_iters+1))
    plt.xlabel("SRTD Iteration")
    plt.ylabel("Residual")
    plt.savefig("results_num_iters/ldc_l1_%.3e_all_20.pdf"%lambda1)
    plt.close()

    # plot residuals
    plt.semilogy(iterations, residuals)
    plt.xticks(iterations)
    plt.xlabel("SRTD Iteration")
    plt.ylabel("Residual")
    plt.savefig("results_num_iters/ldc_l1_%.3e.pdf"%lambda1)
    plt.close()

ldc_file.close()
ldc_residuals_file.close()