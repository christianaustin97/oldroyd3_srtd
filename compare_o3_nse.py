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
    print("Creating ucm_nse_ldc.csv file ...")
    outfile = open('results_compare_ldc/ucm_nse_ldc.csv', 'a')
    writer = csv.writer(outfile)
    writer.writerow(['meshsize', 'eta', 'lambda1', 'l2_diff', 'h1_diff'])
    outfile.close()

if not (os.path.isfile('results_compare_jb/ucm_nse_jb.csv')):
    print("Creating ucm_nse_jb.csv file ...")
    outfile = open('results_compare_jb/ucm_nse_jb.csv', 'a')
    writer = csv.writer(outfile)
    writer.writerow(['meshsize', 'radius', 'eccentricity', 'eta', 'lambda1', 'l2_diff', 'h1_diff'])
    outfile.close()

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

# Differing l1 values to compare to NSE
coeffs = np.arange(1,10,2) # 1,3,5,7,9
pows = np.power(10.0, np.arange(-6,-1)) # 10^-6, 10^-5, ..., 10^-2
l1vals = np.concatenate((coeffs*pows[0], coeffs*pows[1], coeffs*pows[2], coeffs*pows[3])) #1e-6, 3e-6, ..., 9e-6, 1e-5,..., 1e-3,...9e-3, 

# L2 and H1 differences in UCM and 
ldc_l2_differences = np.zeros(len(l1vals))
ldc_h1_differences = np.zeros(len(l1vals))
bearing_l2_differences = np.zeros(len(l1vals))
bearing_h1_differences = np.zeros(len(l1vals))


