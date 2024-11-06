# This program attempts to find the Weissenberg number attainable precisely using a binary search
# just a rough estimate, so 

from fenics import *
import oldroyd_3_SRTD 
import oldroyd_3_EVSS
import time
import matplotlib.pyplot as plt
import numpy as np



h_jb=0.025
h_ldc = h_jb/2.0

eta=1.0
a=1.0 # 1 for UCM, 0 for corotational

speed=1.0

max_srtd_iters = 20
srtd_tol = 1e-8

Wi_tol = 1e-4 # how close to we want our approx of max Wi to be

# JB geometry parameters
rad = 0.5
ecc = 0.25


# SRTD, JB problem
l1_low_jb = 0.1
l1_high_jb = 2.0

l1_low_ldc = 0.05
l1_high_ldc = 1.5
n = np.ceil(np.log2((l1_high_jb-l1_low_jb)/Wi_tol)).astype(int)

for i in range(n):
    # Binary search 
    l1_jb = (l1_high_jb+l1_low_jb)/2.0
    mu1_jb = a*l1_jb
    print("current l1 for jb: %.5f"%l1_jb)

    l1_ldc = (l1_high_ldc+l1_low_ldc)/2.0
    mu1_ldc = a*l1_ldc
    print("current l1 for ldc: %.5f"%l1_ldc)

    soln_jb = oldroyd_3_EVSS.oldroyd_3_JB_EVSS(h_jb, rad, ecc, speed, eta, l1_jb, mu1_jb)
    soln_ldc = oldroyd_3_EVSS.oldroyd_3_LDC_EVSS(h_ldc, speed, eta, l1_ldc, mu1_ldc)

    #soln_jb = oldroyd_3_SRTD.oldroyd_3_JB_SRTD(h_jb, rad, ecc, speed, eta, l1_jb, mu1_jb, max_srtd_iters, srtd_tol)
    #soln_ldc = oldroyd_3_SRTD.oldroyd_3_LDC_SRTD(h_ldc, speed, eta, l1_ldc, mu1_ldc, max_srtd_iters, srtd_tol)
    if(soln_jb.converged):
        l1_low_jb = l1_jb
    else:
        l1_high_jb = l1_jb

    if(soln_ldc.converged):
        l1_low_ldc = l1_ldc
    else:
        l1_high_ldc = l1_ldc

    
    print("Max l1 for EVSS, JB, speed=%.1f, is between %.5f and %.5f"%(speed, l1_low_jb, l1_high_jb))
    print("Max l1 for EVSS, LDC, speed=%.1f, is between %.5f and %.5f"%(speed, l1_low_ldc, l1_high_ldc))
# 
