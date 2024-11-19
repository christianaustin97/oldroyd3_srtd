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

max_srtd_iters = 20
srtd_tol = 1e-8

Wi_tol = 1e-4 # how close to we want our approx of max Wi to be

# JB geometry parameters
rad = 0.5
ecc = 0.25

speed=1.0

# EVSS l1 limits
evss_jb_l1 = [0.05, 0.2]
evss_ldc_l1 = [0.05, 0.15]

# SRTD l1 limits
srtd_jb_l1 = [0.01, 0.1]
srtd_ldc_l1 = [0.01, 0.1]

n1 = np.ceil(np.log2((evss_jb_l1[1]-evss_jb_l1[0])/Wi_tol)).astype(int)
n2 = np.ceil(np.log2((srtd_jb_l1[1]-srtd_jb_l1[0])/Wi_tol)).astype(int)
n=max([n1, n2])

for i in range(n):
    # Binary search 
    current_evss_jb_l1 = (evss_jb_l1[0] + evss_jb_l1[1])/2.0
    current_evss_ldc_l1 = (evss_ldc_l1[0] + evss_ldc_l1[1])/2.0
    current_srtd_jb_l1 = (srtd_jb_l1[0] + srtd_jb_l1[1])/2.0
    current_srtd_ldc_l1 = (srtd_ldc_l1[0] + srtd_ldc_l1[1])/2.0

    a = 0.0 # corotational Maxwell model

    #print("current l1 for evss jb: %.5f"%current_evss_jb_l1)
    print("current l1 for evss ldc: %.5f"%current_evss_ldc_l1)
    #print("current l1 for srtd jb: %.5f"%current_srtd_jb_l1)
    #print("current l1 for srtd ldc: %.5f"%current_srtd_ldc_l1)


    #evss_jb = oldroyd_3_EVSS.oldroyd_3_JB_EVSS(h_jb, rad, ecc, speed, eta, current_evss_jb_l1, a*current_evss_jb_l1)
    evss_ldc = oldroyd_3_EVSS.oldroyd_3_LDC_EVSS(h_ldc, speed, eta, current_evss_ldc_l1, a*current_evss_ldc_l1)

    #srtd_jb = oldroyd_3_SRTD.oldroyd_3_JB_SRTD(h_jb, rad, ecc, speed, eta, current_srtd_jb_l1, a*current_srtd_jb_l1, max_srtd_iters, srtd_tol)"""
    #srtd_ldc = oldroyd_3_SRTD.oldroyd_3_LDC_SRTD(h_ldc, speed, eta, current_srtd_ldc_l1, a*current_srtd_ldc_l1, max_srtd_iters, srtd_tol)

    """if(evss_jb.converged):
        evss_jb_l1[0] = current_evss_jb_l1
    else:
        evss_jb_l1[1] = current_evss_jb_l1"""

    if(evss_ldc.converged):
        evss_ldc_l1[0] = current_evss_ldc_l1
    else:
        evss_ldc_l1[1] = current_evss_ldc_l1

    """if(srtd_jb.converged):
        srtd_jb_l1[0] = current_srtd_jb_l1
    else:
        srtd_jb_l1[1] = current_srtd_jb_l1"""

    """if(srtd_ldc.converged):
        srtd_ldc_l1[0] = current_srtd_ldc_l1
    else:
        srtd_ldc_l1[1] = current_srtd_ldc_l1"""

    
    #print("Max l1 for EVSS, JB, speed=%.1f, is between %.5f and %.5f"%(speed, evss_jb_l1[0], evss_jb_l1[1]))
    print("Max l1 for EVSS, LDC, speed=%.1f, is between %.5f and %.5f"%(speed, evss_ldc_l1[0], evss_ldc_l1[1]))
    #print("Max l1 for SRTD, JB, speed=%.1f, is between %.5f and %.5f"%(speed, srtd_jb_l1[0], srtd_jb_l1[1]))
    #print("Max l1 for SRTD, LDC, speed=%.1f, is between %.5f and %.5f"%(speed, srtd_ldc_l1[0], srtd_ldc_l1[1]))

