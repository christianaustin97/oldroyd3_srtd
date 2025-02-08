#!/bin/bash

# Runs the mesh refinement experiment near the Weissenberg number limit (WNL), JB problem

# Wi limit for EVSS UCM around 2.0 (U=1, l1=1)
# Wi limit for SRTD UCM around 0.1 (U=1, l1=0.05)

# Wi limit for EVSS Corot around 0.22 (U=1, l1=0.11)
# Wi limit for SRTD Corot around 0.11 (U=1, l1=0.05)

# First arg is speed s, 2nd arg is viscosity eta0, 3rd arg is relaxation time lambda1, 4th is slip param a
#   a=1 : ucm
#   a=0 : corotational maxwell
#   a=-1: lcm

python3 mesh_ref_jb_srtd.py 1.0 1.0 6e-2 1.0

python3 mesh_ref_jb_evss.py 1.0 1.0 7.5e-1 0.0
