#!/bin/bash

# Runs the mesh refinement experiment for various Weissenberg numbers/combinations, JB problem

# First arg is speed s, 2nd arg is viscosity eta0, 3rd arg is relaxation time lambda1, 4th is slip param a
#   a=1 : ucm
#   a=0 : corotational maxwell
#   a=-1: lcm

# Weissenberg number limit for SRTD is around 0.06, while for EVSS only slightly higher at 0.09

# Wi = 5e-2, s=1, l1=5e-2
python3 mesh_ref_ldc_srtd.py 1.0 1.0 5e-2 0.0

#Wi = 6e-2, s=1, l1=6e-2
python3 mesh_ref_ldc_evss.py 1.0 1.0 6e-2 0.0
python3 mesh_ref_ldc_srtd.py 1.0 1.0 6e-2 0.0

#Wi = 8e-2, s=1, l1=6e-2
python3 mesh_ref_ldc_evss.py 1.0 1.0 8e-2 0.0





<<comment_flag
echo "this is inside a comment, so should not run"
comment_flag



