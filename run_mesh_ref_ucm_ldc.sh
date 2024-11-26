#!/bin/bash

# Runs the mesh refinement experiment for various Weissenberg numbers/combinations, JB problem

# First arg is speed s, 2nd arg is viscosity eta0, 3rd arg is relaxation time lambda1, 4th is slip param a
#   a=1 : ucm
#   a=0 : corotational maxwell
#   a=-1: lcm

# UCM, l1=0.1
python3 mesh_ref_ldc_evss.py 1.0 1.0 1e-1 1.0  

# UCM, l1=0.01
python3 mesh_ref_ldc_srtd.py 1.0 1.0 1e-2 1.0  

# 3D LDC SRTD, UCM, l1=0.01
#python3 mesh_ref_ldc3d_srtd.py 1.0 1.0 1e-2 1.0  


# Closer to Wi limit
python3 mesh_ref_ldc_evss.py 1.0 1.0 1.5e-1 1.0 
python3 mesh_ref_ldc_evss.py 1.0 1.0 2e-1 1.0 
python3 mesh_ref_ldc_evss.py 1.0 1.0 2.5e-1 1.0 

python3 mesh_ref_ldc_srtd.py 1.0 1.0 5.5e-2 1.0 
python3 mesh_ref_ldc_srtd.py 1.0 1.0 6e-2 1.0 