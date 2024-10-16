#!/bin/bash

# Runs the mesh refinement experiment for various Weissenberg numbers/combinations, JB problem

# First arg is speed s, 2nd arg is viscosity eta0, 3rd arg is relaxation time lambda1, 4th is slip param a
#   a=1 : ucm
#   a=0 : corotational maxwell
#   a=-1: lcm

# Wi = 2e-2: s=1.0, l1=1e-2
python3 mesh_ref_jb_evss.py 1.0 1.0 1e-2 0.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 1e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 1.0 1.0 1e-2 0.0
# Wi = 2e-2 : s=0.5, l1 = 2e-2
python3 mesh_ref_jb_evss.py 0.5 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd.py 0.5 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 0.5 1.0 2e-2 0.0
# Wi = 2e-2 : s=0.25, l1 = 4e-2
python3 mesh_ref_jb_evss.py 0.25 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd.py 0.25 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 0.25 1.0 4e-2 0.0

# Wi = 4e-2 : s=2, l1=1e-2
python3 mesh_ref_jb_evss.py 2.0 1.0 1e-2 0.0
python3 mesh_ref_jb_srtd.py 2.0 1.0 1e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 2.0 1.0 1e-2 0.0

# Wi = 4e-2 : s=1, l1=2e-2
python3 mesh_ref_jb_evss.py 1.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 1.0 1.0 2e-2 0.0


# Wi = 8e-2 : s=2, l1=2e-2
python3 mesh_ref_jb_evss.py 2.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd.py 2.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 2.0 1.0 2e-2 0.0

# Wi = 8e-2 : s=4, l1=1e-2
python3 mesh_ref_jb_evss.py 4.0 1.0 1e-2 0.0
python3 mesh_ref_jb_srtd.py 4.0 1.0 1e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 4.0 1.0 1e-2 0.0

# Wi = 8e-2 : s=1, l1=4e-2
python3 mesh_ref_jb_evss.py 1.0 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 1.0 1.0 4e-2 0.0


#SRTD breaks down around here for UCM, Wi=0.08 for the LDC problem. EVSS still works okay
# Wi = 1.6e-1 : s=4, l1=2e-2
python3 mesh_ref_jb_evss.py 4.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd.py 4.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 4.0 1.0 2e-2 0.0
# Wi = 1.6e-1 : s=2, l1=4e-2
python3 mesh_ref_jb_evss.py 2.0 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd.py 2.0 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 2.0 1.0 4e-2 0.0
# Wi = 1.6e-1 : s=1, l1=8e-2
python3 mesh_ref_jb_evss.py 1.0 1.0 8e-2 0.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 8e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 1.0 1.0 8e-2 0.0
# Wi = 1.6e-1 : s=8, l1=1e-2
python3 mesh_ref_jb_evss.py 8.0 1.0 1e-2 0.0
python3 mesh_ref_jb_srtd.py 8.0 1.0 1e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 8.0 1.0 1e-2 0.0


# Wi = 3.2e-1 : s=2, l1=8e-2
python3 mesh_ref_jb_evss.py 2.0 1.0 8e-2 0.0
python3 mesh_ref_jb_srtd.py 2.0 1.0 8e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 2.0 1.0 8e-2 0.0

# Wi = 3.2e-1 : s=4, l1=4e-2 
python3 mesh_ref_jb_evss.py 4.0 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd.py 4.0 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 4.0 1.0 4e-2 0.0
# Wi = 3.2e-1 : s=8, l1=2e-2
python3 mesh_ref_jb_evss.py 8.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd.py 8.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd_supg.py 8.0 1.0 2e-2 0.0


# Wi = 6.4e-1 : s=4, l1=8e-2
python3 mesh_ref_jb_evss.py 4.0 1.0 8e-2 0.0
# Wi = 6.4e-1 : s=8, l1=4e-2
python3 mesh_ref_jb_evss.py 8.0 1.0 4e-2 0.0
