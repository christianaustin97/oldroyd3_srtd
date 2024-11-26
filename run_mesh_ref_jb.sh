#!/bin/bash

# Runs the mesh refinement experiment for various Weissenberg numbers/combinations, JB problem

# First arg is speed s, 2nd arg is viscosity eta0, 3rd arg is relaxation time lambda1, 4th is slip param a
#   a=1 : ucm
#   a=0 : corotational maxwell
#   a=-1: lcm

# SRTD UCM and corot same approximate Wi limit around 0.1
# SRTD UCM
python3 mesh_ref_jb_srtd.py 1.0 1.0 2e-2 1.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 3e-2 1.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 4e-2 1.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 5e-2 1.0

# SRTD Corot
python3 mesh_ref_jb_srtd.py 1.0 1.0 2e-2 0.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 3e-2 0.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 4e-2 0.0
python3 mesh_ref_jb_srtd.py 1.0 1.0 5e-2 0.0


# EVSS UCM is around 2.0 (l1 = 1.0). EVSS corot around 0.2 (l1=0.1)
# EVSS UCM
python3 mesh_ref_jb_evss.py 1.0 1.0 5e-2 1.0
python3 mesh_ref_jb_evss.py 1.0 1.0 1e-1 1.0
python3 mesh_ref_jb_evss.py 1.0 1.0 5e-1 1.0
python3 mesh_ref_jb_evss.py 1.0 1.0 1.0 1.0

# EVSS Corot
python3 mesh_ref_jb_evss.py 1.0 1.0 1e-2 0.0
python3 mesh_ref_jb_evss.py 1.0 1.0 5e-2 0.0
python3 mesh_ref_jb_evss.py 1.0 1.0 8e-2 1.0
python3 mesh_ref_jb_evss.py 1.0 1.0 1e-1 0.0



