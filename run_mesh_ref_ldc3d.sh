#!/bin/bash

# First arg is speed s, 2nd arg is viscosity eta0, 3rd arg is relaxation time lambda1, 4th is slip param a
#   a=1 : ucm
#   a=0 : corotational maxwell
#   a=-1: lcm

python3 mesh_ref_ldc3d_srtd.py 1.0 1.0 1e-2 1.0 4

python3 mesh_ref_ldc3d_srtd.py 1.0 1.0 2e-2 1.0 4

python3 mesh_ref_ldc3d_srtd.py 1.0 1.0 3e-2 1.0 4
