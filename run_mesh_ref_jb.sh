#!/bin/bash

# First arg is speed s, 2nd arg is viscosity eta0, 3rd arg is relaxation time lambda1

# Wi = 2e-2: s=1.0, l1=1e-2
python3 mesh_ref_jb_ucm_evss.py 1.0 1.0 1e-2
python3 mesh_ref_jb_ucm_srtd.py 1.0 1.0 1e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 1.0 1.0 1e-2
# Wi = 2e-2 : s=0.5, l1 = 2e-2
python3 mesh_ref_jb_ucm_evss.py 0.5 1.0 2e-2
python3 mesh_ref_jb_ucm_srtd.py 0.5 1.0 2e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 0.5 1.0 2e-2
# Wi = 2e-2 : s=0.25, l1 = 4e-2
python3 mesh_ref_jb_ucm_evss.py 0.25 1.0 4e-2
python3 mesh_ref_jb_ucm_srtd.py 0.25 1.0 4e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 0.25 1.0 4e-2

# Wi = 4e-2 : s=2, l1=1e-2
python3 mesh_ref_jb_ucm_evss.py 2.0 1.0 1e-2
python3 mesh_ref_jb_ucm_srtd.py 2.0 1.0 1e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 2.0 1.0 1e-2

# Wi = 4e-2 : s=1, l1=2e-2
python3 mesh_ref_jb_ucm_evss.py 1.0 1.0 2e-2
python3 mesh_ref_jb_ucm_srtd.py 1.0 1.0 2e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 1.0 1.0 2e-2


# Wi = 8e-2 : s=2, l1=2e-2
python3 mesh_ref_jb_ucm_evss.py 2.0 1.0 2e-2
python3 mesh_ref_jb_ucm_srtd.py 2.0 1.0 2e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 2.0 1.0 2e-2

# Wi = 8e-2 : s=4, l1=1e-2
python3 mesh_ref_jb_ucm_evss.py 4.0 1.0 1e-2
python3 mesh_ref_jb_ucm_srtd.py 4.0 1.0 1e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 4.0 1.0 1e-2

# Wi = 8e-2 : s=1, l1=4e-2
python3 mesh_ref_jb_ucm_evss.py 1.0 1.0 4e-2
python3 mesh_ref_jb_ucm_srtd.py 1.0 1.0 4e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 1.0 1.0 4e-2


#SRT appears to break down around here lol, Wi=0.08 for the LDC. Don't know how comparable that is for JB yet
# Wi = 1.6e-1 : s=4, l1=2e-2
python3 mesh_ref_jb_ucm_evss.py 4.0 1.0 2e-2
python3 mesh_ref_jb_ucm_srtd.py 4.0 1.0 2e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 4.0 1.0 2e-2
# Wi = 1.6e-1 : s=2, l1=4e-2
python3 mesh_ref_jb_ucm_evss.py 2.0 1.0 4e-2
python3 mesh_ref_jb_ucm_srtd.py 2.0 1.0 4e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 2.0 1.0 4e-2
# Wi = 1.6e-1 : s=1, l1=8e-2
python3 mesh_ref_jb_ucm_evss.py 1.0 1.0 8e-2
python3 mesh_ref_jb_ucm_srtd.py 1.0 1.0 8e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 1.0 1.0 8e-2
# Wi = 1.6e-1 : s=8, l1=1e-2
python3 mesh_ref_jb_ucm_evss.py 8.0 1.0 1e-2
python3 mesh_ref_jb_ucm_srtd.py 8.0 1.0 1e-2
python3 mesh_ref_jb_ucm_srtd_supg.py 8.0 1.0 1e-2


# Wi = 3.2e-1 : s=2, l1=8e-2
python3 mesh_ref_jb_ucm_evss.py 2.0 1.0 8e-2
# Wi = 3.2e-1 : s=4, l1=4e-2
python3 mesh_ref_jb_ucm_evss.py 4.0 1.0 4e-2
# Wi = 3.2e-1 : s=8, l1=2e-2
python3 mesh_ref_jb_ucm_evss.py 8.0 1.0 2e-2


# Wi = 6.4e-1 : s=4, l1=8e-2
python3 mesh_ref_jb_ucm_evss.py 4.0 1.0 8e-2
# Wi = 6.4e-1 : s=8, l1=4e-2
python3 mesh_ref_jb_ucm_evss.py 8.0 1.0 4e-2

# EVSS appears to break down around here. Good at Wi=0.64. Could probably be improved with a better starting guess for Newton? 

# Wi = 1.28e+0 : s=4, l1=1.6e-1
python3 mesh_ref_jb_ucm_evss.py 4.0 1.0 1.6e-1
# Wi = 1.28e+0 : s=8, l1=8e-2
python3 mesh_ref_jb_ucm_evss.py 8.0 1.0 8e-2
# Wi = 1.28e+0 : s=16, l1=4e-2
python3 mesh_ref_jb_ucm_evss.py 16.0 1.0 4e-2

# Wi = 2.56e+0 : s=8, l1=1.6e-1
python3 mesh_ref_jb_ucm_evss.py 8.0 1.0 1.6e-1
# Wi = 2.56e+0 : s=16, l1=8e-2
python3 mesh_ref_jb_ucm_evss.py 16.0 1.0 8e-2

