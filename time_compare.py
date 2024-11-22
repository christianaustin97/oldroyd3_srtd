from fenics import *
from oldroyd_3_SRTD import *
from oldroyd_3_EVSS import *
import time




start_solve = time.time()
#solution = oldroyd_3_JB_SRTD(h, rad, ecc, speed, eta, lambda1, mu1, max_srtd_iters, srtd_tol)
end_solve = time.time()
solve_time = end_solve - start_solve

