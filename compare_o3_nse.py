"""
Does some comparisons of Oldroyd 3 solution to the NSE/Newtonian solution
"""

from fenics import *
from steady_nse_solver import *
from oldroyd_3_EVSS import *
from oldroyd_3_SRTD import *

import os
import sys
import matplotlib.pyplot as plt
import numpy as np

