from fenics import *
import matplotlib.pyplot as plt
import numpy as np
import csv


with open('results_compare_ldc/ucm_nse_ldc.csv', newline='') as ldc_csv:
    ldc_data = list(csv.reader(ldc_csv))

with open('results_compare_jb/ucm_nse_jb.csv', newline='') as jb_csv:
    jb_data = list(csv.reader(jb_csv))

l1_ldc_vals = [float(row[2]) for row in ldc_data[1:-1]] # ignore header. 2nd column is l1 vals
l2diff_ldc_vals = [float(row[3]) for row in ldc_data[1:-1]]
h1diff_ldc_vals =  [float(row[4]) for row in ldc_data[1:-1]]
l2pdiff_ldc_vals = [float(row[5]) for row in ldc_data[1:-1]]


l1_jb_vals = [float(row[4]) for row in jb_data[1:-1]] # ignore header. 4th column is l1 vals
l2diff_jb_vals = [float(row[5]) for row in jb_data[1:-1]]
h1diff_jb_vals =  [float(row[6]) for row in jb_data[1:-1]]
l2pdiff_jb_vals = [float(row[7]) for row in jb_data[1:-1]]

# LDC plots
p_ldc_l2 = plt.loglog(l1_ldc_vals, l2diff_ldc_vals, marker='o', markerfacecolor='none')
plt.ylabel("$L^{2}$ diff between UCM and Newtonian velocity")
plt.xlabel("$\lambda_{1}$ value")

###########################################################################
# # add slopes using polyfit
# before elbow
elbow = 13

l2_logx1 = np.log(l1_ldc_vals[0:elbow])
l2_logy1 = np.log(l2diff_ldc_vals[0:elbow])
l2_pfit_coeffs1 = np.polyfit(l2_logx1, l2_logy1, 1)
l2_pfit_poly1 = np.poly1d(l2_pfit_coeffs1)

# after elbow
l2_logx2 = np.log(l1_ldc_vals[elbow:19])
l2_logy2 = np.log(l2diff_ldc_vals[elbow:19])
l2_pfit_coeffs2 = np.polyfit(l2_logx2, l2_logy2, 1)
l2_pfit_poly2 = np.poly1d(l2_pfit_coeffs2)

# function to fit linear polynomial
l2_fit_1 = lambda x: np.exp(l2_pfit_poly1(np.log(x)))
l2_fit_2 = lambda x: np.exp(l2_pfit_poly2(np.log(x)))

# add to plots
plt.loglog(l1_ldc_vals[0:elbow], l2_fit_1(l1_ldc_vals[0:elbow]), '--')
plt.loglog(l1_ldc_vals[elbow:19], l2_fit_2(l1_ldc_vals[elbow:19]), '--')

plt.legend(["L2 difference", "Slope=%1.3f"%l2_pfit_coeffs1[0], "Slope=%1.3f"%l2_pfit_coeffs2[0]])
###########################################################################

plt.savefig("results_compare_ldc/ucm_nse_ldc_l2.pdf")
plt.close()

p_ldc_h1 = plt.loglog(l1_ldc_vals, h1diff_ldc_vals, marker='o', markerfacecolor='none')
plt.ylabel("$H^{1}$ diff between UCM and Newtonian velocity")
plt.xlabel("$\lambda_{1}$ value")

###########################################################################
# add slopes using polyfit
# before elbow
elbow = 13

h1_logx1 = np.log(l1_ldc_vals[0:elbow])
h1_logy1 = np.log(h1diff_ldc_vals[0:elbow])
h1_pfit_coeffs1 = np.polyfit(h1_logx1, h1_logy1, 1)
h1_pfit_poly1 = np.poly1d(h1_pfit_coeffs1)

# after elbow
h1_logx2 = np.log(l1_ldc_vals[elbow:19])
h1_logy2 = np.log(h1diff_ldc_vals[elbow:19])
h1_pfit_coeffs2 = np.polyfit(h1_logx2, h1_logy2, 1)
h1_pfit_poly2 = np.poly1d(h1_pfit_coeffs2)

# function to fit linear polynomial
h1_fit_1 = lambda x: np.exp(h1_pfit_poly1(np.log(x)))
h1_fit_2 = lambda x: np.exp(h1_pfit_poly2(np.log(x)))

# add to plots
plt.loglog(l1_ldc_vals[0:elbow], h1_fit_1(l1_ldc_vals[0:elbow]), '--')
plt.loglog(l1_ldc_vals[elbow:19], h1_fit_2(l1_ldc_vals[elbow:19]), '--')

plt.legend(["H1 difference", "Slope=%1.3f"%h1_pfit_coeffs1[0], "Slope=%1.3f"%h1_pfit_coeffs2[0]])
###########################################################################

plt.savefig("results_compare_ldc/ucm_nse_ldc_h1.pdf")
plt.close()


p_ldc_l2p = plt.loglog(l1_ldc_vals, l2pdiff_ldc_vals, marker='o', markerfacecolor='none')
plt.ylabel("$L^{2}$ diff between UCM and Newtonian pressure")
plt.xlabel("$\lambda_{1}$ value")
plt.savefig("results_compare_ldc/ucm_nse_ldc_l2p.pdf")
plt.close()


# JB Plots
p_jb_l2 = plt.loglog(l1_jb_vals, l2diff_jb_vals, marker='o', markerfacecolor='none')
plt.ylabel("$L^{2}$ diff between UCM and Newtonian velocity")
plt.xlabel("$\lambda_{1}$ value")


###########################################################################
# # add slopes using polyfit
# before elbow
elbow = 15

l2_logx1 = np.log(l1_jb_vals[0:elbow])
l2_logy1 = np.log(l2diff_jb_vals[0:elbow])
l2_pfit_coeffs1 = np.polyfit(l2_logx1, l2_logy1, 1)
l2_pfit_poly1 = np.poly1d(l2_pfit_coeffs1)

# after elbow
l2_logx2 = np.log(l1_jb_vals[elbow:19])
l2_logy2 = np.log(l2diff_jb_vals[elbow:19])
l2_pfit_coeffs2 = np.polyfit(l2_logx2, l2_logy2, 1)
l2_pfit_poly2 = np.poly1d(l2_pfit_coeffs2)

# function to fit linear polynomial
l2_fit_1 = lambda x: np.exp(l2_pfit_poly1(np.log(x)))
l2_fit_2 = lambda x: np.exp(l2_pfit_poly2(np.log(x)))

# add to plots
plt.loglog(l1_jb_vals[0:elbow], l2_fit_1(l1_jb_vals[0:elbow]), '--')
plt.loglog(l1_jb_vals[elbow:19], l2_fit_2(l1_jb_vals[elbow:19]), '--')

plt.legend(["L2 difference", "Slope=%1.3f"%l2_pfit_coeffs1[0], "Slope=%1.3f"%l2_pfit_coeffs2[0]])
###########################################################################



plt.savefig("results_compare_jb/ucm_nse_jb_l2.pdf")
plt.close()

p_jb_h1 = plt.loglog(l1_jb_vals, h1diff_jb_vals, marker='o', markerfacecolor='none')
plt.ylabel("$H^{1}$ diff between UCM and Newtonian velocity")
plt.xlabel("$\lambda_{1}$ value")

###########################################################################
# # add slopes using polyfit
# before elbow
elbow = 19

h1_logx1 = np.log(l1_jb_vals[0:elbow])
h1_logy1 = np.log(h1diff_jb_vals[0:elbow])
h1_pfit_coeffs1 = np.polyfit(h1_logx1, h1_logy1, 1)
h1_pfit_poly1 = np.poly1d(h1_pfit_coeffs1)

"""# after elbow"
h1_logx2 = np.log(l1_jb_vals[elbow:19])
h1_logy2 = np.log(h1diff_jb_vals[elbow:19])
h1_pfit_coeffs2 = np.polyfit(h1_logx2, h1_logy2, 1)
h1_pfit_poly2 = np.poly1d(h1_pfit_coeffs2)"""

# function to fit linear polynomial
h1_fit_1 = lambda x: np.exp(h1_pfit_poly1(np.log(x)))
#h1_fit_2 = lambda x: np.exp(h1_pfit_poly2(np.log(x)))

# add to plots
plt.loglog(l1_jb_vals[0:elbow], h1_fit_1(l1_jb_vals[0:elbow]), '--')
#plt.loglog(l1_jb_vals[elbow:19], h1_fit_2(l1_jb_vals[elbow:19]), '--')

plt.legend(["H1 difference", "Slope=%1.3f"%h1_pfit_coeffs1[0], "Slope=%1.3f"%h1_pfit_coeffs2[0]])
###########################################################################


plt.savefig("results_compare_jb/ucm_nse_jb_h1.pdf")
plt.close()

p_jb_l2p = plt.loglog(l1_jb_vals, l2pdiff_jb_vals, marker='o', markerfacecolor='none')
plt.ylabel("$L^{2}$ diff between UCM and Newtonian pressure")
plt.xlabel("$\lambda_{1}$ value")
plt.savefig("results_compare_jb/ucm_nse_jb_l2p.pdf")
plt.close()


