from fenics import *
import matplotlib.pyplot as plt
import numpy as np
import csv


with open('results_compare_ldc/ucm_nse_ldc.csv', newline='') as ldc_csv:
    ldc_data = list(csv.reader(ldc_csv))

with open('results_compare_jb/ucm_nse_jb.csv', newline='') as jb_csv:
    jb_data = list(csv.reader(jb_csv))

l1_ldc_vals = [float(row[2]) for row in ldc_data[1:]] # ignore header. 2nd column is l1 vals
l2diff_ldc_vals = [float(row[3]) for row in ldc_data[1:]]
h1diff_ldc_vals =  [float(row[4]) for row in ldc_data[1:]]
l2pdiff_ldc_vals = [float(row[5]) for row in ldc_data[1:]]


l1_jb_vals = [float(row[4]) for row in jb_data[1:]] # ignore header. 4th column is l1 vals
l2diff_jb_vals = [float(row[5]) for row in jb_data[1:]]
h1diff_jb_vals =  [float(row[6]) for row in jb_data[1:]]
l2pdiff_jb_vals = [float(row[7]) for row in jb_data[1:]]

# LDC plots
p_ldc_l2 = plt.loglog(l1_ldc_vals, l2diff_ldc_vals, marker='o', markerfacecolor='none')
plt.ylabel("$L^{2}$ diff between UCM and Newtonian velocity")
plt.xlabel("$\lambda_{1}$ value")
plt.savefig("results_compare_ldc/ucm_nse_ldc_l2.pdf")
plt.close()

p_ldc_h1 = plt.loglog(l1_ldc_vals, h1diff_ldc_vals, marker='o', markerfacecolor='none')
plt.ylabel("$H^{1}$ diff between UCM and Newtonian velocity")
plt.xlabel("$\lambda_{1}$ value")
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
plt.savefig("results_compare_jb/ucm_nse_jb_l2.pdf")
plt.close()

p_jb_h1 = plt.loglog(l1_jb_vals, h1diff_jb_vals, marker='o', markerfacecolor='none')
plt.ylabel("$H^{1}$ diff between UCM and Newtonian velocity")
plt.xlabel("$\lambda_{1}$ value")
plt.savefig("results_compare_jb/ucm_nse_jb_h1.pdf")
plt.close()

p_jb_l2p = plt.loglog(l1_jb_vals, l2pdiff_jb_vals, marker='o', markerfacecolor='none')
plt.ylabel("$L^{2}$ diff between UCM and Newtonian pressure")
plt.xlabel("$\lambda_{1}$ value")
plt.savefig("results_compare_jb/ucm_nse_jb_l2p.pdf")
plt.close()


