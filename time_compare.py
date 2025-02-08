from fenics import *
from oldroyd_3_SRTD import *
from oldroyd_3_EVSS import *
from steady_nse_solver import *
import time
import csv
import numpy as np


h_array_ldc = np.array([2.5e-2, 1.25e-2, 6.25e-3])
h_array_jb = np.array([5e-2, 2.5e-2, 1.25e-2])

l1_array = np.array([0.04])

eta = 1.0
speed = 1.0
rad = 0.5
ecc = 0.25

tol = 1e-9
max_srtd_iters = 20


# Time comparisons for LDC UCM
problem = 'ldc'
a = 1.0
model = 'ucm'

for l1 in l1_array:
    l1 = float(l1)
    mu1 = a*l1

    # now the problem, model, and l1 are fixed. 
    table_file = open('results_time/' + problem + '_' + model + '_l1=%.3e.csv'%l1, 'w')
    writer = csv.writer(table_file)
    # write header
    writer.writerow(['Formulation', 'h=%.3e 1'%h_array_ldc[0], 'h=%.3e 2'%h_array_ldc[0], 'h=%.3e 3'%h_array_ldc[0], 'h=%.3e avg'%h_array_ldc[0],
                      'h=%.3e 1'%h_array_ldc[1], 'h=%.3e 2'%h_array_ldc[1], 'h=%.3e 3'%h_array_ldc[1], 'h=%.3e avg'%h_array_ldc[1],
                      'h=%.3e 1'%h_array_ldc[2], 'h=%.3e 2'%h_array_ldc[2], 'h=%.3e 3'%h_array_ldc[2], 'h=%.3e avg'%h_array_ldc[2]]) 
    table_file.flush()

    # SRTD tests
    # three different h/meshsize values we want to compare
    outrow = ['SRTD']
    for i in range(3):
        h = h_array_ldc[i]
        # For each h, want to clock 3 times, take the average
        times = np.zeros(3)
        for j in range(3):
            start_solve = time.time()
            solution = oldroyd_3_LDC_SRTD(h, speed, eta, l1, mu1, max_srtd_iters, tol)
            end_solve = time.time()
            times[j] = end_solve - start_solve
        # done with 3 experiments, take the average
        outrow.extend(['%.2f'%val for val in times])
        outrow.extend(['%.2f'%(np.average(times))])

    writer.writerow(outrow)
    table_file.flush()

    # EVSS tests
    # three different h/meshsize values we want to compare
    outrow = ['EVSS']
    for i in range(3):
        h = h_array_ldc[i]
        # For each h, want to clock 3 times, take the average
        times = np.zeros(3)
        for j in range(3):
            start_solve = time.time()
            solution = oldroyd_3_LDC_EVSS(h, speed, eta, l1, mu1)
            end_solve = time.time()
            times[j] = end_solve - start_solve
        # done with 3 experiments, take the average
        outrow.extend(['%.2f'%val for val in times])
        outrow.extend(['%.2f'%(np.average(times))])

    writer.writerow(outrow)
    table_file.flush()
    
    table_file.close()

# LDC Corot
problem = 'ldc'
a = 0.0
model = 'corot'

for l1 in l1_array:
    l1 = float(l1)
    mu1 = a*l1

    # now the problem, model, and l1 are fixed. 
    table_file = open('results_time/' + problem + '_' + model + '_l1=%.3e.csv'%l1, 'w')
    writer = csv.writer(table_file)
    # write header
    writer.writerow(['Formulation', 'h=%.3e 1'%h_array_ldc[0], 'h=%.3e 2'%h_array_ldc[0], 'h=%.3e 3'%h_array_ldc[0], 'h=%.3e avg'%h_array_ldc[0],
                      'h=%.3e 1'%h_array_ldc[1], 'h=%.3e 2'%h_array_ldc[1], 'h=%.3e 3'%h_array_ldc[1], 'h=%.3e avg'%h_array_ldc[1],
                      'h=%.3e 1'%h_array_ldc[2], 'h=%.3e 2'%h_array_ldc[2], 'h=%.3e 3'%h_array_ldc[2], 'h=%.3e avg'%h_array_ldc[2]]) 
    table_file.flush()

    # SRTD tests
    # three different h/meshsize values we want to compare
    outrow = ['SRTD']
    for i in range(3):
        h = h_array_ldc[i]
        # For each h, want to clock 3 times, take the average
        times = np.zeros(3)
        for j in range(3):
            start_solve = time.time()
            solution = oldroyd_3_LDC_SRTD(h, speed, eta, l1, mu1, max_srtd_iters, tol)
            end_solve = time.time()
            times[j] = end_solve - start_solve
        # done with 3 experiments, take the average
        outrow.extend(['%.2f'%val for val in times])
        outrow.extend(['%.2f'%(np.average(times))])

    writer.writerow(outrow)
    table_file.flush()

    # EVSS tests
    # three different h/meshsize values we want to compare
    outrow = ['EVSS']
    for i in range(3):
        h = h_array_ldc[i]
        # For each h, want to clock 3 times, take the average
        times = np.zeros(3)
        for j in range(3):
            start_solve = time.time()
            solution = oldroyd_3_LDC_EVSS(h, speed, eta, l1, mu1)
            end_solve = time.time()
            times[j] = end_solve - start_solve
        # done with 3 experiments, take the average
        outrow.extend(['%.2f'%val for val in times])
        outrow.extend(['%.2f'%(np.average(times))])

    writer.writerow(outrow)
    table_file.flush()
    
    table_file.close()





# Time comparisons for JB UCM
problem = 'jb'
a=1.0
model = 'ucm'

for l1 in l1_array:
    l1 = float(l1)
    mu1 = a*l1

    # now the problem, model, and l1 are fixed. 
    table_file = open('results_time/' + problem + '_' + model + '_l1=%.3e.csv'%l1, 'w')
    writer = csv.writer(table_file)
    # write header
    writer.writerow(['Formulation', 'h=%.3e 1'%h_array_jb[0], 'h=%.3e 2'%h_array_jb[0], 'h=%.3e 3'%h_array_jb[0], 'h=%.3e avg'%h_array_jb[0],
                      'h=%.3e 1'%h_array_jb[1], 'h=%.3e 2'%h_array_jb[1], 'h=%.3e 3'%h_array_jb[1], 'h=%.3e avg'%h_array_jb[1],
                      'h=%.3e 1'%h_array_jb[2], 'h=%.3e 2'%h_array_jb[2], 'h=%.3e 3'%h_array_jb[2], 'h=%.3e avg'%h_array_jb[2]]) 
    table_file.flush()

    # SRTD tests
    # three different h/meshsize values we want to compare
    outrow = ['SRTD']
    for i in range(3):
        h = h_array_jb[i]
        # For each h, want to clock 3 times, take the average
        times = np.zeros(3)
        for j in range(3):
            start_solve = time.time()
            solution = oldroyd_3_JB_SRTD(h, rad, ecc, speed, eta, l1, mu1, max_srtd_iters, tol)
            end_solve = time.time()
            times[j] = end_solve - start_solve
        # done with 3 experiments, take the average
        outrow.extend(['%.2f'%val for val in times])
        outrow.extend(['%.2f'%(np.average(times))])

    writer.writerow(outrow)
    table_file.flush()

    # EVSS tests
    # three different h/meshsize values we want to compare
    outrow = ['EVSS']
    for i in range(3):
        h = h_array_jb[i]
        # For each h, want to clock 3 times, take the average
        times = np.zeros(3)
        for j in range(3):
            start_solve = time.time()
            solution = oldroyd_3_JB_EVSS(h, rad, ecc, speed, eta, l1, mu1)
            end_solve = time.time()
            times[j] = end_solve - start_solve
        # done with 3 experiments, take the average
        outrow.extend(['%.2f'%val for val in times])
        outrow.extend(['%.2f'%(np.average(times))])

    writer.writerow(outrow)
    table_file.flush()
    
    table_file.close()



# Time comparisons for JB Corot

problem = 'jb'
a=0.0
model = 'corot'

for l1 in l1_array:
    l1 = float(l1)
    mu1 = a*l1

    # now the problem, model, and l1 are fixed. 
    table_file = open('results_time/' + problem + '_' + model + '_l1=%.3e.csv'%l1, 'w')
    writer = csv.writer(table_file)
    # write header
    writer.writerow(['Formulation', 'h=%.3e 1'%h_array_jb[0], 'h=%.3e 2'%h_array_jb[0], 'h=%.3e 3'%h_array_jb[0], 'h=%.3e avg'%h_array_jb[0],
                      'h=%.3e 1'%h_array_jb[1], 'h=%.3e 2'%h_array_jb[1], 'h=%.3e 3'%h_array_jb[1], 'h=%.3e avg'%h_array_jb[1],
                      'h=%.3e 1'%h_array_jb[2], 'h=%.3e 2'%h_array_jb[2], 'h=%.3e 3'%h_array_jb[2], 'h=%.3e avg'%h_array_jb[2]]) 
    table_file.flush()

    # SRTD tests
    # three different h/meshsize values we want to compare
    outrow = ['SRTD']
    for i in range(3):
        h = h_array_jb[i]
        # For each h, want to clock 3 times, take the average
        times = np.zeros(3)
        for j in range(3):
            start_solve = time.time()
            solution = oldroyd_3_JB_SRTD(h, rad, ecc, speed, eta, l1, mu1, max_srtd_iters, tol)
            end_solve = time.time()
            times[j] = end_solve - start_solve
        # done with 3 experiments, take the average
        outrow.extend(['%.2f'%val for val in times])
        outrow.extend(['%.2f'%(np.average(times))])

    writer.writerow(outrow)
    table_file.flush()

    # EVSS tests
    # three different h/meshsize values we want to compare
    outrow = ['EVSS']
    for i in range(3):
        h = h_array_jb[i]
        # For each h, want to clock 3 times, take the average
        times = np.zeros(3)
        for j in range(3):
            start_solve = time.time()
            solution = oldroyd_3_JB_EVSS(h, rad, ecc, speed, eta, l1, mu1)
            end_solve = time.time()
            times[j] = end_solve - start_solve
        # done with 3 experiments, take the average
        outrow.extend(['%.2f'%val for val in times])
        outrow.extend(['%.2f'%(np.average(times))])

    writer.writerow(outrow)
    table_file.flush()
    
    table_file.close()




