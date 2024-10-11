import csv
import matplotlib.pyplot as plt


with open('results_num_iters_corot/jb_residuals.csv', newline='') as jb_csv:
    jb_data = list(csv.reader(jb_csv))

l1_vals = [float(row[0]) for row in jb_data[16:21]] # just l1 values/first entry in each row
#print(l1_vals)

jb_plot_data = [[float(i) for i in row[1:]] for row in jb_data[16:21]] # the ones we're interested in are in rows 17-21
#print(jb_plot_data[0])

for i in range(4, -1, -1):
    plt.semilogy(range(1, len(jb_plot_data[i])+1), jb_plot_data[i], label='$\lambda_{1}$=%.3f'%l1_vals[i])

plt.xticks(range(1, 26))
plt.xlabel("SRTD Iteration")
plt.ylabel("Residual")
#plt.legend(loc=4)
plt.legend()
plt.savefig("results_num_iters_corot/jb_corot_num_iters_comparison.pdf")
plt.show()




with open('results_num_iters_corot/ldc_residuals.csv', newline='') as ldc_csv:
    ldc_data = list(csv.reader(ldc_csv))

l1_vals = [float(row[0]) for row in ldc_data[16:21]] # just l1 values/first entry in each row
#print(l1_vals)

ldc_plot_data = [[float(i) for i in row[1:]] for row in ldc_data[16:21]] # the ones we're interested in are in rows 17-21
#print(jb_plot_data[0])

for i in range(4, -1, -1):
    plt.semilogy(range(1, len(ldc_plot_data[i])+1), ldc_plot_data[i], label='$\lambda_{1}$=%.3f'%l1_vals[i])

plt.xticks(range(1, 26))
plt.xlabel("SRTD Iteration")
plt.ylabel("Residual")
#plt.legend(loc=4)
plt.legend()
plt.savefig("results_num_iters_corot/ldc_corot_num_iters_comparison.pdf")
plt.show()