import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

n, it = np.loadtxt('n_it.txt', usecols = (0,1), unpack = True, skiprows = 1)

n_log = np.log10(n)
it_log = np.log10(it)

#print(linregress(n_log, it_log))

plt.rc('axes', titlesize=12)     # fontsize of the axes title
plt.rc('axes', labelsize=12)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels

# plt.plot(n, it)
plt.plot(n, it, '-ok')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.title('Iterations as function of grid size n = 20, 30, ..., 250.', wrap = True)
plt.xlabel('n')
plt.ylabel('Iterations')
plt.savefig('n_it_plot.pdf')
plt.show()


ax = plt.plot(n_log, it_log, '-ok', label = 'Slope = 2.37')
plt.title('Logarithmic plot of iterations as function of grid size n = 20, 30, ..., 250.', wrap = True)
plt.xlabel('n, log10')
plt.ylabel('Iterations, log10')
plt.ylim([2.2,5.5])
plt.xlim([1.2,2.5])
plt.legend()
plt.savefig('n_it_logplot.pdf')

plt.show()

plt.tight_layout()

#Linear regression output:
#LinregressResult(slope=2.3720692189510313, intercept=-0.47772746346677764, rvalue=0.9970661528583901,
#pvalue=4.710792339012157e-26, stderr=0.038824601858391504, intercept_stderr=0.080397025097479)
