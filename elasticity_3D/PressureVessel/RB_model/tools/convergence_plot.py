import matplotlib.pyplot as plt
import numpy as np
from io import StringIO

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


offline_folder_name = "PressureVessel"
folderpath = "../" + offline_folder_name

# Define default LaTeX fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Get the values of maximum error estimators
with open(folderpath + "/post_processing/error_estimator_max.txt") as file:
    rb_data = file.read()
rb_absolute_error_estimator = np.genfromtxt(StringIO(rb_data[1:-1]), delimiter=',')
rb_relative_error_estimator = rb_absolute_error_estimator/rb_absolute_error_estimator[0]

# Get the number of basis functions
with open(folderpath + "/basis/basis.length") as file:
    N_max = int(file.read())
N = np.linspace(0, N_max, num=N_max+1)

# RB convergence plot
fig1, ax1 = plt.subplots()
ax1.plot(N, rb_relative_error_estimator, 'o-')
ax1.set_yscale("log")
ax1.set_ylabel(r'Relative RB error estimator $\Delta_\mathrm{rel}^\mathrm{RB}$', fontsize=14)
ax1.set_xlabel(r'Number of basis functions $N$', fontsize=14)
ax1.grid()
fig1.savefig('rb_convergence_plot.pdf')
# fig1.show()
