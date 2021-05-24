import matplotlib.pyplot as plt
import numpy as np
from io import StringIO

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

# Import text file
with open('error_estimator_max_RB.txt') as file:
    rb_data = file.read()
rb_absolute_error_estimator = np.genfromtxt(StringIO(rb_data[1:-1]), delimiter=',')
rb_relative_error_estimator = rb_absolute_error_estimator/rb_absolute_error_estimator[0]

with open('error_estimator_max_SCM.txt') as file:
    scm_data = file.read()
scm_absolute_error_estimator = np.genfromtxt(StringIO(scm_data[1:-1]), delimiter=',')
scm_relative_error_estimator = scm_absolute_error_estimator/scm_absolute_error_estimator[0]

N = np.linspace(0, 50, num=51)
SCM = np.linspace(0, 12, num=13)

# Define default LaTeX fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# SCM convergence plot
fig1, ax1 = plt.subplots()
ax1.plot(SCM, scm_relative_error_estimator, 'o-')
ax1.set_ylabel(r'Relative SCM error estimator $\Delta_\mathrm{rel}^\mathrm{SCM}$', fontsize=14)
ax1.set_xlabel(r'Number of SCM iterations', fontsize=14)
ax1.grid()
fig1.savefig('scm_convergence_plot.pdf')
# fig1.show()
print(rb_relative_error_estimator)

# RB convergence plot
fig2, ax2 = plt.subplots()
ax2.plot(N, rb_relative_error_estimator, 'o-')
ax2.set_yscale("log")
ax2.set_ylabel(r'Relative RB error estimator $\Delta_\mathrm{rel}^\mathrm{RB}$', fontsize=14)
ax2.set_xlabel(r'Number of basis functions $N$', fontsize=14)
ax2.grid()
fig2.savefig('rb_convergence_plot.pdf')
# fig2.show()