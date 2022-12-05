import numpy as np
from numba import njit,jit
from time import time
import pandas as pd
import matplotlib.pyplot as plt


# 4D data with futlinks creation...
#masses = np.array([1,1.5,2])
#means= [7,16.875,29]
#stds = [2.915,3.06,6.55]
# 4D data with futmatrix...
#masses = np.array([1,1.5,2])
#means2 = [8.8,16,24]
#stds2 = [2.92,3.35,3.79]

# 4D data with futmatrix...
masses = np.array([1,1.5,2,2.5,3])
means= [7.8,15.6,23,41.2,51.2]
stds = [1.7,3.72,5.83,2.04,8.9]




plt.figure()

plt.errorbar(masses**2,means,yerr=stds,capsize=4,ls="-",fmt="o",
    label=r"$3+1D$ Schwarzshild spacetime")
plt.legend()
plt.ylabel(r"Number of links $N$")
plt.xlabel(r"Area$\propto M^2$ [a.u]")
plt.grid(alpha=0.3)
plt.show()