import numpy as np
from numba import njit,jit
from time import time
import pandas as pd
import matplotlib.pyplot as plt

masses = np.array([1,1.5,2])
means= [7,16.875,29]
stds = [2.915,3.06,6.55]

means2 = [8.8,16,24]
stds2 = [2.92,3.35,3.79]

plt.figure()
plt.errorbar(masses**2,means,yerr=stds,capsize=4,ls="-",fmt="x",
    label="3+1D Schwarzshild: futlinks")
plt.errorbar(masses**2,means2,yerr=stds2,capsize=4,ls="-",fmt="x",
    label="3+1D Schwarzshild: cmatrix")
plt.legend()
plt.ylabel("Number of links")
plt.xlabel("Area [a.u]")
plt.grid()
plt.show()