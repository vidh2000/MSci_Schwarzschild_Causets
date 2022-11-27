import numpy as np
from numba import njit,jit
from time import time
import pandas as pd
import matplotlib.pyplot as plt

masses = [1,1.5,2,2.5]
means= [5.45,7.25,9.85,13.15]
stds = [2.60,3.78,3.42,3.15]


plt.figure()
plt.errorbar(masses,means,yerr=stds,capsize=4,ls="-",fmt="x",
    label="2+1D Schwarzshild")
plt.legend()
plt.ylabel("Number of links")
plt.xlabel("Area [a.u]")
plt.grid()
plt.show()