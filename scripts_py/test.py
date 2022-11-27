import numpy as np
from numba import njit,jit
from time import time
import pandas as pd
import matplotlib.pyplot as plt

x = [1,2,3,4]
means= [6.2,13,16.6,21.2]
stds = [1.47,2.36,3.38,3.6]


plt.figure()
plt.errorbar(x,means,yerr=stds,capsize=4,ls="-",fmt="x",label="2+1D Schwarzshild")
plt.legend()
plt.ylabel("Number of links")
plt.xlabel("Area [a.u]")
plt.grid()
plt.show()