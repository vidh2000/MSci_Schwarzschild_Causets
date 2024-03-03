import numpy as np
from numba import njit,jit
from time import time
from scipy.special import gamma
import matplotlib.pyplot as plt

def MM_drelation(d):
    a = gamma(d+1)
    b = gamma(d/2)
    c = 4*gamma(3*d/2)
    return a*b/c - 3.1415926535/8


x = np.arange(0.1,5,0.1)
plt.figure()
plt.plot(x,MM_drelation(x))
plt.show()
