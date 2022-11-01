import numpy as np
from numba import njit,jit
from time import time


def f(N):
    matrix = []
    for i in range(N):
        matrix.append([])
        for j in range(N):
            matrix[i].append(0)
    return matrix

@njit
def f2(matrix,N):
    for i in range(N):
        for j in range(N):
            matrix[i][j] = 1
    return matrix

N=int(1e4)
a = f(N)

s = time()
a = f2(a,N)
print(sum(a[0]))
print("Time:", time()-s)