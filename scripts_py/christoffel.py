# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 20:26:35 2023

@author: Stefano Veroni, Vid Homsak
"""

from sympy import diff, symbols, sin, Matrix, init_printing, sqrt
import numpy as np
import pandas as pd

def Christoffel(g_ij, g_inv, i, j, k):
    return (1/2) \
            * sum(g_inv[i, m]   
            * (  diff(g_ij[k, m], cs[j])
               + diff(g_ij[j, m], cs[k])
               - diff(g_ij[j, k], cs[m])) 
               for m in range(g_ij.shape[0]))

def print4D(g_ij):
    g_ij.simplify()
    g_ij_str = [[str(g_ij[i,j])+"   " for j in range(4)] 
                    for i in range(4)]
    df = pd.DataFrame(g_ij_str)
    print('\n'.join(df.to_string(index = False).split('\n')[1:]))

def printChris(g_ij):
    mus = [0, 1, 2, 3]
    g_inv = g_ij.inv()
    Gamma = np.zeros((4,4,4)).tolist()
    for i in range(4):
        for j in range(4):
            for k in range(4):
                Gamma[i][j][k] = Christoffel(g_ij, g_inv, i, j, k)
                if Gamma[i][j][k]:
                    print(f"G^{mus[i]}_{mus[j]}{mus[k]} = ", Gamma[i][j][k].simplify())
    return Gamma

def trace(g_inv, A):
    S = 0
    for i in range(4):
        for j in range(4):
            S += g_inv[i,j]*A[i,j]
    return S
    

init_printing() 


#%%
print("\n################ Schwartschild Test ##########################\n")
t, r, theta, phi, M = symbols('t r o p M')
cs = [t, r, theta, phi]

g_ij = Matrix([[2*M/r-1,      0     , 0   , 0],
               [0      , 1/(1-2*M/r), 0   , 0],
               [0      ,      0     , r**2, 0],
               [0      ,      0     , 0   , r**2 * (sin(theta))**2]]
              )
g_inv = g_ij.inv()
print("g_munu")
print4D(g_ij)
print("g^munu")
print4D(g_inv)

print("\nThe Scwaharztcshild non-zero Christoffel Symbols are\n")
printChris(g_ij)


#%%
print("\n################ EF Original ##########################\n")
t, r, theta, phi, M = symbols('t r o p M')
cs = [t, r, theta, phi]

g_ij = Matrix([[2*M/r-1,    2*M/r   , 0   , 0],
               [2*M/r  ,   2*M/r+1  , 0   , 0],
               [0      ,      0     , r**2, 0],
               [0      ,      0     , 0   , r**2 * (sin(theta))**2]]
              )
g_inv = g_ij.inv()
print("g_mn")
print4D(g_ij)
print("g^mn")
print4D(g_inv)

print("\n\nThe EF-original non-zero Christoffel Symbols are\n")
printChris(g_ij)

print("\n\nTrace of K_ab\n")
Kab = Matrix([[-sqrt(2)/M  ,    -2*sqrt(2)/M   , 0    , 0],
               [-2*sqrt(2)/M,    -4*sqrt(2)/M   , 0    , 0],
               [     0      ,          0        , -64*M, 0],
               [     0      ,          0        , 0    , -64*M * (sin(theta))**2]]
              ) /32
tr = trace(g_inv.subs(r, 2*M), Kab)
print(tr)
tr = trace(g_inv, Kab)
print(tr.subs(r,2*M))




