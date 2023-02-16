# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 20:26:35 2023

@author: Stefano Veroni, Vid Homsak
"""

from causets_py import relfunctions as rel
from sympy import diff, symbols, sin, Matrix, init_printing, sqrt
import numpy as np
import pandas as pd


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
rel.print4D(g_ij)
print("g^munu")
rel.print4D(g_inv)

print("\nThe Scwaharztcshild non-zero Christoffel Symbols are\n")
rel.printChris(g_ij, cs)


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
rel.print4D(g_ij)
print("g^mn")
rel.print4D(g_inv)

print("\n\nThe EF-original non-zero Christoffel Symbols are\n")
rel.printChris(g_ij, cs)

print("\n\nTrace of K_ab\n")
Kab = Matrix([[-sqrt(2)/M  ,    -2*sqrt(2)/M   , 0    , 0],
               [-2*sqrt(2)/M,    -4*sqrt(2)/M   , 0    , 0],
               [     0      ,          0        , -64*M, 0],
               [     0      ,          0        , 0    , -64*M * (sin(theta))**2]]
              ) /32
tr = rel.trace(g_ij.subs(r, 2*M), Kab)
print(tr)
tr = rel.trace(g_ij, Kab)
print(tr.subs(r,2*M))




