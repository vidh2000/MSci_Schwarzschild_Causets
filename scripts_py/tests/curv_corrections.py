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
print("\n################ EF Original ##########################\n")
t, r, theta, phi, M, rS = symbols('t r o p M rS')
cs = [t, r, theta, phi]

print("\n STEP 1: THE METRIC & SURFACES ################\n")
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

# Define the horizon and hypersurface
rS = 2*M
H = r - rS
Sigma = t
print(f"Horizon : {H} = 0")
print(f"Sigma   : {Sigma} = 0")


print("\n STEP 2: THE VECTORS  ################\n")
# Define the covector n_mu
A = sp.symbols('A', real=True)
n_cov = A*sp.Matrix([sp.diff(Sigma, coord) for coord in (t, r, o, p)])
# Raise the index to find the vector n^mu
g_inv = g.inv()
n_con = g_inv*n_cov
# Impose the condition n^mu n_mu = -1 at r = rS and solve for A
n_s = n_con.subs(r, rS)
eq = sp.Eq(n_s.dot(n_s), -1)
A_val = sp.solve(eq, A)[0]
# Substitute the value of A in n_cov and n_con
n_cov = n_cov.subs(A, A_val)
n_con = n_con.subs(A, A_val)

# Print n_cov and n_con
print("Covector n_mu =")
sp.pprint(n_cov)
print("\nVector n^mu =")
sp.pprint(n_con)


# Define the covector k_mu
lamb = sp.symbols('lamb', real=True)
k_cov = lamb*sp.Matrix([sp.diff(H, coord) for coord in (t, r, o, p)])

# Raise the index to find the vector k^mu
k_con = g_inv*k_cov

# Impose the condition k^mu n_mu = -1/sqrt(2) at H = 0 and Sigma = 0 and solve for lamb
k = k_con.dot(n_cov)
eq1 = sp.Eq(k.subs({H:0, Sigma:0}), -1/sp.sqrt(2))
lamb_val = sp.solve(eq1, lamb)[0]

# Substitute the value of lamb in k_cov and k_con
k_cov = k_cov.subs(lamb, lamb_val)
k_con = k_con.subs(lamb, lamb_val)


# Print n_cov, n_con, k_cov, k_con in general
print("Covector n_mu =")
sp.pprint(n_cov)
print("\nVector n^mu =")
sp.pprint(n_con)
print("\nCovector k_mu =")
sp.pprint(k_cov)
print("\nVector k^mu =")
sp.pprint(k_con)


# Define m^mu
m_con = sp.sqrt(2)*k_con - n_con
# Lower the index to find m_mu
m_cov = g*m_con
# Print m^mu and m_mu in general
print("\nVector m^mu =")
sp.pprint(m_con)
print("\nCovector m_mu =")
sp.pprint(m_cov)

# Define ell^mu
ell_con = sp.symbols('ell^t ell^r ell^o ell^p', real=True)
ell_con[0] = -k_con[0]/k_con[1]
ell_con[1] = 1/k_con[1]
ell_con[2:] = [0, 0]

# Raise the index to find ell^mu
ell_cov = g*ell_con
# Print ell^mu and ell_mu in general
print("\nVector ell^mu =")
sp.pprint(ell_con)
print("\nCovector ell_mu =")
sp.pprint(ell_cov)

# Evaluate the vectors and covectors at H = 0 and Sigma = 0
n_cov_0 = n_cov.subs({H:0, Sigma:0})
n_con_0 = n_con.subs({H:0, Sigma:0})
k_cov_0 = k_cov.subs({H:0, Sigma:0})
k_con_0 = k_con.subs({H:0, Sigma:0})
m_con_0 = m_con.subs({H:0, Sigma:0})
m_cov_0 = m_cov.subs({H:0, Sigma:0})
ell_con_0 = ell_con.subs({H:0, Sigma:0})
ell_cov_0 = ell_cov.subs({H:0, Sigma:0})

# Print n_cov, n_con, k_cov, k_con evaluated at H = 0 and Sigma = 0
print("\nCovector n_mu at H = 0, Sigma = 0:")
sp.pprint(n_cov_0)
print("\nVector n^mu at H = 0, Sigma = 0:")
sp.pprint(n_con_0)
print("\nCovector k_mu at H = 0, Sigma = 0:")
sp.pprint(k_cov_0)
print("\nVector k^mu at H = 0, Sigma = 0:")
sp.pprint(k_con_0)
# Print m^mu and m_mu at H=0 and Sigma=0
print("\nVector m^mu at H=0 and Sigma=0:")
sp.pprint(m_con_0)
print("\nCovector m_mu at H=0 and Sigma=0:")
sp.pprint(m_cov_0)
# Print ell^mu and ell_mu at H=0 and Sigma=0
print("\nVector ell^mu at H=0 and Sigma=0:")
sp.pprint(ell_con_0)
print("\nCovector ell_mu at H=0 and Sigma=0:")
sp.pprint(ell_cov_0)














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




