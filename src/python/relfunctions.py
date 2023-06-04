#!/usr/bin/env python

'''
Created on 16 Feb 2022

@author: Stefano Veroni, Vid Homsak
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sympy import diff, symbols, sin, Matrix, init_printing, sqrt, simplify
from . import causet_helpers as ch


#############################################################################
##### SIMPY BASED FUNCTIONS FOR RELATIVISTIC CALCULATIONS
#############################################################################
def print4D(g_ij):
    g_ij.simplify()
    g_ij_str = [[str(g_ij[i,j])+"   " for j in range(4)] 
                    for i in range(4)]
    df = pd.DataFrame(g_ij_str)
    init_printing() 
    print('\n'.join(df.to_string(index = False).split('\n')[1:]))


def trace(g_ij, A_ij):
    g_inv = g_ij.inv()
    S = 0
    for i in range(4):
        for j in range(4):
            S += g_inv[i,j]*A_ij[i,j]
    return S


def Christoffel_ijk(g_ij, i, j, k, indexes = symbols('t r o p')):
    cs = indexes
    g_inv = g_ij.inv()
    return (1/2) \
            * sum(g_inv[i, m]   
            * (  diff(g_ij[k, m], cs[j])
               + diff(g_ij[j, m], cs[k])
               - diff(g_ij[j, k], cs[m])) 
               for m in range(g_ij.shape[0]))

def printChris(g_ij, indexes = symbols('t r o p')):
    mus = [0, 1, 2, 3]
    g_inv = g_ij.inv()
    cs = indexes
    Gamma = np.zeros((4,4,4)).tolist()
    init_printing() 
    for i in range(4):
        for j in range(4):
            for k in range(4):
                Gamma[i][j][k] = Christoffel_ijk(g_ij, i, j, k, cs)
                if Gamma[i][j][k]:
                    print(f"G^{mus[i]}_{mus[j]}{mus[k]} = ",
                           Gamma[i][j][k].simplify())
    return Gamma

def covDerivative_covector(cov, g_ij, indexes = symbols('t r o p')):
    n_ab = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for a in range(4):
        for b in range(4):
            n_ab[a,b] += diff(cov[a], indexes[b])
            #print(simplify(diff(cov[a], indexes[b])))
            Chris_mu_a_b_Times_cov_mu = 0
            for mu in range(4):
                Chris_mu_a_b = Christoffel_ijk(g_ij, mu, b, a, indexes)
                Chris_mu_a_b_Times_cov_mu += Chris_mu_a_b * cov[mu]
            n_ab[a,b] -= Chris_mu_a_b_Times_cov_mu
    return n_ab
