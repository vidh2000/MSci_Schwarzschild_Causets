#!/usr/bin/env python

'''
Created on 4 Jan 2022

@author: Stefano Veroni, Vid Homsak
'''
import numpy as np
import matplotlib.pyplot as plt
from . import causet_helpers as ch

def dr_dtstar (r, M, eta, pm = -1):
    """
    dr/dtstar (EF-original coords) for a black hole with mass M at position r.
    Eta is E/L and pm_r is the sign of r2-r1 (- for going inside, + outside).
    """
    a = 2*M/r-1
    b = np.sqrt(a + eta**2 * r**2)
    d = b/(- pm * eta * r - 2*M/r*b)
    if np.sign(pm) == np.sign(a*d):
        return a * d
    else:
        return 0

def dphi_dtstar (r, M, eta, pm_r = +1):
    """
    dphi/dtstar (EF-original coords) for a black hole with mass M at position r.
    Eta is E/L and pm_r is the sign of r2-r1 (- for going inside, + outside).
    """
    a =  1 - 2*M/r
    b = np.sqrt(a + eta**2 * r**2)
    d = eta * r**2 + pm_r * 2*M*b
    return a / d

def BH_lightcone (x0, M):
    t0, r0, phi0 = x0
    if r0 > 2*M:
        do stuff
    elif r0 <= 2*M:
       pm = -1
        
    return points




    
