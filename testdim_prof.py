#!/usr/bin/env python
"""
Created on 13 Oct 2022
@authors: Stefano Veroni, Vid Homsak
"""
#%%
from __future__ import annotations
from asyncio.windows_events import CONNECT_PIPE_INIT_DELAY
from typing import List, Tuple 

from causets.causetevent import CausetEvent
from causets.causet import Causet
from causets.sprinkledcauset import SprinkledCauset
from causets.shapes import CoordinateShape
import causets.causetplotting as cplt

from functions import * #profiler
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
warnings.filterwarnings("ignore", category=RuntimeWarning)

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import random


#%% 1. CHECK THAT INTERVAL AND ORDERING FRACTION WORK CORRECTLY

print("\n=========================================================")
print("CHECK THAT INTERVAL AND ORDERING FRACTION WORK CORRECTLY")
print("=========================================================\n")
N = 8
Cm = np.zeros([N,N])
for i in range(N-1):
    for j in range(i+1, N):
        Cm[i,j] = 1
print("TEST1: Using CMatrix (note, causet is a chain)")
print(Cm)
print("\nDifference between retrieved and original Cmatrix, hopefully 0")
C = Causet().FromCausalMatrix(Cm)
print(C.CMatrix(method = "label") - Cm)
print("\nWith Chain Link Matrix (note, causet is a chain)\n", C.LMatrix())
Clist = C.nlist()
print("A1 has ordering Frac = ", C.ord_fr_A(C.Interval(Clist[3], Clist[6],
                                                disjoin = True)))
print("Should have 1\n")
print("A2 has ordering Frac = ", C.ord_fr_ab(Clist[3], Clist[6]))
print("Should have 1\n")
del N, i, j


Cm = np.array([[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
               [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
               [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
               [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
               [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
               [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
C = Causet().FromCausalMatrix(Cm)
print("\nTEST2: Chain Link Matrix\n", C.LMatrix())
Clist = C.nlist(method="label")

A = C.Interval(Clist[3], Clist[8], disjoin = True)
print(f"Cardinality is {len(A)}")
print("Has ordering Frac = ", C.ord_fr_A(A))
print("Has ordering Frac = ", C.ord_fr_ab(Clist[3], Clist[8]))
print("Should have 0.9\n")
A = C.Interval(Clist[2], Clist[8], disjoin = True)
print(f"Cardinality is {len(A)}")
print("Has ordering Frac = ", C.ord_fr_A(A))
print("Has ordering Frac = ", C.ord_fr_ab(Clist[2], Clist[8]))
print("Should have 0.93\n")
A = C.Interval(Clist[3], Clist[9], disjoin = True)
print(f"Cardinality is {len(A)}")
print("Has ordering Frac = ", C.ord_fr_A(A))
print("Has ordering Frac = ", C.ord_fr_ab(Clist[3], Clist[9]))
print("Should have 1\n")
""" 
A = C.Interval(Clist[2], Clist[8], disjoin = True)
print(Clist[2], " preceeds ", Clist[3], " ? ",Clist[2]<=Clist[3])
print(Clist[3], " preceeds ", Clist[4], " ? ",Clist[3]<=Clist[4])
print(Clist[4], " preceeds ", Clist[5], " ? ",Clist[4]<=Clist[5])
print(Clist[5], " preceeds ", Clist[8], " ? ",Clist[5]<=Clist[8])
print(Clist[2], " preceeds ", Clist[8], " ? ",Clist[2]<=Clist[8])
print(f"Cardinality is {len(A)}")
try:
    print("Has ordering Frac = ", C.ord_fr(A))
except ZeroDivisionError:
    print("There was a division by zero in ordering fraction")
"""
del Cm, C, Clist

#%% 2. TEST MM DIM RELATION & FSOLVE OPTIMIZER
print("\n=========================================================")
print("CHECK MMDIM RELATION IS CORRECLTY SOLVED BY SP...FSOLVE")
print("=========================================================\n")
from scipy.optimize import fsolve
from scipy.special import gamma as spgamma

def MM_drelation(d):
            a = spgamma(d+1)
            b = spgamma(d/2)
            c = 4 * spgamma(3*d/2)
            return a*b/c
        
def MM_to_solve(d, ord_fr):
    return MM_drelation(d) - ord_fr/2

print("Ordering Fraction -> d")
print("     1.00 -> ",fsolve(MM_to_solve, 2, 1)[0])
print("     0.50 -> ",fsolve(MM_to_solve, 2, 0.5)[0])
print("     0.40 -> ",fsolve(MM_to_solve, 2, 0.40)[0])
print("     0.35 -> ",fsolve(MM_to_solve, 2, 0.35)[0])
print("     0.30 -> ",fsolve(MM_to_solve, 2, 0.3)[0])
print("     0.23 -> ",fsolve(MM_to_solve, 2, 8/35, full_output=1)[0][0])
print("     0.15 -> ",fsolve(MM_to_solve, 2, 0.15)[0])
print("     0.12 -> ",fsolve(MM_to_solve, 2, 0.12)[0])
print("     0.10 -> ",fsolve(MM_to_solve, 2, 0.1)[0])
print("     0.08 -> ",fsolve(MM_to_solve, 2, 0.08)[0])
print("     0.05 -> ",fsolve(MM_to_solve, 2, 0.05)[0])


#%% 3. CHECK DIMENSION ESTIMATOR IN FLAT SPACETIME FOR ALL COORDINATES
print("\n=========================================================")
print("CHECK MMrdim_est IN FLAT SPACETIME")
print("=========================================================\n")
# For profiling the code
PROFILE = False

from causets.spacetimes import *
st   = [    FlatSpacetime   , deSitterSpacetime, 
       AntideSitterSpacetime, BlackHoleSpacetime]
dims = [  [1,2,3,4],                [2,3,4],            
           [2,3,4],                    [2]      ]
r = 10
dur = 1
##%%% Shapes to be used
ballh_ps = {'name': 'ball',     'radius': r, 'hollow':0.99} 
ball_ps  = {'name': 'ball',     'radius': r, 'hollow':0}
cylh_ps  = {'name': 'cylinder', 'radius': r, 'hollow':0.99, 'duration':dur} 
cyl_ps   = {'name': 'cylinder', 'radius': r, 'hollow':0   , 'duration':dur}
cub_ps   = {'name': 'cube'    , 'edge'  : r} 
diamh_ps = {'name': 'diamond',  'radius': r, 'hollow':0.99} 
diam_ps  = {'name': 'diamond',  'radius': r, 'hollow':0}
##%%% Code running
shapes = [
        # ['ball_hollow'    , ballh_ps ],
        # ['ball'           , ball_ps  ],
        # ['cylinder_hollow',cylh_ps   ],
        # ['cylinder'       , cyl_ps   ],
        # ['cube'           , cub_ps   ],
        # ['diamond_hollow' , diamh_ps ],
          ['diamond'        , diam_ps  ]
         ]

# Define
# cardinalities to test, repetitions to average, cuts for coarse graining
Ns = [1024, 512, 256, 128, 64, 32, 16]
repetitions = 10
cuts = np.array([0]+Ns[:-1])-np.array([0]+Ns[1:])

print(f"\nEstimating MMd, N0={Ns[0]}:")
print(f"-{repetitions} repeats of {len(cuts)} times coarse-grained causets")

if __name__ == '__main__':
 
    for sps in shapes:
        # To collect estimates and stds from methods "random" and "big"
        rdim_est = []; bdim_est = [] 
        rdim_std = []; bdim_std = []
        nskip = 1 #skip d = 1, ..., nskip
        for i in range(len(dims[0])-nskip):
            d = dims[0][i+nskip] 
            rdim_est.append([]); bdim_est.append([])
            rdim_std.append([]); bdim_std.append([])
            for rep in tqdm(range(repetitions),f"{sps[0]} D={d}"):
                rdim_est[i].append([]); bdim_est[i].append([])
                rdim_std[i].append([]); bdim_std[i].append([])
                
                kwargs = {"card"      : Ns[0],
                          "spacetime" : FlatSpacetime(d),
                          "shape"     : CoordinateShape(d, **sps[1])}
                try:
                    if PROFILE:
                        C = profiler(SprinkledCauset, **kwargs)
                    else:
                        C = SprinkledCauset(**kwargs)

                except ZeroDivisionError:
                    print(f"At dimension {d} did NOT use {sps[0]}")
                    C = SprinkledCauset(card=Ns[0],
                                        spacetime=FlatSpacetime(d))
                
                for cut in cuts:
                    if cut != 0:
                        C.coarsegrain(card = cut) #cgrain Ns[i] -> Ns[i+1]
                    MMd = C.MMdim_est(Nsamples = 20, 
                                        method = "random",
                                        #ptime_constr=lambda t:t<2.5*r,
                                        size_min = min(1000, int(len(C)/2)),
                                        #size_max = 50,
                                        full_output = True)
                    rdim_est[i][rep].append(MMd[0]) # add to rth repetition
                    MMd = C.MMdim_est(Nsamples = 20, 
                                        method = "big",
                                        #ptime_constr=lambda t:t<2.5*r,
                                        size_min = min(1000, int(len(C)/2)),
                                        #size_max = 50,
                                        full_output = True)
                    bdim_est[i][rep].append(MMd[0]) # add to rth repetition 
            
            #Average over repetitions
            #print(f"rdim_est:\n{rdim_est}")
            rdim_std[i] = np.nanstd (rdim_est[i], axis = 0,dtype=np.float64)
            rdim_est[i] = np.nanmean(rdim_est[i], axis = 0,dtype=np.float64)
            bdim_std[i] = np.nanstd (bdim_est[i], axis = 0,dtype=np.float64)
            bdim_est[i] = np.nanmean(bdim_est[i], axis = 0,dtype=np.float64)

        #PLOT
          
        # Plot the random case
        fig = plt.figure(f"MMFlatDim ('random') {sps[0]}")
        plt.title(f"Myrheim-Mayers ('random') in {sps[0]} Minkowski")
        plt.xlabel("Cardinality")
        plt.ylabel("Dimension")
        for i in range(len(dims[0])-nskip):
            ests = rdim_est[i]
            stds = rdim_std[i]
            finmask = np.isfinite(ests)
            lbl  = f"Dimension {dims[0][i+nskip]}"
            plt.errorbar(np.array(Ns)[finmask], ests[finmask], 
                        yerr=stds[finmask], fmt=".", capsize = 4, label = lbl)
        plt.legend()
        plt.hlines([2, 3, 4], 0, Ns[0]*1.1, ls = "dashed", color = "r")
        plt.xscale('log')
        #plt.ylim((0.99, 4.5))
        plt.show()

        # Plot the big case
        fig = plt.figure(f"MMFlatDim ('big') {sps[0]}")
        plt.title(f"Myrheim-Mayers ('big') in {sps[0]} Minkowski")
        plt.xlabel("Cardinality")
        plt.ylabel("Dimension")
        for i in range(len(dims[0])-nskip):
            ests = bdim_est[i]
            stds = bdim_std[i]
            finmask = np.isfinite(ests)
            lbl  = f"Dimension {dims[0][i+nskip]}"
            plt.errorbar(np.array(Ns)[finmask], ests[finmask], 
                        yerr=stds[finmask], fmt=".", capsize = 4, label = lbl)
        plt.legend()
        plt.hlines([2, 3, 4], 0, Ns[0]*1.1, ls = "dashed", color = "r")
        plt.xscale('log')
        #plt.ylim((0.99, 4.5))
        plt.show()

del d, sps, i, rep, cut, nskip
del Ns, cuts, repetitions
del st, dims, shapes, r, dur, ballh_ps, ball_ps, cylh_ps, cyl_ps, cub_ps
del ests, stds, lbl, fig
