#!/usr/bin/env python
"""
Created on 13 Oct 2022

@author: Stefano Veroni
"""
#%%
from __future__ import annotations
from typing import List, Tuple 

from causets.causetevent import CausetEvent
from causets.causet import Causet
from causets.sprinkledcauset import SprinkledCauset
from causets.shapes import CoordinateShape
import causets.causetplotting as cplt

import numpy as np
import random
from tqdm import tqdm
import matplotlib.pyplot as plt


#%% CHECK THAT INTERVAL AND ORDRING FRACTION WORK CORRECTLY

N = 8
Cm = np.zeros([N,N])
for i in range(N-1):
    for j in range(i+1, N):
        Cm[i,j] = 1
print("CMatrix")
print(Cm)
print("\nDifference between retrieved and original Cmatrix")
C = Causet().FromCausalMatrix(Cm)
print(C.CMatrix(method = "label") - Cm)
print("\nChain Link Matrix\n", C.LMatrix())
Clist = C.nlist()
print("Has ordering Frac = ", C.ord_fr_A(C.Interval(Clist[3], Clist[6],
                                                disjoin = True)))
print("Should have 1\n")
print("Has ordering Frac = ", C.ord_fr_ab(Clist[3], Clist[6]))
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
print("\nChain Link Matrix\n", C.LMatrix())
Clist = C.nlist(method="label")
#print("Clist:\n",Clist)

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

#%% TEST MM DIM RELATION & OPTIMIZER
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
print("1.00 -> ",fsolve(MM_to_solve, 2, 1))
print("0.50 -> ",fsolve(MM_to_solve, 2, 0.5))
print("0.40 -> ",fsolve(MM_to_solve, 2, 0.40))
print("0.35 -> ",fsolve(MM_to_solve, 2, 0.35))
print("0.30 -> ",fsolve(MM_to_solve, 2, 0.3))
print("0.23 -> ",fsolve(MM_to_solve, 2, 8/35, full_output=1)[0])
print("0.15 -> ",fsolve(MM_to_solve, 2, 0.15))
print("0.12 -> ",fsolve(MM_to_solve, 2, 0.12))
print("0.10 -> ",fsolve(MM_to_solve, 2, 0.1))
print("0.08 -> ",fsolve(MM_to_solve, 2, 0.08))
print("0.05 -> ",fsolve(MM_to_solve, 2, 0.05))

#%% CHECK DIMESNION ESTIMATOR IN FLAT SPACETIME FOR ALL COORDINATES
from causets.spacetimes import *
st   = [    FlatSpacetime   , deSitterSpacetime, 
       AntideSitterSpacetime, BlackHoleSpacetime]
dims = [  [1,2,3,4],                [2,3,4],            
           [2,3,4],                    [2]      ]
r = 2
dur = 1
ballh_ps = {'name': 'ball',     'radius': r, 'hollow':0.99} 
ball_ps  = {'name': 'ball',     'radius': r, 'hollow':0}
cylh_ps  = {'name': 'cylinder', 'radius': r, 'hollow':0.99, 'duration':dur} 
cyl_ps   = {'name': 'cylinder', 'radius': r, 'hollow':0,    'duration':dur}
cub_ps   = {'name': 'cube'    , 'edge'  : r} 
shapes = [['ball_hollow'    , ballh_ps],
          ['ball'           , ball_ps ],
          ['cylinder_hollow',cylh_ps  ],
          ['cylinder'       , cyl_ps  ],
          ['cube'           , cub_ps  ] ]

Ns = [2048, 1024, 512, 384, 256, 192,  128, 64, 32, 16]# cardinalities to test
repetitions = 10                                       #repetitions to average
cuts = np.array([0]+Ns[:-1])-np.array([0]+Ns[1:])      #for coarse graining
for sps in shapes:
    dim_est = []
    dim_std = []
    x = 1 #skip d = 1, ..., x
    for i in range(len(dims[0])-x):
        d = dims[0][i+x] 
        dim_est.append([])
        dim_std.append([])
        for rep in range(repetitions):
            dim_est[i].append([])
            dim_std[i].append([])
            
            S: CoordinateShape = CoordinateShape(d, **sps[1])
            try:
                C: SprinkledCauset = SprinkledCauset(card=Ns[0],
                                                spacetime=FlatSpacetime(d), 
                                                shape=S)
            except ZeroDivisionError:
                print(f"At dimension {d} did NOT use {sps[0]}")
                C: SprinkledCauset = SprinkledCauset(card=Ns[0],
                                                spacetime=FlatSpacetime(d))
            
            for cut in tqdm(cuts, f"{sps[0]} (dim {d})"):
                if cut != 0:
                    C.coarsegrain(card = cut) #cgrain Ns[i]->Ns[i+1]
                MMd = C.MMdim_est(Nsamples = 50, 
                                    ptime_constr=lambda t:t<2.5*r,
                                    size_min = min(50, int(len(C)/4)),
                                    full_output = True)
                dim_est[i][rep].append(MMd[0]) # add to rth repetition 
        
        #Average over repetitions:
        try:
            dim_std[i] = np.nanstd (dim_est[i], axis = 0)
            dim_est[i] = np.nanmean(dim_est[i], axis = 0)
        except (TypeError, ZeroDivisionError):
            try:
                dim_std[i]= np.nanstd (np.array(dim_est[i]).astype(np.float64))
                dim_est[i]= np.nanmean(np.array(dim_est[i]).astype(np.float64))
            except (TypeError, ZeroDivisionError):
                print("SECOND EXCEPT USED")
                beforeerror = dim_est[i]
                dim_std[i] = np.nanstd (np.array(dim_est[i], dtype=np.float64),
                                        axis = 0)
                dim_est[i] = np.nanmean(np.array(dim_est[i], dtype=np.float64),
                                        axis = 0)
    try:
        fig = plt.figure(f"MMFlatDim {sps[0]}")
        #Ns.reverse()
        plt.title(f"Myrheim-Mayers in {sps[0]} Minkowski")
        plt.xlabel("Cardinality")
        plt.ylabel("Dimension")
        for i in range(len(dims[0])-x):
            ests = dim_est[i]#np.flip(dim_est[i])
            stds = dim_std[i]#np.flip(dim_std[i])
            lbl  = f"Dimension {dims[0][i+x]}"
            plt.errorbar(Ns,ests, yerr=stds, fmt=".", capsize = 4, label = lbl)
        plt.legend()
        plt.xscale('log')
        fig.show()
    except TypeError:
        continue

del d, sps, i, rep, cut
del S, Ns, cuts, repetitions
del st, dims, shapes, r, dur, ballh_ps, ball_ps, cylh_ps, cyl_ps, cub_ps
del ests, stds, lbl, fig


# %%
print(np.nanstd(beforeerror, dtype = np.float64, axis = 0))