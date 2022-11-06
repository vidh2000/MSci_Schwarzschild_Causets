#!/usr/bin/env python
'''
Created on 9 Oct 2020

@author: Stefano Veroni, Vid Homsak
'''

# Small change
# Change
#%%
from __future__ import annotations
from typing import List, Tuple  
from causets.sprinkledcauset import SprinkledCauset  
from causets.embeddedcauset import EmbeddedCauset
from causets.spacetimes import FlatSpacetime, deSitterSpacetime  
from causets.shapes import CoordinateShape  
from causets.causetevent import CausetEvent  
import causets.causetplotting as cplt  

import matplotlib.pyplot as plt
import numpy as np
#%%
d = 2
edges = [4] + [1]*(d-2) + [4]
card = 50

# DO PERIODIC
S = CoordinateShape(d, 'cuboid', edges = edges)
period = tuple([0]*(d-2) + [S.Parameter("edges")[d-1]])
M = FlatSpacetime(d, period = period)
C = SprinkledCauset(card = card, dim = d, spacetime = M, shape = S) 

#PLOT PERIODIC
cplt.setDefaultColors('UniYork')  # using University of York brand colours
cplt.figure() 
periodic_dimensions = []
for i in range(d-1):
    if period[i]!=0:
        periodic_dimensions.append(i)
dims = list(range(1, d)) + [0] #choose order of plot dimensions
cplt.plot(C, dims=dims, events={'alpha': 1},
            links={'alpha': 0.5, 'linewidth': 0.5}, labels=False)
ax: cplt.Axes = cplt.gca()
ax.set_xlabel('x' if dims[0] == 1 else 't')
if d == 3:
    ax.set_ylabel('y' if dims[1] == 2 else 't')
    if len(dims) > 2:
        ax.set_zlabel('space' if dims[2] > 0 else 't')
elif d== 2:
    ax.set_ylabel('t')
ax.set_title(f"Periodic in Dimension(s) {periodic_dimensions}")
ax.grid(False)
cplt.show()

#PLOT NON-PERIODIC EQUIVALENT
M2 = FlatSpacetime(d)
C2 = EmbeddedCauset(M2, S, C.get_coords())
cplt.figure() 
dims = list(range(1, d)) + [0] #choose order of plot dimensions
cplt.plot(C2, dims=dims, events={'alpha': 1},
            links={'alpha': 0.5, 'linewidth': 0.5}, labels=False)
ax: cplt.Axes = cplt.gca()
ax.set_xlabel('x' if dims[0] == 1 else 't')
if d == 3:
    ax.set_ylabel('y' if dims[1] == 2 else 't')
    if len(dims) > 2:
        ax.set_zlabel('space' if dims[2] > 0 else 't')
elif d== 2:
    ax.set_ylabel('t')
ax.set_title(f"Non Periodic")
ax.grid(False)
cplt.show()
# %%
