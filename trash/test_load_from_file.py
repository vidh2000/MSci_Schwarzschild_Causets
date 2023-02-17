
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
import os 

########################################################
# Load file
dim = 2
card = 25
use_redge_in_name_file = True
edge = 4.00
mass = 1

isBH = True
isEFv = False
isS = False



C: EmbeddedCauset = EmbeddedCauset()
#S: CoordinateShape = CoordinateShape(3,"cylinder",radius=2,duration=3)

path = os.getcwd() # folder path
print(path)
# os.chdir("../")
# path = os.getcwd()
# print(path)
#file_name = os.path.join(path, 'data/flatspace_bicone_causet.txt')
#file_name = os.path.join(path, 'data/known_causet_from_matrixSetsTest.txt')

if isBH:
    if isEFv and isS:
        file_name = f"data/blackhole_EFvToS_{dim}D_N{card}"
    elif isEFv:
        file_name = f"data/blackhole_EFv_{dim}D_N{card}"
    elif isS:
        file_name = f"data/blackhole_S_{dim}D_N{card}"
    else:
        file_name = f"data/blackhole{dim}D_N{card}"
else:
    file_name = f"data/flat{dim}D_N{card}"
file_name = path + "/"+ file_name
#path = os.chdir("scripts_py")
if use_redge_in_name_file:
    edge_string = str(edge)
    if len(edge_string) == 3: edge_string += "0"
    file_name += "_redge"+ edge_string
file_name += ".txt"

print(file_name)
C.create_EmbeddedCauset_from_file(file_name)


print("Dim:", C.dim)
# print("Coords:\n", C.coords)
# print("pasts:\n", C.pasts)
# print("futures:\n", C.futures)
# print("Fut_links:\n", C.fut_links)
# print("Past_links:\n", C.past_links)

# Plotting setup:
cplt.setDefaultColors('UniYork')  # using University of York brand colours

if C.dim==3:
    dims: List[int] = [1,2,0]  # choose the (order of) plot dimensions
elif C.dim == 2:
    dims: List[int] = [1,0]
cplt.figure(figsize=(15.0, 12.0))
if len(dims) > 2:
    cplt.figure(figsize=(15.0, 12.0))
#S.plot(dims)  # plot the embedding shape
# Add causet plots and show result:
cplt.plot(C, dims=dims, events={'alpha': 1},
          links={'alpha': 0.4, 'linewidth': 0.5}, labels=False)

ax: cplt.Axes = cplt.gca()
ax.set_xlabel('space' if dims[0] > 0 else 'time')
ax.set_ylabel('space' if dims[1] > 0 else 'time')
if len(dims) > 2:
    ax.set_zlabel('space' if dims[2] > 0 else 'time')
    ax.grid(False)

# Plot horizon in 2D
ax.axvline(2,0,20,color="red", ls="--", label = "Horizon")
ylims = ax.get_ylim()

# Plot cones inside horizon crossing a point (t0, r0)
if len(dims) == 2:
    rs = [0.72, 1.35, 2.77]
    ts = [1.15, 2.66, 1.35]
    
    def upper_null(r, t0, r0, mass=1):
        return t0 + r - r0 + 4*mass*np.log( (2*mass-r)/(2*mass-r0) )
    def lower_null(r, t0, r0):
        return t0 + r0 - r

    for i, (r0, t0) in enumerate(zip(rs, ts)):
        rs = np.linspace(0, r0, 1000)
        ax.plot(rs, lower_null(rs, t0, r0), ls = "--", c = "green", alpha = 0.8)
        if r0 == 2*mass:
            continue
        elif r0 > 2*mass:
            rs = np.linspace(r0, edge, 1000)
        if i != len(ts) - 1:
            ax.plot(rs, upper_null(rs, t0, r0), ls = "--", c = "green", alpha = 0.8)
        else:
            ax.plot(rs, upper_null(rs, t0, r0), ls = "--", c = "green", alpha = 0.8,
                    label = "Light Cone")
ax.set_ylim(ylims)
ax.legend()

#ax.axvline(-2,0,20,color="red", ls="--")
cplt.show()

# %%
