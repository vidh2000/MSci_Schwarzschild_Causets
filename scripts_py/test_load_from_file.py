
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
import os 

########################################################
# Load file
dim = 2
card = 500
use_redge_in_name_file = True
edge = 5.000000

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
    file_name += "_redge_"+ str(edge)
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
ax.axvline(2,0,20,color="red", ls="--")
#ax.axvline(-2,0,20,color="red", ls="--")
cplt.show()

# %%
