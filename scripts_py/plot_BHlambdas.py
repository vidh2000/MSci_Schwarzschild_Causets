
#%% 
from __future__ import annotations
from typing import List, Tuple  
from causets.sprinkledcauset import SprinkledCauset
from causets.embeddedcauset import EmbeddedCauset  
from causets.spacetimes import FlatSpacetime, deSitterSpacetime  
from causets.shapes import CoordinateShape  
from causets.causetevent import CausetEvent  
import causets.causetplotting as cplt  

import causets.OURcausetplotting as ocplt

import matplotlib.pyplot as plt
import numpy as np
import os 

########################################################
# Load file
dim = 2
card = 25
use_redge_in_name_file = True
edge = 5.00

isBH = True
isEFv = False
isS = False


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
        file_name = f"data/blackhole_and_lambdas{dim}D_N{card}"
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

ocplt.plot_causet_and_lambdas_and_lambdasdistr(file_name)

# %%
