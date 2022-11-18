import numpy as np
from numba import njit,jit
from time import time
import pandas as pd

file = "D:\Documents/Sola/Imperial College London/Year 4/MSci project/MSci_Schwarzschild_Causets/data/flatspace_bicone_causet.txt"

size = 100

lines = [[v for v in line.split(",")] for line in open(file)]

print("LINES:")
for i,l in enumerate(lines):
    print(l)
    if i>15:
        break

pasts = [[int(v) for v in line if v != "\n"] for line in lines[6:6+size]]
futures = [[int(v) for v in line if v != "\n"] for line in
                                lines[6+size+1:6+2*size+1]]
plinks = [[int(v) for v in line if v != "\n"] for line in
                                lines[6+2*size+2:6+3*size+2]]
flinks = [[int(v) for v in line if v != "\n"] for line in
                                lines[6+3*size+3:6+4*size+3]]
coords = [[float(v) for v in line if v != "\n"] for line in
                                lines[-size:]]
print("PASTS")
print(pasts)
print("\nFUTURES")
print(futures)
print("\nplinks")
print(plinks)
print("\nflinks")
print(flinks)
print("\nCoords")
print(coords)

print("Lengths:", len(pasts), len(futures),len(plinks), len(flinks),len(coords))