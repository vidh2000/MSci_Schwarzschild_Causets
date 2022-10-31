#%%
from pickletools import int4
import numpy as np
from time import time
from tqdm import tqdm
import multiprocessing as mp

def make_coords(N,dim,min=0,max=1,sorted=True):
    coords = []
    for i in range(N):
        pos_i = list(min + (max-min)*np.random.rand(dim))
        coords.append(pos_i)
    if sorted:
        coords.sort(key=lambda x: x[0])
    return coords

def areTimelike(x,y,dim):
    dt2 = (x[0]-y[0])*(x[0]-y[0])
    dspace2 = sum([(x[i]-y[i])*(x[i]-y[i]) for i in range(1,dim)])
    ds2 = dt2-dspace2 
    if ds2>0:
        return True
    else:
        return False



class Causet():
    def __init__(self, coordinates, mode = "cmatrix"):
        self.coords = coordinates
        self.dim = len(coordinates[0])
        self.size = len(coordinates)
        self.pasts = [set() for i in range(self.size)]

        if mode == "cmatrix":
            self.make_cmatrix()
        elif mode=="pasts":
            self.make_past_sets()
        else:
            raise ValueError("Please provide 'cmatrix' or 'pasts' as method\
                                parameter.")
    
    def make_past_sets(self):
        
        for i in tqdm(range(self.size)):
            for j in reversed(range(i)):
                if j in self.pasts[i]:
                    continue
                else:
                    if areTimelike(self.coords[i],self.coords[j],self.dim):
                        self.pasts[i].add(j)
                        self.pasts[i].update(self.pasts[j])
                
    def make_cmatrix(self):
        
        dims = (self.size,self.size)
        cm = np.zeros(dims,dtype=int)
        for i in tqdm(range(self.size)):
            for j in reversed(range(i)):
                if cm[i][j] ==1:
                    continue
                else:
                    if areTimelike(self.coords[i],coords[j],self.dim):
                        cm[i][j] = 1
                        # Inherit the past
                        cm[i,:] = np.bitwise_or(cm[i,:], cm[j,:])
        self.causal_matrix = cm

dim=2
N=int(10**4)

coords = make_coords(N,dim)


# 'pasts'
start =  time()
Causet(coords,mode = "pasts")
duration = round(time()-start,3)
print(f"Time taken for 'pasts', t = {duration}s")


# 'cmatrix'
start =  time()
Causet(coords,mode = "cmatrix")
duration = round(time()-start,3)
print(f"Time taken for 'cmatrix', t = {duration}s")


