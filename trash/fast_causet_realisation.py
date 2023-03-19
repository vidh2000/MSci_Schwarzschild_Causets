#%%

import numpy as np
from time import time
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
from numba import *

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
    if (dt2-dspace2)>0:
        return True
    else:
        return False

@jit
def make_causet_matrix(coords):
    dims, dim, size = (len(coords[0]),len(coords[0])), len(coords[0]),len(coords)
    cm = np.zeros(dims,dtype=int)
    for i in range(size):
        for j in range(i):
            if cm[i][j] ==1:
                continue
            else:
                dt2 = (coords[i][0]-coords[j][0])*(coords[i][0]-coords[j][0])
                dspace2 = sum([(coords[i][p]-coords[j][p])*
                        (coords[i][p]-coords[j][p]) for p in range(1,dim)])
                ds2 = dt2-dspace2 
                if ds2>0:
                    cm[i][j] = 1
                    ## Inherit the past
                    #cm[i,:] = np.bitwise_or(cm[i,:], cm[j,:])
                else:
                    pass
            
    return cm
        
def mp_add_pasts(pasts, coords, dim):
    for i in range(len(pasts)):
        for j in reversed(range(i)):
            if j in pasts[i]:
                continue
            else:
                if areTimelike(coords[i],coords[j],dim):
                    pasts[i].add(j)
                    pasts[i].update(pasts[j])
    return pasts

class Causet():
    def __init__(self, coordinates, mode = "cmatrix", use_transitivity=True):
        self.coords = coordinates
        self.dim = len(coordinates[0])
        self.size = len(coordinates)
        self.pasts = [set() for i in range(self.size)]

        if mode == "cmatrix":
            self.make_cmatrix(use_transitivity)
        elif mode=="pasts":
            self.make_past_sets()
        else:
            raise ValueError("Please provide 'cmatrix' or 'pasts' as method\
                                parameter.")
    

    def make_past_sets(self,multiprocessing=False):
        if multiprocessing:
            with mp.Pool() as pool:
                self.pasts = pool.starmap(mp_add_pasts,
                            zip(self.pasts, self.coords,
                                [self.dim for i in range(self.size)]))
            pool.close()
            pool.join()
        else:
            for i in tqdm(range(self.size)):
                for j in reversed(range(i)):
                    if j in self.pasts[i]:
                        continue
                    else:
                        if areTimelike(self.coords[i],self.coords[j],self.dim):
                            self.pasts[i].add(j)
                            self.pasts[i].update(self.pasts[j])
          
    def make_cmatrix(self,use_transitivity=True,jit=False):
        # Uses jit 

        if jit:
            self.causal_matrix = make_causet_matrix(self.coords)
        else:
            dims = (self.size,self.size)
            cm = np.zeros(dims,dtype=np.int32)
            if use_transitivity:
                print("Making cmatrix with transitivity")
                for i in range(1,self.size):
                    for j in range(i-1,-1,-1):
                        if cm[i][j] == 1:
                            continue
                        else:
                            if areTimelike(self.coords[i],self.coords[j],self.dim):
                                cm[i][j] = 1
                                # Inherit the past
                                cm[i,:] = np.bitwise_or(cm[i,:], cm[j,:])
                            else:
                                pass
            else:
                print("Making cmatrix without transitivity")
                for i in range(self.size):
                    for j in range(i-1,-1,-1):
                        if areTimelike(self.coords[i],self.coords[j],self.dim):
                            cm[i][j] = 1

if __name__=="__main__":
    dim=4
    N=int(5000)

    coords = make_coords(N,dim)

    # 'pasts'
    start =  time()
    c = Causet(coords,mode = "cmatrix",use_transitivity=True)
    duration = round(time()-start,3)
    print(f"Time taken for 'pasts' for N={N} in D={dim}\nt = {duration}s")


    ## 'cmatrix'
    #start =  time()
    #Causet(coords,mode = "cmatrix")
    #duration = round(time()-start,3)
    #print(f"Time taken for 'cmatrix', t = {duration}s")



# %%
