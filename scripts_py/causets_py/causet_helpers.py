#!/usr/bin/env python

'''
Created on 15 Jan 2022

@author: Stefano Veroni, Vid Homsak
'''

import numpy as np
import functools

def get_causet_attrs (lambdasfile_ext):
    """

    Parameters
    ----------
    lambdasfile_ext : str
        Name of file from which import info
    
    Returns
    ---------
    storage_option : str (#0)

    size : int (#1)

    dim : int (#2)

    shapename : str (#3)

    spacetimename : str (#4)

    coords : list<list<float>> (#5)

    fut_links : list<list<int>> (#6)
        Note, some are empty

    r_S : float (#7)
    
    lambdas : list<list<int>> (#7)

    distribution : list<list<int>> (#8)
    """
    with open(lambdasfile_ext, 'r') as fl: 
            f = fl.readlines() 

            storage_option = str(f[0].split(",")[1])
            size           = int(f[1].split(",")[1])
            dim            = int(f[2].split(",")[1])
            shapename      = str(f[3].split(",")[1])
            spacetimename  = str(f[4].split(",")[1])
            if dim > 3:
                raise AttributeError(f"Dim is {dim}: too big!!!")

            distribution = []
            lambdas = []
            coords = []
            r_S = 0

            go = -1
            while go < 0:
                row = f[go].split(",")
                key = row[0]
                if key[0] == "L":
                    lambdas.append([])
                    for label in row[1:]:
                        if label != "" and label != "\n":
                            lambdas[-1].append(int(label))
                    go -= 1
                elif key[0] == "N":
                    distribution.append([int(key[-1]), int(row[1])])
                    go -= 1
                elif key == "r_S" or key == "r_S\n":
                    r_S = float(row[1])
                    go -= 1
                elif key == "Coordinates" or key == "Coordinates\n":
                    break
                elif key == "":
                    go -= 1
                else:
                    coords.insert(0, [])
                    for i in range(dim):
                        coords[0].append(float(row[i]))
                    go -= 1
        
            if storage_option == "cmatrix" or storage_option == "cmatrix\n":
                cmatrix = []
                for i in range(6, 6+size):
                    cmatrix.append(f[i].split(","))
                cmatrix = np.array(cmatrix, dtype = 'int')
                cmatrix2 = np.matmul(cmatrix, cmatrix)
                for i in range(size):
                    fut_links.append[[]]
                    for j in range(size):
                        if (cmatrix[i,j] and cmatrix2[i,j]==0):
                            fut_links[i].append[j]
            else:
                lines= [[v for v in line.split(",")] for line in open(lambdasfile_ext)]
                fut_links = [[int(v) for v in line if (v != "\n" and v != "")] 
                                for line in lines[6+3*size+3:6+4*size+3]]

    return [storage_option, #0
            size,           #1
            dim,            #2
            shapename,      #3
            spacetimename,  #4

            coords,         #5
            fut_links,      #6
            r_S,            #7

            lambdas,        #7
            distribution]   #8

#from from https://stackoverflow.com/questions/15301999/default-arguments-with-args-and-kwargs
def default_kwargs(**defaultKwargs):
    """
    Decorator to make defaultkawargs the default kwargs in a function.
    """
    def actual_decorator(fn):
        @functools.wraps(fn)
        def g(*args, **kwargs):
            defaultKwargs.update(kwargs)
            return fn(*args, **defaultKwargs)
        return g
    return actual_decorator

def std_combine(Ns, mus, stds):
    """
    Function combining the stds of M rounds i of Ns[i[ measurements, each 
    having mean mus[i] and std stds[i].

    Note: M = len(Ns) = len(ms) = len(stds).

    Parameters
    ----------
    Ns : arraylike[int]
        Measurements of each ith round.
    
    mus : arraylike[float]
        Mean of each ith round of measurements.
    
    stds : arraylike[float]
        Std of each ith round of measurements.
    
    Returns
    -------
    std : float
        Std of the combination of all measurements.
    """
    coeffs = []
    M = len(Ns)
    N = sum(Ns)
    for i in range(M):
        Ni = Ns[i]
        mui = mus[i]
        stdi = stds[i]
        term_i = (Ni-1)/(N-1)*stdi**2
        coeffs.append(term_i)
        for j in range(i, M):
            Nj = Ns[j]
            muj = mus[j]
            stdj = stds[j]
            term_mixed = Ni*Nj/(N*(N-1)) * (mui - muj)**2
            coeffs.append(term_mixed)
    return np.sqrt(sum(coeffs))