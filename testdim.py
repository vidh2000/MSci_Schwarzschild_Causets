#!/usr/bin/env python
"""
Created on 13 Oct 2022

@author: Stefano Veroni
"""
#%%
from __future__ import annotations
from typing import List, Tuple  # @UnusedImport

from causets.causet import Causet
import causets.causetplotting as cplt

import numpy as np
x = [0,1]
print(sum(x[1:]))
#%%
N = 8
Cm = np.zeros([N,N])
for i in range(N-1):
    Cm[i, i+1] += 1

C = Causet().FromFutureMatrix(Cm)
print("Causet has cardinality ", len(C))
#cplt.plot(C)

#%%
Clist = list(C)
Clist.sort(key = lambda e : e.PastCard)
for e in Clist:
    print(e.PastCard)
#%%
a = Clist[3]
b = Clist[6]
afuture = a.Future
bpast = b.Past
A = C.Interval(a, b)
Alist = list(A)
Alist.sort(key = lambda e : e.PastCard)

c = Alist[0]
print(c.PastCard)
# %%
