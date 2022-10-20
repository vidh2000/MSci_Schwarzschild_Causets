#%%
import numpy as np 

a = set((1,2,3,4))

b = set((4,5,6))

c=  a & b
print(c)
print(type(c))

c = len(a)>len(b)
print(c)
