import numpy as np
from scipy import special as spsp
import matplotlib.pyplot as plt

ts = np.linspace(0, 4, 4000)
ys = (
        2/np.sqrt(3)*ts**2*np.exp(-np.pi/3*ts**4) 
        + spsp.erfc(np.sqrt(np.pi/3)*ts**2)    
      )*np.sqrt(3)/10

plt.plot(ts, ys*1e4) #1e4 for the area
plt.hlines(1, plt.xlim()[0], plt.xlim()[1], ls = "--", color = "green")
plt.hlines(1/100, plt.xlim()[0], plt.xlim()[1], ls = "--", color = "red")
plt.show()
