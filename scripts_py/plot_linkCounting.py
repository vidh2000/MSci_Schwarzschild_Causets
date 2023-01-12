import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd
from os.path import expanduser
import os


# Home Directory
home = expanduser("~")
plotsDir = home + "/MSci_Schwarzschild_Causets/figures/Nlinks_vs_Area/"
dataDir = home + "/MSci_Schwarzschild_Causets/data/linkcounting_files"

# Variables to select runs you want
#Ms = [1.0,1.2,1.4,2.0,2.4,3.0] 
Ms = None #=None if you want to choose all existinig masses
rho = 1000


# Get data
if Ms == None:
    masses = []
else:
    masses = Ms
Nlinks_avgs     = []
Nlinks_stds     = []
totReps_arr     = []
for i,filename in enumerate(os.listdir(dataDir)):
    
    # Select wanted files
    cont = True
    if f"Rho={rho}" not in filename:
        cont = False
    if Ms == None:
        M = float(filename[2:6])
        masses.append(M)
    else:
        for M in Ms:
            if f"M={M}" not in filename:
                cont = False
    if cont==False:
        continue

    # Extract useful information and calculate avgs/stds
    print(filename)
    filename = dataDir +"/"+filename
    data = pd.DataFrame(pd.read_csv(filename, sep = ",",header=None))
    data = data.iloc[:,[3,4,5]].values
    
    avg = 0
    std = 0
    for gen in data:
        avg += gen[0]*gen[1]
        std += gen[0]*gen[2]

    totReps= sum(data[:,0])
    totReps_arr.append(totReps)

    Nlinks_avgs.append(avg/totReps)
    Nlinks_stds.append(std/totReps)



plt.figure()
plt.errorbar(np.array(masses)**2, Nlinks_avgs, yerr=Nlinks_stds,
            capsize=4,ls="",fmt="o", color="black",
            label=r"$3+1D$ Schwarzshild spacetime")
plt.legend()
plt.ylabel(r"Number of links $N$")
plt.xlabel(r"Area$\propto M^2$ [a.u]")
plt.grid(alpha=0.3)

figname = plotsDir+f"Nlinks_vs_Area_4D_Rho={rho}.png"
plt.savefig(figname)
plt.show()