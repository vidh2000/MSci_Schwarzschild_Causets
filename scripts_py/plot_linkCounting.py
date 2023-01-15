import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd
from os.path import expanduser
import os



# Home Directory
home = expanduser("~")

Poiss="True"
plotsDir = home + "/MSci_Schwarzschild_Causets/figures/Nlinks_vs_Area/"
dataDir = home + f"/MSci_Schwarzschild_Causets/data/linkcounting_files/Poiss={Poiss}"

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
cards_arr       = []
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
    card = int(filename[filename.find("Card")+5:filename.find("_r=")])
    print(filename)
    print(card)
    filename = dataDir +"/"+filename
    data = pd.DataFrame(pd.read_csv(filename, sep = ",",header=None,
                                    skiprows=0))
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
x = 4*np.pi*np.array(masses)**2
y = Nlinks_avgs
plt.errorbar(x, y, yerr=Nlinks_stds,
            capsize=4,ls="",fmt="o", color="black",
            label=r"$3+1D$ Schwarzshild spacetime")

# Linear fit
coef = np.polyfit(x,y,1)
print(f"===============================================================\n\
        Linear fit coefficients: {coef}\n\
        ===============================================================")
poly1d_fn = np.poly1d(coef) 
x = np.linspace(min(x),max(x),100)
plt.plot(x, poly1d_fn(x), '--', color="red",
        label=f"Linear fit")

plt.legend()
plt.ylabel(r"Number of links $N$")
plt.xlabel(r"Area$\propto M^2$ [a.u]")
plt.grid(alpha=0.3)

figname = plotsDir+f"Nlinks_vs_Area_4D_Rho={rho}_Poiss={Poiss}.png"
plt.savefig(figname)
plt.show()