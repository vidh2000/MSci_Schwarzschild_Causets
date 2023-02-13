import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd
from os.path import expanduser
import os
from scipy.optimize import curve_fit



usehome = True
#Poiss="False"

# Home Directory
home = expanduser("~")
# Path
path = os.getcwd()

if usehome:
    plotsDir = home + "/MSci_Schwarzschild_Causets/figures/Nlinks_vs_Area/"
    dataDir = home + f"/MSci_Schwarzschild_Causets/data/links/"
    #dataDir = home + f"/MSci_Schwarzschild_Causets/data/linkcounting_files/Poiss=False/"
else:
    plotsDir = path + "/figures/Nlinks_vs_Area/"
    dataDir = path + f"/data/linkcounting_files/Poiss={Poiss}"

# Variables to select runs you want
#Ms = [1.0,1.2,1.4,2.0,2.4,3.0] 
Ms = None #=None if you want to choose all existinig masses
rho = 20000


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
    if f"Rho={rho}_" not in filename:

        cont = False
    if cont==False:
        continue
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



### Plot the N_vs_Area [l^2]

# Determine the scale
scale = (rho/(4*np.pi/3*26))**(-1/4)
print("Scale = ", scale)
x = 4*np.pi*(4*np.array(masses)**2)
x = x/scale**2
y = Nlinks_avgs

plt.figure()
plt.errorbar(x, y, yerr=Nlinks_stds,
            capsize=4,ls="",fmt="o", color="black",
            label=r"$3+1D$ Schwarzshild spacetime")

### Linear fit

# Fitting function == linear through vertex
def lin_func(x,a):
    return a*x

# fit the function to the data
popt, pcov = curve_fit(lin_func, x, y, sigma=Nlinks_stds,
                        absolute_sigma=True)
unc = np.sqrt(np.diag(pcov))

# print the optimal parameter values
print(f"===============================================================")
print(f"\tGradient factor w.r.t \sqrt(rho) = {round(popt[0],3)} +- {round(unc[0],3)}")
print(f"===============================================================")
 
x = np.linspace(min(x),max(x),100)
plt.plot(x, lin_func(x,*popt), '--', color="red",
        label=f"Linear fit")

plt.legend()
plt.ylabel(r"Number of links $N$")
plt.xlabel(r"Horizon Area [$\ell^2$]")
plt.grid(alpha=0.3)

figname = plotsDir+f"Nlinks_vs_Area_4D_Rho={rho}_final.png"
plt.savefig(figname)
figname = plotsDir+f"Nlinks_vs_Area_4D_Rho={rho}_final.pdf"
plt.savefig(figname)
plt.show()