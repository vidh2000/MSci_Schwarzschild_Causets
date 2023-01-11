import matplotlib.pyplot as plt 
import numpy as np 
from os.path import expanduser
import os as os


# Home Directory
home = expanduser("~")
plotsDir = home + "/MSci_Schwarzcshild_Causets/figures/Nlinks_vs_Area/"
dataDir = home + "/MSci_Schwarzschild_Causets/data/linkcounting_files"
# Get data
for file in os.dirname(dataDir):
    print(file)



# plt.figure()
# plt.errorbar(masses**2,means,yerr=stds,capsize=4,ls="-",fmt="o",
#     label=r"$3+1D$ Schwarzshild spacetime")
# plt.legend()
# plt.ylabel(r"Number of links $N$")
# plt.xlabel(r"Area$\propto M^2$ [a.u]")
# plt.grid(alpha=0.3)

# figname = plotsDir+"Nlinks_vs_Area_4D_Nmult=1000_M=1.0-5.0._100reps.png"
# plt.savefig(figname)
# plt.show()