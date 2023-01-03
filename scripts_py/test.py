import matplotlib.pyplot as plt
import numpy as np 

# 4D data with futlinks creation...
#masses = np.array([1,1.5,2])
#means= [7,16.875,29]
#stds = [2.915,3.06,6.55]
# 4D data with futmatrix...
#masses = np.array([1,1.5,2])
#means2 = [8.8,16,24]
#stds2 = [2.92,3.35,3.79]

# 4D data with futmatrix...
# masses = np.array([1,1.5,2,2.5,3])
# means= [7.8,15.6,23,41.2,51.2]
# stds = [1.7,3.72,5.83,2.04,8.9]

#4D data from cluster runs: [100,100,100,50,20,20] realisations
# masses = np.array([1,1.5,2,2.5,3,3.5,4])
# means= [7.06,15.33,26.63,38.93,54.26,71.75,93.25]
# stds = [2.89,4.47,6.24,6,7.2,9.6,11.1]

#4D data from cluster runs: 100 realisations per M
masses = np.array([1,1.5,2,2.5,3,3.5,4,4.5,5])
means = [11.4,23.83,37.57,60.55,83,111.55,144.98,184.08,221.84]
stds = [3.7,5.8,6.5,8.3,10.5,10.5,12.9,13.7,17.6]




plt.figure()

plt.errorbar(masses**2,means,yerr=stds,capsize=4,ls="-",fmt="o",
    label=r"$3+1D$ Schwarzshild spacetime")
plt.legend()
plt.ylabel(r"Number of links $N$")
plt.xlabel(r"Area$\propto M^2$ [a.u]")
plt.grid(alpha=0.3)
plt.savefig("/rds/general/user/vh119/home/MSci_Schwarzschild_Causets/figures/Schwarzschild/Nlinks_vs_Area_4D_Nmult=1000_M=1.0-5.0._100reps.png")
plt.show()