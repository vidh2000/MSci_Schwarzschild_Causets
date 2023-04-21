import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# #=======================================================================
# # A CPP FILE, scripts_cpp/tests/MMdim_test_forpy.cpp, IS RUN. IT
# # PERFORMS THE Mirheim-Mayer ESTIMATOR FOR CAUSETS IN FLAT SPACETIME,
# # WITH:
# # - dimensions go from 1 to 4
# # - sizes given by sizes variable below
# # - number of repetitions per realisations given by Nrep variable below
# #
# # THE CPP FILE SAVES A TXT FILE, WITH NAME file_saved_by_cpp, WHICH IS
# # THEN LOADED BY THIS PYTHON FILE, SO THAT RESULTS CAN BE PLOTTED. 
# #=======================================================================

#plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {#'text.usetex' : True,
          'font.size' : 18,
          'font.family' : 'Times New Roman',
          'axes.labelsize':18,
          'legend.fontsize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'axes.prop_cycle':plt.cycler(color=
                            plt.rcParams['axes.prop_cycle'].by_key()['color']
                            +['magenta'])
          }
plt.rcParams.update(params)

####################################################################
# 0. SET SETTINGS
####################################################################
Nreps = 20
maxsize = 16384

file_saved_by_cpp =\
    f"data/test_MMdim_forpy/MMdim_Flat_Nreps{Nreps}_UpTo{maxsize}.txt"


#####################################################################
## 1. COMPILE C++ FILE RUNNING MMDIM
#####################################################################
cwd = os.getcwd()

# Path to MSci_Schwarzschild_Causets
file_dir = os.path.dirname(os.path.realpath(__file__))
main_dir = os.path.dirname(file_dir)
#remove scripts_py
main_dir = main_dir[:-10]
print(main_dir)

# Specify which file you want to run: 
runfileDir = f"{main_dir}scripts_cpp/tests/test_causets_general/"
runfilename = "MMdim_test_forpy.cpp"


#####################################################################
## 3. LOAD SAVE FILE AND PLOT
#####################################################################
labels = np.loadtxt(f"{main_dir}{file_saved_by_cpp}", dtype = str, 
                  delimiter = ",", skiprows = 0, max_rows=1)
data = np.loadtxt(f"{main_dir}{file_saved_by_cpp}", 
                  delimiter = ",", skiprows = 1)
sizes = []

plt.figure(figsize = (10,6), tight_layout = True)
for (dim, line) in enumerate(data):
    ests = []
    stds = []
    for i in range(int(len(line)/2)):
        if dim == 0:
            sizes.append(int(labels[2*i][:-3]))
        ests.append(line[2*i])
        stds.append(line[2*i+1])
    
    plt.errorbar(sizes, ests, stds,
                 fmt = '.', capsize = 2, label = f"{dim+1}D")
    plt.hlines(np.arange(2, dim+3), 0, plt.xlim()[1], 
               ls = "--", color = "r")
    plt.xscale("log")
    plt.yticks(np.arange(1.75, 4.01, 0.25))
    plt.grid(alpha = 0.3)
    

plt.xlabel("Cardinality")
plt.ylabel("MM Dimension")
plt.savefig(f"{main_dir}figures/MMd/MMdim_Flat_Nreps{Nreps}_UpTo{max(sizes)}.png")
plt.savefig(f"{main_dir}figures/MMd/MMdim_Flat_Nreps{Nreps}_UpTo{max(sizes)}.pdf")
plt.show()
    