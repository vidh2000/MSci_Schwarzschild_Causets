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


####################################################################
# 0. SET SETTINGS
####################################################################
Nreps = 5
sizes = [64, 128, 12288]#, 512, 1024, 2048, 4096]

file_saved_by_cpp =\
    f"data/test_MMdim_forpy/MMdim_Flat_Nreps{Nreps}_UpTo{max(sizes)}.txt"


#####################################################################
## 1. COMPILE C++ FILE RUNNING MMDIM
#####################################################################
cwd = os.getcwd()

# Path to MSci_Schwarzschild_Causets
file_dir = os.path.dirname(os.path.realpath(__file__))
main_dir = os.path.dirname(file_dir)
print(main_dir)
#remove scripts_py
main_dir = main_dir[:-17]
print(main_dir)

# Specify which file you want to run: 
runfileDir = f"{main_dir}scripts_cpp/tests/test_causets_general/"
runfilename = "MMdim_test_forpy.cpp"

# # Directory path of causets_cpp
# causets_cpp = main_dir + "scripts_cpp/causets_cpp/"

# # Compile commands
# compiler_stuff = [
#     f"'C:\\Program Files\\MinGW-w64\\bin\\g++.exe'",
#     f" -g '{runfileDir+runfilename}'",
#     f"'{causets_cpp}causet.cpp'",
#     f"'{causets_cpp}embeddedcauset.cpp'",
#     f"'{causets_cpp}sprinkledcauset.cpp'",
#     f"'{causets_cpp}shapes.cpp'",
#     f"'{causets_cpp}spacetimes.cpp'",
# ]

# # Create Executable
# ############################################
# exec = [ " ",
#     "-o",
#     f'"{runfileDir}{runfilename[:-4]}.exe"'
# ]

# #Optimisation flags
# ############################################
# optimisations = [ " ",
#     #"-std=c++17",   
#     "-Ofast",
#     "-fopenmp",
#     ]

# # Include paths
# ###########################################
# from os.path import expanduser
# home = expanduser("~")
# includes = [" ",
#     "-I",
#     '"C:\\Program Files\\boost\\boost_1_80_0"',
#     ]

# # Compile
# ############################################
# to_render = " ".join(compiler_stuff) + " ".join(exec) + \
#             " ".join(optimisations) + " ".join(includes)

# print("\nC++ files to compile, flags used, executable produced:")
# print(to_render)
# os.system(to_render)
# print(".............................Compiled............................!\n")



# #####################################################################
# ## 2. RUN C++ EXECUTABLE RUNNING MMDIM
# #####################################################################
# print(f"'{runfileDir}{runfilename[:-4]}.exe' {Nreps} {sizes}")
# os.system(f"'{runfileDir}{runfilename[:-4]}.exe' {Nreps} {sizes}")
# print(".............................Run............................!\n")



#####################################################################
## 3. LOAD SAVE FILE AND PLOT
#####################################################################
info = np.loadtxt(f"{main_dir}{file_saved_by_cpp}", skiprows = 1)

plt.figure((7,7))
for (d,line) in enumerate(info):
    ests = []
    stds = []
    for i in range(len(line)/2):
        ests.append(line[2*i])
        stds.append(line[2*i+1])
    
    plt.errorbar(sizes, ests, stds,
                 fmt = '.', capsize = 2, label = f"{d+1}D")
    plt.grid(alpha = 0.3)

plt.xlabel("Cardinality")
plt.ylabel("MM Dimension")
plt.savefig(f"{main_dir}figures/MMd/MMdim_Flat_Nreps{Nreps}_UpTo{max(sizes)}.png")
plt.savefig(f"{main_dir}figures/MMd/MMdim_Flat_Nreps{Nreps}_UpTo{max(sizes)}.pdf")
plt.show()
    