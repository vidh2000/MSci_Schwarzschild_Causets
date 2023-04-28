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
Nreps = 60
sizes = [128,256,512,1024,2048,4096, 8192,16384,32768]#, 512, 1024, 2048, 4096]

file_saved_by_cpp =\
    f"data/test_MMdim_forpy/MMdim_Flat_Nreps{Nreps}_UpTo{max(sizes)}.txt"


#####################################################################
## 1. COMPILE C++ FILE RUNNING MMDIM
#####################################################################
main_dir = os.getcwd()
main_dir = main_dir+"/"
# # Path to MSci_Schwarzschild_Causets
file_dir = os.path.dirname(os.path.realpath(__file__))
#main_dir = os.path.dirname(file_dir)
# print(main_dir)
#remove scripts_py
# main_dir = main_dir[:-17]
# print(main_dir)

# # Specify which file you want to run: 
# runfileDir = f"{main_dir}scripts_cpp/tests/test_causets_general/"
# runfilename = "MMdim_test_forpy.cpp"
 
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
params = {'text.usetex' : True,
          'font.size' : 25,
          'font.family' : 'Times New Roman',
          'axes.labelsize': 28,
          'legend.fontsize': 20,
          'xtick.labelsize': 22,
          'ytick.labelsize': 22,
          'figure.figsize': [8.5, 6.5], 
          'axes.prop_cycle':plt.cycler(color=
                            plt.rcParams['axes.prop_cycle'].by_key()['color']
                            +['magenta'])
          }
plt.rcParams.update(params)


info = np.loadtxt(f"{main_dir}{file_saved_by_cpp}", skiprows = 1,
                  delimiter=",")

print(info)

plt.figure()
for (d,line) in enumerate(info):
    ests = []
    stds = []
    for i in range(int(len(line)/2)):
        ests.append(line[2*i])
        stds.append(line[2*i+1])
    
    plt.errorbar(sizes, ests, yerr=stds,
                 fmt = '.', capsize = 4, label = f"D={d+2}")
    plt.hlines(d+2,128,32800,linestyles="dashed", color="red")
    #plt.grid(alpha = 0.3)

plt.xlabel(r"Cardinality")
plt.ylabel(r"MM Dimension")
plt.legend(loc="lower right")
plt.xscale("log")
plt.ylim((1,4.1))
plt.tight_layout()
plt.savefig(f"{main_dir}figures/MMd/MMdim_Flat_Nreps{Nreps}_UpTo{max(sizes)}.png")
plt.savefig(f"{main_dir}figures/MMd/MMdim_Flat_Nreps{Nreps}_UpTo{max(sizes)}.pdf")
plt.show()
    