import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os
from os.path import expanduser
from causets_py import causet_helpers as ch
from scipy.optimize import curve_fit
from scipy.stats import chisquare, pearsonr, skew, kurtosis
from scipy.sparse.csgraph import connected_components

###############################################################################
# Open files HRVs_connectivity/M=<>_Rho=<>_Card=<>_r=<>_hollow=<>_dur=<>_n.txt
# each containing the HRVs molecules from the nth run of a causet with 
# specified parameters.
#
# Then save a file with columns being # clusters of HRVs with certain size
# and associated standard deviation (like with lambdas), of title
# HRV_connectivity_compress/M=<>_Rho=<>_Card=<>_r=<>_hollow=<>_dur=<>.txt
###############################################################################


def get_info_of_HRVs(file_ext):
    """

    Parameters
    ----------
    - file_ext : str
        File from which import info

    Returns
    -------
    - size : int (#0)

    - dim : int (#1)

    - shapename : str (#2)

    - spacetimename : str (#3)

    - coords : list<list<float>> (#4)

    - r_S : float (#5)
    
    - molecules : list<list<int>> (#6) /not all files have that,
    so might return empty list if molecules are not saved.

    - Adj_matrix : np.ndarray(int, int)
        Adj_ij = 1 if i < j or j < i
        Irreflexive, i.e 0 diagonal.
    """
    with open(file_ext, 'r') as fl: 
            f = fl.readlines() 

            size           = int(f[1].split(",")[1])
            dim            = int(f[2].split(",")[1])
            shapename      = str(f[3].split(",")[1])
            spacetimename  = str(f[4].split(",")[1])

            mols = []
            coords = []
            r_S = 0

            #start reading file from last line
            go = -1
            while go < 0:
                row = f[go].split(",")
                key = row[0]

                if key == "":
                    go -= 1
                
                elif key[0:3] == "HRV":
                    mols.append([])
                    for label in row[1:]:
                        if label != "" and label != "\n":
                            mols[-1].append(int(label))
                    go =- 1

                elif "r_S" in key:
                    r_S = float(row[1])
                    go -= 1

                elif "Coordinates" in key:
                    break
                else: #the coordinates#
                    coords.insert(0, [])
                    for i in range(dim):
                        coords[0].append(float(row[i]))
                    go -= 1
            # finished reading file from last line

            # make adjacency matrix of connections
            N = np.amax(mols) +1
            Adj_matrix = np.zeros((N,N))
            for mol in mols:
                Adj_matrix[mol[0], mol[1]] = 1
                Adj_matrix[mol[1], mol[0]] = 1
                Adj_matrix[mol[0], mol[2]] = 1
                Adj_matrix[mol[2], mol[0]] = 1

    return [size,           #0
            dim,            #1
            shapename,      #2
            spacetimename,  #3

            coords,         #4
            r_S,            #5

            mols,           #6
            Adj_matrix]     #7
    

##############################################################################
# 0. SET USER VARIABLES
##############################################################################
molecules = "HRVs" #lambdas, HRVs
varying_var = "M"     #variable varying: can be M, Rho, Nmult (M best)
fixed_var = "Rho"     #variable fixed:   can be, M, Rho, Nmult (Rho best)
fixed_val = 5000     #value of fixed_var 
use_selected_masses = True #gives equal spacing

#plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 20,
          'font.family' : 'lmodern',
          'axes.labelsize':34,
          'legend.fontsize': 18,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'figure.figsize': [8.5, 6.5], 
          'axes.prop_cycle':plt.cycler(color=
                            plt.rcParams['axes.prop_cycle'].by_key()['color']
                            +['magenta'])
          }
plt.rcParams.update(params)


##############################################################################
# 1 SET USEFUL STRINGS & OTHER VARIABLES
##############################################################################
home = expanduser("~")
plotsDir = f"{home}/MSci_Schwarzschild_Causets/figures/N{molecules}_vs_Area/"
dataDir = f"{home}/MSci_Schwarzschild_Causets/data/{molecules}_connectivity/"
#dataDir = home + f"/MSci_Schwarzschild_Causets/data/linkcounting_files/Poiss=False/"
if not os.path.exists(dataDir):
    path = os.getcwd()
    plotsDir = path + f"/figures/N{molecules}_vs_Area/"
    dataDir = path + f"/data/{molecules}_connectivity/"
print(dataDir)
# ensure that, if using Mor Rho, it is decimal with exactly 2 dec digits

if fixed_var == "M" or fixed_var == "Rho":
    fixed_string = f"{fixed_var}=" + format(fixed_val,".2f")
else:
    fixed_string = f"{fixed_var}={fixed_val}"


selected_masses = np.array([0.53, 0.75, 0.92, 1.06, 1.19, 
                            1.30, 1.40, 1.50, 
                            1.59, 1.68, 1.76, 1.84, 1.91, 1.98, 
                            2.05, 2.12, 2.19, 2.25, 2.31, 2.37, 2.43,
                            0.65, 0.84, 0.99, 1.13, 1.24, 1.35, 1.45,
                            1.55, 1.63, 1.72, 1.88, 1.95,
                            1.8,
                            2.02, 2.09, 2.15, 2.22, 2.28, 2.34, 
                            2.4, 2.46 ])



##############################################################################
# 2. GET INFO FROM FILES
##############################################################################

# an entry for each M
Nreps = []
varying_values = []
connectivities = []


for root, dirs, files in os.walk(dataDir):
    # for each file file_i
    for i, file_i in enumerate(files):
        if fixed_string in file_i:
            go_on = False
            if varying_var+"=" in file_i:
                go_on = True
                file_i_path = os.path.join(root, file_i)
                file_i_pieces = file_i.split("_")
                # get the value of the varying variable of file_i
                for piece_j in file_i_pieces:
                    if varying_var+"=" in piece_j:
                        varying_var_i = float(piece_j.split("=")[1])
                        if use_selected_masses \
                            and varying_var_i not in selected_masses:
                            go_on = False
                            break #j loop
                        if varying_var_i not in varying_values:
                            print(varying_var+"=", varying_var_i)
                            varying_values.append(varying_var_i)
                            Nreps.append([])
                            connectivities.append([])
                        break
            # if the file is correct
            if go_on: 
                # 2.1 get # of components and info on dstribution of sizes of
                # components
                results = get_info_of_HRVs(file_i_path)
                Aij = results[-1]
                n_components, labels = connected_components(Aij)
                sizes = np.bincount(labels)
                mean = np.mean(sizes)
                std  = np.std(sizes)
                sk   = skew(sizes)
                kur  = kurtosis(sizes)

                # update
                index = varying_values.index(varying_var_i)
                Nreps[index] += 1
                connectivities[index].append(
                    n_components, mean, std, sk, kur)

#############################################################
# INFO ON NREPS
varying_values_copy = np.array(varying_values)
combined = list(zip(varying_values_copy, Nreps)) 
sorted_combined = sorted(combined, key=lambda x: x[0]) # sort by values
vals = [x[0] for x in sorted_combined] 
Nreps = [x[1] for x in sorted_combined]
maxrep = max(Nreps)

vals_not_max = [vals [i] for i in range(len(vals)) if Nreps[i] < maxrep]
reps_not_max = [Nreps[i] for i in range(len(vals)) if Nreps[i] < maxrep]
repstable = pd.DataFrame(
                np.column_stack(
                [vals_not_max, reps_not_max, maxrep-np.array(reps_not_max)]),
                columns = ["M Value", "Current Nreps", "Reps to Max"])
print(f"Maximum number of Nreps is {maxrep}")
print("For the following, the maximum number of rep was not reached")
print(repstable)


#############################################################
# INFO ON CONNECTIVITIES
con_means = [] # array of means of n_components, mean, std, sk, kur per Mass
con_stds  = [] # array of stds of  n_components, mean, std, sk, kur per Mass
for i in range(len(varying_values)):
    con_means.append(np.mean(connectivities[i], axis = 0))
    con_stds .append(np.std(connectivities[i], axis = 0, ddof = 1))



############################################################
# PLOTS OF RESULTS

# Set the right x for the varying-fixed values -> x  is Area in l^2 units
x = np.array(varying_values)
if varying_var=="M":
    x = 4*np.pi*(2*x)**2 #==Area
    if fixed_var == "Nmult":    
        Rho = fixed_val/(4/3*np.pi*26)
        print(f"You're doing Nmult {fixed_val} fixed!") 
    elif fixed_var == "Rho":
        Rho = fixed_val
    x /= Rho**(-1/2)
r_S_norm = np.sqrt(x / 4 / np.pi)   
print(f"Rho = {Rho:.0f}")
fixed_string = rf"Rho = {Rho:.0f}"


plt.figure("n_components per Mass", tight_layout = True)
plt.errorbar(x, con_means[:,0], con_stds[:,0], 
             fmt = '.', capsize = 2, color = "black",
            zorder = 10)
plt.ylabel("Number of Components")
plt.xlabel(r"Horizon Area $[\ell^2]$")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_number.png")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_number.pdf")


plt.figure("components' mean size per Mass", tight_layout = True)
plt.errorbar(x, con_means[:,1], con_stds[:,1], 
             fmt = '.', capsize = 2, color = "black",
            zorder = 10)
plt.ylabel("Mean of Components' Size")
plt.xlabel(r"Horizon Area $[\ell^2]$")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_size_mean.png")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_size_mean.pdf")


plt.figure("components' sdev on size per Mass", tight_layout = True)
plt.errorbar(x, con_means[:,2], con_stds[:,2], 
             fmt = '.', capsize = 2, color = "black",
            zorder = 10)
plt.ylabel("Std of Components' Size")
plt.xlabel(r"Horizon Area $[\ell^2]$")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_size_sdev.png")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_size_sdev.pdf")


plt.figure("components' skew on size per Mass", tight_layout = True)
plt.errorbar(x, con_means[:,3], con_stds[:,3], 
             fmt = '.', capsize = 2, color = "black",
            zorder = 10)
plt.ylabel("Skew of Components' Size")
plt.xlabel(r"Horizon Area $[\ell^2]$")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_size_skew.png")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_size_skew.pdf")


plt.figure("components' kurtosis on size per Mass", tight_layout = True)
plt.errorbar(x, con_means[:,4], con_stds[:,4], 
             fmt = '.', capsize = 2, color = "black",
            zorder = 10)
plt.ylabel("Skew of Components' Size")
plt.xlabel(r"Horizon Area $[\ell^2]$")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_size_thekurtosis.png")
plt.savefig(plotsDir + f"{fixed_string}_HRV_components_size_thekurtosis.pdf")
