import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os
from os.path import expanduser
from scipy.sparse.csgraph import connected_components

###############################################################################
# Open files HRVs_connectivity/M=<>_Rho=<>_Card=<>_r=<>_hollow=<>_dur=<>_n.txt
# each containing the HRVs molecules from the nth run of a causet with 
# specified parameters.
#
# Then save a file with columns being # clusters of HRVs with certain # of points
# and associated standard deviation (like with lambdas), of title
# HRV_connectivity_compress/M=<>_Rho=<>_Card=<>_r=<>_hollow=<>_dur=<>.txt
###############################################################################

a = [1,1,1,2]
print(np.bincount(a))
print(np.bincount(a)[1:])


def get_info_of_HRVs_for_connectivity(file_ext):
    """

    Parameters
    ----------
    - file_ext : str
        File from which import info

    Returns
    -------
    - Adj_matrix : np.ndarray(int, int)
        Adj_ij = 1 if i < j or j < i
        Irreflexive, i.e 0 diagonal.
    """
    with open(file_ext, 'r') as fl: 
            f = fl.readlines() 
            mols = []

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
                    go -= 1

                elif "r_S" in key:
                    r_S = float(row[1])
                    break
                elif "Coordinates" in key:
                    break
                else: #the coordinates#
                    break
            # finished reading file

            # make adjacency matrix of connections
            N = np.amax(mols) +1
            Adj_matrix = np.zeros((N,N))
            for mol in mols:
                Adj_matrix[mol[0], mol[1]] = 1
                Adj_matrix[mol[1], mol[0]] = 1
                Adj_matrix[mol[0], mol[2]] = 1
                Adj_matrix[mol[2], mol[0]] = 1

    return Adj_matrix
    

##############################################################################
# 0. SET USER VARIABLES
##############################################################################
molecules = "HRVs" #lambdas, HRVs
varying_var = "M"     #variable varying: can be M, Rho, Nmult (M best)
fixed_var = "Rho"     #variable fixed:   can be, M, Rho, Nmult (Rho best)
fixed_val = 5000     #value of fixed_var 
use_selected_masses = True #gives equal spacing


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
newDataDir = dataDir[:-1] + "_compress/"
print(dataDir)
# ensure that, if using Mor Rho, it is decimal with exactly 2 dec digits

if fixed_var == "M" or fixed_var == "Rho":
    fixed_string = f"{fixed_var}=" + format(fixed_val,".2f")
else:
    fixed_string = f"{fixed_var}={fixed_val}"


selected_masses = np.array([
                            0.53 ,0.65, 0.75, 0.84, 0.92, 0.99, 1.06, 
                            1.13, 1.19, 1.24, 1.30, 1.35, 1.40, 1.45, 
                            1.50, 1.55, 1.59, 1.63, 1.68, 1.72, 1.76, 
                            1.80, 1.84, 1.88, 1.91, 1.95, 1.98, 2.02, 
                            2.05, 2.09, 2.12, 2.15, 2.19, 2.22, 2.25, 
                            2.28, 2.31, 2.34, 2.37, 2.40, 
                            2.43, 2.46, 2.49,
                            ])



##############################################################################
# 2. GET INFO FROM FILES
##############################################################################

# an entry for each M
file_names = [] #store name of files (one per value of varying_value)
Nreps = []
varying_values = []
sizes_i_n = [] # the ith array corresponds to  a certain varying_value and is
               # and is an Nreps[i]-long array of arrays of sizes of the HRVs
               # clusters in the nth repetition


for root, dirs, files in os.walk(dataDir):
    # for each file file_i
    for i, file_i in enumerate(files):
        if fixed_string in file_i:
            go_on = False
            if varying_var+"=" in file_i:
                go_on = True
                file_i_path = os.path.join(root, file_i)
                file_i_pieces = file_i.split("_")
                # get the value of the varying variable (M) of file_i
                for piece_j in file_i_pieces:
                    if varying_var+"=" in piece_j:
                        varying_var_i = float(piece_j.split("=")[1])
                        if use_selected_masses \
                            and varying_var_i not in selected_masses:
                            go_on = False
                            break #j loop
                        if varying_var_i not in varying_values:
                            print(varying_var+" =", varying_var_i)
                            varying_values.append(varying_var_i)
                            Nreps.append(0)
                            sizes_i_n.append([])
                            # Set name of file in which saving compress info
                            # of HRVs of certain mass
                            filename_to_use = ""
                            for piece_i in file_i_pieces[:-1]:
                                filename_to_use += piece_i + "_"
                            filename_to_use = filename_to_use[:-1]+"txt"
                            file_names.append(filename_to_use)
                        break #piece_j loop 20ish lines above

            # if the file is correct
            if go_on: 
                # 2.1 get # of components and info on dstribution of sizes of
                # components
                Aij = get_info_of_HRVs_for_connectivity(file_i_path)
                #1. Get the number of clusters and label elements by custer
                n_components, labels = connected_components(Aij, directed = False)
                #2. Bin the labels to get the size of each ith cluster 
                #Each entry of size_clusters is the size of the ith cluster
                size_clusters = np.bincount(labels)
                #3. Bin again to get the number of clusters of a certain size
                # excluding size 0, 1 and 2 (which are not possible as size of
                # HRV is 3).
                number_per_size = np.bincount(size_clusters)[3:]

                # update
                index = varying_values.index(varying_var_i)
                Nreps[index] += 1
                sizes_i_n[index].append(number_per_size)


#############################################################
# INFO ON NREPS
varying_values_copy = np.array(varying_values)
combined = list(zip(varying_values_copy, Nreps)) 
sorted_combined = sorted(combined, key=lambda x: x[0]) # sort by values
vals = [x[0] for x in sorted_combined] 
Nreps = [x[1] for x in sorted_combined]
maxrep = max(Nreps)

vals_not_max = [vals [i] for i in range(len(vals)) if Nreps[i] != maxrep]
reps_not_max = [Nreps[i] for i in range(len(vals)) if Nreps[i] != maxrep]
repstable = pd.DataFrame(
                np.column_stack(
                [vals_not_max, reps_not_max, maxrep-np.array(reps_not_max)]),
                columns = ["M Value", "Current Nreps", "Reps to Aim"])
print(f"Maximum number of Nreps is {maxrep}")
print("The following values have a number of Nreps different than the desired")
print(repstable)


#############################################################
# SAVE INFO ON SIZES OF CLUSTERS OF HRVS 
print("\nStarting the saving file procedure")
for i in range(len(varying_values)):
    # sizes_sum is an array, whose ith element is the sum of the # of 
    # HRV clusters of size i among all repetitions. Dividing by Nreps[i] gives
    # avg number of said cluster per repetitions.
    # sizes2_sum is the sum of the squares.
    L = len(sizes_i_n[i][0])
    sizes_sum = sizes_i_n[i][0]*1.
    sizes2_sum= sizes_i_n[i][0]*sizes_i_n[i][0]
    for n in range(1,Nreps[i]):
        try:
            sizes_in = sizes_i_n[i][n]
        except IndexError:
            print(f"We are at rep {n}, within varying value {i}")
            print(f"Nreps[i] is {Nreps[i]}")
            print(f"Length size_i_n[i] is {len(sizes_i_n[i])}")
            sizes_in = sizes_i_n[i][n]
        len_diff = len(sizes_in) - L
        if len_diff > 0:
            temp_list = sizes_sum.tolist()
            temp_list2 = sizes2_sum.tolist()
            for ld in range(len_diff):
                temp_list.append(0)              
                temp_list2.append(0)
                L += 1
            sizes_sum = np.array(temp_list)
            sizes2_sum = np.array(temp_list)
        elif len_diff < 0:
            temp_list = sizes_in.tolist()
            for ld in range(abs(len_diff)):
                temp_list.append(0)
            sizes_in = np.array(temp_list)
        sizes_sum += sizes_in
        sizes2_sum+= sizes_in*sizes_in
    mean_number_per_size = sizes_sum/Nreps[i]
    std_number_per_size = np.sqrt(sizes2_sum/Nreps[i]
                                  - mean_number_per_size*mean_number_per_size)
    table_dict = {"Nreps" : Nreps[i]}
    for s in range(len(sizes_sum)):
        table_dict[f"{s+3}avg"] = round(mean_number_per_size[s],3)
        table_dict[f"{s+3}std"] = round(std_number_per_size[s],3)
    ith_file = file_names[i][:-6] + ".txt"
    ith_table = pd.DataFrame(table_dict, index = [0])
    ith_table.to_csv(newDataDir+ith_file, sep = ",", 
                     header = True, index = False)
    print(f"Saved {newDataDir+ith_file}")