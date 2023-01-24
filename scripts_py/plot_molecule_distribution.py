import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import expanduser
from causets_py import causet_helpers as ch


##############################################################################
# 0. SET USER VARIABLES
##############################################################################
molecules = "lambdas" #lambdas, HRVs
varying_var = "M"     #variable varying: can be M, Rho, Nmult
fixed_var = "Nmult"   #variable fixed: can be, M, Rho, Nmult
fixed_val = 40000     #value of fixed_var

plot_boundaries = True
plot_molecules = False


##############################################################################
# 1 SET USEFUL STRINGS & OTHER VARIABLES
##############################################################################
home = expanduser("~")
plotsDir = f"{home}/MSci_Schwarzschild_Causets/figures/N{molecules}_vs_Area/"
dataDir = f"{home}/MSci_Schwarzschild_Causets/data/{molecules}/"
#dataDir = home + f"/MSci_Schwarzschild_Causets/data/linkcounting_files/Poiss=False/"
if not os.path.exists(dataDir):
    path = os.getcwd()
    plotsDir = path + f"/figures/N{molecules}_vs_Area/"
    dataDir = path + f"/data/{molecules}"


fixed_string = f"{fixed_var}={fixed_val}"
# ensure that, if using Mor Rho, it is decimal with exactly 2 dec digits
if fixed_var == "M" or fixed_var == "Rho":
    if "." not in fixed_string:
            fixed_string += ".00"
    elif fixed_string.index(".") != len(fixed_string)-3:
        while fixed_string.index(".") != len(fixed_string)-3:
            fixed_string.append("0")
            



##############################################################################
# 2. GET INFO FROM FILES
##############################################################################

# an entry for each rho
varying_values = []
Nreps = []
outermosts = []
outermosts_std = []
innermosts = []
innermosts_std = []
mintimes = []
mintimes_std = []
# matrix: a row for each type of molecule 
#           -> molecules_distr[i] = array of occurrences of ith mol for each rho
#         a column for each rho
molecules_distr = []     #each entry is a row of entries
molecules_distr_std = []


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
                        if varying_var_i not in varying_values:
                            print(varying_var+"=", varying_var_i)
                            varying_values.append(varying_var_i)
                        break
            elif want_Nmult and "Nmult" in file_i:
                go_on = True
                file_i_path = os.path.join(root, file_i)
                file_i_pieces = file_i.split("_")
                # get nmult of file_i
                for piece_j in file_i_pieces:
                    if piece_j[0:3] == "Nmult":
                        nmult_i = float(piece_j.split("=")[1])
                        if nmult_i not in nmults:
                            nmults.append(nmult_i)
                        break
            # if the file is correct
            if go_on: 
                # 2.1 get Nreps ############################################
                Nreps_i = np.loadtxt(file_i_path, delimiter = ",", skiprows=1,
                                            usecols = 0)
                # 2.2 get outermost, innermost, mintime info ################
                outermosts_avg_arr_i = np.loadtxt(file_i_path, delimiter = ",", 
                                            skiprows=1,
                                            usecols = 1)
                outermosts_std_arr_i = np.loadtxt(file_i_path, delimiter = ",", 
                                            skiprows=1,
                                            usecols = 2)
                innermosts_avg_arr_i = np.loadtxt(file_i_path, delimiter = ",", 
                                            skiprows=1,
                                            usecols = 3)
                innermosts_std_arr_i = np.loadtxt(file_i_path, delimiter = ",", 
                                            skiprows=1,
                                            usecols = 4)
                mintimes_avg_arr_i = np.loadtxt(file_i_path, delimiter = ",", 
                                            skiprows=1,
                                            usecols = 5)
                mintimes_std_arr_i = np.loadtxt(file_i_path, delimiter = ",", 
                                            skiprows=1,
                                            usecols = 6)
                outermost_i, outermost_std_i = ch.combine_meass(Nreps_i,
                                                        outermosts_avg_arr_i,
                                                        outermosts_std_arr_i,)
                innermost_i, innermost_std_i = ch.combine_meass(Nreps_i,
                                                        outermosts_avg_arr_i,
                                                        outermosts_std_arr_i,)
                mintime_i, mintime_std_i = ch.combine_meass(Nreps_i,
                                                        outermosts_avg_arr_i,
                                                        outermosts_std_arr_i,)
                outermosts    .append(outermost_i)
                outermosts_std.append(outermost_std_i)
                innermosts    .append(innermost_i)
                innermosts_std.append(innermost_std_i)
                mintimes      .append(mintime_i)
                mintimes_std  .append(mintime_std_i)
                # 2.3 get molecules info #####################################
                mol_info_i = np.loadtxt(file_i_path, delimiter = ",", 
                                        skiprows=1, dtype = str)[:,7:-1]\
                                        .astype(float)
                ntypes = mol_info_i.shape[1]/2
                if (ntypes - int(ntypes)): 
                    raise Exception(f"Number of cols {ntypes*2} is odd")
                else:
                    ntypes = int(ntypes)
                molecules_distr.append([])
                molecules_distr_std.append([])
                for n in range(ntypes):
                    avgs_mol_n_i = mol_info_i[:,2*n]
                    stds_mol_n_i = mol_info_i[:,2*n+1]
                    n_mol_n_i, std_mol_n_i = ch.combine_meass(Nreps,
                                                            avgs_mol_n_i,
                                                            stds_mol_n_i)
                    # if the nth molecule had already been found
                    try:
                        molecules_distr[n]    .append(n_mol_n_i)
                        molecules_distr_std[n].append(std_mol_n_i)
                    #if the nth molecule was never found
                    except IndexError: 
                        # make zeros for previous (general for rho and nmult)
                        n_prev_rhos = max(len(varying_values), len(nmults)) -1
                        zeros = np.zeros(n_prev_rhos).tolist()
                        molecules_distr    .append(zeros + [n_mol_n_i]) 
                        molecules_distr_std.append(zeros + [std_mol_n_i])
# Clear a bit
# del root, dirs, files
# del file_i_path, file_i_pieces, piece_j
# del i, file_i
# del outermosts_avg_arr_i, outermosts_std_arr_i, outermost_i, outermost_std_i
# del innermosts_avg_arr_i, innermosts_std_arr_i, innermost_i, innermost_std_i
# del mintimes_avg_arr_i, mintimes_std_arr_i, mintime_i, mintime_std_i 
# del mol_info_i, avgs_mol_n_i, stds_mol_n_i, n_mol_n_i, std_mol_n_i, 
# del n, ntypes, zeros


##############################################################################
# 2. PLOT
##############################################################################
print("loop done")
A = 4*np.pi*(2*mass)**2
x = np.array(varying_values) if want_rho else np.array(nmults)
print(x)
x_A = A/np.sqrt(x)

# 2.1 MOLECULES'S BOUNDARIES ##############################################
if plot_boundaries:
    rc = 3; c = 1; r = 3
    figsize = (c * 4, r * 3)
    plt.figure(f'Boundaries for {fixed_string}',
                figsize = figsize, tight_layout = True)
    # MINTIME #######################################################
    ax = plt.subplot(r, c, 1)
    plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(x, mintimes, mintimes_std, 
                fmt = '.', ls = 'dashed', capsize = 4, 
                zorder = 5)
    props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
                ls = '', alpha=0.5)
    ax.text(0.05, 0.95, fixed_string, transform=ax.transAxes, fontsize=14, 
            va='bottom', ha = 'left', bbox=props)
    plt.xlabel('Rho [a.u.]')
    plt.ylabel("Oldest Molecule's Time")
    plt.grid(alpha = 0.4) 

    # INNERMOST #######################################################
    ax = plt.subplot(r, c, 2)
    plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(x, innermosts, innermosts_std, 
                fmt = '.', ls = 'dashed', capsize = 4, 
                zorder = 5)
    props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
                ls = '', alpha=0.5)
    ax.text(0.05, 0.95, fixed_string, transform=ax.transAxes, fontsize=14, 
            va='bottom', ha = 'left', bbox=props)
    plt.xlabel('Rho [a.u.]')
    plt.ylabel("Innermost Molecule's r")
    plt.grid(alpha = 0.4) 

    # OUTERMOST #######################################################
    ax = plt.subplot(r, c, 3)
    plt.annotate ("c)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(x, outermosts, outermosts_std, 
                fmt = '.', ls = 'dashed', capsize = 4, 
                zorder = 5)
    props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
                ls = '', alpha=0.5)
    ax.text(0.05, 0.95, fixed_string, transform=ax.transAxes, fontsize=14, 
            va='bottom', ha = 'left', bbox=props)
    plt.xlabel('Rho [a.u.]')
    plt.ylabel("Outermost Molecule's r")
    plt.grid(alpha = 0.4) 

    plt.savefig(plotsDir + f"{fixed_string}_Boundaries.png")
    plt.show()


# 2.1 MOLECULES'S DISTRIBUTION ##############################################
if plot_molecules:
    fig, ax = plt.figure(f"Molecules for {fixed_string}")
    for n in range(len(molecules_distr)):
        y    = molecules_distr[n]
        yerr = molecules_distr_std[n]
        prefix = ("Open " if n==0 else "Close ") if molecules == "HRVs" \
                else f"{n}-"
        plt.errorbar(x_A, y, yerr, 
                    fmt = '.', ls = 'dashed', capsize = 4, 
                    label = f'{prefix}{molecules}')
    props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
                ls = '', alpha=0.5)
    ax.text(0.05, 0.95, fixed_string, transform=ax.transAxes, fontsize=14, 
            va='bottom', ha = 'left', bbox=props)
    plt.legend()
    plt.xlabel('Horizon Area [\ell^2]')
    plt.ylabel("Occurrences")
    plt.grid(alpha = 0.2) 




                    
            
           