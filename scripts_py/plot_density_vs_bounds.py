import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import expanduser
from causets_py import causet_helpers as ch
import pandas as pd
import re

##############################################################################
# 0. SET USER VARIABLES
##############################################################################
usehome = True
want_rho = True

mass = 2.0

plot_boundaries = True
plot_molecules = True

want_Nmult = not want_rho



##############################################################################
# 1 SET USEFUL STRINGS
##############################################################################
# Home Directory
home = expanduser("~")
# Path
path = os.getcwd()
if usehome:
    plotsDir = f"{home}/MSci_Schwarzschild_Causets/figures/density_vs_bounds/"
    dataDir = f"{home}/MSci_Schwarzschild_Causets/data/test_boundary_vs_density/"
else:
    plotsDir = path + f"/figures/density_vs_bounds/"
    dataDir = path + f"/data/test_boundary_vs_density/"

# mass_string must have 6 characters
mass = round(mass, 2)
mass_string = f"M={round(mass,2)}"




##############################################################################
# 2. GET INFO FROM FILES
##############################################################################

# an entry for each rho
rhos = []
nmults = []
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
        

        print(file_i)
        if f"M={round(mass,2)}" not in file_i:
            continue

        # Store the density for each file in array    
        rhos.append(float(file_i.split("_M=")[0][4:]))

        file_i = root + file_i

        # Get data from the file_i
        with open(file_i, "r") as f:
            lines = list(f)[1:] #skip first row
            data = []
            for i, line in enumerate(lines):
                data.append([float(x) for x in line.split()])
        data = pd.DataFrame(data)

        # 2.1 get Nreps ############################################
        Nreps_i = [int(x) for x in data.iloc[:,0].tolist()] 
        # 2.2 get outermost, innermost, mintime info ################
        outermosts_avg_arr_i = data.iloc[:,1].tolist()
        outermosts_std_arr_i = data.iloc[:,2].tolist()
        innermosts_avg_arr_i = data.iloc[:,3].tolist()
        innermosts_std_arr_i = data.iloc[:,4].tolist()
        mintimes_avg_arr_i = data.iloc[:,5].tolist()
        mintimes_std_arr_i = data.iloc[:,6].tolist()

        print(Nreps_i)
        print(outermosts_avg_arr_i)
        print(outermosts_std_arr_i)

        outermost_i, outermost_std_i = ch.std_combine(Nreps_i,
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
        mol_info_i = np.loadtxt(file_i, delimiter = ",", 
                                skiprows=1) [:,7:]
        ntypes = mol_info_i.shape[1]/2
        if (ntypes - int(ntypes)): 
            raise Exception(f"Number of cols {ntypes*2} is odd")
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
                n_prev_rhos = max(len(rhos), len(nmults)) -1
                zeros = np.zeros(n_prev_rhos).tolist()
                molecules_distr    .append(zeros + [n_mol_n_i]) 
                molecules_distr_std.append(zeros + [std_mol_n_i])


print(rhos)
##############################################################################
# 2. PLOT
##############################################################################
A = 4*np.pi*(2*mass)**2
x = np.array(rhos) if want_rho else np.array(nmults)
x_A = A/np.sqrt(x)

# 2.1 MOLECULES'S BOUNDARIES ##############################################
if plot_boundaries:
    rc = 3; c = 1; r = 3
    figsize = (c * 4, r * 3)
    plt.figure(f'Boundaries for {mass_string}',
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
    ax.text(0.05, 0.95, mass_string, transform=ax.transAxes, fontsize=14, 
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
    ax.text(0.05, 0.95, mass_string, transform=ax.transAxes, fontsize=14, 
            va='bottom', ha = 'left', bbox=props)
    plt.xlabel('Rho [a.u.]')
    plt.ylabel("Innermost Molecule's r")
    plt.grid(alpha = 0.4) 

    # OUTERMOST #######################################################
    ax = plt.subplot(r, c, 2)
    plt.annotate ("c)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(x, outermosts, outermosts_std, 
                fmt = '.', ls = 'dashed', capsize = 4, 
                zorder = 5)
    props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
                ls = '', alpha=0.5)
    ax.text(0.05, 0.95, mass_string, transform=ax.transAxes, fontsize=14, 
            va='bottom', ha = 'left', bbox=props)
    plt.xlabel('Rho [a.u.]')
    plt.ylabel("Outermost Molecule's r")
    plt.grid(alpha = 0.4) 

    plt.savefig(plotsDir + f"{mass_string}_Boundaries")
    plt.show()


# 2.1 MOLECULES'S DISTRIBUTION ##############################################
"""
if plot_molecules:
    fig, ax = plt.figure(f"Molecules for {mass_string}")
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
    ax.text(0.05, 0.95, mass_string, transform=ax.transAxes, fontsize=14, 
            va='bottom', ha = 'left', bbox=props)
    plt.legend()
    plt.xlabel('Horizon Area [\ell^2]')
    plt.ylabel("Occurrences")
    plt.grid(alpha = 0.2) 
"""



                    
            
           