import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import expanduser
from causets_py import causet_helpers as ch
import pandas as pd
import re
from matplotlib.ticker import FuncFormatter,ScalarFormatter

##############################################################################
# 0. SET USER VARIABLES
##############################################################################
usehome = True

mass = 2.0


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
mass_string = "M="+format(mass,".2f")


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

        outermost_i, outermost_std_i = ch.combine_meass(Nreps_i,
                                                outermosts_avg_arr_i,
                                                outermosts_std_arr_i,)
        innermost_i, innermost_std_i = ch.combine_meass(Nreps_i,
                                                innermosts_avg_arr_i,
                                                innermosts_std_arr_i,)
        mintime_i, mintime_std_i = ch.combine_meass(Nreps_i,
                                                mintimes_avg_arr_i,
                                                mintimes_std_arr_i,)
        outermosts    .append(outermost_i)
        outermosts_std.append(outermost_std_i)
        innermosts    .append(innermost_i)
        innermosts_std.append(innermost_std_i)
        mintimes      .append(mintime_i)
        mintimes_std  .append(mintime_std_i)


##############################################################################
# 2. PLOT
##############################################################################


x = np.array(rhos)

# 2.1 MOLECULES'S BOUNDARIES ##############################################

rc = 3; c = 1; r = 3
figsize = (c * 4, r * 3)
plt.figure(f'Boundaries for {mass_string}',
            figsize = figsize, tight_layout = True)
# MINTIME #######################################################
ax = plt.subplot(r, c, 1)
plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
                va='bottom', ha = 'left')
plt.errorbar(x, -np.array(mintimes), mintimes_std, 
            fmt = '.', capsize = 4, #ls = 'dashed', 
            zorder = 5)
props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
            ls = '', alpha=0.5)
ax.text(0.95, 0.95, mass_string, transform=ax.transAxes, fontsize=14, 
        va='top', ha = 'right', bbox=props)
plt.xlabel(r'Number density $\rho$ [a.u.]')
plt.ylabel(r"$-t_{min}$ [a.u]")
plt.grid(alpha = 0.4) 
ax.set_xscale("log")

# INNERMOST #######################################################
ax = plt.subplot(r, c, 2)
plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction",  
                va='bottom', ha = 'left')
plt.errorbar(x, innermosts, innermosts_std, 
            fmt = '.', capsize = 4, 
            zorder = 5)
props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
            ls = '', alpha=0.5)
ax.text(0.95, 0.05, mass_string, transform=ax.transAxes, fontsize=14, 
        va='bottom', ha = 'right', bbox=props)
plt.xlabel(r'Number density $\rho$ [a.u.]')
plt.ylabel(r"$r_{min}$ [a.u]")
plt.grid(alpha = 0.4) 
ax.set_xscale("log")

# OUTERMOST #######################################################
ax = plt.subplot(r, c, 3)
ax.annotate ("c)", (-0.05, 1.05), xycoords = "axes fraction", 
                va='bottom', ha = 'left')
ax.errorbar(x, outermosts, outermosts_std, 
            fmt = '.', capsize = 4, 
            zorder = 5)
props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
            ls = '', alpha=0.5)
ax.text(0.95, 0.95, mass_string, transform=ax.transAxes, fontsize=14, 
        va='top', ha = 'right', bbox=props)
ax.set_xlabel(r'Number density $\rho$ [a.u.]')
ax.set_ylabel(r"$r_{max}$ [a.u]")
ax.grid(alpha = 0.4) 
ax.set_xscale("log")

plt.savefig(plotsDir + f"{mass_string}_Boundaries.png")
plt.show()


                    
            
           