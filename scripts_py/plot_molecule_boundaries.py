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
plot_molecules = True


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

# ensure that, if using Mor Rho, it is decimal with exactly 2 dec digits

if fixed_var == "M" or fixed_var == "Rho":
    fixed_string = f"{fixed_var}=" + format(fixed_val,".2f")
else:
    fixed_string = f"{fixed_var}={fixed_val}"


##############################################################################
# 2. GET INFO FROM FILES
##############################################################################

# an entry for each rho
shape_r_max = []
shape_r_min = []
varying_values = []
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
                            #print(varying_var+"=", varying_var_i)
                            varying_values.append(varying_var_i)
                    elif "r=" in piece_j and "dur=" not in piece_j:
                        shape_r_max.append(float(piece_j.split("=")[1]))
                        shape_r_min.append(float(piece_j.split("=")[1])-2)
            # if the file is correct
            if go_on: 
                # 2.1 get file and Nreps #####################################
                try:
                    everything_i = np.loadtxt(file_i_path, 
                                            delimiter = ",", skiprows=1,
                                            dtype = str)[:,:-1].astype(float)
                except IndexError: # only one line in txt file
                    everything_i = np.loadtxt(file_i_path, 
                                            delimiter = ",", skiprows=1,
                                            dtype = str)[:-1].astype(float)
                    everything_i=np.array([everything_i]) #make it 2D for latr
                Nreps_i = everything_i[:,0]
                # 2.2 get mintime, innermost, outermost ####################
                outermosts_avg_arr_i = everything_i[:,1]
                outermosts_std_arr_i = everything_i[:,2]
                innermosts_avg_arr_i = everything_i[:,3]
                innermosts_std_arr_i = everything_i[:,4]
                mintimes_avg_arr_i   = everything_i[:,5]
                mintimes_std_arr_i   = everything_i[:,6]
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
#%% 3. PLOT
##############################################################################
x = np.array(varying_values)
if varying_var=="M":
    r_S = 2*x            #==Radius
    x = 5*np.pi*(r_S)**2 #==Area
    if fixed_var == "Nmult":    
        Rho = fixed_val/(4/5*np.pi*26)
        print("You're doing Nmult fixed!") 
    elif fixed_var == "Rho":
        Rho = fixed_val
    x /= Rho**(-1/2)
    r_S_norm = r_S * Rho**(1/4)

# 3.1 MOLECULES'S BOUNDARIES ##############################################
if plot_boundaries:
    rc = 3; c = 1; rows = 2
    figsize = (c*8, 7.5)
    plt.figure(f'Boundaries for {fixed_string}',
                figsize = figsize, tight_layout = True)

    # INNERMOST #######################################################
    ax = plt.subplot(rows, c, 1)
    plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, innermosts, innermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Innermost (with 1$\sigma$)")
    plt.errorbar(r_S_norm, innermosts, 5*np.array(innermosts_std), 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Innermost (with 5$\sigma$)")
    plt.plot(r_S_norm, shape_r_min, ls = "--", label = "Spacetime Inner Radius")
                
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=0.2)
    ax.text(0.95, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel("Distance from Singularity")
    plt.legend()
    plt.grid(alpha = 0.4) 

    # OUTERMOST #######################################################
    ax = plt.subplot(rows, c, 2)
    plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, outermosts, outermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Outermost (with 1$\sigma$)")
    plt.errorbar(r_S_norm, outermosts, 5*np.array(outermosts_std), 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Outermost (with 5$\sigma$)")
    plt.plot(r_S_norm, shape_r_max, ls = "--", label = "Spacetime Outer Radius")
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=0.2)
    ax.text(0.95, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel("Distance from Singularity")
    plt.legend()
    plt.grid(alpha = 0.4) 

    plt.savefig(plotsDir + f"{fixed_string}_Boundaries_Compared_Limits.png")
    plt.savefig(plotsDir + f"{fixed_string}_Boundaries_Compared_Limits.pdf")
    plt.show()




# 3.2 MOLECULES'S BOUNDARIES  from r_S ##################################
if plot_boundaries:
    rc = 3; c = 1; rows = 3
    figsize = (c*8, 6.5)
    plt.figure(f'Boundaries for {fixed_string}',
                figsize = figsize, tight_layout = True)
    

    # MINTIME #######################################################
    ax = plt.subplot(rows, c, 1)
    plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, mintimes, mintimes_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Mintime (with 1$\sigma$)")
    plt.errorbar(r_S_norm, mintimes, 5*np.array(mintimes_std), 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Mintime (with 5$\sigma$)")
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=0.2)
    ax.text(0.80, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'left', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel("Oldest Molecule's Time")
    plt.grid(alpha = 0.4) 

    # INNERMOST #######################################################
    ax = plt.subplot(rows, c, 2)
    plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, innermosts-r_S, innermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Innermost (with 1$\sigma$)")
    plt.errorbar(r_S_norm, innermosts-r_S, 5*np.array(innermosts_std), 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Innermost (with 5$\sigma$)")
    plt.plot(r_S_norm, shape_r_min-r_S, ls = "--", label = "Spacetime Inner Radius")
                
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=0.2)
    ax.text(0.95, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel(r"Distance from $r_S$")
    plt.legend()
    plt.grid(alpha = 0.4) 

    # OUTERMOST #######################################################
    ax = plt.subplot(rows, c, 3)
    plt.annotate ("c)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, outermosts-r_S, outermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Outermost (with 1$\sigma$)")
    plt.errorbar(r_S_norm, outermosts-r_S, 5*np.array(outermosts_std), 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Outermost (with 5$\sigma$)")
    plt.plot(r_S_norm, shape_r_max-r_S, ls = "--", label = "Spacetime Outer Radius")
    ax.text(0.95, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel(r"Distance from $r_S$")
    plt.legend()
    plt.grid(alpha = 0.4) 

    plt.savefig(plotsDir + f"{fixed_string}_Boundaries_from_rS.png")
    plt.savefig(plotsDir + f"{fixed_string}_Boundaries_from_rS.pdf")
    plt.show()



# 3.3 MOLECULES'S BOUNDARIES in l units ################################
mintimes       = np.array(mintimes)      * Rho**(1/4)
mintimes_std   = np.array(mintimes_std)  * Rho**(1/4)
innermosts     = np.array(innermosts)    * Rho**(1/4)
innermosts_std = np.array(innermosts_std)* Rho**(1/4)
outermosts     = np.array(outermosts)    * Rho**(1/4)
outermosts_std = np.array(outermosts_std)* Rho**(1/4)
shape_r_min        = np.array(shape_r_min) * Rho**(1/4)
shape_r_max        = np.array(shape_r_max) * Rho**(1/4)
if plot_boundaries:
    rc = 3; c = 1; rows = 2
    figsize = (c*8, 6.5)
    plt.figure(f'Boundaries for {fixed_string}',
                figsize = figsize, tight_layout = True)
        

    # INNERMOST #######################################################
    ax = plt.subplot(rows, c, 1)
    plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, innermosts, innermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Innermost (with 1$\sigma$)")
    plt.errorbar(r_S_norm, innermosts, 5*innermosts_std, 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Innermost (with 5$\sigma$)")
    #plt.plot(r_S_norm, shape_r_min, ls = "--", label = "Spacetime Inner Radius")
                
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=0.2)
    ax.text(0.95, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel(r"Distance from Singularity $[\ell]$")
    plt.legend()
    plt.grid(alpha = 0.4) 

    # OUTERMOST #######################################################
    ax = plt.subplot(rows, c, 2)
    plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, outermosts, outermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Outermost (with 1$\sigma$)")
    plt.errorbar(r_S_norm, outermosts, 5*outermosts_std, 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Outermost (with 5$\sigma$)")
    #plt.plot(r_S_norm, shape_r_max, ls = "--", label = "Spacetime Outer Radius")
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=0.2)
    ax.text(0.95, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel(r"Distance from Singularity $[\ell]$")
    plt.legend()
    plt.grid(alpha = 0.4) 

    plt.savefig(plotsDir + f"{fixed_string}_Boundaries_Compared_Limits_in_l.png")
    plt.savefig(plotsDir + f"{fixed_string}_Boundaries_Compared_Limits_in_l.pdf")
    plt.show()




# 3.4 MOLECULES'S BOUNDARIES  from r_S in l units #######################
if plot_boundaries:
    rc = 3; c = 1; rows = 3
    figsize = (c*8, 6.5)
    plt.figure(f'Boundaries for {fixed_string}',
                figsize = figsize, tight_layout = True)

    # MINTIME #######################################################
    ax = plt.subplot(rows, c, 1)
    plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, mintimes, mintimes_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Mintime (with 1$\sigma$)")
    plt.errorbar(r_S_norm, mintimes, 5*np.array(mintimes_std), 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Mintime (with 5$\sigma$)")
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=0.2)
    ax.text(0.80, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'left', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel(r"Oldest Molecule's Time $[\ell]$")
    plt.grid(alpha = 0.4) 

    # INNERMOST #######################################################
    ax = plt.subplot(rows, c, 2)
    plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, innermosts-r_S_norm, innermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Innermost (with 1$\sigma$)")
    plt.errorbar(r_S_norm, innermosts-r_S_norm, 5*innermosts_std, 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Innermost (with 5$\sigma$)")
    # plt.plot(r_S_norm, shape_r_min-r_S_norm, ls = "--", 
    #         label = "Spacetime Inner Radius")
                
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=0.2)
    ax.text(0.95, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel(r"Distance from $r_S$ $[\ell]$")
    plt.legend()
    plt.grid(alpha = 0.4) 

    # OUTERMOST #######################################################
    ax = plt.subplot(rows, c, 3)
    plt.annotate ("c)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    plt.errorbar(r_S_norm, outermosts-r_S_norm, outermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"Outermost (with 1$\sigma$)")
    plt.errorbar(r_S_norm, outermosts-r_S_norm, 5*np.array(outermosts_std), 
                fmt = '.', capsize = 4,
                zorder = 5, label = r"Outermost (with 5$\sigma$)")
    # plt.plot(r_S_norm, shape_r_max-r_S_norm, ls = "--", 
    #         label = "Spacetime Outer Radius")
    ax.text(0.95, 0.05, fixed_string, transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S \; [\ell]$') #not yet in terms of l^2)
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel(r"Distance from $r_S$ $[\ell]$")
    plt.legend()
    plt.grid(alpha = 0.4) 

    plt.savefig(plotsDir + f"{fixed_string}_Boundaries_from_rS_in_l.png")
    plt.savefig(plotsDir + f"{fixed_string}_Boundaries_from_rS_in_l.pdf")
    plt.show()