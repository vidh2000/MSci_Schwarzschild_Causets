import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import expanduser
from causets_py import causet_helpers as ch
from scipy.optimize import curve_fit

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
# 0. SET USER VARIABLES
##############################################################################
molecules = "HRVs" #lambdas, HRVs
varying_var = "M"     #variable varying: can be M, Rho, Nmult
fixed_var = "Rho"   #variable fixed: can be, M, Rho, Nmult
fixed_val = 5000     #value of fixed_var

plot_histogram_Nreps = True
plot_boundaries = 1
plot_molecules = 1


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

# an entry for each M
Nreps = []
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
                            print(varying_var+"=", varying_var_i)
                            varying_values.append(varying_var_i)
                        break
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
                tot_Nreps_i = sum(Nreps_i)
                Nreps.append(tot_Nreps_i)
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
                # 2.3 get molecules info #####################################
                mol_info_i = everything_i[:,7:]
                ntypes_i = mol_info_i.shape[1]/2
                # should be even
                if (ntypes_i - int(ntypes_i)): 
                    raise Exception(f"Number of cols {ntypes_i*2} is odd")
                else:
                    ntypes_i = int(ntypes_i)
                #molecules_distr.append([])
                #molecules_distr_std.append([])
                for n in range(ntypes_i):
                    avgs_mol_n_i = mol_info_i[:,2*n]
                    stds_mol_n_i = mol_info_i[:,2*n+1]
                    n_mol_n_i, std_mol_n_i = ch.combine_meass(Nreps_i,
                                                            avgs_mol_n_i,
                                                            stds_mol_n_i)
                    # if the nth molecule had already been found
                    try:
                        molecules_distr[n]    .append(n_mol_n_i)
                        molecules_distr_std[n].append(std_mol_n_i)
                    #if the nth molecule was never found
                    except IndexError: 
                        # make zeros for previous (general for rho and nmult)
                        n_prev_values = len(varying_values) -1
                        zeros = np.zeros(n_prev_values).tolist()
                        molecules_distr    .append(zeros + [n_mol_n_i]) 
                        molecules_distr_std.append(zeros + [std_mol_n_i])


n_max = max([len(l) for l in molecules_distr])
for l,ls in zip(molecules_distr, molecules_distr_std):
    while len(l) < n_max:
        l.append(0)
        ls.append(0)

##############################################################################
# 2. PLOT
##############################################################################

# Fitting function == linear through vertex
def lin_func(x,a):
    " x is Area in l units"
    return a*x

# Fitting function2 == corrected-linear through vertex
def corrected_lin_func(x,a, coeff = - 0.0544):
    " x is A in l units"
    M = np.sqrt(np.array(x)/(16*np.pi))
    return (a + coeff/M)*x

# Set the right x for the varying-fixed values -> Area in l^2 units
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
    # x_A = x/divide by fund units

print(f"Rho = {Rho:.0f}")
fixed_string = rf"Rho = {Rho:.0f}"


###########################################################################
# 2.0 Nreps ##############################################
if plot_histogram_Nreps:
    varying_values_copy = np.array(varying_values)
    combined = list(zip(varying_values_copy, Nreps)) 
    sorted_combined = sorted(combined, key=lambda x: x[0]) # sort by values
    vals = [x[0] for x in sorted_combined] 
    Nreps = [x[1] for x in sorted_combined]
    print(vals)
    print(Nreps)
    plt.figure("Histogram Nreps", (20, 6))
    plt.plot(vals, Nreps, "x", markersize = 5, ls ="")
    plt.xticks(vals, fontsize = 5)
    plt.vlines(vals, 0, 200, alpha = 0.1, ls ="--")
    plt.ylabel("Nreps")
    plt.xlabel("M")
    plt.tight_layout()
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}Nreps.png")
    


###########################################################################
# 2.1 MOLECULES'S BOUNDARIES ##############################################
if plot_boundaries:
    mintimes       = np.array(mintimes)      * Rho**(1/4)
    mintimes_std   = np.array(mintimes_std)  * Rho**(1/4)
    innermosts     = np.array(innermosts)    * Rho**(1/4)
    innermosts_std = np.array(innermosts_std)* Rho**(1/4)
    outermosts     = np.array(outermosts)    * Rho**(1/4)
    outermosts_std = np.array(outermosts_std)* Rho**(1/4)

    rc = 2; c = 1; r = 2
    figsize = (6 * c, 14 / r)
    plt.figure(f'Boundaries for {fixed_string}',
                figsize = (10, 6), tight_layout = True)

    # MINTIME #######################################################
    ax = plt.subplot(r, c, 1)
    # plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
                    # va='bottom', ha = 'left')
    plt.errorbar(x, mintimes, mintimes_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"1$\sigma$")
    plt.errorbar(x, mintimes, 3*np.array(mintimes_std), 
                fmt = '.', capsize = 4, color = "C0",
                zorder = 5, label = r"3$\sigma$")
    ax.set_xticklabels([])
    plt.ylabel(r"$t_{\mathrm{min}}$ $[\ell]$")
    plt.legend()
    plt.grid(alpha = 0.4) 

    # INNERMOST & OUTMOST ####################################################
    ax = plt.subplot(r, c, 2)
    # plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction", 
    #                 va='bottom', ha = 'left')
    plt.errorbar(x, innermosts-r_S_norm, innermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"1$\sigma$")
    plt.errorbar(x, innermosts-r_S_norm, 3*innermosts_std, 
                fmt = '.', capsize = 4, color = "C0",
                zorder = 5, label = r"3$\sigma$")
    plt.errorbar(x, outermosts-r_S_norm, outermosts_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10)#, label = r"Outermost (with 1$\sigma$)")
    plt.errorbar(x, outermosts-r_S_norm, 5*np.array(outermosts_std), 
                fmt = '.', capsize = 4, color = "C0",
                zorder = 5)#, label = r"Outermost (with 5$\sigma$)")

    # props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
    #             ls = '', alpha=0.2)
    # ax.text(0.95, 0.05, rf"$\rho \: = \: {Rho:.0f}$", transform=ax.transAxes, 
    # fontsize=10, va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'Horizon Area $[\ell^2]$')
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel(r"$\Delta r_{\mathrm{max}}$ $[\ell]$")
    from matplotlib.ticker import FormatStrFormatter
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.legend()
    plt.grid(alpha = 0.4) 
    plt.savefig(plotsDir + f"{fixed_string}__{molecules}Boundaries_to_rs_in_l.png")
    plt.savefig(plotsDir + f"{fixed_string}__{molecules}Boundaries_to_rs_in_l.pdf")
    



#############################################################################
# 2.1 MOLECULES'S DISTRIBUTION ##############################################
print("")
if plot_molecules:

    plt.figure(f"Molecules for {fixed_string}", figsize=figsize, 
            tight_layout = True)

    gradients = []
    gradients_unc = []
    for n in range(len(molecules_distr)):
        y    = molecules_distr[n]
        yerr = molecules_distr_std[n]
        label = ("Open HRV" if n==0 else "Closed HRV ") if molecules == "HRVs"\
                else str(n+1) + r"$\mathbf{-}\Lambda$"
        plt.errorbar(x, y, yerr, 
                    fmt = '.', capsize = 4, 
                    label = label)

        # print(f"mean of # of {n+1}-Lambda: {[round(yi,8) for yi in y]}")
        # print(f"stds of # of {n+1}-Lambda: {[round(yerri,8) for yerri in yerr]}")
        # n-Lambda Linear fit to line through vertex
        try:
            popt, pcov = curve_fit(lin_func, x, y, sigma=yerr,
                                    absolute_sigma=True)
            unc = np.sqrt(np.diag(pcov))
        except RuntimeError:
            non0_indices = tuple([np.where(np.array(y) != 0)])
            xs2      = np.array(x)[non0_indices][0]
            ys2      = np.array(y)[non0_indices][0]
            if len(ys2) > 0:
                yerrs2   = np.array(yerr)[non0_indices][0]
                popt, pcov = curve_fit(lin_func, xs2, ys2, sigma=yerrs2,
                                        absolute_sigma=True)
                unc = np.sqrt(np.diag(pcov))
            else:
                popt = [0]
                unc = [0]
        print(f"Gradient factor for {label} = {round(popt[0],5)} +- {round(unc[0],5)}")
        
        xfit = np.linspace(min(x),max(x),100)
        #plt.plot(xfit, lin_func(xfit,*popt), '--', color="red")
        gradients.append(popt[0])
        gradients_unc.append(unc[0])
    plt.legend()
    plt.xlabel(r'Horizon Area $[\ell^2]$') #not yet in terms of l^2
    plt.ylabel(f"Number of HRVs")
    plt.grid(alpha = 0.2)
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}.png") 
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}.pdf") 
    


    ##########################################################################
    ##########################################################################
    ### DO THE NUMERICAL ANALYSIS
    ##########################################################################
    ##########################################################################
    noninf_indices = tuple([np.where(np.isfinite(gradients_unc))[0]])
    gradients      = np.array(gradients)[noninf_indices]
    gradients_unc  = np.array(gradients_unc)[noninf_indices]


    ##########################################################################
    ### Find Entropy and C_hv and Discretness Length Scale
    ##########################################################################
    if gradients[1] != 0:
        grad_sum = sum(gradients)
        grad_sum_unc = np.sqrt(sum([g**2 for g in gradients_unc]))

        hrv_probs = gradients/sum(gradients)
        hrv_probs_uncs = [0,0]
        hrv_probs_uncs[0] = 1/grad_sum**2 \
                            * np.sqrt(
                                (gradients[0]*gradients_unc[1])**2\
                               +(gradients[1]*gradients_unc[0])**2
                            )

        hrv_probs_uncs[1] = 1/grad_sum**2 \
                            * np.sqrt(
                                (gradients[0]*gradients_unc[1])**2\
                               +(gradients[1]*gradients_unc[0])**2
                            )
        
        
        print(" \n###PROBABILITIES ###")
        print(f"p_op = {round(hrv_probs[0],4)}+-{round(hrv_probs_uncs[0],4)}")
        print(f"p_cl = {round(hrv_probs[1],4)}+-{round(hrv_probs_uncs[1],4)}")


        print(" \n#### FINAL RESULTS ####")
        # Find overall proportionality and discreteness scale
        C_hv = sum(gradients*np.log(grad_sum/gradients))
        dC_da_i = np.log(grad_sum/gradients)
        C_hv_unc = np.sqrt( sum(gradients_unc**2 * dC_da_i**2)) 

    else:
        print(" AAAAAAAA GRADIENTS[1] == 0!!!!!")
        C_hv = gradients[0] * np.log(2)
        C_hv_unc = gradients_unc[0] * np.log(2)

    print(f"Overall Proportion C_hv = {round(C_hv,5)} +- {round(C_hv_unc,5)}")

    # Discreteness length in terms of Planckian length
    l = 2*np.sqrt(C_hv)
    l_unc = C_hv_unc/np.sqrt(C_hv)
    print(f"Discreteness scale      = {round(l,5)} +- {round(l_unc,5)} l_p")
 

plt.show()