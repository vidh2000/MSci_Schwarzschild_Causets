import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import expanduser
from causets_py import causet_helpers as ch
from scipy.optimize import curve_fit


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
    return a*x


x = np.array(varying_values)
if varying_var=="M":
    x = 4*np.pi*(2*x)**2 #==Area
    if fixed_var == "Nmult":    
        Rho = fixed_val/(4/3*np.pi*26)
        print("You're doing Nmult fixed!") 
    elif fixed_var == "Rho":
        Rho = fixed_val
    x /= Rho**(-1/2)
    
    # x_A = x/divide by fund units

print(f"Rho = {Rho:.0f}")
fixed_string = rf"Rho = {Rho:.0f}"

# 2.1 MOLECULES'S BOUNDARIES ##############################################
if plot_boundaries:
    rc = 3; c = 1; r = 3
    figsize = (c * 6, r * 2.3)
    plt.figure(f'Boundaries for {fixed_string}',
                figsize = figsize, tight_layout = True)
    # MINTIME #######################################################
    ax = plt.subplot(r, c, 1)
    plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    if varying_var == "M":
        plt.errorbar(2*np.array(varying_values), mintimes,
                mintimes_std, fmt = '.', capsize = 4, 
                zorder = 5)
    else:
        plt.errorbar(np.array(varying_values), mintimes, mintimes_std, 
                fmt = '.', capsize = 4, 
                zorder = 5)
    props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
                ls = '', alpha=0.2)
    ax.text(0.80, 0.05, rf"$\rho \: = \: {Rho:.0f}$", transform=ax.transAxes, fontsize=10,  #text==fixed_string
            va='bottom', ha = 'left', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S$ [a.u.]')
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel("Oldest Molecule's Time")
    plt.grid(alpha = 0.4) 

    # INNERMOST #######################################################
    ax = plt.subplot(r, c, 2)
    plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    if varying_var == "M":
        plt.errorbar(2*np.array(varying_values), innermosts,
                innermosts_std, fmt = '.', capsize = 4, 
                zorder = 5)
    else:
        plt.errorbar(np.array(varying_values), innermosts, innermosts_std, 
                fmt = '.', capsize = 4, 
                zorder = 5)

    props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
                ls = '', alpha=0.2)
    ax.text(0.95, 0.05, rf"$\rho \: = \: {Rho:.0f}$", transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S$ [a.u.]')
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel("Innermost Molecule's r")
    plt.grid(alpha = 0.4) 

    # OUTERMOST #######################################################
    ax = plt.subplot(r, c, 3)
    plt.annotate ("c)", (-0.05, 1.05), xycoords = "axes fraction", 
                    va='bottom', ha = 'left')
    if varying_var == "M":
        plt.errorbar(2*np.array(varying_values), outermosts,
                outermosts_std, fmt = '.', capsize = 4, 
                zorder = 5)
    else:
        plt.errorbar(np.array(varying_values), outermosts, outermosts_std, 
                fmt = '.', capsize = 4, 
                zorder = 5)
    props = dict(boxstyle='round', facecolor='wheat', edgecolor = 'grey', 
                ls = '', alpha=0.2)
    ax.text(0.95, 0.05, rf"$\rho \: = \: {Rho:.0f}$", transform=ax.transAxes, fontsize=10, 
            va='bottom', ha = 'right', bbox=props)
    if varying_var == "M":
        plt.xlabel(r'$r_S$ [a.u.]')
    else:
        plt.xlabel(f'{varying_var} [a.u.]')
    plt.ylabel("Outermost Molecule's r")
    plt.grid(alpha = 0.4) 

    plt.savefig(plotsDir + f"{fixed_string}_Boundaries.png")
    plt.savefig(plotsDir + f"{fixed_string}_Boundaries.pdf")
    plt.show()

# 2.1 MOLECULES'S DISTRIBUTION ##############################################
print("")
if plot_molecules:

    plt.figure(f"Molecules for {fixed_string}")

    gradients = []
    gradients_unc = []
    for n in range(len(molecules_distr)):
        y    = molecules_distr[n]
        yerr = molecules_distr_std[n]
        label = ("Open HRV" if n==0 else "Closed HRV ") if molecules == "HRVs" \
                else str(n+1) + r"$\mathbf{-}\Lambda$"
        plt.errorbar(x, y, yerr, 
                    fmt = '.', capsize = 4, 
                    label = label)
        # Linear fit
        # Fitting function == linear through vertex

        # fit the function to the data
        popt, pcov = curve_fit(lin_func, x, y, sigma=yerr,
                                absolute_sigma=True)
        unc = np.sqrt(np.diag(pcov))
        print(f"Gradient factor for {n+1}-lambda = {round(popt[0],3)} +- {round(unc[0],3)}")
        
        xfit = np.linspace(min(x),max(x),100)
        #plt.plot(xfit, lin_func(xfit,*popt), '--', color="red")
        gradients.append(popt[0])
        gradients_unc.append(unc[0])
    
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=1)
    plt.annotate(rf"$\rho \: = \: {Rho:.0f}$", (0.95, 0.5), xycoords = "axes fraction",
            fontsize=12, va='center', ha = 'right', bbox=props)
    plt.legend()
    plt.xlabel(r'Horizon Area $[\ell^2]$') #not yet in terms of l^2
    plt.ylabel(r"Number of $n\mathbf{-}\Lambda$")
    plt.grid(alpha = 0.2)
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}.png") 
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}.pdf") 
    plt.show()


    plt.figure(f"Large molecules for {fixed_string}")
    for n in range(3, len(molecules_distr)):
        y    = molecules_distr[n]
        yerr = molecules_distr_std[n]
        label = ("Open HRV" if n==0 else "Closed HRV ") if molecules == "HRVs" \
                else str(n+1) + r"$\mathbf{-}\Lambda$"
        plt.errorbar(x, y, yerr, 
                    fmt = '.', capsize = 4, 
                    label = label)
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=1)
    plt.annotate(rf"$\rho \: = \: {Rho:.0f}$", (0.95, 0.5), xycoords = "axes fraction",
            fontsize=12, va='center', ha = 'right', bbox=props)
    plt.legend()
    plt.xlabel(r'Horizon Area $[\ell^2]$') #not yet in terms of l^2
    plt.ylabel(r"Number of $n\mathbf{-}\Lambda$")
    plt.grid(alpha = 0.2)
    plt.savefig(plotsDir + f"{fixed_string}_large_{molecules}.png")
    plt.savefig(plotsDir + f"{fixed_string}_large_{molecules}.pdf") 
    plt.show()

    ### Find entropy and C_hv and discretness length scale

    noninf_indices = tuple([np.where(np.isfinite(gradients_unc))[0]])
    gradients = np.array(gradients)[noninf_indices]
    gradients_unc = np.array(gradients_unc)[noninf_indices]



    # Uncertainty propagation
    coefsum = sum([(n+1)*gradients[n] for n in range(len(gradients))])
    coefsum_unc = np.sqrt(
        sum([((n+1)*gradients_unc[n])**2 for n in range(len(gradients))]))
    print(f"\nGradients weighted sum = {round(coefsum,5)} +- {round(coefsum_unc,5)}")

    # Find probability distribution of n-lambdas
    grad_sum = sum(gradients)
    grad_sum_unc = np.sqrt(sum([g**2 for g in gradients_unc]))

    lambd_probs = gradients/sum(gradients)
    lambd_probs_uncs = np.sqrt(gradients_unc**2/grad_sum**2 +
                               gradients**2*grad_sum_unc**2/grad_sum**4)
    

    plt.figure("n-lambda probability distribution")
    plt.bar(np.arange(1,len(lambd_probs)+1,1), lambd_probs,
            label = r"$n\mathbf{-}\Lambda$ probability distribution")
    plt.errorbar(np.arange(1,len(lambd_probs)+1,1), lambd_probs,
            yerr=lambd_probs_uncs,capsize=7,fmt="",ls="",ecolor="red")
    plt.xlabel(r"$n$")
    plt.ylabel("Probability")
    plt.legend(loc="upper right")
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=1)
    #plt.annotate(rf"$\rho \: = \: {Rho:.0f}$", (0.95, 0.95), xycoords = "axes fraction",
    #        fontsize=12, va='top', ha = 'right', bbox=props)
    plt.grid(alpha=0.2)
    plt.savefig(plotsDir + "n_lambda_probability_distribution.png")
    plt.savefig(plotsDir + "n_lambda_probability_distribution.pdf")
    plt.show()
    
    print("Distribution of n-lambdas:\n",
            [round(l,7) for l in lambd_probs])

    # Find discreteness scale
    P = sum([p*np.log(p) for p in lambd_probs])
    P_unc = np.sqrt(sum((lambd_probs+1)**2 * lambd_probs_uncs**2))
    C_hv = - sum(gradients)*sum([p*np.log(p) for p in lambd_probs])
    C_hv_unc = np.sqrt(P**2 * grad_sum_unc**2 +
                       grad_sum**2 * P_unc**2)     
    # Discreteness length in terms of Planckian length
    l = 2*np.sqrt(C_hv)
    l_unc = C_hv_unc/np.sqrt(C_hv)
    print(f"\nDiscreteness scale = {round(l,5)} +- {round(l_unc,5)} l_p")
 