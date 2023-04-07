import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os
from os.path import expanduser
from causets_py import causet_helpers as ch
from scipy.optimize import curve_fit
from scipy.stats import chisquare, pearsonr


##############################################################################
# 0. SET USER VARIABLES
##############################################################################
molecules = "lambdas" #lambdas, HRVs
varying_var = "M"     #variable varying: can be M, Rho, Nmult (M best)
fixed_var = "Rho"     #variable fixed:   can be, M, Rho, Nmult (Rho best)
fixed_val = 5000     #value of fixed_var 
use_selected_masses = True #gives equal spacing

#stef_txt_in_file = 1 #used when _stef.txt was added from Stef's jobs


plot_histogram_Nreps = True
plot_boundaries = 0
plot_molecules = 0
do_also_not_main_plots = 0 #those NOT for poster

#plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {#'text.usetex' : True,
          'font.size' : 20,
          #'font.family' : 'lmodern',
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
dataDir = f"{home}/MSci_Schwarzschild_Causets/data/{molecules}/"
#dataDir = home + f"/MSci_Schwarzschild_Causets/data/linkcounting_files/Poiss=False/"
if not os.path.exists(dataDir):
    path = os.getcwd()
    plotsDir = path + f"/figures/N{molecules}_vs_Area/"
    dataDir = path + f"/data/{molecules}/"
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
                        if use_selected_masses \
                            and varying_var_i not in selected_masses:
                            go_on = False
                            break #j loop
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
                except IndexError: # only one line of data in txt file
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
# 3. FIT & PLOT
##############################################################################

# Fitting function == linear through vertex
def lin_func(x,a):
    " x is Area in l units"
    return a*x

# Fitting function2 == corrected-linear through vertex
def corrected_lin_func(x,a, coeff = - 0.0464):
    " x is A in l units"
    M = np.sqrt(np.array(x)/(16*np.pi))
    return (a + coeff/M)*x

def i_exp(n, I):
    """ (1-I) * I^(n-1) """
    return (1. - I) * I**(n-1)

def i_exp_on_n(n, I):
    """ - [ I / ln(1-I) ] * I^(n-1) """
    A = - I / np.log(1-I)
    return A * I**(n-1) /n

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
# 3.0 Nreps ##############################################
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
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}Nreps.png")
    #plt.show()

    vals_not_200 = [vals [i] for i in range(len(vals)) if Nreps[i] != 200]
    reps_not_200 = [Nreps[i] for i in range(len(vals)) if Nreps[i] != 200]
    repstable = pd.DataFrame(
                  np.column_stack(
                    [vals_not_200, reps_not_200, 200-np.array(reps_not_200)]),
                  columns = ["M Value", "Current Nreps", "Reps to 200"])
    print(repstable)


###########################################################################
###########################################################################
###########################################################################
###########################################################################
# 3.1 MOLECULES'S BOUNDARIES ##############################################
if plot_boundaries:
    mintimes       = np.array(mintimes)      * Rho**(1/4)
    mintimes_std   = np.array(mintimes_std)  * Rho**(1/4)
    innermosts     = np.array(innermosts)    * Rho**(1/4)
    innermosts_std = np.array(innermosts_std)* Rho**(1/4)
    outermosts     = np.array(outermosts)    * Rho**(1/4)
    outermosts_std = np.array(outermosts_std)* Rho**(1/4)

    rc = 2; c = 1; r = 2
    plt.figure(f'Boundaries for {fixed_string}',
                tight_layout = True)

    # MINTIME #######################################################
    ax = plt.subplot(r, c, 1)
    #plt.annotate ("a)", (-0.05, 1.05), xycoords = "axes fraction", 
    #                va='bottom', ha = 'left')
    plt.errorbar(x, mintimes, mintimes_std, 
                fmt = '.', capsize = 2, color = "black",
                zorder = 10, label = r"1$\sigma$")
    plt.errorbar(x, mintimes, 3*np.array(mintimes_std), 
                fmt = '.', capsize = 4, color = "C0",
                zorder = 5, label = r"3$\sigma$")
    ax.set_xticklabels([])
    plt.ylabel(r'$t_{\mathrm{min}} [\ell]$')
    plt.legend()
    plt.grid(alpha = 0.4) 

    # INNERMOST & OUTMOST ####################################################
    ax = plt.subplot(r, c, 2)
    #plt.annotate ("b)", (-0.05, 1.05), xycoords = "axes fraction", 
    #                va='bottom', ha = 'left')
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
    plt.grid(alpha = 0.4) 
    from matplotlib.ticker import FormatStrFormatter
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}Boundaries_to_rs_in_l.png")
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}Boundaries_to_rs_in_l.pdf")
    #plt.show()



###########################################################################
###########################################################################
###########################################################################
#############################################################################
# 2.1 MOLECULES'S DISTRIBUTION ##############################################
print("")
if plot_molecules:

    plt.figure(f"Molecules for {fixed_string}")

    gradients     = [] # list of gradients of linear fit
    gradients_unc = [] # uncertainties of linear fit
    unsafe_start  = len(molecules_distr)
                       # first lambda that returns exception in curve_fit,
                       # implying not enough statistics is there
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
        Chi2, pvalue = None, None
        try:
            popt, pcov = curve_fit(lin_func, x, y, sigma=yerr,
                                    absolute_sigma=True)
            unc = np.sqrt(np.diag(pcov))
            if (unsafe_start == len(molecules_distr) and np.isnan(unc[0])):
                unsafe_start = n 
            expected = lin_func(x, *popt)
            #Chi2, pvalue = chisquare(y, expected, len(popt)) #not ideal
        except RuntimeError:
            if (unsafe_start == len(molecules_distr)):
                unsafe_start = n 
            non0_indices = tuple([np.where(np.array(y) != 0)])
            xs2      = np.array(x)[non0_indices][0]
            ys2      = np.array(y)[non0_indices][0]
            yerrs2   = np.array(yerr)[non0_indices][0]
            popt, pcov = curve_fit(lin_func, xs2, ys2, sigma=yerrs2,
                                    absolute_sigma=True)
            unc = np.sqrt(np.diag(pcov))
        print(f"Gradient factor for {n+1}-lambda"+
                f" = {round(popt[0],5)} +- {round(unc[0],5)}")
        r, pvalue = pearsonr(x, y)
        print(f"Linearity of {n+1}-lambda data:"+ 
                f"Pearson r = {round(r,3)}, p-value = {round(1-pvalue,3)}")
        #print(f"The associated Chi2 = {Chi2} and p-value = {pvalue}")
        
        
        xfit = np.linspace(min(x),max(x),100)
        #plt.plot(xfit, lin_func(xfit,*popt), '--', color="red")
        gradients.append(popt[0])
        gradients_unc.append(unc[0])
    
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'black', 
                ls = '-', alpha=1)
    #plt.annotate(rf"$\rho \: = \: {Rho:.0f}$", (0.95, 0.5), xycoords = "axes fraction",
    #        fontsize=12, va='center', ha = 'right', bbox=props)
    plt.legend(ncol = 2)
    plt.xlabel(r'Horizon Area $[\ell^2]$')
    plt.ylabel(r"Number of $\Lambda_{n}$")
    plt.grid(alpha = 0.2)
    plt.tight_layout()
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}.png") 
    plt.savefig(plotsDir + f"{fixed_string}_{molecules}.pdf") 
    #plt.show()

    if do_also_not_main_plots:
        plt.figure(f"Medium molecules for {fixed_string}")
        for n in range(3, min(6, len(molecules_distr))):
            y    = molecules_distr[n]
            yerr = molecules_distr_std[n]
            label = ("Open HRV" if n==0 else "Closed HRV ") if molecules == "HRVs" \
                    else str(n+1) + r"$\mathbf{-}\Lambda$"
            plt.errorbar(x, y, yerr, 
                        fmt = '.', capsize = 4, 
                        label = label)
        plt.legend()
        plt.xlabel(r'Horizon Area $[\ell^2]$') #not yet in terms of l^2
        plt.ylabel(r"Number of $n\mathbf{-}\Lambda$")
        plt.grid(alpha = 0.2)
        plt.tight_layout()
        plt.savefig(plotsDir + f"{fixed_string}_large_{molecules}.png")
        plt.savefig(plotsDir + f"{fixed_string}_large_{molecules}.pdf") 
        #plt.show()


        plt.figure(f"Large molecules for {fixed_string}")
        for n in range(min(6, len(molecules_distr)), len(molecules_distr)):
            y    = molecules_distr[n]
            yerr = molecules_distr_std[n]
            label = ("Open HRV" if n==0 else "Closed HRV ") if molecules == "HRVs" \
                    else str(n+1) + r"$\mathbf{-}\Lambda$"
            plt.errorbar(x, y, yerr, 
                        fmt = '.', capsize = 4, 
                        label = label)
        plt.legend()
        plt.xlabel(r'Horizon Area $[\ell^2]$') #not yet in terms of l^2
        plt.ylabel(r"Number of $n\mathbf{-}\Lambda$")
        plt.grid(alpha = 0.2)
        plt.tight_layout()
        plt.savefig(plotsDir + f"{fixed_string}_XXL_{molecules}.png")
        plt.savefig(plotsDir + f"{fixed_string}_XXL_{molecules}.pdf") 
        #plt.show()



        plt.figure(f"Small molecules uncertainty for {fixed_string}")
        for n in range(0, 4):
            yerr = molecules_distr_std[n]
            label = ("Open HRV" if n==0 else "Closed HRV ") if molecules == "HRVs" \
                    else str(n+1) + r"$\mathbf{-}\Lambda$"
            plt.errorbar(x, yerr, 
                        fmt = '.', capsize = 4, 
                        label = label)
        plt.legend()
        plt.xlabel(r'Horizon Area $[\ell^2]$') #not yet in terms of l^2
        plt.ylabel(r"Uncertainty of $n\mathbf{-}\Lambda$")
        plt.grid(alpha = 0.2)
        plt.tight_layout()
        plt.savefig(plotsDir+f"{fixed_string}_uncertainty_small_{molecules}.png")
        plt.savefig(plotsDir+f"{fixed_string}_uncertainty_small_{molecules}.pdf") 
        #plt.show()


    ##########################################################################
    ##########################################################################
    ### DO THE NUMERICAL ANALYSIS
    ##########################################################################
    ##########################################################################
    noninf_indices = tuple([np.where(np.isfinite(gradients_unc))[0]])
    gradients      = np.array(gradients)[noninf_indices]
    gradients_unc  = np.array(gradients_unc)[noninf_indices]


    ##########################################################################
    ### Find Link proportionality (Barton et al)
    ##########################################################################
    print(" \n#### LINKS ANALYSIS (for Barton et al a = 0.17321) ####")
    # Link Counting from weighted sum of gradients (to compare to Barton et al.)
    coefsum = sum([(n+1)*gradients[n] for n in range(len(gradients))])
    coefsum_unc = np.sqrt(
        sum([((n+1)*gradients_unc[n])**2 for n in range(len(gradients))]))
    print(f"Gradients Weighted Sum = {round(coefsum,4)} +- {round(coefsum_unc,4)}")

    # Find total number of links for each x
    links      = np.zeros(len(molecules_distr[0]))
    links_std2 = np.zeros(len(molecules_distr[0]))
    for n in range(len(molecules_distr)):
        links      += (n+1)    *np.array(molecules_distr[n])
        links_std2 += (n+1)**2 * np.array(molecules_distr_std[n])**2
    links_std = np.sqrt(links_std2)
    
    # Fit to linear fit
    popt, pcov = curve_fit(lin_func, x, links, sigma=links_std,
                            absolute_sigma=True)
    unc = np.sqrt(np.diag(pcov))
    print(f"Gradient of Links      = {round(popt[0],4)} +- {round(unc[0],4)}")

    # Fit to corrected linear fit
    def schwarz_lin_func(x, a):
        return corrected_lin_func(x, a, coeff = +0.0464)
    popt, pcov = curve_fit(schwarz_lin_func, x, links, sigma=links_std,
                            absolute_sigma=True)
    unc = np.sqrt(np.diag(pcov))
    print(f"Curv-Correct Gradient  = {round(popt[0],4)} +- {round(unc[0],4)}")

    r, pvalue = pearsonr(x, links)
    print(f"Linearity links: Pearson r = {round(r,3)}, p-val = {round(1-pvalue,3)}")


    
    ##########################################################################
    ### Find Probability Distribution of Lambda
    ##########################################################################
    grad_sum = sum(gradients)
    grad_sum_unc = np.sqrt(sum([g**2 for g in gradients_unc]))
    
    lambd_occurs = [sum(molecules_distr[n]) 
                        for n in range(len(molecules_distr))]
    lambd_probs = gradients/sum(gradients)
    lambd_probs_uncs = np.sqrt(gradients_unc**2/grad_sum**2 +
                               gradients**2*grad_sum_unc**2/grad_sum**4)

    ###################################################################
    # Fit to exponential (1-I) * I**n

    unsafe_start = len(noninf_indices[0]) #gets rid of data with bad statistics
    ns = np.arange(1, unsafe_start+1)
    popt, pcov = curve_fit(i_exp, ns, lambd_probs[:unsafe_start], 
                            p0 = 1-lambd_probs[0],
                            sigma=lambd_probs_uncs[:unsafe_start],
                            absolute_sigma=True)
    unc = np.sqrt(np.diag(pcov))
    #expected = i_exp(ns, *popt)
    #Chi2, pvalue = chisquare(lambd_occurs[:unsafe_start], expected, len(popt))
    print(" \n#### N-LAMBDAS DISTRIBUTION ####")
    print("Distribution of n-lambdas:\n",
            [round(pn,7) for pn in lambd_probs])
    print("Ratio of n+1_to_n-lambdas:\n",
            [round(lambd_probs[n+1]/lambd_probs[n],5) 
            for n in range(len(lambd_probs)-1)])
    print("Ratio of n_to_n+1-lambdas:\n",
            [round(lambd_probs[n]/lambd_probs[n+1],5) 
            for n in range(len(lambd_probs)-1)])
    print("\nFIT RESULTS")
    print("Exponential fit (1-I) I^(n-1) has:")
    print(f" - I = {round(popt[0],4)} +- {round(unc[0],4)}")
    print(f"With e^-x rather than I, it has:")
    print(f" - x = {-round(np.log(popt[0]),4)} +- {round(unc[0]/popt[0],4)}")
    #print(f"The associated Chi2 = {Chi2} and p-value = {pvalue}")

    ###################################################################
    # Do numerical calculation rather than fit for exponential 

    r1 = grad_sum/coefsum
    r1unc = grad_sum_unc/coefsum
    I = 1 - r1
    Iunc = r1unc
    chi = - np.log(I)
    chiunc = Iunc/I
    chi_ord = 0
    chiunc_copy = chiunc*1.
    while abs(chiunc_copy) < 1:
        chiunc_copy *= 10
        chi_ord += 1
    print("\nPLAIN NUMERICAL RESULTS")
    print(f"Exponential model (1-I) I^(n-1) has:")
    print(f" - I = {round(I,4)} +- {round(Iunc,4)}")
    print(f"With e^-x rather than I, it has:")
    print(f" - x = {round(chi,4)} +- {round(chiunc,4)}")
    print(f"p1 = {round(lambd_probs[0],4)} == gradsum/a_links = {round(r1,4)}?")


    ###################################################################
    # Do numerical calculation rather than fit for exponential - linear fall

    r0 = gradients[0]/coefsum
    r0unc = gradients_unc[0]/coefsum
    I2 = 1 - r0
    I2unc = r0unc
    chi2 = - np.log(I2)
    chi2unc = I2unc/I2
    print("\nPLAIN NUMERICAL RESULTS")
    print(f"Exponential-LinearFall model 1/ln[1/(1-I)] I^(n)/n has:")
    print(f" - I = {round(I2,4)} +- {round(I2unc,4)}")
    print(f"With e^-x rather than I, it has:")
    print(f" - x = {round(chi2,4)} +- {round(chi2unc,4)}")

    

    #################################################################
    # PLot Distribution (all, all in logscale, small in logscale)
    if do_also_not_main_plots:
        plt.figure("n-lambda probability distribution")
        plt.errorbar(np.arange(1,len(lambd_probs)+1,1), lambd_probs,
                yerr=lambd_probs_uncs,capsize=7,fmt="",ls="",ecolor="red",
                label = r"$n\mathbf{-}\Lambda$ probability distribution")
        xs = np.linspace(1, len(lambd_probs)+1,100)
        plt.plot(xs, i_exp(xs, *popt), ls = "--", color = "green",
                label = r"$(1-\mathcal{I})$ $\mathcal{I}^{(n-1)}$"+ 
                f", I = {round(I,3)}+-{round(Iunc,3)}")
        plt.xlabel(r"$n$")
        plt.ylabel("Probability")
        plt.legend(loc="upper right")
        plt.grid(alpha=0.2)
        plt.tight_layout()
        plt.savefig(plotsDir + "n_lambda_probability_distribution.png")
        plt.savefig(plotsDir + "n_lambda_probability_distribution.pdf")
        #plt.show()


        plt.figure("n-lambda probability distribution (logscale)")
        plt.errorbar(np.arange(1,len(lambd_probs)+1,1), lambd_probs,
                yerr=lambd_probs_uncs,capsize=7,fmt="",ls="",ecolor="red",
                label = r"$n\mathbf{-}\Lambda$ distribution")
        xs = np.linspace(1, len(lambd_probs)+1,100)
        plt.plot(xs, i_exp(xs, *popt), ls = "--", color = "gold",
                label = r"$(1-\mathcal{I})$ $\mathcal{I}^{(n-1)}$"+ 
                f", I = {round(I,3)}+-{round(Iunc,3)}")
        # plt.plot(xs, i_exp_on_n(xs, I2), ls = "--", color = "green",
        #          label = r"$\frac{-I}{1-I} \frac{I^{n-1}}{n}$"+
        #          f" I = {round(I2,3)}+-{round(I2unc,3)}")
        plt.xlabel(r"$n$")
        plt.ylabel("Probability")
        plt.yscale("log")
        plt.legend(loc="upper right")
        plt.grid(alpha=0.2)
        plt.savefig(plotsDir + "n_lambda_probability_distribution_I_logy.png")
        plt.savefig(plotsDir + "n_lambda_probability_distribution_I_logy.pdf")
        #plt.show()

        plt.figure("n-lambda probability distribution (safe)")
        plt.errorbar(np.arange(1,len(lambd_probs[:unsafe_start])+1,1), 
                lambd_probs[:unsafe_start],
                yerr=lambd_probs_uncs[:unsafe_start],
                capsize=7,fmt="",ls="",ecolor="red")
        xs = np.linspace(1, len(lambd_probs[:unsafe_start])+1,50)
        plt.plot(xs, i_exp(xs, *popt), ls = "--", color = "gold",
                label = r"(1-$e^{- \chi (n-1)}$) $e^{- \chi (n-1)}$"+ 
                r", $\chi$"+ 
                f" = {round(chi,3)}+-{round(chiunc,3)}")
        plt.xlabel(r"$n$")
        plt.ylabel(r"Probability $p_n$")
        plt.legend(loc="upper right")
        plt.grid(alpha=0.2)
        plt.tight_layout()
        plt.savefig(plotsDir + "n_lambda_probability_distribution_small.png")
        plt.savefig(plotsDir + "n_lambda_probability_distribution_small.pdf")
        #plt.show()

    x = 8
    plt.figure("n-lambda exp probability distribution (logscale)")
    plt.errorbar(np.arange(1,len(lambd_probs)+1,1), lambd_probs,
            yerr=lambd_probs_uncs,capsize=7,fmt="",ls="",ecolor="red",
            label = r"$\Lambda_n$ distribution")
    xs = np.linspace(1, len(lambd_probs)+0.2,100)
    plt.plot(xs, i_exp(xs, *popt), ls = "--", color = "darkorange",
            label = r"$(e^{\chi}-1)$ $e^{- n \chi}$, $\chi$"+ 
            f" = {round(chi,chi_ord)}"+
            f"({int(round(chiunc,chi_ord)*10**chi_ord)})")
    plt.xlabel(r"$n$")
    plt.ylabel("Probability")
    plt.yscale("log")
    plt.legend(loc="upper right", fontsize = 24)
    plt.grid(alpha=0.2)
    plt.xticks(np.arange(1,len(lambd_probs)+1,1))
    plt.tight_layout()
    plt.savefig(plotsDir + "n_lambda_probability_distribution_expx_logy.png")
    plt.savefig(plotsDir + "n_lambda_probability_distribution_expx_logy.pdf")
    


    


    ##########################################################################
    ### Find Entropy and C_hv and Discretness Length Scale
    ##########################################################################
    print(" \n#### FINAL RESULTS ####")
    # Find overall proportionality and discreteness scale
    C_hv = sum(gradients*np.log(grad_sum/gradients))
    dC_da_i = np.log(grad_sum/gradients)
    C_hv_unc = np.sqrt( sum(gradients_unc**2 * dC_da_i**2))  
    print(f"Overall Proportion C_hv = {round(C_hv,5)} +- {round(C_hv_unc,5)}")

    # Discreteness length in terms of Planckian length
    l = 2*np.sqrt(C_hv)
    l_unc = C_hv_unc/np.sqrt(C_hv)
    print(f"Discreteness scale      = {round(l,5)} +- {round(l_unc,5)} l_p\n")

plt.show()
 