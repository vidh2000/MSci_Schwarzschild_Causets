from causets_py import causetplotting as cplt
import matplotlib.pyplot as plt
import numpy as np
import os 

#%% 
#############################################################################
## SET PARAMETERS FOR PLOTTING
#############################################################################
# Load file
dim  = 4
card = 86000
edge = 2.36 #or radius
h    = 0.70
want_save_file = 1 #see line 31 for saving name

use_redge_in_name_file = 1 #the size is given by _redge<number> in filename
use_h                  = 1 #there is the bit _ha.bd before .txt
centre_cube_in_horizon = 0

#choices: [0->causet, 1->molecules, 2->molecules in horizon, 3->causet+molecul]
plot_choice = 2
molecule = "lambda"

# Only plotting the <which_phi_interval>th interval of size phi0
# If plot_choice = 0 and 3D and project, project data in interval in 2D r-t
phi0 = 6.29/1#64#
which_phi_interval = 0
projection = 1

##############################################################################
##############################################################################
plot_choice_string = ["causet","molecules","horizon","all"]
if want_save_file:
    savefile_ext = "data/BHplotting/"\
                 f"BH_{molecule}_plot{plot_choice_string[plot_choice]}_"
    if plot_choice == 0 and projection:
        savefile_ext += f"proj{round(phi0,2)}_"
    savefile_ext += f"{dim}D_N{card}_redge{edge}_h{h}.pdf"
else:
    savefile_ext = 0

ps = {#"text.usetex": True,
      "font.size" : 16,
      "font.family" : "Times New Roman",
      "axes.labelsize": 16,
      "legend.fontsize": 14,
      "xtick.labelsize": 14,
      "ytick.labelsize": 14,
      "figure.figsize": [7.5, 6],
      "mathtext.default": "default"
       }
plt.rcParams.update(ps)
del ps


###################################################################
#%%### SET PHI LIMITS
###################################################################
phi_limits = (0, 6.29)
if phi0 < 6.28:
    phi_limits = np.array([0, phi0]) + which_phi_interval*phi0


###################################################################
#%%### GET FILE
###################################################################
path = os.getcwd() # folder path
if "scripts_py" in path:
    path = path[:-11]
print("PATH")
print(path)
file_name = "data/data_for_plotting/"
file_name += "blackhole_and_"
file_name += "lambdas" if (molecule == "lambda") else "HRVs"
file_name += f"{dim}D"
file_name += f"_N{card}"
file_name = path + "/"+ file_name

# redge size
if use_redge_in_name_file:
    edge_string = str(round(edge, 2))
    if len(edge_string) == 1:
        edge_string += "."
    while len(edge_string) < 4:
        edge_string += "0"
    file_name += "_redge"+ edge_string

# hollowness h
if use_h:
    hstring = str(round(h, 2))
    if len(hstring) == 1:
        hstring += "."
    while len(hstring) < 4:
        hstring += "0"
    file_name += "_h"+hstring

# centred in horizon?
if (centre_cube_in_horizon):
        file_name += "_horizon_centred"

file_name += ".txt"
print("FILE NAME")
print(file_name)

print("SAVEFILE")
if savefile_ext:
    savefile_ext = path+"/"+savefile_ext
    print(savefile_ext)
else:
    print("savefile is 0")


###################################################################
#%%### PLOT
###################################################################
#choices: [0->causet, 1->molecules, 2->molecules in horizon, 3->causet+molecul]
if plot_choice == 0:
    if projection:
        print(f"Plotting 0 with projection of deltaphi = {round(phi0,2)}")
    else:
        print(f"Plotting 0")
    ax = cplt.plot_causet(file_name, savefile_ext = savefile_ext,
                        phi_limits = phi_limits, 
                        projection=projection)
elif plot_choice == 1:
    print("Plotting 1")
    ax = cplt.plot_lambdas(file_name, savefile_ext = savefile_ext,
                            phi_limits = phi_limits)
elif plot_choice == 2:
    print("Plotting 2")
    ax = cplt.plot_lambdas_horizon(file_name, savefile_ext = savefile_ext,
                                    phi_limits = phi_limits,
                                    figsize = (7,7))
elif plot_choice == 3:
    print("Plotting 3")
    ax = cplt.plot_causet_and_lambdas(file_name,savefile_ext = savefile_ext,
                                    phi_limits = phi_limits)
# ax = cplt.plot_causet_and_lambdas(file_name, 
#                                   phi_limits = phi_limits)


# #Plot cones inside horizon crossing a point (t0, r0)
# if dim == 2 or projection == True:
#     print("adding cones")
#     rs = [0.72, 1.35, 2.77]
#     ts = [1.15, 2.66, 1.35]
    
#     def upper_null(r, t0, r0, mass=1):
#         return t0 + r - r0 + 4*mass*np.log( (2*mass-r)/(2*mass-r0) )
#     def lower_null(r, t0, r0):
#         return t0 + r0 - r

#     for i, (r0, t0) in enumerate(zip(rs, ts)):
#         rs = np.linspace(0, r0, 1000)
#         ax.plot(rs, lower_null(rs, t0, r0), ls = "--", c = "green", alpha = 0.8)
#         if r0 == 2:
#             continue
#         elif r0 > 2:
#             rs = np.linspace(r0, edge, 1000)
#         if i != len(ts) - 1:
#             ax.plot(rs, upper_null(rs, t0, r0), ls = "--", c = "green", 
#                     alpha = 0.8)
#         else:
#             ax.plot(rs, upper_null(rs, t0, r0), ls = "--", c = "green", 
#                     alpha = 0.8,
#                     label = "Light Cone")
plt.show()
# %%
