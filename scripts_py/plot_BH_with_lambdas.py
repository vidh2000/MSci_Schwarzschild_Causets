from causets_py import causetplotting as cplt
import matplotlib.pyplot as plt
import numpy as np
import os 

#%% 
#############################################################################
## SET PARAMETERS FOR PLOTTING
#############################################################################
# Load file
dim  = 3
card = 3000
edge = 4.00 #or radius
h    = 0.00

use_redge_in_name_file = 1 #the size is given by _redge<number> in filename
use_h                  = 0
centre_cube_in_horizon = 0

#choices: [0->causet, 1->molecules, 2->molecules in horizon, 3->causet+molecul]
plot_choice = 0
molecule = "lambda"

phi0 = 6.29/24#4
which_phi_interval = 0
projection = 1
##############################################################################
##############################################################################


ps = {"text.usetex": True,
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
print(file_name)


###################################################################
#%%### PLOT
###################################################################
#choices: [0->causet, 1->molecules, 2->molecules in horizon, 3->causet+molecul]
if plot_choice == 0:
    ax = cplt.plot_causet(file_name, phi_limits = phi_limits, 
                        projection=projection)
elif plot_choice == 1:
    ax = cplt.plot_lambdas(file_name, phi_limits = phi_limits)
elif plot_choice == 2:
    ax = cplt.plot_lambdas_horizon(file_name, phi_limits = phi_limits,
                                    figsize = (7,7))
elif plot_choice == 3:
    ax = cplt.plot_causet_and_lambdas(file_name, phi_limits = phi_limits)
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
