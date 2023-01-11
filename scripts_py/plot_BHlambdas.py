
#%% 
import causetplotting as ocplt

import matplotlib.pyplot as plt
import numpy as np
import os 

########################################################
# Load file
dim = 3
card = 100
edge = 4

phi0 = 6.29/1.
phi_limits = (0, 6.29)
if phi0 < 6.28:
    phi_limits = np.array([0, phi0]) + 0*phi0
projection = 0
centre_cube_in_horizon = False

use_redge_in_name_file = True
isBH = True
isEFv = False
isS = False


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



path = os.getcwd() # folder path
print(path)
# os.chdir("../")
# path = os.getcwd()
# print(path)
#file_name = os.path.join(path, 'data/flatspace_bicone_causet.txt')
#file_name = os.path.join(path, 'data/known_causet_from_matrixSetsTest.txt')

if isBH:
    if isEFv and isS:
        file_name = f"data/blackhole_EFvToS_{dim}D_N{card}"
    elif isEFv:
        file_name = f"data/blackhole_EFv_{dim}D_N{card}"
    elif isS:
        file_name = f"data/blackhole_S_{dim}D_N{card}"
    else:
        file_name = f"data/blackhole_and_lambdas{dim}D_N{card}"
else:
    file_name = f"data/flat{dim}D_N{card}"

file_name = path + "/"+ file_name
#path = os.chdir("scripts_py")

if use_redge_in_name_file:
    edge_string = str(round(edge, 2))
    if len(edge_string) == 1:
        edge_string += "."
    while len(edge_string) < 4:
        edge_string += "0"
    file_name += "_redge"+ edge_string

if (centre_cube_in_horizon):
        file_name += "_horizon_centred"

file_name += ".txt"
print(file_name)

ocplt.plot_causet_and_lambdas(file_name, phi_limits = phi_limits, projection = projection)



# #Plot cones inside horizon crossing a point (t0, r0)
# if dim == 2 or projection == True:
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
#             ax.plot(rs, upper_null(rs, t0, r0), ls = "--", c = "green", alpha = 0.8)
#         else:
#             ax.plot(rs, upper_null(rs, t0, r0), ls = "--", c = "green", alpha = 0.8,
#                     label = "Light Cone")
# %%
