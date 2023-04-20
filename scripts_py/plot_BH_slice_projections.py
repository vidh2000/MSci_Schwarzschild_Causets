from causets_py import causetplotting as cplt
import matplotlib.pyplot as plt
import numpy as np
import os 

ps = {"text.usetex": True,
      "font.size" : 16,
      "font.family" : "lmodern",
      "axes.labelsize": 22,
      "legend.fontsize": 20,
      "xtick.labelsize": 20,
      "ytick.labelsize": 20,
      "figure.figsize": [7.5, 7.5],
      "mathtext.default": "default"
       }
plt.rcParams.update(ps)
del ps


#%% 
#############################################################################
## SET PARAMETERS FOR PLOTTING
#############################################################################

# Load correct file
dim  = 3
card = 100
edge = 4.00 #or radius
phi0 = 0.25
h    = 0.00

use_redge_in_name_file = 1 #the size is given by _redge<number> in filename
use_h                  = 0 #there is the bit _ha.bd before .txt
centre_cube_in_horizon = 0

want_save_file = 1 #save file? (see line 100 for saving name)
n_save_file_max = 50




#1 SET NAME OF FILE TO GET INFO FROM
###################################################################
path = os.getcwd() # folder path
if "scripts_py" in path:
    path = path[:-11]
print("PATH")
print(path)
file_name = "data/data_for_plotting/"
file_name += "blackhole_slice"
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


#2. SET SAVING FILE
###################################################################
print("SAVEFILE")
if want_save_file:
    savefile_ext = "data/BHplotting/"\
                f"BH_plot_slice{round(phi0,2)}_"
    savefile_ext += f"{dim}D_N{card}_redge{edge}_nolcones.pdf"
else:
    savefile_ext = 0
if savefile_ext:
    savefile_ext = path+"/"+savefile_ext
    print(savefile_ext)
else:
    print("savefile is 0")


#3. PLOT
###################################################################
ax = cplt.plot_causet(file_name, savefile_ext = savefile_ext,
                    projection=1,
                    link_alpha = 0.42,
                    figsize = (7,7))


#Plot cones inside horizon crossing a point (t0, r0)
if 1:
    print("adding cones")
    rs = [ 2.9403,  1.210,  1.038,  2.212]
    ts = [-3.1434, -3.228, -0.294, -1.534]
    
    def upper_null(r, t0, r0, mass=1):
        return t0 + r - r0 + 4*mass*np.log( (2*mass-r)/(2*mass-r0) )
    def lower_null(r, t0, r0):
        return t0 + r0 - r

    for i, (r0, t0) in enumerate(zip(rs, ts)):
        rs = np.linspace(0, r0, 1000)
        ax.plot(rs, lower_null(rs, t0, r0), ls = "--", c = "green", 
                  alpha = 0.8)
        if r0 == 2:
            continue
        elif r0 > 2:
            rs = np.linspace(r0, edge, 1000)
        if i != len(ts) - 1:
            ax.plot(rs, upper_null(rs, t0, r0), ls = "--", 
                    c = "green", alpha = 0.8)
        else:
            ax.plot(rs, upper_null(rs, t0, r0), ls = "--", 
                    c = "green", alpha = 0.8, label = "Light Cone")
            
ys = ax.set_ylim()
r_S = 2
ax.vlines(r_S, ys[0], ys[1], ls = "--", color = "red", label = "Horizon")
ax.set_ylim(ys)
plt.legend()
plt.show()
# %%
