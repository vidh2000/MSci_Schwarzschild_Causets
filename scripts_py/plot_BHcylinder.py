from causets_py import causetplotting as cplt
import matplotlib.pyplot as plt
import numpy as np
import os 

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

#%% 
#############################################################################
## SET PARAMETERS FOR PLOTTING
#############################################################################
# Load file
dim  = 3
card = 200
edge = 1.5 #or radius
h    = 0.5
T    = 5 # height in time
want_save_file = 0 #see line 31 for saving name

use_redge_in_name_file = 1 #the size is given by _redge<number> in filename
use_h                  = 1 #there is the bit _h<> before .txt
centre_cube_in_horizon = 0 #there is the bit _horizon_centred before .txt

#choices: [0->causet, 1->molecules, 2->molecules in horizon, 3->causet+molecul]
plot_choice = 0
molecule = "none" # none, or lambda or HRV

# Only plotting the <which_phi_interval>th interval of size phi0
# If plot_choice = 0 and 3D and project, project data in interval in 2D r-t
phi0 = 6.29/1#64#
which_phi_interval = 0
projection = 0

# Plotting kwargs
link_lw = 0

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
file_name += "blackhole"
if molecule == "lambda":
    file_name += "_and_lambdas"
elif molecule == "HRV":
    file_name += "_and_HRVs"
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
                        projection=projection,
                        link_lw = link_lw)
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


# Add cylinder
# from https://stackoverflow.com/questions/26989131/add-cylinder-to-plot 
def cylinder_along_z(center_x,center_y,center_z,radius,height_z):
    z = np.linspace(-(height_z/2), (height_z/2), 100) + center_z
    theta = np.linspace(0, 2*np.pi, 100)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

def circle(center_x, center_y, radius, z):
    theta = np.linspace(0, 2*np.pi, 100)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid[0],y_grid[0],z_grid[0]
    
Xc,Yc,Zc = cylinder_along_z(0, 0, 0, edge, T)
ax.plot_surface(Xc, Yc, Zc, alpha=0.3, 
                color = "red")
for i in range(2):
    t = -T/2 + T*i/1
    ax.plot(*circle(0,0,edge,t), color = "blue", lw = 0.5)

if use_h:
    Xc_in,Yc_in,Zc_in = cylinder_along_z(0, 0, 0, h*edge, T)
    ax.plot_surface(Xc_in, Yc_in, Zc_in, alpha=0.7, 
                color = "red", lw = 0.1)
    for i in range(2):
        t = -T/2 + T*i/1
        ax.plot(*circle(0,0,h*edge,t), color = "blue", lw = 0.5)
    

ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

plt.show()
# %%
