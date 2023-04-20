from causets_py import causetplotting as cplt
import matplotlib.pyplot as plt
import numpy as np
import os 

#%% 
#############################################################################
## SET PARAMETERS FOR PLOTTING
#############################################################################
# Load file
dim  = 2
card = 50
edge = 4.00 #or radius
r_S = 2
two_sided = False

use_redge_in_name_file = 1 #the size is given by _redge<number> in filename

##############################################################################
##############################################################################


ps = {"text.usetex": True,
      "font.size" : 20,
      "font.family" : "lmodern",
      "axes.labelsize": 24,
      "legend.fontsize": 20,
      "xtick.labelsize": 18,
      "ytick.labelsize": 18,
      "figure.figsize": [6, 12],
      "mathtext.default": "default"
       }
plt.rcParams.update(ps)
del ps



###################################################################
#%%### GET FILE
###################################################################
savefile_ext = "data/BHplotting/"
savefile_ext += f"{dim}D_N{card}_redge{edge}"

path = os.getcwd() # folder path
print(path)
file_name = "data/data_for_plotting/"
file_name += f"blackhole{dim}D_N{card}"
file_name = path + "/"+ file_name

# redge size
if use_redge_in_name_file:
    edge_string = str(round(edge, 2))
    if len(edge_string) == 1:
        edge_string += "."
    while len(edge_string) < 4:
        edge_string += "0"
    file_name += "_redge"+ edge_string

file_name += ".txt"
print(file_name)


###################################################################
#%%### PLOT
###################################################################
ax = cplt.plot_causet(file_name,# savefile_ext=savefile_ext+".pdf",
                      link_alpha = 0.5,
                        figsize = (7,7))



#Plot cones inside horizon crossing a point (t0, r0)
if dim == 2:
    print("adding horizon")
    ys = plt.ylim()
    plt.vlines(r_S, ys[0], ys[1], lw = 2.5, ls = "--", color = "red",
               label = "Horizon")
    if two_sided:
        plt.vlines(-r_S, ys[0], ys[1], lw = 2.5, ls = "--",color="r")
    plt.ylim(ys)

    print("adding cones")
    rs = [0.819, 1.195, 3.118]
    ts = [0.292, 2.087, 1.643]
    def upper_null(r, t0, r0, mass=1):
        return t0 + r - r0 + 4*mass*np.log( (2*mass-r)/(2*mass-r0) )
    def lower_null(r, t0, r0):
        return t0 + r0 - r
    for i, (r0, t0) in enumerate(zip(rs, ts)):
        rs = np.linspace(0, r0, 1000)
        ax.plot(rs, lower_null(rs, t0, r0), lw = 2.5, ls = "--", c = "green", 
                alpha = 1)
        if r0 == 2:
            continue
        elif r0 > 2:
            rs = np.linspace(r0, edge, 1000)
        if i != len(ts) - 1:
            ax.plot(rs, upper_null(rs, t0, r0), lw = 2.5, ls = "--", c = "green", 
                    alpha = 1)
        else:
            ax.plot(rs, upper_null(rs, t0, r0), lw = 2.5, ls = "--", c = "green", 
                    alpha = 1,
                    label = "Light cone")
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_yticks([0, 1, 2, 3, 4])
plt.legend(loc = "lower right")
plt.savefig(savefile_ext+"_lcones.pdf")
plt.savefig(savefile_ext+"_lcones.png")
plt.show()