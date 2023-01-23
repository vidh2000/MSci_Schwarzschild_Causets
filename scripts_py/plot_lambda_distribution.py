import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import expanduser

usehome = True
molecules = "lambdas" #links, HRVs
want_Nmult = False
want_rho = True
mass = round(0.25,2)

# Home Directory
home = expanduser("~")
# Path
path = os.getcwd()
if usehome:
    plotsDir = f"{home}/MSci_Schwarzschild_Causets/figures/N{molecules}_vs_Area/"
    dataDir = f"{home}/MSci_Schwarzschild_Causets/data/{molecules}/"
    #dataDir = home + f"/MSci_Schwarzschild_Causets/data/linkcounting_files/Poiss=False/"
else:
    plotsDir = path + f"/figures/N{molecules}_vs_Area/"
    dataDir = path + f"/data/{molecules}"


mass_string = f"M={mass}"

rhos = []
Nreps = []
outermosts = []
outermosts_std = []
innermosts = []
innermosts_std = []
mintimes = []
mintimes_std = []
molecules_distr = []
molecules_distr_std = []


for root, dirs, files in os.walk(dataDir):
    # for each file file_i
    for i, file_i in enumerate(files):
        if mass_string in file_i:
            if want_rho and "Rho" in file_i:
                file_path = os.path.join(root, file_i)
                file_pieces = file_i.split("_")
                # get rho of file_i
                for piece_j in file_pieces:
                    if piece_j[0:3] == "Rho":
                        rho_i = piece_j.split("=")[1]
                #DO SOMETHING
            elif want_Nmult and "Nmult" in file_i:
                file_path = os.path.join(root, file_i)
                file_pieces = file_i.split("_")
                # get rho of file_i
                for piece_j in file_pieces:
                    if piece_j[0:3] == "Nmult":
                        nmult_i = piece_j.split("=")[1]
           