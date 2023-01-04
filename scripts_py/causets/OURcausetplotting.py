#!/usr/bin/env python
'''
Created on 4 Jan 2022

@author: Stefano Veroni, Vid Homsak
'''
import numpy as np
import matplotlib.pyplot as plt

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



def plot_lambdas_and_lambdasdistr(lambdasfile_ext, savefile_ext = 0):
    """
    PLot the lambdas of a causet with Horizon and their distribution.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to

    Returns
    -------
    plt.figure
    plt.axs
    """
    data = np.genfromtxt(lambdasfile_ext, delimiter=",")

    storage_option = data[0,0]
    size           = data[1,1]
    dim            = data[2,1]
    shapename      = data[3,1]
    spacetimename  = data[4,1]
    if dim > 3:
        raise AttributeError(f"Dim is {dim}: too big!!!")

    distribution = []
    lambdas = []
    coords = []
    go = -1
    while go < 0:
        row = data[go]
        key = row[0]
        if key[0] == "L":
            lambdas.append(row[1:])
            go -= 1
        elif key[0] == "N":
            distribution.append([key[-1], row[1]])
            go -= 1
        elif key == "r_S":
            r_S = row[1]
            go -= 1
        elif key == "Coordinates":
            break
        else:
            coords.insert[0, row]


    rc = 2
    r  = 2
    c  = 1
    figsize = (c * 4, r * 2.2)
    fig, axs  = plt.figure(figsize = (10, 15), tight_layout = True)
    plt.suptitle( "Lambdas Information" )
    lambdas_colors = ["green", "yellow", "orange", "red", "purple"]
    Ncolors = len(lambdas_colors)
    
    plt.subplot(r, c, 1)
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            plt.plot(xup, tup, marker="o", 
                        markersize=20, 
                        markeredgecolor="red", 
                        markerfacecolor="cs:black",
                        zorder = 2)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot(x, t, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [tup, t], 
                        markerfacecolor=facecolor,
                        zorder = 3)
    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            yup = coords[uplabel][2]
            plt.plot(xup, tup, marker="o", 
                        markersize=20, 
                        markeredgecolor="red", 
                        markerfacecolor="cs:black",
                        zorder = 2)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                y = coords[label][2]
                plt.plot(x, y, t, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [yup, y], [tup, t], 
                        markerfacecolor=facecolor,
                        zorder = 3)

    
    plt.subplot(r, c, 2)
    plt.xlabel("Lambdas' size (number of links)")
    plt.ylabel("Number of Occurrences")
    plt.step(distribution[:,0], distribution[:,1])
    props = dict(boxstyle='round', facecolor='wheat', 
                edgecolor = 'black', ls = '', alpha=0.5)
    plt.annotate (f"Tot = {sum(distribution[:,1])}", 
                  (0.8, 0.8), xycoords = "axes fraction", 
                  va='bottom', ha = 'left', bbox=props) 

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig, axs




def plot_lambdas(lambdasfile_ext, savefile_ext = 0):
    """
    PLot the lambdas of a causet with Horizon.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to

    Returns
    -------
    plt.figure
    plt.axs
    """
    data = np.genfromtxt(lambdasfile_ext, delimiter=",")

    storage_option = data[0,0]
    size           = data[1,1]
    dim            = data[2,1]
    shapename      = data[3,1]
    spacetimename  = data[4,1]
    if dim > 3:
        raise AttributeError(f"Dim is {dim}: too big!!!")

    lambdas = []
    coords = []
    go = -1
    while go < 0:
        row = data[go]
        key = row[0]
        if key[0] == "L":
            lambdas.append(row[1:])
            go -= 1
        elif key[0] == "N":
            go -= 1
        elif key == "r_S":
            r_S = row[1]
            go -= 1
        elif key == "Coordinates":
            break
        else:
            coords.insert[0, row]


    #PLOT
    fig, axs  = plt.figure(figsize = (15, 15), tight_layout = True)
    lambdas_colors = ["green", "yellow", "orange", "red", "purple"]
    Ncolors = len(lambdas_colors)
    
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            plt.plot(xup, tup, marker="o", 
                        markersize=20, 
                        markeredgecolor="red", 
                        markerfacecolor="cs:black",
                        zorder = 2)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot(x, t, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [tup, t], 
                        markerfacecolor=facecolor,
                        zorder = 3)
        plt.vlines(r_S)

    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            yup = coords[uplabel][2]
            plt.plot(xup, tup, marker="o", 
                        markersize=20, 
                        markeredgecolor="red", 
                        markerfacecolor="cs:black",
                        zorder = 2)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                y = coords[label][2]
                plt.plot(x, y, t, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [yup, y], [tup, t], 
                        markerfacecolor=facecolor,
                        zorder = 3)

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig, axs





def plot_causet_and_lambdas_and_lambdasdistr(lambdasfile_ext, 
                                            savefile_ext = 0):
    """
    PLot the whole causet highlighting the lambdas and their distribution.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to

    Returns
    -------
    plt.figure
    plt.axs
    """
    data = np.genfromtxt(lambdasfile_ext, delimiter=",")

    storage_option = data[0,0]
    size           = data[1,1]
    dim            = data[2,1]
    shapename      = data[3,1]
    spacetimename  = data[4,1]
    if dim > 3:
        raise AttributeError(f"Dim is {dim}: too big!!!")

    distribution = []
    lambdas = []
    coords = []
    go = -1
    while go < 0:
        row = data[go]
        key = row[0]
        if key[0] == "L":
            lambdas.append(row[1:])
            go -= 1
        elif key[0] == "N":
            distribution.append([key[-1], row[1]])
            go -= 1
        elif key == "r_S":
            r_S = row[1]
            go -= 1
        elif key == "Coordinates":
            break
        else:
            coords.insert[0, row]
    
    if data[5,0] == "Matrix":
        cmatrix = data[6:6+size]
        cmatrix2 = np.matmul(cmatrix, cmatrix)
        for i in range(size):
            fut_links.append[[]]
            for j in range(size):
                if (cmatrix[i,j] and cmatrix2[i,j]==0):
                    fut_links[i].append[j]
    else:
        lines= [[v for v in line.split(",")] for line in open(lambdasfile_ext)]
        fut_links = [[int(v) for v in line if v != "\n"] 
                        for line in lines[6+3*size+3:6+4*size+3]]
    
    rc = 2
    r  = 2
    c  = 1
    fig, axs  = plt.figure(figsize = (10, 15), tight_layout = True)
    plt.suptitle( "Lambdas Information" )
    lambdas_colors = ["green", "yellow", "orange", "red", "purple"]
    Ncolors = len(lambdas_colors)
    
    plt.subplot(r, c, 1)
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")        
        
        #START WITH LAMBDAS
        lambdas_labels = []
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            lambdas_labels.append(uplabel)
            plt.plot(xup, tup, marker="o", 
                        markersize=20, 
                        markeredgecolor="red", 
                        markerfacecolor="cs:black",
                        zorder = 2)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                lambdas_labels.append(label)
                plt.plot(x, t, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [tup, t], 
                        markerfacecolor=facecolor,
                        zorder = 3)
        plt.vlines(r_S)

        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i]
                xi = coords[i]
                plt.plot(xi, ti, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor="cyan", zorder = 4)
            for j in fut_links[i]:
                tj = coords[j]
                xj = coords[j]
                plt.plot([xi, xj], [ti, tj], 
                        markerfacecolor="cyan",
                        zorder = 5)

    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")

        #START WITH LAMBDAS
        lambdas_labels = []
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            yup = coords[uplabel][2]
            lambdas_labels.append(uplabel)
            plt.plot(xup, tup, marker="o", 
                        markersize=20, 
                        markeredgecolor="red", 
                        markerfacecolor="cs:black",
                        zorder = 2)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                y = coords[label][2]
                lambdas_labels.append(label)
                plt.plot(x, y, t, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [yup, y], [tup, t], 
                        markerfacecolor=facecolor,
                        zorder = 3)
        
        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i]
                xi = coords[i]
                yi = coords[i]
                plt.plot(xi, ti, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor="cyan", zorder = 4)
            for j in fut_links[i]:
                tj = coords[j]
                xj = coords[j]
                yj = coords[j]
                plt.plot([xi, xj], [yi, yj], [ti, tj], 
                        markerfacecolor="cyan",
                        zorder = 5)

    
    plt.subplot(r, c, 2)
    plt.xlabel("Lambdas' size (number of links)")
    plt.ylabel("Number of Occurrences")
    plt.step(distribution[:,0], distribution[:,1])
    props = dict(boxstyle='round', facecolor='wheat', 
                edgecolor = 'black', ls = '', alpha=0.5)
    plt.annotate (f"Tot = {sum(distribution[:,1])}", 
                  (0.8, 0.8), xycoords = "axes fraction", 
                  va='bottom', ha = 'left', bbox=props) 

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig, axs




def plot_causet_and_lambdas(lambdasfile_ext, savefile_ext = 0):
    """
    PLot the whole causet highlighting the lambdas.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to

    Returns
    -------
    plt.figure
    plt.axs
    """
    data = np.genfromtxt(lambdasfile_ext, delimiter=",")

    storage_option = data[0,0]
    size           = data[1,1]
    dim            = data[2,1]
    shapename      = data[3,1]
    spacetimename  = data[4,1]
    if dim > 3:
        raise AttributeError(f"Dim is {dim}: too big!!!")

    distribution = []
    lambdas = []
    coords = []
    go = -1
    while go < 0:
        row = data[go]
        key = row[0]
        if key[0] == "L":
            lambdas.append(row[1:])
            go -= 1
        elif key[0] == "N":
            distribution.append([key[-1], row[1]])
            go -= 1
        elif key == "r_S":
            r_S = row[1]
            go -= 1
        elif key == "Coordinates":
            break
        else:
            coords.insert[0, row]
    
    if data[5,0] == "Matrix":
        cmatrix = data[6:6+size]
        cmatrix2 = np.matmul(cmatrix, cmatrix)
        for i in range(size):
            fut_links.append[[]]
            for j in range(size):
                if (cmatrix[i,j] and cmatrix2[i,j]==0):
                    fut_links[i].append[j]
    else:
        lines= [[v for v in line.split(",")] for line in open(lambdasfile_ext)]
        fut_links = [[int(v) for v in line if v != "\n"] 
                        for line in lines[6+3*size+3:6+4*size+3]]
    

    fig, axs  = plt.figure(figsize = (10, 15), tight_layout = True)
    lambdas_colors = ["green", "yellow", "orange", "red", "purple"]
    Ncolors = len(lambdas_colors)
    
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")        
        
        #START WITH LAMBDAS
        lambdas_labels = []
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            lambdas_labels.append(uplabel)
            plt.plot(xup, tup, marker="o", 
                        markersize=20, 
                        markeredgecolor="red", 
                        markerfacecolor="cs:black",
                        zorder = 2)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                lambdas_labels.append(label)
                plt.plot(x, t, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [tup, t], 
                        markerfacecolor=facecolor,
                        zorder = 3)
        plt.vlines(r_S)

        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i]
                xi = coords[i]
                plt.plot(xi, ti, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor="cyan", zorder = 4)
            for j in fut_links[i]:
                tj = coords[j]
                xj = coords[j]
                plt.plot([xi, xj], [ti, tj], 
                        markerfacecolor="cyan",
                        zorder = 5)

    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")

        #START WITH LAMBDAS
        lambdas_labels = []
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            yup = coords[uplabel][2]
            lambdas_labels.append(uplabel)
            plt.plot(xup, tup, marker="o", 
                        markersize=20, 
                        markeredgecolor="red", 
                        markerfacecolor="cs:black",
                        zorder = 2)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                y = coords[label][2]
                lambdas_labels.append(label)
                plt.plot(x, y, t, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [yup, y], [tup, t], 
                        markerfacecolor=facecolor,
                        zorder = 3)
        
        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i]
                xi = coords[i]
                yi = coords[i]
                plt.plot(xi, ti, marker="o", 
                        markersize=20, 
                        markeredgecolor="black", 
                        markerfacecolor="cyan", zorder = 4)
            for j in fut_links[i]:
                tj = coords[j]
                xj = coords[j]
                yj = coords[j]
                plt.plot([xi, xj], [yi, yj], [ti, tj], 
                        markerfacecolor="cyan",
                        zorder = 5)

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig, axs
    

    