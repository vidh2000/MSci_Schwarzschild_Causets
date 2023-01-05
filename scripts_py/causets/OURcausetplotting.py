#!/usr/bin/env python
'''
Created on 4 Jan 2022

@author: Stefano Veroni, Vid Homsak
'''
import numpy as np
import matplotlib.pyplot as plt
import functools

#from from https://stackoverflow.com/questions/15301999/default-arguments-with-args-and-kwargs
def default_kwargs(**defaultKwargs):
    def actual_decorator(fn):
        @functools.wraps(fn)
        def g(*args, **kwargs):
            defaultKwargs.update(kwargs)
            return fn(*args, **defaultKwargs)
        return g
    return actual_decorator

my_kwargs = {"figsize"   :(7.5, 15), 
            "link_alpha":0.5, 
            "link_lw"   :0.5,
            "markersize":5,
            "std_color" : "#002147",
            "lambdas_colors": ["green", "gold", "darkorange", "red", "purple"]}

ColorSchemes = {
    'matplotlib':            {'core':        'tab:blue',
                              'black':       'black',
                              'gray':        'tab:gray',
                              'grey':        'tab:gray',
                              'white':       'snow',
                              'purple':      'tab:purple',
                              'blue':        'tab:blue',
                              'cyan':        'tab:cyan',
                              'green':       'tab:green',
                              'lime':        'limegreen',
                              'yellow':      'gold',
                              'orange':      'tab:orange',
                              'red':         'tab:red',
                              'pink':        'tab:pink'},
    # ---------------------------------------------------------------
    # Imperial College London, UK. Brand style added on 19/08/2020
    # http://www.imperial.ac.uk/brand-style-guide/visual-identity/brand-colours/
    'ImperialLondon':        {'core':        '#002147',
                              'black':       '#002147',
                              'gray':        '#EBEEEE',
                              'grey':        '#EBEEEE',
                              'coolgrey':    '#9D9D9D',
                              'white':       '#D4EFFC',
                              'violet':      '#960078',
                              'iris':        '#751E66',
                              'purple':      '#653098',
                              'plum':        '#321E6D',
                              'navy':        '#002147',
                              'darkblue':    '#003E74',
                              'blue':        '#006EAF',
                              'cyan':        '#009CBC',
                              'green':       '#02893B',
                              'kermitgreen': '#66A40A',
                              'lime':        '#BBCEOO',
                              'yellow':      '#FFDD00',
                              'tangerine':   '#EC7300',
                              'orange':      '#D24000',
                              'cherry':      '#E40043',
                              'red':         '#DD2501',
                              'brick':       '#A51900',
                              'pink':        '#C81E78',
                              'raspberry':   '#9F004E'}
}


@default_kwargs(**my_kwargs)
def plot_causet(file_ext, savefile_ext = 0, **plot_kwargs):
    """
    PLot the causet. Highlight horizon if BlackHole.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    file_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to
    
    **plot_kwargs: These include
    - figsize. Default (7.5, 15).
    - link_alpha. Default 0.5.
    - link_lw. Default 0.5.
    - markersize. Default 5.
    - std_color. Default "#002147" (kind of dark blue).

    Returns
    -------
    plt.figure
    """
    #GET KWARGS
    figsize    = plot_kwargs['figsize']
    link_alpha = plot_kwargs['link_alpha']
    link_lw    = plot_kwargs['link_lw']
    markersize = plot_kwargs['markersize']
    std_color  = plot_kwargs['std_color']


    with open(file_ext, 'r') as fl: 
        f = fl.readlines() 

        storage_option = str(f[0].split(",")[1])
        size           = int(f[1].split(",")[1])
        dim            = int(f[2].split(",")[1])
        shapename      = str(f[3].split(",")[1])
        spacetimename  = str(f[4].split(",")[1])
        if dim > 3:
            raise AttributeError(f"Dim is {dim}: too big!!!")

        coords = []
        r_S = 0
        go = -1
        while go < 0:
            row = f[go].split(",")
            key = row[0]
            if key[0] == "L":
                go -= 1
            elif key[0] == "N":
                go -= 1
            elif key == "r_S" or key == "r_S\n":
                r_S += float(row[1])
                go -= 1
            elif key == "Coordinates" or key == "Coordinates\n":
                break
            elif key == "":
                go -= 1
            else:
                coords.insert(0, [])
                for i in range(dim):
                    coords[0].append(float(row[i]))
                go -= 1
        if not r_S:
            r_S += 2
    
        if storage_option == "cmatrix" or storage_option == "cmatrix\n":
            cmatrix = []
            for i in range(6, 6+size):
                cmatrix.append(f[i].split(","))
            cmatrix = np.array(cmatrix, dtype = 'int')
            cmatrix2 = np.matmul(cmatrix, cmatrix)
            for i in range(size):
                fut_links.append[[]]
                for j in range(size):
                    if (cmatrix[i,j] and cmatrix2[i,j]==0):
                        fut_links[i].append[j]
        else:
            lines= [[v for v in line.split(",")] for line in open(file_ext)]
            fut_links = [[int(v) for v in line if (v != "\n" and v != "")] 
                            for line in lines[6+3*size+3:6+4*size+3]]
    

    fig = plt.figure(figsize = figsize, tight_layout = True)
    
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")        

        for i in range(size):
            ti = coords[i][0]
            xi = coords[i][1]
            for j in [i] + fut_links[i]:
                tj = coords[j][0]
                xj = coords[j][1]
                plt.plot([xi, xj], [ti, tj], 
                            marker = "o", markersize = markersize, 
                            markeredgecolor = "black", 
                            markerfacecolor = std_color,
                            ls = "solid", color = std_color, 
                            alpha = link_alpha, lw = link_lw,
                            zorder = 5)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()
        plt.vlines(r_S, ys[0], ys[1], ls = "--", color = "red")
        
        if min(np.array(coords)[:,1]) < 0:
            plt.vlines(-r_S, ys[0], ys[1], ls = "--",color="r")
        plt.ylim(ys)

    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")
        
        for i in range(size):
            ti = coords[i][0]
            xi = coords[i][1]
            yi = coords[i][2]
            for j in [i] + fut_links[i]:
                tj = coords[j][0]
                xj = coords[j][1]
                yj = coords[j][2]
                plt.plot([xi, xj], [yi, yj], [ti, tj], 
                            marker = "o", markersize = markersize, 
                            markeredgecolor = "black", 
                            markerfacecolor = std_color,
                            ls = "solid", color = std_color, 
                            alpha = link_alpha, lw = link_lw,
                            zorder = 5)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()

        h = ys[1] - ys[0]
        cz = (ys[0]+ys[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        plt.plot_surface(Xc, Yc, Zc, alpha=0.2)

        plt.ylim(ys)

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig



@default_kwargs(**my_kwargs)
def plot_lambdas_and_lambdasdistr(lambdasfile_ext, savefile_ext = 0, 
                                    **plot_kwargs):
    """
    PLot the lambdas of a causet with Horizon and their distribution.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to
    
    **plot_kwargs: These include
    - figsize. Default (7.5, 15).
    - link_alpha. Default 0.5.
    - link_lw. Default 0.5.
    - markersize. Default 5.
    - std_color. Default "#002147" (kind of dark blue).

    Returns
    -------
    plt.figure
    """
    #GET KWARGS
    figsize    = plot_kwargs['figsize']
    link_alpha = plot_kwargs['link_alpha']
    link_lw    = plot_kwargs['link_lw']
    markersize = plot_kwargs['markersize']
    std_color  = plot_kwargs['std_color']
    lambdas_colors = plot_kwargs['lambdas_colors']


    with open(lambdasfile_ext, 'r') as fl: 
        f = fl.readlines() 

        storage_option = str(f[0].split(",")[1])
        size           = int(f[1].split(",")[1])
        dim            = int(f[2].split(",")[1])
        shapename      = str(f[3].split(",")[1])
        spacetimename  = str(f[4].split(",")[1])
        if dim > 3:
            raise AttributeError(f"Dim is {dim}: too big!!!")

        distribution = []
        lambdas = []
        coords = []
        go = -1
        while go < 0:
            row = f[go].split(",")
            key = row[0]
            if key[0] == "L":
                lambdas.append([])
                for label in row[1:]:
                    if label != "" and label != "\n":
                        lambdas[-1].append(int(label))
                go -= 1
            elif key[0] == "N":
                distribution.append([int(key[-1]), int(row[1])])
                go -= 1
            elif key == "r_S" or key == "r_S\n":
                r_S = float(row[1])
                go -= 1
            elif key == "Coordinates" or key == "Coordinates\n":
                break
            elif key == "":
                go -= 1
            else:
                coords.insert(0, [])
                for i in range(dim):
                    coords[0].append(float(row[i]))
                go -= 1


    rc = 2
    r  = 2
    c  = 1
    figsize = (c * 4, r * 2.2)
    fig  = plt.figure(figsize = (10, 15), tight_layout = True)
    Ncolors = len(lambdas_colors)
    
    plt.subplot(r, c, 1)
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")        
        
        #START WITH LAMBDAS
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size-1]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            plt.plot(xup, tup, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor,
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot(x, t, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [tup, t], ls = "solid", 
                        color=facecolor, alpha = link_alpha, lw = link_lw,
                        zorder = 5)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()
        plt.vlines(r_S, ys[0], ys[1], ls = "--", color = "red")
        
        if min(np.array(coords)[:,1]) < 0:
            plt.vlines(-r_S, ys[0], ys[1], ls = "--",color="r")
        plt.ylim(ys)

    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")

        #START WITH LAMBDAS
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
                        markersize = markersize, 
                        markeredgecolor="red", 
                        markerfacecolor="black",
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                y = coords[label][2]
                plt.plot(x, y, t, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 10)
                plt.plot([xup, x], [yup, y], [tup, t], ls = "solid", 
                        markerfacecolor=facecolor,
                        zorder = 8)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()
        
        h = ys[1] - ys[0]
        cz = (ys[0]+ys[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        plt.plot_surface(Xc, Yc, Zc, alpha=0.2)

        plt.ylim(ys)

    
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

    return fig

@default_kwargs(**my_kwargs)
def plot_lambdas(lambdasfile_ext, savefile_ext = 0, **plot_kwargs):
    """
    PLot the lambdas of a causet with Horizon.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to

    **plot_kwargs: These include
    - figsize. Default (7.5, 15).
    - link_alpha. Default 0.5.
    - link_lw. Default 0.5.
    - markersize. Default 5.
    - std_color. Default "#002147" (kind of dark blue).

    Returns
    -------
    plt.figure
    """
    #GET KWARGS
    figsize    = plot_kwargs['figsize']
    link_alpha = plot_kwargs['link_alpha']
    link_lw    = plot_kwargs['link_lw']
    markersize = plot_kwargs['markersize']
    std_color  = plot_kwargs['std_color']
    lambdas_colors = plot_kwargs['lambdas_colors']

    with open(lambdasfile_ext, 'r') as fl: 
        f = fl.readlines() 

        storage_option = str(f[0].split(",")[1])
        size           = int(f[1].split(",")[1])
        dim            = int(f[2].split(",")[1])
        shapename      = str(f[3].split(",")[1])
        spacetimename  = str(f[4].split(",")[1])
        if dim > 3:
            raise AttributeError(f"Dim is {dim}: too big!!!")

        distribution = []
        lambdas = []
        coords = []
        go = -1
        while go < 0:
            row = f[go].split(",")
            key = row[0]
            if key[0] == "L":
                lambdas.append([])
                for label in row[1:]:
                    if label != "" and label != "\n":
                        lambdas[-1].append(int(label))
                go -= 1
            elif key[0] == "N":
                distribution.append([int(key[-1]), int(row[1])])
                go -= 1
            elif key == "r_S" or key == "r_S\n":
                r_S = float(row[1])
                go -= 1
            elif key == "Coordinates" or key == "Coordinates\n":
                break
            elif key == "":
                go -= 1
            else:
                coords.insert(0, [])
                for i in range(dim):
                    coords[0].append(float(row[i]))
                go -= 1
    

    fig = plt.figure(figsize = figsize, tight_layout = True)
    Ncolors = len(lambdas_colors)
    
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")        
        
        #START WITH LAMBDAS
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size-1]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            plt.plot(xup, tup, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor,
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot(x, t, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [tup, t], ls = "solid", 
                        color=facecolor, alpha = link_alpha, lw = link_lw,
                        zorder = 5)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()
        plt.vlines(r_S, ys[0], ys[1], ls = "--", color = "red")
        
        if min(np.array(coords)[:,1]) < 0:
            plt.vlines(-r_S, ys[0], ys[1], ls = "--",color="r")
        plt.ylim(ys)

    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")

        #START WITH LAMBDAS
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
                        markersize = markersize, 
                        markeredgecolor="red", 
                        markerfacecolor="black",
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                y = coords[label][2]
                plt.plot(x, y, t, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 10)
                plt.plot([xup, x], [yup, y], [tup, t], ls = "solid", 
                        markerfacecolor=facecolor,
                        zorder = 8)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()
        
        h = ys[1] - ys[0]
        cz = (ys[0]+ys[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        plt.plot_surface(Xc, Yc, Zc, alpha=0.2)

        plt.ylim(ys)

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig



@default_kwargs(**my_kwargs)
def plot_causet_and_lambdas_and_lambdasdistr(lambdasfile_ext, 
                                            savefile_ext = 0,
                                            **plot_kwargs):
    """
    PLot the whole causet highlighting the lambdas and their distribution.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to

    **plot_kwargs: These include
    - figsize. Default (7.5, 15).
    - link_alpha. Default 0.5.
    - link_lw. Default 0.5.
    - markersize. Default 5.
    - std_color. Default "#002147" (kind of dark blue).

    Returns
    -------
    plt.figure
    """
    #GET KWARGS
    figsize    = plot_kwargs['figsize']
    link_alpha = plot_kwargs['link_alpha']
    link_lw    = plot_kwargs['link_lw']
    markersize = plot_kwargs['markersize']
    std_color  = plot_kwargs['std_color']
    lambdas_colors = plot_kwargs['lambdas_colors']

    with open(lambdasfile_ext, 'r') as fl: 
        f = fl.readlines() 

        storage_option = str(f[0].split(",")[1])
        size           = int(f[1].split(",")[1])
        dim            = int(f[2].split(",")[1])
        shapename      = str(f[3].split(",")[1])
        spacetimename  = str(f[4].split(",")[1])
        if dim > 3:
            raise AttributeError(f"Dim is {dim}: too big!!!")

        distribution = []
        lambdas = []
        coords = []
        go = -1
        while go < 0:
            row = f[go].split(",")
            key = row[0]
            if key[0] == "L":
                lambdas.append([])
                for label in row[1:]:
                    if label != "" and label != "\n":
                        lambdas[-1].append(int(label))
                go -= 1
            elif key[0] == "N":
                distribution.append([int(key[-1]), int(row[1])])
                go -= 1
            elif key == "r_S" or key == "r_S\n":
                r_S = float(row[1])
                go -= 1
            elif key == "Coordinates" or key == "Coordinates\n":
                break
            elif key == "":
                go -= 1
            else:
                coords.insert(0, [])
                for i in range(dim):
                    coords[0].append(float(row[i]))
                go -= 1
    
        if storage_option == "cmatrix" or storage_option == "cmatrix\n":
            cmatrix = []
            for i in range(6, 6+size):
                cmatrix.append(f[i].split(","))
            cmatrix = np.array(cmatrix, dtype = 'int')
            cmatrix2 = np.matmul(cmatrix, cmatrix)
            for i in range(size):
                fut_links.append[[]]
                for j in range(size):
                    if (cmatrix[i,j] and cmatrix2[i,j]==0):
                        fut_links[i].append[j]
        else:
            lines= [[v for v in line.split(",")] for line in open(lambdasfile_ext)]
            fut_links = [[int(v) for v in line if (v != "\n" and v != "")] 
                            for line in lines[6+3*size+3:6+4*size+3]]
    
    rc = 2
    r  = 2
    c  = 1
    fig = plt.figure(figsize = figsize, tight_layout = True)
    Ncolors = len(lambdas_colors)
    
    plt.subplot(r, c, 1)
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")        
        
        #START WITH LAMBDAS
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size-1]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            plt.plot(xup, tup, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor,
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot(x, t, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [tup, t], ls = "solid", 
                        color=facecolor, alpha = link_alpha, lw = link_lw,
                        zorder = 5)

        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        lambdas_labels = np.array(lambdas).flatten()
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i][0]
                xi = coords[i][1]
                for j in [i] + fut_links[i]:
                    tj = coords[j][0]
                    xj = coords[j][1]
                    plt.plot([xi, xj], [ti, tj], 
                                marker = "o", markersize = markersize, 
                                markeredgecolor = "black", 
                                markerfacecolor = std_color,
                                ls = "solid", color = std_color, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 5)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()
        plt.vlines(r_S, ys[0], ys[1], ls = "--", color = "red")
        
        if min(np.array(coords)[:,1]) < 0:
            plt.vlines(-r_S, ys[0], ys[1], ls = "--",color="r")
        plt.ylim(ys)

    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")

        #START WITH LAMBDAS
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
                        markersize = markersize, 
                        markeredgecolor="red", 
                        markerfacecolor="black",
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                y = coords[label][2]
                plt.plot(x, y, t, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 10)
                plt.plot([xup, x], [yup, y], [tup, t], ls = "solid", 
                        markerfacecolor=facecolor,
                        zorder = 8)
        
        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        lambdas_labels = np.array(lambdas).flatten()
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i][0]
                xi = coords[i][1]
                yi = coords[i][2]
                for j in [i] + fut_links[i]:
                    tj = coords[j][0]
                    xj = coords[j][1]
                    yj = coords[j][2]
                    plt.plot([xi, xj], [yi, yj], [ti, tj], 
                                marker = "o", markersize = markersize, 
                                markeredgecolor = "black", 
                                markerfacecolor = std_color,
                                ls = "solid", color = std_color, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 5)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()

        h = ys[1] - ys[0]
        cz = (ys[0]+ys[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        plt.plot_surface(Xc, Yc, Zc, alpha=0.2)

        plt.ylim(ys)

    
    plt.subplot(r, c, 2)
    plt.xlabel("Lambdas' size (number of links)")
    plt.ylabel("Number of Occurrences")
    distribution = np.array(distribution)
    plt.step(distribution[:,0], distribution[:,1])
    props = dict(boxstyle='round', facecolor='wheat', 
                edgecolor = 'black', ls = '', alpha=0.5)
    plt.annotate (f"Tot = {sum(distribution[:,1])}", 
                  (0.8, 0.8), xycoords = "axes fraction", 
                  va='bottom', ha = 'left', bbox=props) 

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig

@default_kwargs(**my_kwargs)
def plot_causet_and_lambdas(lambdasfile_ext, savefile_ext = 0, **plot_kwargs):
    """
    PLot the whole causet highlighting the lambdas.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to

    **plot_kwargs: These include
    - figsize. Default (7.5, 15).
    - link_alpha. Default 0.5.
    - link_lw. Default 0.5.
    - markersize. Default 5.
    - std_color. Default "#002147" (kind of dark blue).

    Returns
    -------
    plt.figure
    """
    #GET KWARGS
    figsize    = plot_kwargs['figsize']
    link_alpha = plot_kwargs['link_alpha']
    link_lw    = plot_kwargs['link_lw']
    markersize = plot_kwargs['markersize']
    std_color  = plot_kwargs['std_color']
    lambdas_colors = plot_kwargs['lambdas_colors']

    with open(lambdasfile_ext, 'r') as fl: 
        f = fl.readlines() 

        storage_option = str(f[0].split(",")[1])
        size           = int(f[1].split(",")[1])
        dim            = int(f[2].split(",")[1])
        shapename      = str(f[3].split(",")[1])
        spacetimename  = str(f[4].split(",")[1])
        if dim > 3:
            raise AttributeError(f"Dim is {dim}: too big!!!")

        distribution = []
        lambdas = []
        coords = []
        go = -1
        while go < 0:
            row = f[go].split(",")
            key = row[0]
            if key[0] == "L":
                lambdas.append([])
                for label in row[1:]:
                    if label != "" and label != "\n":
                        lambdas[-1].append(int(label))
                go -= 1
            elif key[0] == "N":
                distribution.append([int(key[-1]), int(row[1])])
                go -= 1
            elif key == "r_S" or key == "r_S\n":
                r_S = float(row[1])
                go -= 1
            elif key == "Coordinates" or key == "Coordinates\n":
                break
            elif key == "":
                go -= 1
            else:
                coords.insert(0, [])
                for i in range(dim):
                    coords[0].append(float(row[i]))
                go -= 1
    
        if storage_option == "cmatrix" or storage_option == "cmatrix\n":
            cmatrix = []
            for i in range(6, 6+size):
                cmatrix.append(f[i].split(","))
            cmatrix = np.array(cmatrix, dtype = 'int')
            cmatrix2 = np.matmul(cmatrix, cmatrix)
            for i in range(size):
                fut_links.append[[]]
                for j in range(size):
                    if (cmatrix[i,j] and cmatrix2[i,j]==0):
                        fut_links[i].append[j]
        else:
            lines= [[v for v in line.split(",")] for line in open(lambdasfile_ext)]
            fut_links = [[int(v) for v in line if (v != "\n" and v != "")] 
                            for line in lines[6+3*size+3:6+4*size+3]]
    

    fig = plt.figure(figsize = figsize, tight_layout = True)
    Ncolors = len(lambdas_colors)
    
    if dim == 2:
        plt.ylabel("t*")
        plt.xlabel("r")        
        
        #START WITH LAMBDAS
        for lambda_i in lambdas:
            #Set Color based on Size
            lambdas_size = len(lambda_i)
            if lambdas_size >= Ncolors:
                lambdas_size = Ncolors-1
            facecolor = lambdas_colors[lambdas_size-1]
            #Plot
            uplabel = lambda_i[0]
            tup = coords[uplabel][0]
            xup = coords[uplabel][1]
            plt.plot(xup, tup, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor,
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot(x, t, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 2)
                plt.plot([xup, x], [tup, t], ls = "solid", 
                        color=facecolor, alpha = link_alpha, lw = link_lw,
                        zorder = 5)

        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        lambdas_labels = np.array(lambdas).flatten()
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i][0]
                xi = coords[i][1]
                for j in [i] + fut_links[i]:
                    tj = coords[j][0]
                    xj = coords[j][1]
                    plt.plot([xi, xj], [ti, tj], 
                                marker = "o", markersize = markersize, 
                                markeredgecolor = "black", 
                                markerfacecolor = std_color,
                                ls = "solid", color = std_color, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 5)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()
        plt.vlines(r_S, ys[0], ys[1], ls = "--", color = "red")
        
        if min(np.array(coords)[:,1]) < 0:
            plt.vlines(-r_S, ys[0], ys[1], ls = "--",color="r")
        plt.ylim(ys)

    if dim == 3:
        plt.zlabel("t*")
        plt.xlabel("x")
        plt.ylabel("y")

        #START WITH LAMBDAS
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
                        markersize = markersize, 
                        markeredgecolor="red", 
                        markerfacecolor="black",
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                y = coords[label][2]
                plt.plot(x, y, t, marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor, zorder = 10)
                plt.plot([xup, x], [yup, y], [tup, t], ls = "solid", 
                        markerfacecolor=facecolor,
                        zorder = 8)
        
        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        lambdas_labels = np.array(lambdas).flatten()
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i][0]
                xi = coords[i][1]
                yi = coords[i][2]
                for j in [i] + fut_links[i]:
                    tj = coords[j][0]
                    xj = coords[j][1]
                    yj = coords[j][2]
                    plt.plot([xi, xj], [yi, yj], [ti, tj], 
                                marker = "o", markersize = markersize, 
                                markeredgecolor = "black", 
                                markerfacecolor = std_color,
                                ls = "solid", color = std_color, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 5)
        
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()

        h = ys[1] - ys[0]
        cz = (ys[0]+ys[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        plt.plot_surface(Xc, Yc, Zc, alpha=0.2)

        plt.ylim(ys)

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig






##################################################################
# HELPERS ########################################################
##################################################################


# from https://stackoverflow.com/questions/26989131/add-cylinder-to-plot 
def cylinder_along_z(center_x,center_y,center_z,radius,height_z):
    z = np.linspace(0, height_z, 50) + center_z
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid



    

    