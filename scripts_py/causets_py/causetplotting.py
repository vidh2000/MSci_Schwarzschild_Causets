#!/usr/bin/env python

'''
Created on 4 Jan 2022

@author: Stefano Veroni, Vid Homsak
'''
import numpy as np
import matplotlib.pyplot as plt
from . import causet_helpers as ch


my_kwargs = {
            "figsize"   :(7.5, 15), 
            "lambda_link_alpha": 1,
            "link_alpha":0.5, 
            "link_lw"   :0.5,
            "markersize":5,
            "std_color" : "#002147",
            "lambdas_colors": ["red", 
                                "orange", 
                                "gold", 
                                "limegreen", 
                                "dodgerblue", 
                                "violet"]
            }

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


@ch.default_kwargs(**my_kwargs)
def plot_causet(file_ext, savefile_ext = 0, 
                phi_limits = (0, 2*3.1416), projection = False,
                **plot_kwargs):
    """
    PLot the causet. Highlight horizon if BlackHole.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    file_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to
    
    phi_limits : tuple (default is (0, 2\pi))
        Sets limits of phi values plotted for 3D plot.
    
    projection : bool (default false).
        Draw 3D plot as projection on 2D (t,r) space. Useful when phi_limits
        are narrow.

    
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


    #############################################################
    #GET CAUSET INFO
    storage_option, \
    size,           \
    dim,            \
    shapename,      \
    spacetimename,  \
    coords,         \
    fut_links,      \
    r_S, lambdas,   \
    distribution =  \
    ch.get_causet_attrs(file_ext)
    

    fig = plt.figure(figsize = figsize, tight_layout = True)
    ax = plt.axes()
    
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
    
    if dim == 3 and projection:
        coords, fut_links = take_projection(coords, fut_links, phi_limits)
        ax.set_xlabel("x")
        ax.set_ylabel("t*")

        for i in range(len(coords)):
            ti   = coords[i][0]
            ri   = coords[i][1]
            phii = coords[i][2]
            if phi_limits[0] < phii and phii < phi_limits[1]:
                xi = ri * np.cos(phii)
                yi = ri * np.sin(phii)
                for j in [i] + fut_links[i]:
                    tj = coords[j][0]
                    rj = coords[j][1]
                    phij = coords[j][2]
                    if phi_limits[0] < phij and phij < phi_limits[1]:
                        ax.plot([ri, rj], [ti, tj], 
                                    marker = "o", markersize = markersize, 
                                    markeredgecolor = "black", 
                                    markerfacecolor = std_color,
                                    ls = "solid", color = std_color, 
                                    alpha = link_alpha, lw = link_lw,
                                    zorder = 5)
        #FINALLY MARK THE HORIZON
        ys = ax.set_ylim()
        ax.vlines(r_S, ys[0], ys[1], ls = "--", color = "red")
        if min(np.array(coords)[:,1]) < 0:
            ax.vlines(-r_S, ys[0], ys[1], ls = "--",color="r")
        ax.set_ylim(ys)

    elif dim == 3:
        ax = plt.axes(projection = "3d")
        ax.set_zlabel("t*")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        
        for i in range(size):
            ti   = coords[i][0]
            ri   = coords[i][1]
            phii = coords[i][2]
            if phi_limits[0] < phii and phii < phi_limits[1]:
                xi = ri * np.cos(phii)
                yi = ri * np.sin(phii)
                meci = "red" if ri < r_S else "black"
                ax.plot([xi], [yi], [ti], 
                        marker = "o", markersize = markersize, 
                        markeredgecolor = meci, 
                        markerfacecolor = std_color,
                        zorder = 5)
                for j in [i] + fut_links[i]:
                    tj = coords[j][0]
                    rj = coords[j][1]
                    phij = coords[j][2]
                    if phi_limits[0] < phij and phij < phi_limits[1]:
                        xj = rj * np.cos(phij)
                        yj = rj * np.sin(phij)
                        mecj = "red" if rj < r_S else "black"
                        ax.plot([xj], [yj], [tj], 
                                marker = "o", markersize = markersize, 
                                markeredgecolor = mecj, 
                                markerfacecolor = std_color,
                                zorder = 5)
                        ax.plot([xi, xj], [yi, yj], [ti, tj], markersize = 0, 
                                ls = "solid", color = std_color, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 4)
        #FINALLY MARK THE HORIZON
        zs = ax.set_zlim()
        h = abs(zs[1] - zs[0])
        cz = (zs[0]+zs[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        ax.plot_surface(Xc, Yc, Zc, alpha=0.05, 
                        color = "#9F004E", label = "Horizon")
        ax.set_zlim(zs)
        xymax = max(list(ax.set_ylim())+list(ax.set_xlim()))
        ax.set_ylim(-xymax, xymax)
        ax.set_xlim(-xymax, xymax)

    if savefile_ext:
        plt.savefig(savefile_ext)
    #plt.show()

    return ax



@ch.default_kwargs(**my_kwargs)
def plot_lambdas_and_lambdasdistr(lambdasfile_ext, savefile_ext = 0, 
                                  phi_limits=(0, 2*3.1416), **plot_kwargs):
    """
    PLot the lambdas of a causet with Horizon and their distribution.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to
    
    phi_limits : tuple (default is (0, 2\pi))
        Sets limits of phi values plotted for 3D plot.
    
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


    #############################################################
    #GET CAUSET INFO
    storage_option, \
    size,           \
    dim,            \
    shapename,      \
    spacetimename,  \
    coords,         \
    fut_links,      \
    r_S, lambdas,   \
    distribution =  \
    ch.get_causet_attrs(lambdasfile_ext)


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
        ax = plt.axes(projection = "3d")
        ax.set_zlabel("t*")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

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
            rup = coords[uplabel][1]
            phiup = coords[uplabel][2]
            xup = rup * np.cos(phiup)
            yup = rup * np.sin(phiup)
            for label in lambda_i[1:]:
                t = coords[label][0]
                r = coords[label][1]
                phi = coords[label][2]
                if phi_limits[0] < phi and phi < phi_limits[1]:
                    x = r * np.cos(phi)
                    y = r * np.sin(phi)
                    if phi_limits[0] < phi and phi < phi_limits[1]:
                        ax.plot([xup, x], [yup, y], [tup, t], marker="o", 
                                markersize = markersize, 
                                markeredgecolor="black", 
                                markerfacecolor=facecolor,  
                                ls = "solid", color = std_color, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 5)
                    else:
                        ax.plot([x], [y], [t], marker="o", 
                                markersize = markersize, 
                                markeredgecolor="black", 
                                markerfacecolor=facecolor,
                                zorder = 5)
        
        #FINALLY MARK THE HORIZON
        zs = ax.set_zlim()
        h = abs(zs[1] - zs[0])
        cz = (zs[0]+zs[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        ax.plot_surface(Xc, Yc, Zc, alpha=0.05, 
                        color = "#9F004E", label = "Horizon")
        ax.set_zlim(zs)
        xymax = max(list(ax.set_ylim())+list(ax.set_xlim()))
        ax.set_ylim(-xymax, xymax)
        ax.set_xlim(-xymax, xymax)

    
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

@ch.default_kwargs(**my_kwargs)
def plot_lambdas(lambdasfile_ext, savefile_ext = 0, phi_limits = (0, 2*3.1416),
**plot_kwargs):
    """
    PLot the lambdas of a causet with Horizon.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to

    phi_limits : tuple (default is (0, 2\pi))

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

    #############################################################
    #GET CAUSET INFO
    storage_option, \
    size,           \
    dim,            \
    shapename,      \
    spacetimename,  \
    coords,         \
    fut_links,      \
    r_S, lambdas,   \
    distribution =  \
    ch.get_causet_attrs(lambdasfile_ext)
    

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
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot([xup, x], [tup, t], marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor,ls = "solid", 
                        color=facecolor, alpha = link_alpha, lw = link_lw,
                        zorder = 5)
        #FINALLY MARK THE HORIZON
        ys = plt.ylim()
        plt.vlines(r_S, ys[0], ys[1], ls = "--", color = "red")
        if min(np.array(coords)[:,1]) < 0:
            plt.vlines(-r_S, ys[0], ys[1], ls = "--",color="r")
        plt.ylim(ys)

    if dim == 3:
        ax = plt.axes(projection = "3d")
        ax.set_zlabel("t*")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

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
            rup = coords[uplabel][1]
            phiup = coords[uplabel][2]
            xup = rup * np.cos(phiup)
            yup = rup * np.sin(phiup)
            for label in lambda_i[1:]:
                t = coords[label][0]
                r = coords[label][1]
                phi = coords[label][2]
                if phi_limits[0] < phi and phi < phi_limits[1]:
                    x = r * np.cos(phi)
                    y = r * np.sin(phi)
                    if phi_limits[0] < phi and phi < phi_limits[1]:
                        ax.plot([xup, x], [yup, y], [tup, t], marker="o", 
                                markersize = markersize, 
                                markeredgecolor="black", 
                                markerfacecolor=facecolor,  
                                ls = "solid", color = std_color, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 5)
                    else:
                        ax.plot([x], [y], [t], marker="o", 
                                markersize = markersize, 
                                markeredgecolor="black", 
                                markerfacecolor=facecolor,
                                zorder = 5)
        #FINALLY MARK THE HORIZON
        zs = ax.set_zlim()
        h = abs(zs[1] - zs[0])
        cz = (zs[0]+zs[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        ax.plot_surface(Xc, Yc, Zc, alpha=0.05, 
                        color = "#9F004E", label = "Horizon")
        ax.set_zlim(zs)
        xymax = max(list(ax.set_ylim())+list(ax.set_xlim()))
        ax.set_ylim(-xymax, xymax)
        ax.set_xlim(-xymax, xymax)

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return fig



@ch.default_kwargs(**my_kwargs)
def plot_causet_and_lambdas_and_lambdasdistr(lambdasfile_ext, 
                                            savefile_ext = 0,
                                            phi_limits = (0, 2*3.1416),
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
    
    phi_limits : tuple (default is (0, 2*\pi))

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

    #############################################################
    #GET CAUSET INFO
    storage_option, \
    size,           \
    dim,            \
    shapename,      \
    spacetimename,  \
    coords,         \
    fut_links,      \
    r_S, lambdas,   \
    distribution =  \
    ch.get_causet_attrs(lambdasfile_ext)
    
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
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot([xup, x], [tup, t], marker="o", 
                        markersize = markersize, 
                        markeredgecolor="black", 
                        markerfacecolor=facecolor,ls = "solid", 
                        color=facecolor, alpha = link_alpha, lw = link_lw,
                        zorder = 5)

        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        lambdas_labels = np.array(lambdas).flatten()
        for i in range(size):
            if not (i in lambdas_labels):
                ti   = coords[i][0]
                ri   = coords[i][1]
                for j in [i] + fut_links[i]:
                    tj = coords[j][0]
                    rj = coords[j][1]
                    if not (j in lambdas_labels):
                        plt.plot([ri, rj], [ti, tj], 
                                    marker = "o", markersize = markersize, 
                                    markeredgecolor = "black", 
                                    markerfacecolor = std_color,
                                    ls = "solid", color = std_color, 
                                    alpha = link_alpha, lw = link_lw,
                                    zorder = 5)
                    else:
                        plt.plot([ri, rj], [ti, tj], 
                                    marker = "o", markersize = 0,
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
        ax = plt.axes(projection = "3d")
        ax.set_zlabel("t*")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

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
            rup = coords[uplabel][1]
            phiup = coords[uplabel][2]
            xup = rup * np.cos(phiup)
            yup = rup * np.sin(phiup)
            for label in lambda_i[1:]:
                t = coords[label][0]
                r = coords[label][1]
                phi = coords[label][2]
                if phi_limits[0] < phi and phi < phi_limits[1]:
                    x = r * np.cos(phi)
                    y = r * np.sin(phi)
                    if phi_limits[0] < phi and phi < phi_limits[1]:
                        ax.plot([xup, x], [yup, y], [tup, t], marker="o", 
                                markersize = markersize, 
                                markeredgecolor="black", 
                                markerfacecolor=facecolor,  
                                ls = "solid", color = std_color, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 5)
                    else:
                        ax.plot([x], [y], [t], marker="o", 
                                markersize = markersize, 
                                markeredgecolor="black", 
                                markerfacecolor=facecolor,
                                zorder = 5)
        
        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        lambdas_labels = np.array(lambdas).flatten()
        for i in range(size):
            if not (i in lambdas_labels):
                ti   = coords[i][0]
                ri   = coords[i][1]
                phii = coords[i][2]
                if phi_limits[0] < phii and phii < phi_limits[1]:
                    xi = ri * np.cos(phii)
                    yi = ri * np.sin(phii)
                    for j in [i] + fut_links[i]:
                        tj = coords[j][0]
                        rj = coords[j][1]
                        phij = coords[label][2]
                        if phi_limits[0] < phij and phij < phi_limits[1]:
                            xj = rj * np.cos(phij)
                            yj = rj * np.sin(phij)
                            if not (j in lambdas_labels):
                                ax.plot([xi, xj], [yi, yj], [ti, tj], 
                                            marker = "o", markersize = markersize, 
                                            markeredgecolor = "black", 
                                            markerfacecolor = std_color,
                                            ls = "solid", color = std_color, 
                                            alpha = link_alpha, lw = link_lw,
                                            zorder = 5)
                            else:
                                ax.plot([xi, xj], [yi, yj], [ti, tj], 
                                            marker = "o", markersize = 0,
                                            ls = "solid", color = std_color, 
                                            alpha = link_alpha, lw = link_lw,
                                            zorder = 5)
        
        #FINALLY MARK THE HORIZON
        zs = ax.set_zlim()
        h = abs(zs[1] - zs[0])
        cz = (zs[0]+zs[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        ax.plot_surface(Xc, Yc, Zc, alpha=0.05, 
                        color = "#9F004E", label = "Horizon")
        ax.set_zlim(zs)
        xymax = max(list(ax.set_ylim())+list(ax.set_xlim()))
        ax.set_ylim(-xymax, xymax)
        ax.set_xlim(-xymax, xymax)

    
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

@ch.default_kwargs(**my_kwargs)
def plot_causet_and_lambdas(lambdasfile_ext, savefile_ext = 0, 
                            phi_limits = (0, 2*3.1416),
                            projection = False,
                            **plot_kwargs):
    """
    PLot the whole causet highlighting the lambdas.
    Time on vertical axis. The spatial dimensions can be 2 or 3.

    Parameters
    ----------
    lambdasfile_ext : string
        Path + filename + ext

    savefile_ext : string (default is int 0)
        Path + filename + ext to save file to
    
    phi_limits : tuple (default is (0, 2*\pi))
    
    projection: bool (default is false)
        Plot the 3D case as if it was a 2D (t,r) plot. Useful when phi_limits
        pick a narrow interval.

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
    #############################################################
    #GET KWARGS
    figsize    = plot_kwargs['figsize']
    lambda_link_alpha = plot_kwargs['lambda_link_alpha']
    link_alpha = plot_kwargs['link_alpha']
    link_lw    = plot_kwargs['link_lw']
    markersize = plot_kwargs['markersize']
    std_color  = plot_kwargs['std_color']
    lambdas_colors = plot_kwargs['lambdas_colors']

    #############################################################
    #GET CAUSET INFO
    storage_option, \
    size,           \
    dim,            \
    shapename,      \
    spacetimename,  \
    coords,         \
    fut_links,      \
    r_S, lambdas,   \
    distribution =  \
    ch.get_causet_attrs(lambdasfile_ext)
    
    ###############################################################
    #PLOTTING
    fig = plt.figure(figsize = figsize, tight_layout = True)
    ax = plt.axes()
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
            plt.plot([xup], [tup], marker="o", 
                        markersize = markersize, 
                        markeredgecolor=facecolor, 
                        markerfacecolor=facecolor,
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                x = coords[label][1]
                plt.plot([x], [t], marker="o", 
                        markersize = markersize, 
                        markeredgecolor=facecolor, 
                        markerfacecolor=facecolor,
                        zorder = 10)
                plt.plot([xup, x], [tup, t], markersize = 0,
                        ls = "solid", color=facecolor, 
                        alpha = lambda_link_alpha, lw = link_lw,
                        zorder = 10)

        ##NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        lambdas_labels = [lbl for lambda_i in lambdas for lbl in lambda_i]
        for i in range(size):
            if not (i in lambdas_labels):
                ti = coords[i][0]
                xi = coords[i][1]
                for j in [i] + fut_links[i]:
                    tj = coords[j][0]
                    xj = coords[j][1]
                    if not (j in lambdas_labels):
                        plt.plot([xi, xj], [ti, tj], 
                                    marker = "o", markersize = markersize, 
                                    markeredgecolor = "black", 
                                    markerfacecolor = std_color,
                                    ls = "solid", color = std_color, 
                                    alpha = link_alpha, lw = link_lw,
                                    zorder = 5)
                    else:
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
    

    # if dim == 3 and projection:
    #     ax.set_xlabel("x")
    #     ax.set_ylabel("y")

    #     START WITH LAMBDAS
    #     for lambda_i in lambdas:
    #         Set Color based on Size
    #         lambdas_size = len(lambda_i)
    #         if lambdas_size >= Ncolors:
    #             lambdas_size = Ncolors-1
    #         facecolor = lambdas_colors[lambdas_size]
    #         Plot
    #         uplabel = lambda_i[0]
    #         tup = coords[uplabel][0]
    #         rup = coords[uplabel][1]
    #         phiup = coords[uplabel][2]
    #         xup = rup * np.cos(phiup)
    #         yup = rup * np.sin(phiup)
    #         if phi_limits[0] < phiup and phiup < phi_limits[1]:
    #             ax.plot([rup], [tup], marker="o", 
    #                     markersize = markersize, 
    #                     markeredgecolor="black", 
    #                     markerfacecolor=facecolor,
    #                     zorder = 5)
    #         for label in lambda_i[1:]:
    #             t = coords[label][0]
    #             r = coords[label][1]
    #             phi = coords[label][2]
    #             if phi_limits[0] < phi and phi < phi_limits[1]:
    #                 x = r * np.cos(phi)
    #                 y = r * np.sin(phi)
    #                 if phi_limits[0] < phiup and phiup < phi_limits[1]:
    #                     mecup = "red" if rup < r_S else "black"
    #                     mec = "red" if r < r_S else "black"
    #                     ax.plot([rup, r], [tup, t], marker="o", 
    #                             markersize = markersize, 
    #                             markeredgecolor=[mecup, mec], 
    #                             markerfacecolor=facecolor,  
    #                             ls = "solid", color = std_color, 
    #                             alpha = link_alpha, lw = link_lw,
    #                             zorder = 5)
    #                 else:
    #                     ax.plot([r], [t], marker="o", 
    #                             markersize = markersize, 
    #                             markeredgecolor=mec, 
    #                             markerfacecolor=facecolor,
    #                             zorder = 5)
        
    #     NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
    #     lambdas_labels = [lbl for lambda_i in lambdas for lbl in lambda_i]
    #     for i in range(size):
    #         if not (i in lambdas_labels):
    #             ti   = coords[i][0]
    #             ri   = coords[i][1]
    #             phii = coords[i][2]
    #             if phi_limits[0] < phii and phii < phi_limits[1]:
    #                 xi = ri * np.cos(phii)
    #                 yi = ri * np.sin(phii)
    #                 for j in [i] + fut_links[i]:
    #                     tj = coords[j][0]
    #                     rj = coords[j][1]
    #                     phij = coords[j][2]
    #                     if phi_limits[0] < phij and phij < phi_limits[1]:
    #                         xj = rj * np.cos(phij)
    #                         yj = rj * np.sin(phij)
    #                         if not (j in lambdas_labels):
    #                             ax.plot([ri, rj], [ti, tj], 
    #                                         marker = "o", markersize = markersize, 
    #                                         markeredgecolor = "black", 
    #                                         markerfacecolor = std_color,
    #                                         ls = "solid", color = std_color, 
    #                                         alpha = link_alpha, lw = link_lw,
    #                                         zorder = 5)
    #                         else:
    #                             ax.plot([ri, rj], [ti, tj], 
    #                                         marker = "o", markersize = 0,
    #                                         ls = "solid", color = std_color, 
    #                                         alpha = link_alpha, lw = link_lw,
    #                                         zorder = 5)
    #     FINALLY MARK THE HORIZON
    #     ys = plt.ylim()
    #     plt.vlines(r_S, ys[0], ys[1], ls = "--", color = "red")
    #     if min(np.array(coords)[:,1]) < 0:
    #         plt.vlines(-r_S, ys[0], ys[1], ls = "--",color="r")
    #     plt.ylim(ys)


    elif dim == 3:
        ax = plt.axes(projection = "3d")
        ax.set_zlabel("t*")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

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
            rup = coords[uplabel][1]
            phiup = coords[uplabel][2]
            xup = rup * np.cos(phiup)
            yup = rup * np.sin(phiup)
            if phi_limits[0] < phiup and phiup < phi_limits[1]:
                plt.plot([xup], [yup], [tup], marker="o", 
                        markersize = markersize, 
                        markeredgecolor="red", 
                        markerfacecolor=facecolor,
                        zorder = 10)
            for label in lambda_i[1:]:
                t = coords[label][0]
                r = coords[label][1]
                phi = coords[label][2]
                if phi_limits[0] < phi and phi < phi_limits[1]:
                    x = r * np.cos(phi)
                    y = r * np.sin(phi)
                    ax.plot([x], [y], [t], marker="o", 
                                markersize = markersize, 
                                markeredgecolor="black", 
                                markerfacecolor=facecolor,
                                zorder = 10)
                    if phi_limits[0] < phiup and phiup < phi_limits[1]:
                        ax.plot([xup, x], [yup, y], [tup, t], 
                                markersize = 0, 
                                ls = "solid", color = facecolor, 
                                alpha = link_alpha, lw = link_lw,
                                zorder = 9)
                        
        
        #NOW DO OTHER POINTS (SKIPPING THOSE IN LAMBDAS)
        lambdas_labels = [lbl for lambda_i in lambdas for lbl in lambda_i]
        for i in range(size):
            if not (i in lambdas_labels):
                ti   = coords[i][0]
                ri   = coords[i][1]
                phii = coords[i][2]
                if phi_limits[0] < phii and phii < phi_limits[1]:
                    xi = ri * np.cos(phii)
                    yi = ri * np.sin(phii)
                    for j in [i] + fut_links[i]:
                        tj = coords[j][0]
                        rj = coords[j][1]
                        phij = coords[j][2]
                        if phi_limits[0] < phij and phij < phi_limits[1]:
                            xj = rj * np.cos(phij)
                            yj = rj * np.sin(phij)
                            ax.plot([xi, xj], [yi, yj], [ti, tj], 
                                    markersize = 0, 
                                    ls = "solid", color = std_color, 
                                    alpha = link_alpha, lw = link_lw,
                                    zorder = 4)
                            if not (j in lambdas_labels):
                                mecj = "red" if rj < r_S else "black"
                                ax.plot([xj], [yj], [tj], 
                                        marker = "o", markersize = markersize, 
                                        markeredgecolor = mecj, 
                                        markerfacecolor = std_color,
                                        zorder = 5)
        #FINALLY MARK THE HORIZON
        zs = ax.set_zlim()
        h = abs(zs[1] - zs[0])
        cz = (zs[0]+zs[1])/2
        Xc,Yc,Zc = cylinder_along_z(0, 0, cz, r_S, h)
        ax.plot_surface(Xc, Yc, Zc, alpha=0.05, 
                        color = "#9F004E", label = "Horizon")
        ax.set_zlim(zs)
        xymax = max(list(ax.set_ylim())+list(ax.set_xlim()))
        ax.set_ylim(-xymax, xymax)
        ax.set_xlim(-xymax, xymax)

    if savefile_ext:
        plt.savefig(savefile_ext)
    plt.show()

    return ax






##################################################################
# HELPERS ########################################################
##################################################################


# from https://stackoverflow.com/questions/26989131/add-cylinder-to-plot 
def cylinder_along_z(center_x,center_y,center_z,radius,height_z):
    z = np.linspace(-(height_z/2), (height_z/2), 100) + center_z
    theta = np.linspace(0, 2*np.pi, 100)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

def point_satisfies_phi_limits(coord_i, phi_limits):
    phi = coord_i[2]
    return phi_limits[0] < phi and phi < phi_limits[1]

def take_projection(coords, futlinks, phi_limits):
    N = len(coords)
    k = 0
    count = 0
    while k < N:
        #print(k, ", ",N)
        coord_k = coords[k]
        if point_satisfies_phi_limits(coord_k, phi_limits):
            k += 1
        else:
            count += 1
            coords = np.delete(coords, k, axis = 0)
            #add future links of k to future links of i if i<k is link
            for i in range(k):
                for j in futlinks[i]:
                    if j==k: #if k is in future of i
                        for j_k in futlinks[k]:
                            if j_k not in futlinks[i]:
                                futlinks[i].append(j_k)
                        futlinks[i].remove(k)       
            #remove kth set of future links
            futlinks = np.delete(futlinks, k, axis = 0)
            #reduce of one all labels above k
            for i in range(N-1):
                for j_index in range(len(futlinks[i])):
                    if futlinks[i][j_index] > k:
                        futlinks[i][j_index] -= 1
            N -= 1
    return coords, futlinks



    

    