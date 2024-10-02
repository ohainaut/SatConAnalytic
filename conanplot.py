#!/usr/bin/env python3
'''SatConAnalytic - Satellite Constellation Analytic simulations

ploting functions supporting conAn 
'''
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import conan as ca
import telescopes

from matplotlib.colors import ListedColormap, LinearSegmentedColormap
colors = ["black", "lawngreen", "yellow", "orange", "red", "darkred"]

# color map for losses
gyrd = LinearSegmentedColormap.from_list("mycmap", colors)

# colors for airmass
colors = ["darkred","red","orange","yellow","green","lawngreen"]
cairmap = LinearSegmentedColormap.from_list("mycmap", colors)

# colors for sun altitude
colors = ['darkblue','b','deepskyblue','paleturquoise']
csunmap = LinearSegmentedColormap.from_list("mycmap", colors)


#------------------------------------------------------------------------------

def getTelescope(myargs):
    '''return the Telescope object requested by myargs.code
    
    overloaded with other myargs specifications.

    myargs is an argparse object.
    '''
   
    if myargs.code is None: 
        myargs.code = 'DEFAULT'
        
    myTel = findTelescope(myargs.code)

    
    for what in [
        'expt',
        'fovl',
        'magbloom',
        'maglim',
        'resol',
        'trailf',
        'lat',
        'telescope',
        'instrument'
    ]:
        if myargs.__dict__[what] is not None:
            myTel.__dict__[what] = float(myargs.__dict__[what])

    if myargs.fovw is not None:
        myTel.fovw = myargs.fovw
    else:
        myTel.fovw = myTel.fovw *1.

    return myTel

def findTelescope(telinslabel):
    allTel = telescopes.readTelescopes()
    try:
        return allTel.byCode[telinslabel]
    except KeyError:
        raise KeyError(f'{telinslabel} not in {allTel.list}')





#------------------------------------------------------------------------------

def initPolPlot(ax):
    '''initialize a polar plot'''
    
    plt.rcParams.update({'font.size': 15})
    ax.set_theta_zero_location('N')  
    ax.set_rticks([30,60,70,80,90])
    ax.set_yticklabels([])
    ax.set_xticks(np.radians(range(0,360,90)))
    ax.set_xticklabels(['N','E','S','W'])
    ax.set_rlim(0.,90)
    #ax.set_facecolor("k")

#------------------------------------------------------------------------------

def plotOneDens(ax, AzEl, dens, label):
    '''plot a density function
    IN
    - ax: axis object where to plot
    - AzEl: array with the azimut and elevations
    - dens: the density array to be plotted
    - label: label...
    OUT: updated plot in ax.
    '''

    cmap = 'magma'
    clab = 'cyan'
    ccon = 'cyan'

    _ = initPolPlot(ax)
            
    vmin = np.amin(dens)
    vmax = np.amax(dens)

    if (vmax-vmin) >0:
        ccd = ax.contour(np.radians(AzEl[0]), 90.-AzEl[1], dens,
             colors=ccon , linewidths=0.5,
             locator=ticker.LogLocator(subs=(.2,.5,1.)))
        ax.clabel(ccd, colors=clab, fmt='%.2g')
    
        cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1], dens , 
            locator=ticker.LogLocator(subs=10.**np.arange(.1,1.0,.025)),
            cmap=cmap)
    
    ax.set_title(label)
    return

#------------------------------------------------------------------------------

def plotDens(AzEl, nplot, densities, labels, plotLabel):
    '''plot a series of density distributions'''
    
    if nplot == 1:
        fig = plt.figure(figsize=(8,8))
        ax =  fig.subplots(1,nplot,subplot_kw={'projection': 'polar'}) 
        dens = densities
        label = labels
        plotOneDens(ax, AzEl, dens, label)
    else:
        fig = plt.figure(figsize=(5*nplot,6))
        axs =  fig.subplots(1,nplot,subplot_kw={'projection': 'polar'}) 
        for ax, dens, label in zip(axs, densities, labels):
            plotOneDens(ax, AzEl, dens, label)
  
    #        cbar0 = fig.colorbar(cfd, ax=axDens)
    #        cbar0.set_label('Trail/exp.')
    
    fig.suptitle(plotLabel,  size=20)
    
    return

#------------------------------------------------------------------------------

def azlab(ax,x,y,lab,size=12,alpha=0):
    '''write labels at (x,y) on a polar plot in ax'''

    rlab = np.sqrt(x*x + y*y)*90.
    azlab = np.arctan2(-x,y)
       
    if x < 0 :
        halig = 'left'
    elif x == 0:
        halig = 'center'
    else:
        halig = 'right'
    
    t = ax.text(azlab,rlab,lab, horizontalalignment=halig,size=size)
    t.set_bbox(dict(alpha=alpha,facecolor='white',edgecolor='none'))
    return azlab,rlab

#------------------------------------------------------------------------------

def drawHADec(lat):
    '''draw the HA and Dec lines in an altaz plot'''

    corl = 0.3
    corc = 'r'

    #draw RA
    HAdef = np.arange(-12.,12.01)*15.
    Decdef = np.arange(-90.,91.,1.)
    for i in np.arange(0,len(HAdef)):
        azs,els = ca.radec2azel(HAdef[i], Decdef, lat)
        az = azs[els>0]
        el = els[els>0]
        plt.plot(np.radians(az), 90.- el, lw=corl, c=corc)
    
    #label RA
    Dec = 0.
    HAdef = np.arange(-12.,12.01,3.)
    az,el = ca.radec2azel(HAdef*15., Dec, lat)
    for i in np.arange(0,len(HAdef)):
        if el[i] > 0. :
            plt.text(np.radians(az[i]),90. -el[i],'{:.0f}'.format(HAdef[i]), color=corc, fontsize=10)

    #draw Dec
    HAdef = np.arange(-12.,12.01,0.1)*15.
    Decdef = np.arange(-90.,91.,10.)
    for i in np.arange(0,len(Decdef)):
        azs,els = ca.radec2azel(HAdef, Decdef[i], lat)
        az = azs[els>0]
        el = els[els>0]
        if Decdef[i] == 0.:
            lw = 1.
        else:
            lw = corl
        plt.plot(np.radians(az), 90.- el, lw=lw, c=corc)
    
    #label Dec
    HA = 0.
    Decdef = np.arange(-60.,61.,30.)
    az,el = ca.radec2azel(HA, Decdef, lat)
    for i in np.arange(0,len(Decdef)):
        if el[i] > 0. :
            plt.text(np.radians(az[i]),90.-el[i],'{:.0f}$^\circ$'.format(Decdef[i]), color=corc, fontsize=10)



