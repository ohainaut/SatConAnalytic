#!/usr/bin/env python3
# CatAn
# Constellation Analytic simulations
#
# plots functions
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

import ConAn.conan as ca

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

def findPreset(telinslabel):
    # load parameters for a preset telinslabel telescope/instrument
    
    
    maglimoff = 1.75 # 1.75: 1/5sig ; 1.: 2/5sig
    magbloomoff = -2.5

    fovw = -99
    
    # -------IMAGING
    if  telinslabel == "FlyEye":
        telescope   = 'ESA FlyEye'
        lat         = -29.25 #latitude
        instrument  = " "
        expt        = 40. #s
        fovl        = 6.5 #deg
        resol       = 1./3600. # deg
        maglim      = maglimoff + 14.7   # detection limit for expt
        magbloom    = magbloomoff + 4.3  # threshold magnitude for blinding instruments
        trailf = 5./(fovl*3600.)         #" wide trail

    if  telinslabel == "ZTF":
        telescope   = '1.2m Oschin'
        lat         = 33.3 #latitude
        instrument  = "ZTF "
        expt        = 30. #s
        fovl        = 7.7 #deg
        resol       = 1./3600. # deg
        maglim      = maglimoff + 20.5   # detection limit for expt
        magbloom    = magbloomoff + 4.3  # threshold magnitude for blinding instruments
        trailf      = 10./(fovl*3600.)         #" wide trail

        
    elif telinslabel == "Photo":
        telescope   = ' '
        lat         = 50 #latitude
        instrument  = "Camera"
        expt        = 30. #s
        fovl        = 75 #deg
        fovw        = 55 #deg
        resol       = 60./3600. # deg
        maglim      = maglimoff + 10.   # detection limit for expt
        magbloom    = magbloomoff + maglim -15.  # threshold magnitude for blinding instruments
        trailf = 5./(fovl*3600.) #" wide trail
        
    elif telinslabel == "Cat703":
        telescope   = '0.7m Catalina'
        lat         = 30.41 #latitude
        instrument  = "CCD"
        expt        = 30. #s
        fovl        = 4.4 #deg
        resol       = 3./3600. # deg
        maglim      = maglimoff + 19.8   # detection limit for expt
        magbloom    = magbloomoff + maglim -10.  # threshold magnitude for blinding instruments
        trailf = 5./(fovl*3600.) #" wide trail

    elif telinslabel == "CatG96":
        telescope   = '1.52m Catalina'
        lat         = 32.43 #latitude
        instrument  = "CCD"
        expt        = 30. #s
        fovl        = 2.2 #deg
        resol       = 2./3600. # deg
        maglim      = maglimoff + 21.4   # detection limit for expt
        magbloom    = magbloomoff + maglim -10.  # threshold magnitude for blinding instruments
        trailf = 5./(fovl*3600.) #" wide trail
        
    elif telinslabel == "WFI":
        telescope   = 'MPE/ESO 2.2m'
        lat         = -29.25 #latitude
        instrument  = "WFI"
        expt        = 300. #s
        fovl        = 0.5 #deg
        resol       = 1./3600. # deg
        maglim      = maglimoff + 23.3   # detection limit for expt
        magbloom    = magbloomoff + 14.3  # threshold magnitude for blinding instruments
        trailf = 5./(fovl*3600.) #" wide trail
        
    elif telinslabel == "VST":
        telescope   = 'VST'
        lat         = -24.62 #latitude
        instrument  = "OmegaCam"
        expt        = 300. #s
        fovl        = 1. #deg
        resol       = 0.8/3600. # deg
        maglim      = maglimoff + 23.9   # detection limit for expt
        magbloom    = magbloomoff + 14.8  # threshold magnitude for blinding instruments
        trailf = 5./(fovl*3600.) #" wide trail
    
    elif telinslabel == "EFOSC":
        telescope   = 'NTT'
        lat         = -29.25 #latitude
        instrument  = "EFOSC2"
        expt        = 300. #s
        fovl        = 4./60. #deg
        resol       = 1./3600. # deg
        maglim      = maglimoff + 24.2   # detection limit for expt
        magbloom    = magbloomoff + 15.2  # threshold magnitude for blinding instruments
        trailf = 5./(fovl*3600.) #" wide trail

    elif telinslabel == "FORSimg":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "FORS2 Img"
        expt        = 300. #s
        fovl        = 6./60. #deg
        resol       = 1/3600. # deg
        maglim      = maglimoff + 25.2 # detection limit for expt
        magbloom    = 18.5 # threshold magnitude for blinding instruments
        trailf      = 10./(fovl*3600.) #" wide trail


    elif telinslabel == "Binoc":
        telescope   = 'Brussels'
        lat         =  50 #latitude
        instrument  = 'Binocular 10x70'
        expt        = 0.  #s
        fovl        = 5. #deg
        resol       = 10./60. # 10arcmin
        maglim      = 6. # detection limit for expt
        magbloom    = 99 # threshold magnitude for blinding instruments
        trailf      = 1. #" wide trail

    
    elif telinslabel == "LSST":
        telescope   = 'VRO'
        lat         = -30.24 #latitude
        instrument  = "LSST"
        expt        = 15 #s
        fovl        = 3. #deg
        resol       = 0.8/3600. 
        maglim      = maglimoff + 24.6
        trailf      = 10./(fovl*3600.) #10" wide trail
        magbloom    = 18.3 # = 7 -2.5log10( 0.8"/ (30*60 "/s) / 15s)
        trailf      = 10./(fovl*3600.) #10" wide trail

    elif telinslabel == "HAWKI":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "HAWKI"
        expt        = 60. #s
        fovl        = 7.5/60. #deg
        resol       = 0.6/3600. # deg
        maglim      = maglimoff + 21.4 # detection limit for expt
        magbloom    = magbloomoff + 13. # threshold magnitude for blinding instruments
        trailf      = 5./(fovl*3600.) #" wide trail


    elif telinslabel == "MICADO":
        telescope   = 'ELT'
        lat         = -24.6 #latitude
        instrument  = "MICADO"
        expt        = 60. #s
        fovl        = 50./3600. #deg
        resol       = 0.015/3600. # deg
        maglim      = maglimoff + 24.9 # detection limit for expt
        magbloom    = 24.9 # threshold magnitude for blinding instruments
        trailf      = 5./(fovl*3600.) #" wide trail

        #--- SPECTRO

    elif telinslabel == "FORSspec":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "FORS2 Spec"
        expt        = 1200. #s
        fovl        = 6./60. #deg
        fovw        = 1./3600. #deg
        resol       = 1./3600. # deg
        maglim      = maglimoff + 21.95 # detection limit for expt
        magbloom    = magbloomoff + maglim - 8.5 # threshold magnitude for blinding instruments
        trailf      = 5./(fovl*3600.) #" wide trail

    elif telinslabel == "UVES":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "UVES"
        expt        = 1200. #s
        fovl        = 10./3600. #deg
        fovw        = 1./3600. #deg
        resol       = 1./3600. # deg
        maglim      = maglimoff + 17. # detection limit for expt
        magbloom    = magbloomoff + 8.5 # threshold magnitude for blinding instruments
        trailf      = 5./(fovl*3600.) #" wide trail


    elif telinslabel == "XSHOOTER":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "X-SHOOTER"
        expt        = 1000. #s
        fovl        = 10./3600. #deg
        fovw        = 1./3600. #deg
        resol       = 1./3600. # deg
        maglim      = maglimoff + 20. # detection limit for expt
        magbloom    = magbloomoff + 8.5 # threshold magnitude for blinding instruments
        trailf      = 5./(fovl*3600.) #" wide trail

    elif telinslabel == "UVES1h":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "UVES"
        expt        = 3600. #s
        fovl        = 10./3600. #deg
        fovw        = 1./3600. #deg
        resol       = 1./3600. # deg
        maglim      = maglimoff + 18.2 # detection limit for expt
        magbloom    = magbloomoff + 10. # threshold magnitude for blinding instruments
        trailf      = 5./(fovl*3600.) #" wide trail
        
    elif telinslabel == "4MOST": 
        telescope   = 'VISTA'
        lat         = -24.62 #latitude
        instrument  = "4MOST"
        expt        = 1200. #s
        fovl        = 2.3   #deg FoV of the focal plane
        resol       = 3/3600. # deg Fibre = 3"
        maglim      = maglimoff + 19.75 # 18.75#=from Genv simul # 20.51 # detection limit for expt for LRS
        maglim      = maglimoff +  20.51 # detection limit for expt for LRS
        magbloom    = magbloomoff + 11. # threshold magnitude for blinding instruments
        trailf      = 1.3/2436 # 1.3 fibre hit per trail

    elif telinslabel == "4MOSTH": 
        telescope   = 'VISTA'
        lat         = -24.62 #latitude
        instrument  = "4MOST"
        expt        = 1200. #s
        fovl        = 2.3   #deg FoV of the focal plane
        resol       = 3/3600. # deg Fibre = 3"
        maglim      = maglimoff + 18.6 # detection limit for expt for LRS
        magbloom    = magbloomoff + 11. # threshold magnitude for blinding instruments
        trailf      = 1.3/2436 # 1.3 fibre hit per trail

    elif telinslabel == "LAMOST": 
        telescope   = 'LAMOST'
        lat         = 40.40 #latitude
        instrument  = "MRS"
        expt        = 1200. #s
        fovl        = 5.   #deg FoV of the focal plane
        resol       = 3.3/3600. # deg Fibre = 3"
        maglim      = maglimoff + 99. # detection limit for expt for MRS
        magbloom    = magbloomoff + 11. # threshold magnitude for blinding instruments
        trailf      = 0.43 # 1.3 fibre hit per trail, /4000 for fraction

    elif telinslabel == "MOONS-LR": 
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "MOONS-LR"
        expt        = 150. #s
        fovl        = .210   #deg FoV of the focal plane
        resol       = 0.5/3600. # deg Fibre = 3"
        maglim      = maglimoff + 21.3 # detection limit for expt for MRS
        magbloom    = magbloomoff -99. # threshold magnitude for blinding instruments
        trailf      = 2.1/2000 # fibre hit per trail, /n_fibre for fraction
    elif telinslabel == "MOONS-HR": 
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "MOONS-LR"
        expt        = 600. #s
        fovl        = .210   #deg FoV of the focal plane
        resol       = 0.5/3600. # deg Fibre = 3"
        maglim      = maglimoff + 18.73 # detection limit for expt for MRS
        magbloom    = magbloomoff -99. # threshold magnitude for blinding instruments
        trailf      = 2.1/2000 # fibre hit per trail, /n_fibre for fraction

    elif telinslabel == "ESPRESSO":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "ESPRESSO"
        expt        = 1200. #s
        fovl        = 2./3600. #deg
        resol       = 2./3600. # deg
        maglim      = maglimoff + 15.8 # detection limit for expt
        magbloom    = magbloomoff + 6. # threshold magnitude for blinding instruments
        trailf      = 1. # fibre lost

        #--- Thermal IR
    elif telinslabel == "VISIR":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "VISIR"
        expt        = 10. #s
        fovl        = 60./3600. #deg
        resol       = 0.2/3600. # deg
        trailf      =  5./(fovl*3600.)  #" wide trail
        maglim      = 99
        magbloom = -99. # threshold magnitude for blinding instruments
        elevs       = 0.

    elif telinslabel == "VISIR1":
        telescope   = 'VLT'
        lat         = -24.62 #latitude
        instrument  = "VISIR"
        expt        = 0.02 #s
        fovl        = 60./3600. #deg
        resol       = 0.2/3600. # deg
        trailf      =  5./(fovl*3600.)  #" wide trail
        maglim      = 99
        magbloom = -99. # threshold magnitude for blinding instruments

        elevs       = 0.


    elif telinslabel == "ELTHrm":
        telescope = 'ELT'
        instrument = "HARMONY"
        expt = 600. #s
        fovl = 8./3600. #deg
        lat = -24.6 #latitude
        trailf = 1
        magbloom = -999.

    elif telinslabel == "ELTMetImg":
        telescope = 'ELT'
        lat         = -24.62 #latitude
        instrument = "METIS Imaging"
        expt = 10. #s
        fovl = 10./3600. #deg
        resol = 0.03/3600.
        maglim = 99.
        trailf = 1
        magbloom = -999.


    elif telinslabel == "ELTMetImg":
        telescope = 'ELT'
        instrument = "METIS Imaging"
        expt = 60. #s
        fovl = 10./3600. #deg
        lat = -24.6 #latitude
        trailf = 1
        magbloom = -999.

    elif telinslabel == "ELTMetLng":
        telescope = 'ELT'
        instrument = "METIS Long Slit"
        expt = 300. #s
        fovl = 10./3600. #deg
        lat = -24.6 #latitude
        trailf = 1
        magbloom = -999. 
    

    elif telinslabel == "ELTMicImg":
        telescope = 'ELT'
        instrument = "MICADO Imaging"
        expt = 300. #s
        fovl = 50./3600. #deg
        lat = -24.6 #latitude
        trailf = 1
        magbloom = -999.     

    elif telinslabel == "ELTMicLng":
        telescope = 'ELT'
        instrument = "MICADO Long Slit "
        expt = 600. #s
        fovl = 3./3600. #deg
        lat = -24.6 #latitude
        trailf = 1
        magbloom = -999.     

    
    elif telinslabel == "TrailDens":
        telescope = ' '
        lat = -24.6 #latitude
        instrument = "Density"
        expt = 1. #s
        fovl = 1. #deg
        resol       = 1./3600. # deg
        trailf = 1. # trailf factor for faint satellites
        maglim      = 99 # 5sig detection limit for expt
        magbloom = -999.# not relevant: dest=1

    elif telinslabel == "SatDens":
        telescope = ' '
        lat = -24.6 #latitude
        instrument = "Density"
        telescope   = ' '
        expt        = 0. #s
        fovl        = 1. #deg
        resol       = 1./3600. # deg
        maglim      = 99 # 5sig detection limit for expt
        magbloom    = -99. # threshold magnitude for blinding instruments
        trailf = 1.
        
    
    elif telinslabel == "skyMag":
        telescope   = ' '
        lat = -24.6 #latitude
        instrument  = "Sky brightness"
        expt        = 1. #s
        fovl        = 1. #deg
        resol       = 0.8/3600. # deg
        maglim      = 99 # 5sig detection limit for expt
        magbloom    = -99. # threshold magnitude for blinding instruments
        trailf = 1.
        
    elif telinslabel == "skyFrac":
        telescope   = ' '
        lat = -24.6 #latitude
        instrument  = "Sky brightness"
        expt        = 1. #s
        fovl        = 1. #deg
        resol       = 0.8/3600. # deg
        maglim      = 99 # 5sig detection limit for expt
        magbloom    = -99. # threshold magnitude for blinding instruments
        trailf = 1.
        

    else:
        print("Telescope "+telinslabel + " unknown")
        exit()


    if fovw <0:
        fovw = fovl*1.
    return telescope, lat, instrument, expt, fovl, fovw, resol, maglim, magbloom, trailf



#------------------------------------------------------------------------------

def findConstellations(constellationsll):
    # Assemble a set of constellations for preset constellationll

    from ConAn.constellations import SL1, SL2, hSL1, hSL2, OW2, OW2r, GW, AK, yesturday, constoday, SLtoday, csl1a1, ESP

    constDir = {
        "YESTURDAY": {  "constellations" : yesturday,
                        "constellationsl" : "Pre-const LEOs" },
        "TODAY": {      "constellations" : np.concatenate((yesturday,constoday)),
                        "constellationsl" : "LEOs today" },
        "SLtoday": {    "constellations" : SLtoday,
                        "constellationsl" : "SL today"},
        "SL2OW2": {     "constellations" : np.concatenate((SL2,OW2)),
                        "constellationsl" : "SL2 & OW2"},
        "csl1a1": {     "constellations" : [csl1a1],
                        "constellationsl" : "SL1 A1"},
        "SLOW2": {      "constellations" : np.concatenate((SL1, SL2, OW2)),
                        "constellationsl" : "SL1+2 & OW2"},
        "OW2r": {       "constellations" : OW2r,
                        "constellationsl" : "OW2 reduced"},
        "OW2": {        "constellations" : OW2,
                        "constellationsl" : "OW2 original"},
        "SL1": {        "constellations" : (SL1),
                        "constellationsl" : "SL1"},
        "SL": {         "constellations" : np.concatenate((SL1, SL2)),
                        "constellationsl" : "SL1+2"},
        "SLhigh": {     "constellations" : np.concatenate((hSL1, hSL2)),
                        "constellationsl" : "SL1+2 HIGH"},
        "SLOWGW": {     "constellations" : np.concatenate((SL1, SL2, OW2r, GW)),
                        "constellationsl" : "SL1+2, OW2r, GW1+2"},
        "SLOWGWAK": {   "constellations" : np.concatenate((SL1, SL2, OW2r, GW, AK)),
                        "constellationsl" : "SL1+2, OW2r, GW1+2, AK"},
        "SLOWr": {      "constellations" : np.concatenate((SL1, SL2, OW2r)),
                        "constellationsl" :  "SL1+2 & OW2 red."},
        "ALL": {        "constellations" : np.concatenate((SL1, SL2, OW2r, GW, AK, ESP)),
                        "constellationsl" : "SL1+2, OW2r, GW1+2, AK, E-Sp"},
        "ESPACE": {     "constellations" : (ESP),
                        "constellationsl" :  "E-Space"}
        }

    if constellationsll in constDir:
        return constDir[constellationsll]['constellationsl'], constDir[constellationsll]['constellations']
    else:
        print('List of available constellations:')
        for c in constDir:
            print(c ,":  ", constDir[c]['constellationsl'])
        exit(1)
    

#------------------------------------------------------------------------------

def initPolPlot(ax):
    # initialize a polar plot
    
    plt.rcParams.update({'font.size': 15})
    _ = ax.set_theta_zero_location('N')  
    _ = ax.set_rticks([30,60,70,80,90])
    _ = ax.set_yticklabels([])
    _ = ax.set_xticks(np.radians(range(0,360,90)))
    _ = ax.set_xticklabels(['N','E','S','W'])
    _ = ax.set_rlim(0.,90)
    return

#------------------------------------------------------------------------------

def plotOneDens(ax, AzEl, dens, label):
    # plot a density function
    # ax: axis object where to plot
    # AzEl: array with the azimut and elevations
    # dens: the density array to be plotted
    # label: label...
    
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
    # plot a series of density distributions
    
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
    # write labels at (x,y) on a polar plot in ax

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
    # draw the HA and Dec lines in an altaz plot

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



