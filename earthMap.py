#!/usr/bin/env python3
''' SatConAnalytic 

Plot a map of the Earth, with the number of observable satellites above  a given elevation.

Overplot a possible realization of the satellites position.

'''

##import matplotlib
##matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from astropy.table import  vstack


import argparse

# import ConAn routines
import conan as caLib
import conanplot as cpLib
import constants as caCst
import satDots


#===== CONFIG ===================================================================
# elevation limits [deg], must be decreasing
elLim = np.array([60.,30.,20., 10.,0.])

# resolution for the scanning of the Earth [deg]
#  1: ideal;   10: fast; 30: debug
groundStep = 3. # [deg]

# resolution for the scanning of the sky [deg]
# 1: production, slow;  3 ok-ish debug; 5-30: quick-and-dirty debug
skyStep = 2.

#===== PARAMETERS ===================================================================
parser = argparse.ArgumentParser(description=
'''Satellite constellations: Map of the Earth showing the number of satellites visible from each point''')


parser.add_argument('-d','--deltaSun', default=-20.,
                    help="Sun: Declination of the Sun [deg]")

parser.add_argument('-C','--constellations', 
                    default='SL1',
#                    default='SLOWGWAK',
                    help="ID of the constellation group; list for a list")


choices=[f'{e:.0f}' for e in elLim]


# observatory
parser.add_argument('-e','--elevCut', 
                    default="30",
                    choices=choices,
                    help="Observatory: Elevation cut-off")

#parser.add_argument('-r','--resol', 
#                    help="Observatory: Resolution of the instrument [deg]")
#parser.add_argument('-t','--expt', 
#                    help="Observatory: Exposure time [sec]")
#parser.add_argument('-f','--fovl',
#                    help="Observatory: Field of view of the instrument. Length or diametre [deg]")
#parser.add_argument('-w','--fovw', 
#                    help="Observatory: Field of view of the instrument. Width. Equal to Length if omitted [deg]")
#parser.add_argument('-m','--maglim', 
#                    help="Observatory: Detection limit magnitude of the instrument [5sigma Mag during expTime]")
#parser.add_argument('--magbloom',
#                    help="Observatory: Magnitude over which the instrument saturates [Mag for expTime]. Default: -99 (no blooming)")
#parser.add_argument('-k','--trailf', 
#                    help="Observatory: Trail filling fraction (width of the trail as fraction of FoV)")
#parser.add_argument('-s','--telescope', 
#                    help="Observatory: Name of the telescope")
#parser.add_argument('-i','--instrument', 
#                    help="Observatory: Name of the instrument")
#parser.add_argument('-l','--lat', 
#                    help="Observatory: latitude")
#parser.add_argument('-T','--code', 
#                    help='''Observatory: Predefined telescope/instrument with
#                         extptime, FoVl, FoVw, maglim, magbloom, trailf,  
#                         telescope, instrument, resolution, latitude. 
#                         Use individual options to overwrite presets.
#                         SPECIAL CODES: SatDens for satellite density; 
#                         TrailogDensity for trail density;
#                         skyMag for sky surface brightness [mag];
#                         skyFrac for sky surface brightness as a fraction.
#''')
#parser.add_argument('-M','--magSelect', default="OBS", 
#                    choices=['ALL', 'OBS', 'BRIGHT', 'FAINT', 'EFFECT'],
#                    help="selection on magnitude; default=OBServable")


myargs = parser.parse_args()





# #
# # OBSERVATORY TELESCOPE INSTRUMENT
# #
# print('=====TELESCOPE/INSTRUMENT SETUP=======================================')
# myTel = cpLib.getTelescope(myargs)
# print(myTel)



# Sun
sunAlpha = 0. # for computation, the sun is always on RA=0.
sunDelta = float(myargs.deltaSun)


# create a grid on the Earth
Lat        = np.arange( (-90. +groundStep), 90., groundStep)
Long       = np.arange( 0., 180.5, groundStep)
           # ONLY half the longitudes - this is symmetric

fillLat, fillLong = np.meshgrid(Lat, Long)
fillSatCount      = np.zeros(  (len(Long),len(Lat),len(elLim) ) )
fillSatObsCount   = np.zeros(  (len(Long),len(Lat),len(elLim) ) )
fillSunElev       = np.zeros(  (len(Long),len(Lat) ) )

print(f'Len Lat: {len(Lat)}, Long {len(Long)}')
print(f'Shape: {fillLat.shape}, {fillLong.shape}')
print(f'storage: {fillSatCount.shape}')
print(f'elev:     {fillSunElev.shape}')


# sky at the considered location:
AzEl = caLib.fill_AzEl(skyStep)        

# expand the constellation Id into a list of real constellations
CONSTELLATIONS = caLib.findConstellations(myargs.constellations)
print('=====CONSTELLATIONS===================================================')
print( CONSTELLATIONS.ToC )
print()

# assemble folder name and make sure it exists 
dirName = f'MAP_{CONSTELLATIONS.name}_{int(sunDelta)}'
Path(dirName).mkdir(parents=True, exist_ok=True)

# Scan the Earth
try: # do we already have the files?
    fillLong        = np.reshape(np.fromfile(dirName+'/fillLong.dat'),     (len(Long), len(Lat)) )
    fillLat         = np.reshape(np.fromfile(dirName+'/fillLat.dat'),      (len(Long), len(Lat)) )
    fillSunElev     = np.reshape(np.fromfile(dirName+'/fillSunElev.dat'),  (len(Long), len(Lat)) )
    fillSatCount    = np.reshape(np.fromfile(dirName+'/fillSatCount.dat'), (len(Long), len(Lat), len(elLim)) )
    #TBF# fillSatObsCount = np.reshape(np.fromfile(dirName+'/fillSatObs.dat'),   (len(Long), len(Lat), len(elLim)) )



except: # no... recompute everything

    sunElMax = -999. # init the maximum sun elevation 
                     # for which no satellites are visible.
                     # We can skip calculation below that limit

    idxLg = np.indices( fillLong.shape )[0].flatten()
    idxLt = np.indices( fillLong.shape )[1].flatten()

    long0 = 999.

    for iii in range(len(idxLg)):
        myLong = fillLong[idxLg[iii],idxLt[iii]]
        myLat  = fillLat[idxLg[iii],idxLt[iii]]

        # console output
        if myLong != long0: # new line
            long0 = myLong
            print(f'\n{myLong:03.0f}:', end="")


        # sun elevation at this observatory
        sunHA = (myLong - sunAlpha)%360.
        sunAz,sunEl = caLib.radec2azel(sunHA, sunDelta, myLat)
        fillSunElev[idxLg[iii],idxLt[iii]]  = sunEl
    

        # init arrays with the various results; same array geometry as AzEl 
        densSatAll = np.zeros_like(AzEl[0]) # density of satellites     (all sat)
        #TBD# densSatObs = np.zeros_like(AzEl[0]) # density of satellites, only observable
        #TBD# densVelObs = np.zeros_like(AzEl[0]) # velocity density of satellites, only observable
        elCount    = np.zeros_like(elLim )


        if sunEl < 0. and sunEl >= sunElMax: 
            # it's worth computing the full satellite distribution

            # Scan the constellation shells
            for myShell in CONSTELLATIONS.shells:                 
                # model the shell
                densSi, veli, magi =  myShell.modelOneShell(AzEl,myLat, sunHA,sunDelta )
            
                # all sat:
                densSatAll += densSi

                #TBD# TBD TBD: account for magnitude, trailing, observability
                #TBD#
                #TBD# # effective magnitude and extremes
                #TBD# # does the satellite trail more than 1 resolution element during expT:
                #TBD# trailing =  veli*myTel.expt >= myTel.resol  # bolean for non-zero trailing
                #TBD# mageffi = magi*1.  # init effective mag for non-trailing
                #TBD# mageffi[trailing] = magi[trailing]\
                #TBD#     - 2.5*np.log10(myTel.resol/ (veli[trailing]*myTel.expt)  )
                #TBD#     # correct for trailing sat
                #TBD#
                #TBD# # only observable ones
                #TBD# densSatObsi = np.copy(densSi)
                #TBD# densSatObsi[ mageffi > myTel.maglim] = 0.
                #TBD# densSatObs += densSatObsi
                #TBD# densVelObs += densSatObsi * veli
                #TBD#
                #TBD# here add the 
                #TBD# - mag selection
                #TBD# - effect on data

            # integrate density -> satellite counts
            elCount = caLib.integrateSat(elLim,AzEl,densSatAll)
            fillSatCount[idxLg[iii],idxLt[iii]] = elCount

            #TBD# TBD TBD
            #TBD# # integrate density -> observable satellite counts
            #TBD# elCount = caLib.integrateSat(elLim,AzEl,densSatObs)
            #TBD# fillSatObsCount[idxLg[iii],idxLt[iii]] = elCount
            #TBD# # observable trail density
            #TBD# densTrailObs =  ( densSatObs * myTel.fovl*myTel.fovw + 
            #TBD#                   densVelObs * myTel.fovl * myTel.expt )

            maxElCount = max(elCount)
            if maxElCount == 0.:  # not a single satellite
                sunElMax = sunEl    # this is the new limit for sunEl

                print('.', end="")
            elif maxElCount <= 10.:
                print('+', end="")
            elif maxElCount <= 100.:
                print('o', end="")
            elif maxElCount <= 1000.:
                print('O', end="")
            else:
                print('X', end="")

            #<- finished scan of the Earth
        else:
            print('-', end="")

    print()
    # save the results to file.
    fillLong.tofile(       dirName+'/fillLong.dat')
    fillLat.tofile(        dirName+'/fillLat.dat')
    fillSatCount.tofile(   dirName+'/fillSatCount.dat')
    #TBD# fillSatObsCount.tofile(dirName+'/fillSatObs.dat')
    fillSunElev.tofile(    dirName+'/fillSunElev.dat')
    print('files are saved')
    #<- finished full re-calculation




## PLOT ##

# reshape output for map:
# duplicate each point of the half-map to get full map coverage
plotLat    = Lat
plotLong   = np.append(- Long[::-1], Long )
fillPlotLat, fillPlotLong = np.meshgrid(plotLat, plotLong)
fillPlotSunElev           = np.concat( [fillSunElev, fillSunElev[::-1,:]] )
fillPlotSatCount          = np.concat( [fillSatCount,  fillSatCount[::-1,:,:] ])
#TBD# fillPlotSatObsCount       = np.concat( [fillSatObsCount,  fillSatObsCount[::-1,:,:]] )


# which elevation cutoff to plot
elevCut = list(elLim).index(  float(myargs.elevCut))

# select the observable satellite count above elevation cutoff
#TBD# plotIt = fillPlotSatObsCount[:,:,elevCut]
plotIt = fillPlotSatCount[:,:,elevCut]

# Flag points in day-time 
plotIt[ fillPlotSunElev > 0.] = -1.

# PLOT
m = Basemap(projection='cyl',
            llcrnrlat=-90,urcrnrlat=90,
            llcrnrlon=-180,urcrnrlon=180,
            resolution='c')



# convert long,lat to x,y; not critical for "cyl"
xpt,ypt = m(fillPlotLong,fillPlotLat) 

# Plot satellite counts
if 1:

    mymax = np.amax(plotIt)
    levels = np.linspace(-1.,mymax,100)
    cmap = 'inferno'
    cs = plt.contourf(xpt,ypt,plotIt, 
                      levels=levels,
                      cmap=cmap,
                      extend='neither')


    # beautify color bar
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    clevels = np.arange(0., 
                        10.**(int(np.log10( mymax )))*9.,
                        10.**(int(np.log10( mymax )-1.)))
    clevels = clevels[ clevels < mymax ]
    while len(clevels) > 20:
        clevels = clevels[::2]

    clabel = [ f'{l:.0f}' for l in clevels]
    cbar.set_ticks(clevels)
    cbar.set_ticklabels(clabel)
    cbar.set_label(f'Number of Sat above {elLim[elevCut]:.0f}'+u'$^o$ elevation')


# Sun
if 1:
    for myElevMin, myAlpha  in zip([-18., -12., -6., 0.], [.3,.3,.5,1.]):
        m.contourf( xpt, ypt, fillPlotSunElev, 
                     levels=[myElevMin,90.],
                     colors='royalblue',
                     alpha=myAlpha)
    csun = m.contour(xpt, ypt,fillPlotSunElev, 
                 levels=[-18.,-12.,-6.,0],
                 linewidths=[1.,1.,1.,2.], 
                 linestyles='solid', 
                 colors='k')
    lsun = plt.clabel(csun, fmt='%.0f$^o$')



# Discrete satellites
if 1:
    SAT =  vstack(
        [ myShell.orbitalElementsTable() for myShell in CONSTELLATIONS.shells]
    )
    SAT["anomr"] = satDots.propagateTimeAnom(SAT["anom0"], SAT["omega"], 1.234) # 1.234: randomish time
    SAT["noder"] = satDots.propagateTimeNode(SAT["node0"], 1.234)
    SAT["latr"],SAT["longr"] = satDots.elementsToLongLat(SAT["noder"], SAT["inc"], SAT["anomr"])
    SAT["xg"],SAT["yg"],SAT["zg"]   = satDots.LongLatToGeoXYZ(SAT["alt"], SAT['latr'], SAT["longr"])
    SAT["bIlluminated"] = satDots.illuminatedSat(SAT["xg"],SAT["yg"],SAT["zg"], sunAlpha, sunDelta)

    SAT["longr"] = SAT["longr"] -np.pi
    SAT["xPt"], SAT["yPt"]  = m(np.degrees(SAT["longr"]), np.degrees(SAT["latr"] ))
    m.scatter(SAT["xPt"], SAT["yPt"], c=SAT["bIlluminated"] , 
              cmap='YlOrRd_r',marker='.', s=1, alpha=0.25)


m.drawcoastlines(linewidth=0.5, color='g')# Sun and twilights
m.drawparallels( np.arange(-80,81,20), color='g', linewidth=0.5)

title = (f'Constellation: {CONSTELLATIONS.name}, {CONSTELLATIONS.totSat} sat. '
         +u' - Sun: $\delta=$'+f'{sunDelta:.0f}'+u'$^o$')
plt.title(title)
plt.savefig( dirName + f'/satObs_{myargs.elevCut}.pdf')
plt.savefig( 'w.png')
print(f'Figure in {dirName}/satObs_{myargs.elevCut}.pdf')
#plt.show()