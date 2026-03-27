#!/usr/bin/env python3
''' SatConAnalytic 

Plot a map of the Earth, with the number of observable satellites above  a given elevation.

Overplot a possible realization of the satellites position.

'''

from html import parser
import logging
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt



from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from astropy.table import  vstack

import SatConAnalytic.conan as caLib
import SatConAnalytic.conanplot as cpLib
import SatConAnalytic.telescopes as telescopeLib
import SatConAnalytic.utils as utLib
import SatConAnalytic.satDots as satDots
import SatConAnalytic.constellations as constLib

log = logging.getLogger('satcon')
log.setLevel(logging.DEBUG)  
for logger_name in ("matplotlib", "matplotlib.pyplot", "matplotlib.font_manager", "matplotlib.backends"):
    logging.getLogger(logger_name).setLevel(logging.WARNING)


#---------------------------------------------------------------------------
def create_argument_parser():
    """Create and return the argument parser"""
    #--- command line arguments
    parser = argparse.ArgumentParser(description='''Satellite constellations: Map of the Earth showing the number of satellites visible from each point''')
    
    elLim = np.array([60.,30.,20., 10.,0.])
    choices=[f'{e:.0f}' for e in elLim]

    parser.add_argument('-a','--sunAlpha', default=0.,
                    help="Sun: longitude of the Sun [deg]")


    parser.add_argument('-d','--deltaSun', default=-20.,
                    help="Sun: Declination of the Sun [deg]")

    parser.add_argument('-C','--constellations', 
                    default='SLOWGWAK',
                        help="Constellation code (or meta-code); list for a list.")
    parser.add_argument('--constFile', default="constellations.json",
                        help="Constellation: Which constellation file to use (default: constellations.json)")
    
    parser.add_argument('--log', action='store_true',
                        help="log scale")

    parser.add_argument('--dots', action='store_true',
                        help="Plot: Plot the satellite dots")


    parser.add_argument('-e','--elevCut', 
                    default="30",
                    choices=choices,
                    help="Observatory: Elevation cut-off")
    parser.add_argument('--pdf', action='store_true',
                        help="OUTPUT: output file in pdf (default is png)")



    return parser

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
def main(args=None):
    """Main function to run the Earth map plotting"""

    #===== CONFIG ===========================================================
    # elevation limits [deg], must be decreasing
    elLim = np.array([60.,30.,20., 10.,0.])

    # resolution for the scanning of the Earth [deg]
    #  1: ideal;   10: fast; 30: debug
    groundStep = 2. # [deg]

    # resolution for the scanning of the sky [deg]
    # 1: production, slow;  3 ok-ish debug; 5-30: quick-and-dirty debug
    skyStep = 2.

    #----- Parse arguments
    parser = create_argument_parser()
    if args is None:
        myargs = parser.parse_args()
    else:
        myargs = parser.parse_args(args)

    if myargs.pdf :
        myargs.outputformat = ".pdf"
    else:
        myargs.outputformat = ".png"
        


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

    log.debug(f'Len Lat: {len(Lat)}, Long {len(Long)}')
    log.debug(f'Shape: {fillLat.shape}, {fillLong.shape}')
    log.debug(f'storage: {fillSatCount.shape}')
    log.debug(f'elev:     {fillSunElev.shape}')


    # sky at the considered location:
    AzEl = caLib.fill_AzEl(skyStep)        


    #===========================================================================
    # CONSTELLATIONS
    #===========================================================================
    log.info('=====CONSTELLATIONS==================================')
    CONSTELLATIONS = constLib.findConstellations(myargs.constellations, constFile=myargs.constFile)  # constellationS object

    log.debug(f'Found constellations:\n{CONSTELLATIONS.ToC}')
    log.info(f'Constellations: {CONSTELLATIONS.totSat} satellites in {CONSTELLATIONS.totShells} shells.')

    # assemble folder name and make sure it exists 
    dirName = f'MAP_{myargs.constellations}_{int(sunDelta)}'
    Path(dirName).mkdir(parents=True, exist_ok=True)

    # Scan the Earth
    try: # do we already have the files?
        fillLong        = np.reshape(np.fromfile(dirName+'/fillLong.dat'),     (len(Long), len(Lat)) )
        fillLat         = np.reshape(np.fromfile(dirName+'/fillLat.dat'),      (len(Long), len(Lat)) )
        fillSunElev     = np.reshape(np.fromfile(dirName+'/fillSunElev.dat'),  (len(Long), len(Lat)) )
        fillSatCount    = np.reshape(np.fromfile(dirName+'/fillSatCount.dat'), (len(Long), len(Lat), len(elLim)) )
        #TBF# fillSatObsCount = np.reshape(np.fromfile(dirName+'/fillSatObs.dat'),   (len(Long), len(Lat), len(elLim)) )

        log.info(f'Found cached results in {dirName}/fill*.dat; reusing them to plot the map. Remove the directory {dirName} and its content to reset the cache and re-run the calculation.')   


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
            elCount    = np.zeros_like(elLim )


            if sunEl < 0. and sunEl >= sunElMax: 
                # it's worth computing the full satellite distribution

                # Scan the constellation shells
                for myShell in CONSTELLATIONS.shells:                 
                    # model the shell
                    densSi, _, _ =  myShell.modelOneShell(AzEl,myLat, sunHA,sunDelta )
                
                    # all sat:
                    densSatAll += densSi

                # integrate density -> satellite counts
                elCount = caLib.integrateSat(elLim,AzEl,densSatAll)
                fillSatCount[idxLg[iii],idxLt[iii]] = elCount

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
        fillSunElev.tofile(    dirName+'/fillSunElev.dat')


        log.info(f'Finished scan of the Earth. Results saved in {dirName}/fill*.dat')
        log.info(f'If you want to re-run the calculation, remove the directory {dirName} and its content to reset the cache.')

        #<- finished full re-calculation




    log.info('====CONVERSIONS================================')

    # reshape output for map:
    # duplicate each point of the half-map to get full map coverage
    plotLat    = Lat
    plotLong   = np.append(- Long[::-1], Long )
    fillPlotLat, fillPlotLong = np.meshgrid(plotLat, plotLong)
    fillPlotSunElev           = np.concat( [fillSunElev, fillSunElev[::-1,:]] )
    fillPlotSatCount          = np.concat( [fillSatCount,  fillSatCount[::-1,:,:] ])


    # adjust the position of the sun according to the input sunAlpha
    sunAlphaNew = (sunAlpha + float(myargs.sunAlpha))%360. # [deg], solar time UT
    sunAlphaStep = int(sunAlphaNew / groundStep + 0.5)
    print(f'sunAlpha: {sunAlpha:.1f} -> {sunAlphaNew:.1f} (step {sunAlphaStep})')   

    # shift the maps by the corresponding number of steps
    fillPlotSunElev = np.roll(fillPlotSunElev, sunAlphaStep, axis=0)
    fillPlotSatCount = np.roll(fillPlotSatCount, sunAlphaStep, axis=0)


    # which elevation cutoff to plot
    elevCut = list(elLim).index(  float(myargs.elevCut))

    # select the observable satellite count above elevation cutoff
    plotIt = fillPlotSatCount[:,:,elevCut]


    # log scale
    if myargs.log:
        plotIt = np.log10( plotIt + 0.1) # +1 to avoid log(0)

    # Flag points in day-time 
    plotIt[ fillPlotSunElev > 0.] = -1.

    # PLOT
    log.info('====PLOT================================')

    m = Basemap(projection='cyl',
                llcrnrlat=-90,urcrnrlat=90,
                llcrnrlon=-180,urcrnrlon=180,
                resolution='c')

    # convert long,lat to x,y; not critical for "cyl"
    xpt,ypt = m(fillPlotLong,fillPlotLat) 

    # Plot satellite counts
    if 1:
        log.info('Plotting satellite contours...')
        mymax = np.amax(plotIt)
        levels = np.linspace(-1.,mymax,100)
        cmap = 'inferno'
        cs = plt.contourf(xpt,ypt,plotIt, 
                        levels=levels,
                        cmap=cmap,
                        extend='neither')


        # beautify color bar
        cbar = m.colorbar(cs,location='bottom',pad="5%")
        if myargs.log:
            clevels = np.arange(0., 
                            int( mymax + 1. ),
                            int(mymax -1.)/5)
            clevels = clevels[ clevels < mymax ]
            clabel = [ f'{10.**l:.0f}' for l in clevels]

        else:
            clevels = np.arange(0., 
                            10.**(int(np.log10( mymax )))*9.,
                            10.**(int(np.log10( mymax )-1.)))
            clevels = clevels[ clevels < mymax ]
            clabel = [ f'{l:.0f}' for l in clevels]

        while len(clevels) > 20:
            clevels = clevels[::2]
            clabel = clabel[::2]

        cbar.set_ticks(clevels)
        cbar.set_ticklabels(clabel)
        cbar.set_label(f'Number of Sat above {elLim[elevCut]:.0f}'+u'$^o$ elevation')


    # Sun
    if 1:
        for myElevMin, myAlpha  in zip([-18., -12., -6., 0.], [.2,.2,.2,.2]):
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
        plt.scatter(180., sunDelta, marker='*', s=200, color='yellow', edgecolor='k', zorder=10)
        plt.scatter(-180., sunDelta, marker='*', s=200, color='yellow', edgecolor='k', zorder=10, label='Sun')


    # Discrete satellites
    if myargs.dots:
        log.info('Plotting satellite dots...')
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

    title = (f'Constellation: {"".join(CONSTELLATIONS.name.split("\n"))} '+
             f'({CONSTELLATIONS.totSat:.0f} sat.)'+
             u' - Sun: $\delta=$'+f'{sunDelta:.0f}'+u'$^o$')
    plt.title(title)



    plt.savefig( dirName + f'/satObs_{myargs.elevCut}{myargs.outputformat}')
    plt.savefig( 'w.png')
    print(f'Figure in {dirName}/satObs_{myargs.elevCut}{myargs.outputformat}')



if __name__ == "__main__":

    utLib.init_logger(log)
    log.info('===EarthMap===')


    main()
