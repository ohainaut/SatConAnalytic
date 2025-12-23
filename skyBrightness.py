#! /usr/bin/env python3
# skyBrightness.py

# Script to compute and plot the sky brightness due to satellite constellations
# Validation: 
# -C FullMoon places the full moon at the anti-sun position
# -C Stars loads a set of bright stars to simulate their contribution to sky brightness

#
## ToBeDone:
# - parameter to adjust the longitude of ascending node
# - set the MagOffset in constellation.json


import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import argparse

import SatConAnalytic.conanplot as cpLib
import SatConAnalytic.skyBrightnessLib as skyLib
import SatConAnalytic.conan as caLib
import SatConAnalytic.satDots as satDots

def create_argument_parser():
    """Create and return the argument parser"""
    parser = argparse.ArgumentParser(description=
    '''Satellite constellations: sky map of satellite trail density or related information.
    Define the position of the observatory, the position of the Sun, the constellation(s),
    the instrument and its characteristics, and the output required.''')
    parser.add_argument('-d','--deltaSun', default=0.,
                        help="Sun: Declination of the Sun [deg]")
    parser.add_argument('-a','--alphaSun',
                        help="Sun: Hour Angle of the Sun [deg]. If present, overwrites elevSun")
    parser.add_argument('-e','--elevSun', default=24.,
                        help="Sun: Elevation of the Sun BELOW the horizon. Should probably be >0 in most cases [deg]")
    parser.add_argument('-C','--constellations', default='SLOWGWAK',
                        help="ID of the constellation group; list for a list")
    
    parser.add_argument('-m','--magOffset', default=0.,
                        help="magnitude offset to apply to all satellites (e.g., to simulate different satellite brightness) [mag]")


    parser.add_argument('-l','--lat', default=-24.6,
                        help="Observatory: Latitude of the observatory [deg]")
    parser.add_argument('-s','--telescope',
                        help="Observatory: Name of the telescope")
    parser.add_argument('-i','--instrument',
                        help="Observatory: Name of the instrument")
    parser.add_argument('-T','--code',
                        help='''Observatory: Predefined telescope/instrument with
                             extptime, FoVl, FoVw, maglim, magbloom, trailf,
                             telescope, instrument, resolution, latitude.
                             Use individual options to overwrite presets.
                             SPECIAL CODES: SatDens for satellite density;
                             TrailogDensity for trail density;
                             skyMag for sky surface brightness [mag];
                             skyFrac for sky surface brightness as a fraction.
    ''')

    parser.add_argument('--noplot', action='store_false',
                        help="Plot: Don't generate the plot (mostly debug)")

    parser.add_argument('--noscalebar', action='store_false',
                        help="Plot: Don't include scalebar")
    parser.add_argument('--noalmuc', action='store_false',
                        help="Plot: Don't write sat count on the almucantars")
    parser.add_argument('--nolabel', action='store_false',
                        help="Plot: Don't label the plot")
    parser.add_argument('--noDots', action='store_false',
                        help="Plot: Don't plot the satellite dots")
    parser.add_argument('--pdf', action='store_true',
                        help="Plot: output file in pdf (default is png)")
    return parser

def main(args=None):

    #----- config
    step = 1. #deg ~1. Smaller values take forever

    # Parse arguments
    parser = create_argument_parser()
    if args is None:
        myargs = parser.parse_args()
    else:
        myargs = parser.parse_args(args)
    myargs.fovw = None  # Force fovw to be computed from fovl

    #===========================================================================
    # OBSERVATORY TELESCOPE INSTRUMENT
    #===========================================================================

    print('=====TELESCOPE/INSTRUMENT SETUP=======================================')
    myTel = cpLib.getTelescope(myargs)
    print(myTel)

    #===========================================================================
    # SUN
    #===========================================================================
    sunAlpha, sunDelta, sunElev = caLib.consolidate_sun(
        myargs.alphaSun, myargs.deltaSun, myargs.elevSun, myTel.lat)

    print('SUN:')
    print(f'\tLocal time: {((180+sunAlpha)/15.)%24:.2f}h')
    print(f'\tHA = {sunAlpha:.1f}deg  = {(sunAlpha/15.)%24:.2f}h, Dec = {sunDelta:.1f}d')
    print(f'\tElevation: {sunElev:.2f}d')

    sunAz,sunEl = caLib.radec2azel(sunAlpha, sunDelta, myTel.lat)
    print(f'Validation: az= {sunAz:.2f}, el= {sunEl:.2f}d\n')



    # Azimuth-Elevation mesh:
    AzEl = caLib.fill_AzEl(step)               # mesh of Azimut-Elevation
    surfBrightness = np.zeros_like(AzEl[0]) # initialize surface brightness array


    if myargs.constellations == "FullMoon":
        # Validation: Full Moon at anti-sun position         
        print('=====FULLMOON===================================================')

        Sb = Table()
        Sb["mag"] = [-12.74]
        Sb["dot"] = [50.]
        Sb["Azr"] = [np.radians(sunAz +180.)]  # opposite to the sun
        Sb["ZD"]  = [90. + sunEl]
        Sb["bIlluminated"] = [True]
        Sb["bVis"] = [True]

        CONST_name = "FullMoon"
        CONST_totsat = 1
        CONST_illum = 1
        CONST_range = 1

        print(f'FullMoon: {Sb["mag"][0]:.2f} mag at Az: {np.degrees(Sb["Azr"][0]):.2f}, ZD: {Sb["ZD"][0]:.2f}')


    elif myargs.constellations == "Stars":
        # Validation: Bright stars
        print('=====STARS===================================================')
        Sb = skyLib.loadBrightStars(obsLat=myTel.lat, sidTime=6.)

        CONST_name = "Stars"
        CONST_totsat = len(Sb)
        CONST_illum = len(Sb)
        CONST_range = len(Sb)


    else:
        # Normal case: satellites
        print('=====SATELLITE CONSTELLATION===================================')

        CONSTELLATIONS = caLib.findConstellations(myargs.constellations)
        print( CONSTELLATIONS.ToC )
        print()
        CONST_name = CONSTELLATIONS.name
        CONST_totsat = CONSTELLATIONS.totSat




        # get the satellite positions and brightnesses
        Sv = satDots.makeConstellationStatTable(CONSTELLATIONS, 
                                        sunAlpha, sunDelta, 
                                        myTel.lat, 
                                        sunAlpha*240)#(180+sunAlpha)*360.) ## slowed 10x
        # Sv returned with 
        #  bIlluminated          bVis       Delta              
        #  Azr                ZD                mag                 dot       

        CONST_range = len(Sv)

        Sv["mag"] = Sv["mag"] - float(myargs.magOffset)  # apply mag offset
        
        Sb = Sv[  Sv["bIlluminated"] ]
        CONST_illum = len(Sb)

        print(f'Total satellites illuminated: {len(Sb)}')
        print(f'Average mag: {np.mean(Sb["mag"]):.2f}, min: {np.min(Sb["mag"]):.2f}, max: {np.max(Sb["mag"]):.2f}')
    


    #===========================================================================
    for i in range(len(Sb)):
        theta_grid = skyLib.angular_distance(AzEl[1], AzEl[0], 90.-Sb["ZD"][i], np.degrees(Sb["Azr"][i]))
        surfBrightness += skyLib.Scattered_brightness_nanoLambert(Sb["mag"][i], theta_grid, Sb["ZD"][i])

    if myargs.constellations == "Stars":
         surfBrightness *= 2.2 # flux in cataloge to mag=8 is 2.2x lower than flux in all stars down to mag=21

    #skyBrightness = skyLib.nanoLambert_to_magarc(surfBrightness)
    skyBrightness = surfBrightness / 54. # fraction of natural sky brightness


    print('=====PLOT=====================================================')

    fig = plt.figure(figsize=(8,8))
    ax =  fig.subplots(1,1,subplot_kw={'projection': 'polar'})
    cpLib.initPolPlot(ax)
    ax.set_facecolor("k")
    cmap = "magma"  # _r for magnitude

    # clean up and find min, max
    skyValid = skyBrightness[ np.isfinite(skyBrightness) ]
    zmin = 0.
    zmax = np.percentile(skyValid, 99)*1.1
    print(f'Sky brightness range: {np.min(skyValid):.2e} to {np.max(skyValid):.2e}  (99th pct: {zmax/1.1:.2e})' )


    #plot the sky brightness map
    cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1], skyBrightness,
                        levels=np.linspace(zmin,zmax,100),  # NUMBER OF LEVELS
                        vmin=zmin, vmax=zmax ,
                        extend='both',
                        cmap=cmap)
    


    ax.contour(np.radians(AzEl[0]), 90.-AzEl[1], skyBrightness, linewidths=0.5,colors='gray')
    

    cbar = fig.colorbar(cfd)
    if zmax > 1.:
            print('>1')
            bMin = 0. 
            bMax =  int(zmax    *3. +1)/3.
            print(f'bMax: {bMax}')

            barLen = 0.
            barStep = 3.
            while barLen < 5:
                barStep = barStep / 3.
                barTicks = np.arange(bMin, bMax +.01 , barStep)
                barLen = len(barTicks)
            while barLen > 15:
                barStep = barStep * 3.
                barTicks = np.arange(bMin, bMax +.01 , barStep)
                barLen = len(barTicks)


            print(f'barStep: {barStep}, barLen: {barLen}')
            barTickLabels = [ f'{x:.1f}' for x in barTicks]
            barLabel = 'Sky brightness [unit of natural sky brightness]'

    elif zmax > .009:        
            print('>0.01')    
            bMin = 0. 
            bMax =  int(zmax*100    *3. +1)/3.
            print(f'bMax: {bMax}')

            barLen = 0.
            barStep = 3.
            while barLen < 5:
                barStep = barStep / 3.
                barTicks = np.arange(bMin, bMax +.01 , barStep)/100.
                barLen = len(barTicks)
            while barLen > 15:
                barStep = barStep * 3.
                barTicks = np.arange(bMin, bMax +.01 , barStep)/100.
                barLen = len(barTicks)


            print(f'barStep: {barStep}, barLen: {barLen}')
            barTickLabels = [ f'{x:.1f}' for x in barTicks*100]
            barLabel = 'Sky brightness [% of natural sky brightness]'

    else:
            print('else', zmax)
            bMin = 0.
            logMax = int(np.log10(zmax*1.1)  )
            bMax = 10.**(logMax)        
            barTicks = np.linspace(zmin, zmax, 5)
            barTickLabels = [ f'{x/bMax:.1f}' for x in barTicks]
            barLabel = f'Sky brightness [10$^{{{logMax}}}$  of natural sky brightness]'

    cbar.set_ticks( barTicks )
    cbar.set_ticklabels( barTickLabels)
    cbar.set_label( barLabel )

    # plot the dots
    ax.scatter(Sb["Azr"],Sb["ZD"], s=Sb["dot"], c="white")



    #----------------------------------------------------------
    #All the labels

    if 1:

        #Sun symbol on horizon
        plt.text(np.radians(sunAz), 93.,r'$\odot$', va="center", ha='center') # raw string for LaTeX
        

        #top left corner
        x = -1.
        y = 1.2
        dy = 0.08
        
        cpLib.azlab(ax,x,y,'Observatory: {} Lat.: {:.1f}$^o$'.format(myTel.telescope, myTel.lat))
        y -= dy
        

                
        # bottom left
        x= -1.
        y= -1.08
        cpLib.azlab(ax,x,y,r'$\odot$ Sun:',14)

        y -= dy
        loct = (sunAlpha/15.+12.)%24
        
        loch = int(loct)
        locm = int( (loct-loch)*60.)
        cpLib.azlab(ax,x,y,f'Loc.time: {loch:02d}:{locm:02d}')
        y -= dy
        cpLib.azlab(ax,x,y,r'$\delta: '+f'{sunDelta:.2f}^o$, Elev: {sunElev:.2f}$^o$')
        y -= dy


        # top right
        x=1.
        y=1.2
        cpLib.azlab(ax,x,y,'Constellation:',14)
        y -= dy
        cpLib.azlab(ax,x,y,CONST_name)
        
        y -= dy
        cpLib.azlab(ax,x,y, f'Total: {CONST_totsat:.0f}')

        y -= dy
        cpLib.azlab(ax,x,y, f'In range: {CONST_range:.0f}')

        y -= dy
        cpLib.azlab(ax,x,y, f'Illum.: {CONST_illum:.0f}')




        #bottom right
        x = 1.
        y = -1.08 -dy
        if 1:
            lab = f'Mag. V in [{np.min(Sb["mag"]):.1f},  {np.max(Sb["mag"]):.1f}]'
            cpLib.azlab(ax,x,y,lab)
            y -= dy

            lab = f'Average mag: {np.mean(Sb["mag"]):.1f}'
            cpLib.azlab(ax,x,y,lab)
            y -= dy

            lab = f'Zenith: {skyBrightness[-1,0]:.2g} x natural sky '
            cpLib.azlab(ax,x,y,lab)
            y -= dy


    plt.show()



if __name__ == "__main__":
    main()
