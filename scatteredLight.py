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

import sys
import logging
import argparse
import json

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from astropy.table import Table
import numpy as np

import SatConAnalytic.conanplot as cpLib
import SatConAnalytic.skyBrightnessLib as skyLib
import SatConAnalytic.conan as caLib
import SatConAnalytic.satDots as satDots
import SatConAnalytic.constellations as constLib


log = logging.getLogger('satcon')
for logger_name in ("matplotlib", "matplotlib.pyplot", "matplotlib.font_manager", "matplotlib.backends"):
    logging.getLogger(logger_name).setLevel(logging.WARNING)




def create_argument_parser():
    """Create and return the argument parser"""
    parser = argparse.ArgumentParser(description=
    '''Sky map of the light scattered by the satellites.
    Define the position of the observatory, the position of the Sun, the constellation(s), and the output required.''')
    parser.add_argument('-d','--deltaSun', default=0.,
                        help="Sun: Declination of the Sun [deg]")
    parser.add_argument('-a','--alphaSun',
                        help="Sun: Hour Angle of the Sun [deg]. If present, overwrites elevSun")
    parser.add_argument('-e','--elevSun', default=24.,
                        help="Sun: Elevation of the Sun BELOW the horizon. Should probably be >0 in most cases [deg]")
    parser.add_argument('-C','--constellations', default='SLOWGWAK',
                        help="ID of the constellation group; list for a list, or 'FullMoon' or 'Stars' for validation cases")
    
    parser.add_argument('--constFile', default="constellations.json",
                        help="Which constellation file to use (default: constellations.json)")
    parser.add_argument(     '--anom0',    help='''Initial anomaly [deg]''')

    
    parser.add_argument('-m','--magOffset', default=0.,
                        help="magnitude offset to apply to all satellites (e.g., to simulate different satellite brightness) [mag]")

    parser.add_argument('-T','--code',
                        help='''Observatory: Predefined telescope/instrument with telescope, instrument, latitude.
                        Use individual options to overwrite presets.
    ''')
    parser.add_argument('-l','--lat', default=-24.6,
                        help="Observatory: Latitude of the observatory [deg]")
    parser.add_argument('-s','--telescope',
                        help="Observatory: Name of the telescope")
    parser.add_argument('-i','--instrument',
                        help="Observatory: Name of the instrument")

    parser.add_argument('-u','--unit', default='frac',
                        choices=['nanoLambert', 'frac', 'magarcsec2', 'microCandelaPerM2'],
                        help="Unit for sky brightness (default: frac)")

    parser.add_argument('--pdf', action='store_true',
                        help="Plot: output file in pdf (default is png)")

    return parser

#---------------------------------------------------------------------
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

    if myargs.anom0 is not None:
        anom0 = float(myargs.anom0)
    else:
        anom0 = 0.

    #===========================================================================
    # OBSERVATORY TELESCOPE INSTRUMENT
    #===========================================================================

    log.info('=====TELESCOPE/INSTRUMENT SETUP================================')

    myargs.expt = None
    myargs.resol = None
    myargs.fovl = None
    myargs.fovw = None

    myTel = cpLib.getTelescope(myargs)
    log.debug(myTel)

    #===========================================================================
    # SUN
    #===========================================================================
    log.info('=====SUN SETUP================================')    
    sunAlpha, sunDelta, sunElev = caLib.consolidate_sun(
        myargs.alphaSun, myargs.deltaSun, myargs.elevSun, myTel.lat)
    
    log.info('Sun position:')
    log.info(f'\tLocal time: {((180+sunAlpha)/15.)%24:.2f}h')
    log.info(f'\tHA = {sunAlpha:.1f}deg  = {(sunAlpha/15.)%24:.2f}h, Dec = {sunDelta:.1f}d')
    log.info(f'\tElevation: {sunElev:.2f}d')

    # get Azimuth (and check Elevation)    
    sunAz,_ = caLib.radec2azel(sunAlpha, sunDelta, myTel.lat)
    log.info(f'\tAzimuth: {sunAz:.2f}d, Elevation (validation): {_:.2f}d\n')
    
    


    # Azimuth-Elevation mesh:
    AzEl = caLib.fill_AzEl(step)               # mesh of Azimut-Elevation
    surfBrightness_nL = np.zeros_like(AzEl[0]) # initialize surface brightness array


    #===========================================================================
    # OBJECTS
    #===========================================================================

    if myargs.constellations == "FullMoon":
        # Validation: Full Moon at anti-sun position         
        log.info('=====FULLMOON===============================================')

        Sb = Table()
        Sb["mag"] = [-12.74]
        Sb["dot"] = [50.]
        Sb["Azr"] = [np.radians(sunAz +180.)]  # opposite to the sun
        Sb["ZD"]  = [90. + sunElev]
        Sb["bIlluminated"] = [True]
        Sb["bVis"] = [True]

        CONST_name = "FullMoon"
        CONST_totsat = 1
        CONST_illum = 1
        CONST_range = 1

        log.info(f'FullMoon: {Sb["mag"][0]:.2f} mag at Az: {np.degrees(Sb["Azr"][0]):.2f}, ZD: {Sb["ZD"][0]:.2f}')

        totalFlux = np.sum(skyLib.magnitude_to_flux(Sb["mag"]))
    elif myargs.constellations == "Stars":
        # Validation: Bright stars
        log.info('=====STARS=================================================')
        Sb = skyLib.loadBrightStars(obsLat=myTel.lat, sidTime=10.)

        CONST_name = "Stars"
        CONST_totsat = len(Sb)
        CONST_illum = len(Sb)
        CONST_range = len(Sb)
        log.info(f'Loaded {len(Sb)} bright stars from the catalog')


    else:
        # Normal case: satellites
        log.info('=====SATELLITE CONSTELLATION===============================')

        CONSTELLATIONS = constLib.findConstellations(myargs.constellations, constFile=myargs.constFile)
        log.debug( CONSTELLATIONS.ToC )
        log.debug(" ")
        log.debug( CONSTELLATIONS)
        log.debug('-------------------------------------------')
        CONST_name = CONSTELLATIONS.name
        CONST_totsat = CONSTELLATIONS.totSat

        # get the satellite positions and brightnesses
        Sv = satDots.makeConstellationStatTable(CONSTELLATIONS, 
                                        sunAlpha, sunDelta, 
                                        myTel.lat, 
                                        sunAlpha,
                                        anom0=anom0)
        # Sv returned with:
        #  bIlluminated, bVis, Delta, Azr, ZD, mag, dot     
        # 

        CONST_range = len(Sv)
        Sv["mag"] = Sv["mag"] - float(myargs.magOffset)  # apply mag offset
        
        Sb = Sv[  Sv["bIlluminated"] ]
        CONST_illum = len(Sb)



        log.info(f'Total satellites illuminated: {len(Sb)}')
        if len(Sb)>0:
             log.info(f'Average mag: {np.mean(Sb["mag"]):.2f}, min: {np.min(Sb["mag"]):.2f}, max: {np.max(Sb["mag"]):.2f}')
    

    totalFlux = np.sum(skyLib.magnitude_to_flux(Sb["mag"]))
    log.info(f'Total flux: {totalFlux:.2e} Wm-2')

    #===========================================================================
    # SKY BRIGHTNESS CALCULATION
    #===========================================================================
    log.info('=====SKY BRIGHTNESS CALCULATION===============================')
    for i in range(len(Sb)): # scan the sources (satellites or stars)

        # compute the angular distance between the source and each point in the AzEl grid
        theta_grid = skyLib.angular_distance(AzEl[1], AzEl[0], 90.-Sb["ZD"][i], np.degrees(Sb["Azr"][i]))

        # compute the brightness contribution of this source to each point in the grid, and add it to the total surface brightness
        surfBrightness_nL += skyLib.Scattered_brightness_nanoLambert(Sb["mag"][i], theta_grid, Sb["ZD"][i])


    # find indices of Az=270, El=30
    idx_sunAz = np.argmin( np.abs( (AzEl[0][0,:] - sunAz)%360 ) )
    idx_30 = np.argmin( np.abs( AzEl[1][:,0] -30. ) )

    log.debug(f'Indices for Az=270, El=45: {idx_sunAz}, {idx_30}')
    log.debug(f'AzEl: {AzEl[:,idx_30,idx_sunAz]}')
    log.debug(f'Sky brightness at Az={sunAz}, El=30: {surfBrightness_nL[idx_30, idx_sunAz]:.2e} nanoLambert')

    if myargs.constellations == "Stars":
         surfBrightness_nL *= 2.2 # flux in cataloge to mag=8 is 2.2x lower than flux in all stars down to mag=21


    # Sky brightness converted in the requested unit

    # surfBrightness comes in nanoLambert [nL]
    # Unit conversion:
    if  myargs.unit == 'magarcsec2':
        log.info("...SkyMagArcsec2")
        myUnit = "MpSA"
        skyBrightness = skyLib.lambert_to_magarc2(surfBrightness_nL* 1e-9) 
        logDensity = -skyBrightness 
        barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_negMag(logDensity)
    else:
        if myargs.unit == 'frac':
            log.info("...fraction of sky brightness")
            myUnit = 'fraction of dark sky'
            skyBrightness = skyLib.lambert_to_skyFraction(surfBrightness_nL* 1e-9)  # convert from nanoLambert to fraction of sky brightness
        
            logDensity = cpLib.log10Sky( skyBrightness )  

            barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)                            

            print(f'logDensity range: {logMinValue:.2f} to {logMaxValue:.2f}')

        else:
    
            if myargs.unit == 'microCandelaPerM2':
                log.info("...microCandela/m2")
                myUnit =  r'$\mu$cd/m$^2$'       # raw string for LaTeX

                skyBrightness =  skyLib.lambert_to_luminance(surfBrightness_nL)/1000.  #  1e-9L x 1e6 muCdm-2 -> ./1000.

            elif myargs.unit == 'nanoLambert':
                log.info("...nanoLambert")
                myUnit = 'nL'
                skyBrightness = surfBrightness_nL
            else:
                log.error(f'Unknown unit: {myargs.unit}')
                return

            logDensity = cpLib.log10Sky( skyBrightness )  
            barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)                            

    barLabel = 'Surface brightness ['+ myUnit + ']'
    pick_zenith = skyBrightness[-1,0]
    pick_sun30  = skyBrightness[idx_30, idx_sunAz]

    log.info(f'At zenith: {pick_zenith:.2e} {myUnit}')
    log.info(f'At Az={sunAz:.1f}d, El=30d: {pick_sun30:.2e} {myUnit}')


    log.info('=====PLOT=====================================================')
    fig = plt.figure(figsize=(10,8))
    ax =  fig.subplots(1,1,subplot_kw={'projection': 'polar'})
    cpLib.initPolPlot(ax)
    ax.set_facecolor("k")
    cmap = "magma"  # _r for magnitude

    #plot the sky brightness map

    skyBrightness = cpLib.log10Sky(skyBrightness)

    cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1], skyBrightness,
                        levels=np.linspace(logMinValue, logMaxValue, 100),  
                        vmin=logMinValue, vmax=logMaxValue, extend='both',
                        cmap=cmap)
    


    cbar = fig.colorbar(cfd, ticks=ticker.MaxNLocator(nbins=7))
    cbar.set_label( f'Sky brightness [{myUnit}]' )

    cbar.set_ticks( barTicks )
    cbar.set_ticklabels(barTickLabels)

    # # plot the dots
    ax.scatter(Sb["Azr"],Sb["ZD"], s=Sb["dot"], c="white")


    # mark the position of the Zenith and of the 30deg elevation on the plot
    ax.scatter(np.radians(AzEl[0, idx_30, idx_sunAz]), 90.-AzEl[1, idx_30, 
    idx_sunAz], s=100, c="red", marker='*', edgecolors='k', linewidths=.5, label='Sun')
    ax.scatter(0., 0., s=100, c="red", marker='*', edgecolors='k', linewidths=.5, label='30$^o$')


    #----------------------------------------------------------
    # LABELS 

    #Sun symbol on horizon
    plt.text(np.radians(sunAz), 93.,r'$\odot$', va="center", ha='center') # raw string for LaTeX
    
    if 1:

        #top left corner
        x = -1.
        y = 1.2
        dy = 0.08
        
        y = cpLib.put_labels(ax,x,y,dy, [f'Observatory: {myTel.telescope}'], fontsize=14, weight='bold')
        y = cpLib.put_labels(ax,x,y,dy, [f'Latitude: {myTel.lat:.1f}$^o$'])
            
            
        # top right
        x=1.1
        y=1.2
        if myargs.constellations == "FullMoon":
            lab = 'Full Moon'
            cpLib.azlab(ax,x,y,lab,size=14)
        elif myargs.constellations == "Stars":
            lab = 'Bright Stars'
            cpLib.azlab(ax,x,y,lab,size=14)
        else:
            cpLib.azlab(ax,x,y,'Constellation:',fontsize=14, weight='bold')
            y -= dy
            y = cpLib.put_labels(ax,x,y,dy,
                [f'Total: {CONST_totsat:.0f}',
                f'In range: {CONST_range:.0f}',
                f'Illum.: {CONST_illum:.0f}']
                )

            y = cpLib.put_labels(ax,x,y,dy,
                [f'Total {CONSTELLATIONS.totSat:.0f} sat.',
                CONSTELLATIONS.name,
                ]
                )
                
        # bottom left
        x= -1.
        y= -1.08
        cpLib.azlab(ax,x,y,r'$\odot$ Sun:',fontsize=14, weight='bold')
        y -= dy
        loct = (sunAlpha/15.+12.)%24
        
        loch = int(loct)
        locm = int( (loct-loch)*60.)
        lab = f'Loc.time: {loch:02d}:{locm:02d}'
        cpLib.azlab(ax,x,y,lab)
        y -= dy


        lab = r'$\delta: '+f'{sunDelta:.2f}^o$, Elev: {sunElev:.2f}$^o$'
        cpLib.azlab(ax,x,y,lab)
        y -= dy

        lab= f'Observatory: {myTel.telescope} Lat.: {myTel.lat:.1f}$^o$'
        cpLib.azlab(ax,x,y,lab)
        y -= dy
  

    if 1:

        #bottom right
        x = .95
        y = -1.08
        cpLib.azlab(ax,x,y,'Scattered light:',fontsize=14, weight='bold')
        y -= dy
        y = cpLib.put_labels(ax,x,y,dy,[
                  f'Z: {cpLib.format_tick(skyBrightness[-1,0])} {myUnit}',
                  f'30$^0$: {cpLib.format_tick(skyBrightness[idx_30, idx_sunAz])} {myUnit}'
        ])

    if 1:
        plt.tight_layout()

        outfileroot  = f'skyBright_{myargs.code}_{myargs.constellations}_'
        outfileroot += f'{int(myTel.lat):02d}_{int(sunAlpha*10):04d}'

        plt.savefig(f'{outfileroot}.{"pdf" if myargs.pdf else "png"}', dpi=300)
        plt.savefig('w.png', dpi=300)
        #plt.show()



    results = {
        'CONST_name': CONST_name,
        'CONST_totsat': CONST_totsat,
        'CONST_range': CONST_range,
        'CONST_illum': CONST_illum,
        'sunAlpha': float(sunAlpha),
        'sunDelta': float(sunDelta),
        'sunElev': float(sunElev),
        'skyBrightness_zenith': float(skyBrightness[-1,0]),
        'skyBrightness_30': float(skyBrightness[idx_30, idx_sunAz]),
        'unitLabel': myUnit,
    }
    return results

#=================================================================
#=================================================================
#=================================================================
if __name__ == "__main__":

    # Configure the root logger to DEBUG
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # Format
    log_format = logging.Formatter('[%(levelname)-8s %(name)s/%(funcName)s] %(message)s')

    # FileHandler for DEBUG and above
    file_handler = logging.FileHandler('satCon.log')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(log_format)

    # create a StreamHandler for INFO and above
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(log_format)

    # Attach handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    
    log.info('===skyBrightness===')

    results = main()
    with open("scatteredLight.json", "w") as f:
        json.dump(results, f, indent=4, cls=cpLib.NumpyEncoder) 

