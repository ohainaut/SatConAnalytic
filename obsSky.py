#!/usr/bin/env python3
'''SatConAnalytic - Satellite Constellation Analytic simulations

Plot satellite density over the map of the sky.

Alternatively, the plot can show
- velocity density of the satellites
- the number of satellites / detectable satellites / saturating satellites / 
  non-saturating satellites in the field of view of an instrument
- the effect of the satellites on the observations (%loss)
- the sky brighness increase caused by satellites.

Optionally, a discrete realization of the satellites can be overplotted
(i.e. positions of the satellites as dots)


Interactive:  option -h for detailed usage
Batch: call   obsSky.main(args) with args a list of arguments
See obsSky_batch.py for examples.


'''

import json
import logging
import argparse

import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np


# import ConAn routines
import SatConAnalytic.conan as caLib
import SatConAnalytic.conanplot as cpLib
import SatConAnalytic.satDots as satDots
import SatConAnalytic.skyBrightnessLib as skyLib
import SatConAnalytic.constellations as constLib
from   SatConAnalytic.conanplot import gyrd

log = logging.getLogger('satcon')
log.setLevel(logging.DEBUG)  
for logger_name in ("matplotlib", "matplotlib.pyplot", "matplotlib.font_manager", "matplotlib.backends"):
    logging.getLogger(logger_name).setLevel(logging.WARNING)

    


# some bright stars for debug; RA in hours, Dec in deg; converted to deg later
myStars = np.array([ [17., 20.], [18., 10],
    [19, -10], [20,12], [21, -43],[22, -4],[23, -64],[0,-14],[1,-32],[2,-13.],
    [20.44,-14], [22.1,-20],[1,4],[2,20],
    [23.1,-60],
    [21.39,-11.], [2., -61.]
    ])
myStars = myStars* np.array([-15.,1.])  # convert hours to deg

#---------------------------------------------------------------------------
def create_argument_parser():
    """Create and return the argument parser"""

    parser = argparse.ArgumentParser(description=
    '''Satellite constellations: sky map of satellite trail density or related information. 
    Define the position of the observatory, the position of the Sun, the constellation(s),
    the instrument and its characteristics, and the output required.''')
    # Sun
    parser.add_argument('-d','--deltaSun', default=0.,
                        help="Sun: Declination of the Sun [deg]")
    parser.add_argument('-a','--alphaSun', 
                        help="Sun: Hour Angle of the Sun [deg]. If present, overwrites elevSun")
    parser.add_argument('-e','--elevSun', default=24.,
                        help="Sun: Elevation of the Sun BELOW the horizon. Should probably be >0 in most cases [deg]")
    
    # Constellation
    parser.add_argument('-C','--constellations', default='SLOWGWAK',
                        help="Constellation: ID of the constellation group; 'list' for a list")
    parser.add_argument('--constFile', default="constellations.json",
                        help="Constellation: Which constellation file to use (default: constellations.json)")
    parser.add_argument('--mag550', 
                        help="Constellation: overwrite absolute mag for all satellites [Mag at 550km]")

    # Observatory and instrument parameters
    parser.add_argument('-T','--code', 
                        help='''*Observatory: Predefined telescope/instrument with
                             extptime, FoVl, FoVw, maglim, magbloom, trailf,  
                             telescope, instrument, resolution, latitude. 
                             Use individual options to overwrite presets.''')

    parser.add_argument('-l','--lat', default=-24.6,
                        help="Observatory: Latitude of the observatory [deg]")
    parser.add_argument('-r','--resol', 
                        help="Observatory: Resolution of the instrument (pixel or seeing) [deg]")
    parser.add_argument('-t','--expt', 
                        help="Observatory: Exposure time [s]")
    parser.add_argument('-f','--fovl',
                        help="Observatory: Field of view of the instrument. Length or diametre [deg]")
    parser.add_argument('-w','--fovw', 
                        help="Observatory: Field of view of the instrument. Width. Equal to Length if omitted [deg]")
    parser.add_argument('-m','--maglim', 
                        help="Observatory: Detection limit magnitude of the instrument [5sigma Mag during expTime]")
    parser.add_argument('--magbloom',
                        help="Observatory: Magnitude over which the instrument saturates [Mag for expTime]. Default: -99 (no blooming)")
    parser.add_argument('-k','--trailf', 
                        help="Observatory: Trail filling fraction (width of the trail as fraction of FoV)")
    parser.add_argument('-s','--telescope', 
                        help="Observatory: Name of the telescope")
    parser.add_argument('-i','--instrument', 
                        help="Observatory: Name of the instrument")
    parser.add_argument('-M','--magSelect', default="all", 
                        choices=['all', 'detected', 'oversaturated', 'notDetected'],
                        help="selection on magnitude; default=all")
    
    # output type
    parser.add_argument('-O','--output', default="trails",
                        choices=['trails', 'losses', 'satDens', 'TrailDens', 'skyDiffuse', 'skyScattered'],
                        help="""Output type: 
                        >trails: number of trail per exposure [default];  
                        >sat: number of satellites per exp. (instantaneous, no trailing effect);  
                        >losses: fraction of FoV lost;  
                        >satDens: number of satellites per square degree (instantaneous);  
                        >TrailDens: number of satellite trails per square degree and per second (i.e. satDens * velocity; mostly debug);  
                        >skyDiffuse: surface brightness increase caused by the satellites; consider -M FAINT; see units below;  
                        >skyScattered: surface brightness increase caused by satLight scattering in the atmosphere; use -M ALL; see units below;
                        """  )

    parser.add_argument('-u','--unit', default='frac',
                        choices=['nanoLambert', 'frac', 'magarcsec2', 'microCandelaPerM2'],
                        help="Unit for sky brightness (default: frac); used only for output types skyDiffuse and skyScattered.")
                        



    # plot parameters
    parser.add_argument('--noPlot', action='store_false',
                        help="Plot: Don't generate the plot (debug/batch)")
    
    parser.add_argument('--noconan', action='store_true',
                        help="Plot: Don't plot the conAn simulation; only the sat dots")

    parser.add_argument('--noshade', action='store_true',
                        help="Plot: Don't shade low elevations")
    parser.add_argument('--noRADecGrid', action='store_true',
                        help="Plot: Don't draw RA/Dec grid")
    parser.add_argument('--noscalebar', action='store_false',
                        help="Plot: Don't include scalebar")
    parser.add_argument('--almuc', action='store_true',
                        help="Plot: Write sat count on the almucantars")
    parser.add_argument('--nolabel', action='store_false',
                        help="Plot: Don't label the plot")
    parser.add_argument('--dots', action='store_true',
                        help="Plot: Plot the satellite dots")
    parser.add_argument('--plotStars', action='store_true',
                        help="Plot: Do plot the stars")
    parser.add_argument('--minmax', nargs=2, 
                        help="Plot: min, max value for the colorscale")

    parser.add_argument('--pdf', action='store_true',
                        help="OUTPUT: output file in pdf (default is png)")
    
    log.debug(f'Argument parser created with {len(parser._actions)} actions.')
    return parser




#---------------------------------------------------------------------------
def main(args=None):
    """Main function that can be called with arguments or from command line
    
    Args:
        args: List of arguments (if None, will parse from command line)
    
    Returns:
        dict: Dictionary containing results and output information
    """
    
    log.info('>>>SatConAnalytic: obsSky<<<')
    
    #----- config
    step = 1. #deg. Use 1 in operation;
              # conversions may not work for other values.
              # Smaller values take forever; 10 for debug
    #outpath = "/home/ohainaut/public_html/outsideWorld/"
    outpath = "./"
    
    #----- Parse arguments
    parser = create_argument_parser()
    if args is None:
        myargs = parser.parse_args()
    else:
        myargs = parser.parse_args(args)


    
    # flags
    myargs.plotflag      = myargs.noPlot
    myargs.shadeflag     = myargs.noshade
    myargs.scalebarflag  = myargs.noscalebar
    myargs.labelplotflag = myargs.nolabel
    myargs.almucantar    = myargs.almuc
    if myargs.pdf :
        myargs.outputformat = ".pdf"
    else:
        myargs.outputformat = ".png"
        
    
    # hack for skyFlux
    if myargs.code == "skyFlux":
        myargs.code = "skyMag"
        skyFluxFlag = True
    else:
        skyFluxFlag = False

    #===========================================================================
    #===========================================================================
    # INITIALIZATION
    #===========================================================================
    #===========================================================================
    
    #===========================================================================
    # OBSERVATORY TELESCOPE INSTRUMENT
    #===========================================================================

    log.info('=====TELESCOPE/INSTRUMENT SETUP================================')
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
    
    
    #===========================================================================
    # objects
    #===========================================================================

    if myargs.plotStars:
        log.info('=====STARS SETUP================================')
        log.info(f'Processing {len(myStars)} stars...')


        starsAz, starsEl = caLib.radec2azel(myStars[:,0] +( sunAlpha - 180.), myStars[:,1], myTel.lat)
        for i in range(len(myStars)):
            log.debug(f'  RA={myStars[i,0]:6.2f}d Dec={myStars[i,1]:6.2f}d -> Az={starsAz[i]:6.2f}d El={starsEl[i]:6.2f}d') 
    else:
        log.debug('no Stars:  No stars will be plotted.')


    #===========================================================================
    # CONSTELLATIONS
    #===========================================================================
    log.info('=====CONSTELLATIONS==================================')

    # expand the constellation Id into a list of real constellations
    CONSTELLATIONS = constLib.findConstellations(myargs.constellations, constFile=myargs.constFile)  # constellationS object

    log.debug(f'Found constellations:\n{CONSTELLATIONS.ToC}')
    log.info(f'Constellations: {CONSTELLATIONS.totSat} satellites in {CONSTELLATIONS.totShells} shells.')

    # output file
    outfileroot  = f'{myargs.code}_{myargs.constellations}_{myargs.output}_'
    outfileroot += f'{myargs.magSelect}_{int(myTel.lat):02d}_{int(sunAlpha*10):04d}'
    

    if myargs.mag550 is not None:
        mag550 = float(myargs.mag550)
        for myShell in CONSTELLATIONS.shells:     
            myShell.mag550 = mag550

    #===========================================================================
    #===========================================================================
    # CALCULATIONS
    #===========================================================================
    #===========================================================================
    
    log.info('=====CALCULATIONS===============================================')


    # Azimuth-Elevation mesh:
    AzEl = caLib.fill_AzEl(step)            # mesh of Azimut-Elevation
    elementArea  = caLib.surface_El( AzEl[1], step) 
                                            # area of each element in sq.deg

    # Arrays with the various results; same array geometry as AzEl 
    densSatAll = np.zeros_like(AzEl[0])    # density of satellites     (all sat)
    densVelAll = np.zeros_like(AzEl[0])    # density of satVelocities (all sat)
    densSatObs = np.zeros_like(AzEl[0])    # density of Observable satellites (brighter than limiting mag)
    densVelObs = np.zeros_like(AzEl[0])    #                       satVel
    densSatBloom = np.zeros_like(AzEl[0])  # density of blooming satellites  (brighter than bloom limit)
    densVelBloom = np.zeros_like(AzEl[0])  #                     satVel

    luminanceUnresolved_tot = np.zeros_like(AzEl[0])  # total flux (for all sat)  
    luminanceUnresolved_totNotDet = np.zeros_like(AzEl[0])  # total flux (for faint sat)  

    totMag_element = np.zeros_like(AzEl[0])  # total mag of all sat in element (for debug)
    totFlux_element = np.zeros_like(AzEl[0])  # total flux of all sat in element (for debug)

    nSat_element = np.zeros_like(AzEl[0])  # number of satellites in the element; for debug

    mag_max    = -99.
    mag_min    =  99.
    mageff_max = -99.
    mageff_min =  99.

    mag550_min = 9999.
    mag550_max = -9999.

    amag=0.

    # Scan the constellation shells
    for myShell in CONSTELLATIONS.shells:     
        # keep the brightest and faintest satellite in the shells 
        mag550_min = min(mag550_min, myShell.mag550)
        mag550_max = max(mag550_max, myShell.mag550)

        # model the shell i
        densSi, veli, magi =  myShell.modelOneShell(AzEl,myTel.lat, sunAlpha,sunDelta )

        # out: 
        # -ensSi:  density of illuminated satellites in the element (Nsat/sq.deg)
        # -veli:   velocity of satellites in the element (deg/sec)
        # -magi:   magnitude of satellites in the element (as if illuminated even if not illuminated; for debug)

        # process and integrate the shell:

        # extreme magnitudes
        mag_max = max(mag_max,np.amax(magi))
        mag_min = min(mag_min,np.amin(magi))
        

        # Effective magnitude and extremes
        # does the satellite trail more than 1 resolution element during expT:
        trailing =  veli*myTel.expt >= myTel.resol  
                 # bolean for non-zero trailing
        mageffi = magi.copy()  # init effective mag for non-trailing
        mageffi[trailing] = magi[trailing]   - 2.5*np.log10(myTel.resol/ (veli[trailing]*myTel.expt)  )
                 # correct for trailing sat
        mageff_max = max(mageff_max,np.amax(mageffi))
        mageff_min = min(mageff_min,np.amin(mageffi))


        # all sat:
        densSatAll += densSi
        densVelAll += densSi * veli

        # observable Sat:
        densSobsi = np.copy(densSi)
        densSobsi[ mageffi > myTel.maglim] = 0.
        densSatObs += densSobsi
        densVelObs += densSobsi * veli
        
        # super bright bloomer satellites:
        densSbloomi = np.copy(densSi)
        densSbloomi[ mageffi > myTel.magbloom ] = 0.
        densSatBloom += densSbloomi
        densVelBloom += densSbloomi * veli

        #luminance from the satellites in the Element:
        #   f = lumZP*10^(-0.4*mag)   
        #        # flux from one satellite converted in Cd/m^2
        #        # Correct for 1sat/sq.arcsec
        #     * densSi / 3600^2    # n sat per sq.arc 
        lumini = skyLib.luminanceZP *10.**(-0.4*magi) * densSi /3600.**2 
        luminNotDeti = lumini.copy()
        luminNotDeti[ mageffi <= myTel.maglim] = 0.  # only the faint ones
        luminanceUnresolved_tot       += lumini
        luminanceUnresolved_totNotDet += luminNotDeti

        # total flux of all sat in the element:
        # (essentially the same info as luminance_tot;
        #  Validated: 
        #    difference of -5log10(3600acsec/deg)= 17.78 mag between the two: 
        #    OK)
        # Nsat = densSi [n/sq.deg] * elementArea [sq.deg]
        # totFlux: in W/m2, Vband, total flux in the element.
        totFlux_elementi = skyLib.magnitude_to_flux( magi ) * densSi * elementArea        
        totFlux_element += totFlux_elementi
        
        nSat_elementi =  densSi * elementArea        
        nSat_element += nSat_elementi

        amag += np.average(magi[ densSi>0 ]) 

    totalFlux = np.sum(totFlux_element)
    log.info(f'Total flux: {totalFlux:.2e} Wm-2')
    totalArea = np.sum(elementArea)
    log.info(f'Total area: {totalArea:.2f} sq.deg')
    totalSat = np.sum(nSat_element)
    log.info(f'Total number of satellites: {totalSat:.2f}')
    

    # FINISHED SCANNING THE SHELLS:
    # Now, we have the following arrays defined:
    #   densSatAll: dens in Nsat/sq.deg
    #   densVelAll: density if trails in Ntrail/deg/sec

    #   densSatObs: same, only for sat with mag < myTel.maglim
    #   densVelObs:
    
    #   densSatBloom: same, only for sat with mag < myTel.magbloom
    #   densVelBloom
    
    #   luminanceUnresolved_tot: luminance from the unresolved satellites in element, Cd/m^2. 
    #   luminanceUnresolved_totNotDet: luminance from the unresolved satellites not detected in element, Cd/m^2. 
    #   totFlux: in W/m2, Vband, total flux in the element.


    log.info(f'Constellation density calculation done.')
    log.debug(f'Satellite magnitudes in [{mag_max:.2f},{mag_min:.2f}]')
    log.debug(f'Satellite eff. mag.  in [{mageff_max:.2f},{mageff_min:.2f}]')
    log.debug(f'Following output for zenith: (AzAlt={AzEl[:,-1,-1]})')
    log.debug(f'- Sat density [n/sq.dg] {round( densSatAll[-1,-1],  4)}')
    log.debug(f'- Sat velocity (for last constellation) [deg/s] {round(veli[-1,-1],4)}')  
    log.debug(f'- Diffuse mag               [mag/sq/arcsec]: {skyLib.luminance_to_magarc2(luminanceUnresolved_tot[-1,-1]):.2f}' )
    log.debug(f'- Diffuse mag of notDet sat [mag/sq/arcsec]: {skyLib.luminance_to_magarc2(luminanceUnresolved_totNotDet[-1,-1]):.2f}' )
    


    #===========================================================================
    #===========================================================================
    # PREPARE THE PLOT
    #===========================================================================
    #===========================================================================


    # default colormap; overwritten later for special cases
    cmap = "magma"


    # sat count for almucantars
    almucantars = [60.,30.,20., 10.,0.] #  elevation [deg]
    almucantars = np.arange(85.,-1.,-5.) #  elevation [deg]
    almucantarCounts = caLib.integrateSat(almucantars,AzEl,densSatAll)
            # almucantarCounts: number of sat higher than almucantar

    # initialize integrated values for various elev.
    effect_elev = [0., 20., 30.] # limits at which the effects are computed
    effect_surfTot  = np.zeros_like(effect_elev)
    effect_tot      = np.zeros_like(effect_elev)

    # find indices of Az=270, El=30
    idx_sunAz = np.argmin( np.abs( (AzEl[0][0,:] - sunAz)%360 ) )
    idx_30 = np.argmin( np.abs( AzEl[1][:,0] -30. ) )







    # Calculate local time
    locTime_h = (sunAlpha/15.+12.)%24

    # ==magSelect==
    # Select effective densities depending on the requested output:
    #   ds: density of satellites;
    #   dv: density of satellite trails (i.e. density of satellites * velocity)
    ds, dv, selectionLab =cpLib.magSelect(myargs.magSelect, 
                    densSatAll, densVelAll, 
                    densSatObs, densVelObs, 
                    densSatBloom, densVelBloom,
                    myTel)



    log.info('=====DENSITY CONVERSION========================================')

    # ==output==
    # Compute what to plot (--> logDensity)
    # Set the info for the colorbar.

    #--------------------------------------------------------------------------
    if 0:
        log.info('DEBUG')
        logDensity = -skyLib.flux_to_magnitude( totFlux_element)
        barLabel = 'DEBUG: MAG'

        barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_negMag(logDensity)


    elif myargs.output == "TrailDens":
        log.info("-> Trail Density")
        barLabel = "Trail/deg/sec."
        codeTitle = ["Trail","density"]
        myUnit = "Trail/deg/s"
        logDensity = cpLib.log10Sky( dv )
        barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.LsetBarLim_standardLog(logDensity)

        pick_zenith = dv[-1,0]
        pick_sun30 = dv[idx_30, idx_sunAz]
        log.info(f'At zenith: {pick_zenith:.2e} {barLabel}')
        log.info(f'At Az={sunAz:.1f}d, El=30d: {pick_sun30:.2e} {barLabel}')

    #--------------------------------------------------------------------------
    elif myargs.output == "satDens":
        log.info("-> Sat Density")
        barLabel = "Number of sat./sq.deg."
        codeTitle = ["Satellite","density"]
        myUnit = "sat./sq.deg."

        logDensity = cpLib.log10Sky( ds )
        barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)

        pick_zenith = ds[-1,0]
        pick_sun30 = ds[idx_30, idx_sunAz]
        log.info(f'At zenith: {pick_zenith:.2e} {barLabel}')
        log.info(f'At Az={sunAz:.1f}d, El=30d: {pick_sun30:.2e} {barLabel}')

    #--------------------------------------------------------------------------
    elif myargs.output == "skyScattered":
        log.info("-> skyBrightness for scattered light...")
        codeTitle = ["Scattered","brightness"]

        # total mag  in each element
        magV_element = skyLib.flux_to_magnitude( totFlux_element)


        # Validation: FullMoon, one element with magV=-12.7 at Az=90, El=45, and the rest with magV=99 (i.e. no contribution)
        # if 0:
        #     magV_element = np.full_like(elementArea, 99.)
        #     magV_element[ 45, 90] = -12.7  # for debug: set a bright element at Az=90, El=45

        
        surfBrightness = np.zeros_like(AzEl[0]) # initialize surface brightness array

        # Scattering

        # scan the elements:
        for i_elementEl, elementEl in enumerate(AzEl[1,:,1]): 
                                            # scan elevation rings

            print('.',end='', flush=True)  # progress indicator
            for i_elementAz, elementAz in enumerate(AzEl[0,i_elementEl,:]): 
                                            # scan azimuths

                #print(i_elementEl, elementEl, i_elementAz, elementAz)
                theta_grid = skyLib.angular_distance(
                    AzEl[1], AzEl[0], 
                    elementEl, elementAz)
                surfBrightness += skyLib.Scattered_brightness_nanoLambert(magV_element[i_elementEl, i_elementAz], theta_grid, 90. - elementEl)  
        
        print()  # new line after progress indicator



        # surfBrightness comes in nanoLambert [nL]
        # Unit conversion:
        if  myargs.unit == 'magarcsec2':
            log.info("...SkyMagArcsec2")
            myUnit = "MpSA"
            skyBrightness = skyLib.lambert_to_magarc2(surfBrightness* 1e-9) 
            logDensity = -skyBrightness 
            barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_negMag(logDensity)

        else:
            if myargs.unit == 'frac':
                log.info("...fraction of sky brightness")
                myUnit = 'fraction of dark sky'
                skyBrightness = skyLib.lambert_to_skyFraction(surfBrightness* 1e-9)  # convert from nanoLambert to fraction of sky brightness
            
                logDensity = cpLib.log10Sky( skyBrightness )  

                barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)                            


                print(f'logDensity range: {logMinValue:.2f} to {logMaxValue:.2f}')

                if 0:# logMaxValue > -3: # green-red colormap for the fraction of sky brightness
                    cmap = gyrd    
                    logMinValue = -4.2  # log limits for the colour scale
                    logMaxValue = max(0.1, np.max(logDensity))  # 2.2

                    barMin = int(logMinValue*3.)/3. 
                    barMax = int(logMaxValue*3. -1)/3.
                    barTicks = np.arange(barMin, barMax , .333333)
                    barTickLabels = [ f'{x:.1g}' for x in 10.**barTicks]
                    print(f'logDensity range: {logMinValue:.2f} to {logMaxValue:.2f}')

            else:
        
                if myargs.unit == 'microCandelaPerM2':
                    log.info("...microCandela/m2")
                    myUnit =  r'$\mu$cd/m$^2$'       # raw string for LaTeX

                    skyBrightness =  skyLib.lambert_to_luminance(surfBrightness) /1000.
                    #  1e-9L x 1e6 muCdm-2 -> ./1000.

                elif myargs.unit == 'nanoLambert':
                    log.info("...nanoLambert")
                    myUnit = 'nL'
                    skyBrightness = surfBrightness


                logDensity = cpLib.log10Sky( skyBrightness )  
                barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)                            

        barLabel = 'Surface brightness ['+ myUnit + ']'
        pick_zenith = skyBrightness[-1,0]
        pick_sun30  = skyBrightness[idx_30, idx_sunAz]

        log.info(f'At zenith: {pick_zenith:.2e} {myUnit}')
        log.info(f'At Az={sunAz:.1f}d, El=30d: {pick_sun30:.2e} {myUnit}')



    #--------------------------------------------------------------------------
    elif myargs.output == "skyDiffuse": 
        log.info("-> skyBrightness for unresolved satellites...")
        codeTitle = ["Diffuse","brightness"]

        # luminanceUnresolved_totNotDet in Cd/m^2.
        # Unit conversion:
        if myargs.unit == 'magarcsec2':
            log.info("...SkyMagArcsec2")
            myUnit = "MpSA"
            skyBrightness = skyLib.luminance_to_magarc2( luminanceUnresolved_totNotDet ) 

            logDensity = - skyBrightness
            barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_negMag(logDensity)

        else:
            if myargs.unit == 'frac':
                log.info("...fraction of sky brightness")
                myUnit = 'fraction of dark sky'
                skyBrightness = skyLib.luminance_to_skyFraction( luminanceUnresolved_totNotDet ) 

                logDensity = cpLib.log10Sky( skyBrightness ) 
                
                cmap = gyrd    

                logMinValue = -4.2  # log limits for the colour scale
                logMaxValue = 0.1  # 2.2
                barMin = int(logMinValue*3.)/3. 
                barMax = int(logMaxValue*3. -1)/3.
                barTicks = np.arange(barMin, barMax , .333333)
                barTickLabels = [ f'{x:.1g}' for x in 10.**barTicks]


            else:
                if myargs.unit == 'microCandelaPerM2':
                    log.info("...microCandela/m2")
                    myUnit =  r'$\mu$cd/m$^2$'       # raw string for LaTeX
                    skyBrightness = luminanceUnresolved_totNotDet * 1e6
                        # convert to microCandela/m2 

                elif myargs.unit == 'nanoLambert':
                    log.info("...nanoLambert")
                    myUnit = 'nL'
                    skyBrightness = skyLib.luminance_to_lambert( luminanceUnresolved_totNotDet ) * 1e9
                        # convert to nanoLambert           

                logDensity = cpLib.log10Sky( skyBrightness ) 
                barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)


        barLabel = 'Surface brightness ['+ myUnit + ']'
        pick_zenith = skyBrightness[-1,0]
        pick_sun30  = skyBrightness[idx_30, idx_sunAz]

        log.info(f'At zenith: {pick_zenith:.2e} {myUnit}')
        log.info(f'At Az={sunAz:.1f}d, El=30d: {pick_sun30:.2e} {myUnit}')




    #--------------------------------------------------------------------------
    elif myargs.output == "losses":
        ''' LOSSES accounts for the width of the trail 
        (as a fraction of the FoV),
        For blooming satellites: trailF = 1, full FoV is affected; 
                (1-trailF) bc trailF is already accounted for.
        '''
        log.info("-> fraction lost...")
        log.info("...overrides magSelect with losses")
        codeTitle = ["Fraction","of FoV","lost"]

        ds = myTel.trailf * densSatObs + (1.-myTel.trailf)* densSatBloom
        dv = myTel.trailf * densVelObs + (1.-myTel.trailf)* densVelBloom
        selectionLab  = 'Selection: all satellites, scaled for losses. '
        selectionLab += f'Detected: V$_{{eff}}$ < {myTel.maglim:.1f} '
        selectionLab += f'Bleeding: V$_{{eff}}$ < {myTel.magbloom:.1f}'

        # density: number of satellites in the FoV, 
        #          plus number of trails crossing FoV during the exposure
        #          accounting for trailF and blooming:
        dens    =  ds* myTel.fovl*myTel.fovw + dv * myTel.fovl * myTel.expt


        # deal with empty sky
        logDensity = cpLib.log10Sky(dens)


        # compute average and total losses
        for i in np.arange(0,len( AzEl[1,:,1]) ): # scan elevation rings

            # Compute total losses in ring:
            # Need to multiply each element's density by its surface area
            # At given elevation, all elements have the same surface area
            ringElement_area = caLib.surface_El( AzEl[1,i,0], step)
            ring_losses  = np.sum( dens[i,:]    * ringElement_area )

            # integrate above effect limits:
            for ieffect in np.arange(0,len(effect_elev)):
                if AzEl[1,i,1] >= effect_elev[ieffect]:
                    effect_tot[ieffect]  += ring_losses
                    effect_surfTot[ieffect] += ringElement_area * len(AzEl[0,i,:])

        # compute average losses above effect limits:
        for ieffect in np.arange(0,len(effect_elev)):
            effect_tot[ieffect]   = effect_tot[ieffect]  /effect_surfTot[ieffect]

            log.info(f'losses on exposures (aver above {effect_elev[ieffect]}): Loss fraction: {effect_tot[ieffect]:.5g}/1.')
            log.debug(f'SurfTot: {effect_surfTot[ieffect]:.2f} sq.deg ')


        log.info(f'losses on exposures (at Zenith      ): Loss fraction: {dens[-1,-1]:.5g}/1.')
        
        # bar info:
        barLabel = "Fraction of FoV lost"
        cmap = gyrd    
        
        logMinValue = -3.5  # log limits for the colour scale
        logMaxValue = 0.5  # 2.2
        barMin = int(logMinValue*3.)/3. 
        barMax = int(logMaxValue*3. -1)/3.
        barTicks = np.arange(barMin, barMax , .333333)
        barTickLabels = [ f'{x:.1g}' for x in 10.**barTicks]

        pick_zenith = dens[-1,0]
        pick_sun30  = dens[idx_30, idx_sunAz]
        myUnit = 'fraction lost'



    #--------------------------------------------------------------------------
    elif myargs.output == "trails":
        log.info("-> number of trails/exp...")
        codeTitle = ["Trails per","exposure"]
    
        # density: number of satellites in the FoV, 
        #          plus number of trails crossing FoV during the exposure.
        #          Magnitude selection already applied in ds and dv.
        dens    =  ds* myTel.fovl*myTel.fovw + dv * myTel.fovl * myTel.expt


        # deal with empty sky
        logDensity = cpLib.log10Sky(dens)

        # compute average and total trails
        for i in np.arange(0,len( AzEl[1,:,1]) ): # scan elevation rings

            # Compute total trails in ring:
            # Need to multiply each element's density by its surface area
            # At given elevation, all elements have the same surface area
            ringElement_area  = caLib.surface_El(AzEl[1,i,0], step)
            ring_trails   = np.sum( dens[i,:]    * ringElement_area )

            # integrate above effect limits:
            for ieffect in np.arange(0,len(effect_elev)):
                if AzEl[1,i,1] >= effect_elev[ieffect]:
                    effect_tot[ieffect] += ring_trails
                    effect_surfTot[ieffect]  += ringElement_area * len(AzEl[0,i,:])

        # compute average trails above effect limits:
        for ieffect in np.arange(0,len(effect_elev)):
            effect_tot[ieffect]   = effect_tot[ieffect]  /effect_surfTot[ieffect]

            log.info(f'Trails on exposures (aver above {effect_elev[ieffect]}): Trails: {effect_tot[ieffect]:.5g} trails/exp')
            log.debug(f'SurfTot: {effect_surfTot[ieffect]:.2f} sq.deg ')


        log.info(f'Trails on exposures (at Zenith):  Trails: {dens[-1,-1]:.5g} trails/exp')
        

        # bar info:
        barLabel = "Number of trails per exp."
        barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)

        pick_zenith = dens[-1,0]
        pick_sun30  = dens[idx_30, idx_sunAz]
        myUnit = 'trails/exp'
        log.info(f'At zenith: {pick_zenith:.2e} {myUnit}')
        log.info(f'At Az={sunAz:.1f}d, El=30d: {pick_sun30:.2e} {myUnit}')


    #--------------------------------------------------------------------------
    else:
        log.error(f'Invalid output type {myargs.output}')
        exit(1)

    #======================================================================
    #======================================================================
    #======================================================================
    # PLOT
    #======================================================================
    #======================================================================
    #======================================================================

    # here, we have the following arrays defined:
    #   logDensity: the quantity to plot, in log10 scale;
    #   barLabel: label for the colorbar;
    #   barTicks, barTickLabels, logMinValue, logMaxValue: info for the colorbar scale.





    if myargs.plotflag:
        log.info('=====PLOT==================================================')
        fig = plt.figure(figsize=(10,8))
        ax =  fig.subplots(1,1,subplot_kw={'projection': 'polar'}) 
        cpLib.initPolPlot(ax)
        ax.set_facecolor("k")
        skyColor = cpLib.get_skyColor(sunElev)
        
        #----------------------------------------------------------------------
        # conAn density plot

        if myargs.noconan:
            logDensity = np.zeros_like(AzEl[0]) - 1000.  # empty sky
        

        cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1], logDensity , 
                        levels=np.linspace(logMinValue,logMaxValue,100),  # NUMBER OF LEVELS
                        vmin=logMinValue, vmax=logMaxValue ,
                        extend='both',
                        cmap=cmap)
    
        if cmap == "magma": 
            cfd.cmap = cpLib.cmap_sky(skyColor, cmap_name=cmap)
        cfd.cmap.set_under(skyColor) # below minimum -> black


        #----------------------------------------------------------------------
        # Stars
        if myargs.plotStars:
            ax.scatter( np.radians(starsAz), 90.-starsEl,
                    s=150, c='r', edgecolors='k', marker='*', zorder=5)

        #----------------------------------------------------------------------
        # Scalebar
        if myargs.scalebarflag:
            cbar = fig.colorbar(cfd)
            cbar.set_ticks( barTicks )
            cbar.set_ticklabels( barTickLabels)
            cbar.set_label(barLabel)

            # Title
            x = 1.6
            y = 1.2
            y = cpLib.put_labels(ax,x,y,.1,codeTitle, fontsize=14, weight='bold', va='top', ha='right')


        #----------------------------------------------------------------------
        # airmass shade
        if myargs.shadeflag:
            cpLib.draw_airmassShade(ax)


        #-----------------------------------------------------------------------
        # Draw RA,Dec lines
        if not myargs.noRADecGrid:
            cpLib.drawHADec(myTel.lat)


        #----------------------------------------------------------
        #All the labels
        if myargs.labelplotflag:

            #Sun symbol on horizon
            plt.text(np.radians(sunAz), 93.,r'$\odot$', va="center", ha='center') # raw string for LaTeX
            

            #top left corner
            x = -1.
            y = 1.2
            dy = 0.08
            
            y = cpLib.put_labels(ax,x,y,dy, [f'Observatory: {myTel.telescope}'], fontsize=14, weight='bold')
            y = cpLib.put_labels(ax,x,y,dy, [f'Latitude: {myTel.lat:.1f}$^o$'])
            
            if myargs.output in ['losses', "trails"]:
                instrument_labels = cpLib.prepareLabel_instrument(myTel)
                y = cpLib.put_labels(ax,x,y,dy,instrument_labels)

            if myargs.output in ['skyDiffuse']:
                if myTel.maglim > -9.:
                    y = cpLib.put_labels(ax,x,y,dy,
                        [f'Limiting magnitude: {myTel.maglim:.1f}'])
                                     



                    
            # bottom left
            loch = int(locTime_h)
            locm = int( (locTime_h-loch)*60.)

            x= -1.
            y= -1.
            cpLib.azlab(ax,x,y,r'$\odot$ Sun:',size=14)
            y -= dy

            y = cpLib.put_labels(ax,x,y,dy,
                [f'Loc.time: {loch:02d}:{locm:02d}',
                r'$\delta: '+f'{sunDelta:.2f}^o$, Elev: {sunElev:.2f}$^o$'
                ])




            # top right
            x=1.
            y=1.2
            cpLib.azlab(ax,x,y,'Constellation:',size=14)
            y -= dy
            y = cpLib.put_labels(ax,x,y,dy,
                [f'Total {CONSTELLATIONS.totSat:.0f} sat.',
                CONSTELLATIONS.name,
                ]
                )



            #bottom right
            x = 1.
            y = -1.08 
            


            labels  = ["Sat. magnitudes: V$_{550km}=$" +\
                    f'{mag550_min:3.1f}'] 

            if myargs.output in ['losses', "trails", 'skyDiffuse']:
                labels.append(
                    "  V$_{sat}$"+f" in [{mag_max:.1f}, {mag_min:.1f}]" +
                    "  V$_{eff}$"+ f" in [{mageff_max:.1f}, {mageff_min:.1f}]")
            else:    
                labels.append("  V$_{sat}$"+f" in [{mag_max:.1f}, {mag_min:.1f}]")


            if myargs.output == "effect":
                labels.append("Selection: all satellites, scaled for effect")
                labels.append(
                    f"Detected: V$_{{eff}}$ < {myTel.maglim:.1f} "+
                    f"   Bleeding: V$_{{eff}}$ < {myTel.magbloom:.1f}")
            elif myargs.magSelect in ["oversaturated", "detected", "notDetected"]  : # ie, not "all"
                labels.append(selectionLab)


            labels.append(f'Value at zenith: {cpLib.format_tick(pick_zenith)};' +
                          f' at 30$^0$: {cpLib.format_tick(pick_sun30)} [{myUnit}]' )
            
    
            y = cpLib.put_labels(ax,x,y,dy,labels)

            




        # sat count on almucantars
        if myargs.almucantar:
            for we, wi in zip(almucantars, almucantarCounts):
                cpLib.azlab(ax,
                            -0,(90.-we-5)/90.,
                            '{:.0f} sat.>{:.0f}$^o$:'.format(wi,we),
                            size=9,alpha=0.5*(1.-we/100.))
        

        # Discrete satellites as dots on the plot
        if myargs.dots:
            log.info('=====DOTS===============================================')

            if myargs.output == "effect":
                myargs.magSelect = "all"
            log.info(f'Plotting dots, selection: {myargs.magSelect}')


            Sv = satDots.makeConstellationStatTable(CONSTELLATIONS, 
                                                    sunAlpha, sunDelta, 
                                                    myTel.lat, 
                                                    sunAlpha*240)#(180+sunAlpha)*360.) ## slowed 10x
            



            if myargs.magSelect == "all":
                # not illuminated in grey
                Si = Sv[  ~Sv["bIlluminated"] ]
                ax.scatter(Si["Azr"],Si["ZD"], s=Si["dot"], c="grey", alpha=.3)

            Si = Sv[  Sv["bIlluminated"] ]
            if myargs.magSelect in ["all",  "detected"] :
                # not-bright in yellow
                Sb = Si[ Si["mag"] >= 7 ]

                log.info(f'Plotting {len(Sb)} non-bright satellites (mag >= 7) in yellow')
                ax.scatter(Sb["Azr"],Sb["ZD"], s=Sb["dot"], c="yellow", alpha=0.4)


            if myargs.magSelect in ["all", "detected", "oversaturated"] :
                # bright in red
                Sb = Si[ Si["mag"] < 7 ]
                log.info(f'Plotting {len(Sb)} bright satellites (mag <7 ) in orange')
                ax.scatter(Sb["Azr"],Sb["ZD"], s=Sb["dot"], c="orange")

                Sb = Si[ Si["mag"] < 6 ]
                log.info(f'Plotting {len(Sb)} bright satellites (mag < 6) in red')
                ax.scatter(Sb["Azr"],Sb["ZD"], s=Sb["dot"], c="red")


        if not myargs.noconan and not myargs.dots:
            # mark the position of the Zenith and of the 30deg elevation on the plot
            ax.scatter(np.radians(AzEl[0, idx_30, idx_sunAz]), 90.-AzEl[1, idx_30, 
            idx_sunAz], s=100, c="red", marker='*', edgecolors='k', linewidths=.5, label='Sun')
            ax.scatter(0., 0., s=100, c="red", marker='*', edgecolors='k', linewidths=.5, label='30$^o$')



        # save plot
        fig.tight_layout()
        plt.savefig(outpath+'w.png')
        filename = outpath+outfileroot+myargs.outputformat

        plt.savefig(filename)
        log.info(f'Plot saved as {filename}')
        #-- end plot
    else:
        log.debug('noPlot: No plot will be generated.')

 
    results = {
        'outfileroot': outfileroot,
        'meta': {
            'constellations': CONSTELLATIONS.list,
            'const_stats': CONSTELLATIONS.totSat,
            'const_shells' : CONSTELLATIONS.totShells,
            'const_nConst': CONSTELLATIONS.totConst,
            'const_v550km': mag550_min,
            'telescope': myTel.ToC,
            'step': step
        },
        'sunAlpha': sunAlpha,
        'sunDelta': sunDelta,
        'sunElev': sunElev,
        'local_time': locTime_h,
        'value' :{
            'zenith': pick_zenith,
            'sun30': pick_sun30,
            'unit': myUnit,
            'type': myargs.output},
        'sat': {
            'elev': almucantars,
            'count': almucantarCounts,
            'unit': 'integrated number of satellites above each elev'}
    }
    
    if myargs.output in ['losses', "trails"]:
        results.update({    
            'effect': {
                'elev': effect_elev,
                'value_elev': effect_tot,
                'value_zenith': dens[-1,-1],
                'surfTot': effect_surfTot,
                'doc': 'Effect on exposures, averaged above each effect limit',
                'type' : f'Effect type: {myargs.output}'
            }
        })
    
    log.info("...done.")
    log.debug(f'Results: {results}')


    return results



if __name__ == "__main__":

    cpLib.init_logger(log)
    log.info('===obsSky===')

    results = main()
    with open("obsSky.json", "w") as f:
        json.dump(results, f, indent=4, cls=cpLib.NumpyEncoder) 
