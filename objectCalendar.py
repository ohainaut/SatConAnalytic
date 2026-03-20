#!/usr/bin/env python3
# SatConAnalytic - Satellite Constellation Analytic simulations
#
# plot the satellite density over a calendar for a given object

import logging
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt

import SatConAnalytic.conan as caLib
import SatConAnalytic.conanplot as cpLib
import SatConAnalytic.constellations as constLib
import SatConAnalytic.skyBrightnessLib as skyLib

from   SatConAnalytic.conanplot import gyrd

log = logging.getLogger('satcon')
log.setLevel(logging.DEBUG)  
for logger_name in ("matplotlib", "matplotlib.pyplot", "matplotlib.font_manager", "matplotlib.backends"):
    logging.getLogger(logger_name).setLevel(logging.WARNING)

    
#outpath = "/home/ohainaut/public_html/outsideWorld/"
outpath = "./"

#---------------------------------------------------------------------------
def create_argument_parser():
    """Create and return the argument parser"""
    #--- command line arguments
    parser = argparse.ArgumentParser(description='Constellation density calendar for an object')

    # Object  position
    parser.add_argument('-a','--RA',         default=0.,
                        help='''Right Ascension [deg]''')
    parser.add_argument('-d','--DEC',        default=0.,
                        help='''Declination [deg]''')
    
    parser.add_argument('-n','--objlabel',  default="",
                        help='''Name of the object for label''')
    
    # Constellation
    parser.add_argument('-C','--constellations', default='SLOWGWAK',
                        help="Constellation code (or meta-code); list for a list.")
    parser.add_argument('--constFile', default="constellations.json",
                        help="Constellation: Which constellation file to use (default: constellations.json)")
    parser.add_argument('--mag550', 
                        help="Constellation: overwrite absolute mag for all satellites [Mag at 550km]")


    # observatory
    parser.add_argument('-T','--code',  default="FORSimg",
                        help='''Observatory: Telescope/Instrument code''')
    parser.add_argument('-l','--lat',   
                        help='''Observatory: Latitude of the observatory [deg] (OVERWRITE preset)''')
    parser.add_argument('-t','--expt',  
                        help='''Observatory: Individual exposure time [s] (OVERWRITE preset)''')
    parser.add_argument('-r','--resol', 
                        help='''Observatory: Resolution element (seeing, pixel) [deg] (OVERWRITE preset)''')
    parser.add_argument('-f','--fovl',  
                        help='''Observatory: Length of the field-of-view [deg] (OVERWRITE preset)''')
    parser.add_argument('-w','--fovw',  help='''Width of the field-of-view [deg] (Default=Fovl; OVERWRITE preset)''')
    parser.add_argument('-k','--trailf',
                        help='''Observatory: Fraction of the exposure destroyed by a trail (1=full) (OVERWRITE preset)''')
    parser.add_argument('-m','--maglim',     
                        help='''Observatory: Limiting magnitude [mag] (detection limit for expTime) (OVERWRITE preset)''')
    parser.add_argument('--magbloom',   
                        help='''Observatory: Saturation magnitude [mag]. Brighter object destroy the full exposure: their trailf=1 (OVERWRITE preset)''')
    parser.add_argument(     '--instrument', 
                        help='''Observatory: Name of the instrument for label (OVERWRITE preset)''')
    parser.add_argument(     '--telescope',  help='''Observatory: Name of the telescope for label (OVERWRITE preset)''')


    parser.add_argument('-M','--magSelect', default="all", 
                        choices=['all', 'detected', 'oversaturated', 'notDetected'],
                        help="selection on magnitude; default=all")

    parser.add_argument('-O','--output', default="trails",
                        choices=['trails', 'losses', 'satDens', 'TrailDens', 'skyDiffuse'],
                        help="""Output type: 
                        >trails: number of trail per exposure [default];  
                        >sat: number of satellites per exp. (instantaneous, no trailing effect);  
                        >losses: fraction of FoV lost;  
                        >satDens: number of satellites per square degree (instantaneous);  
                        >TrailDens: number of satellite trails per square degree and per second (i.e. satDens * velocity; mostly debug);  
                        >skyDiffuse: surface brightness increase caused by the satellites; consider -M FAINT; see units below;  
                        >skyScattered: NOT SUPPORTED
                        """  )
    parser.add_argument('-u','--unit', default='frac',
                        choices=['nanoLambert', 'frac', 'magarcsec2', 'microCandelaPerM2'],
                        help="Unit for sky brightness (default: frac); used only for output types skyDiffuse and skyScattered.")
                        


    parser.add_argument('--pdf', action='store_true',
                        help="OUTPUT: output file in pdf (default is png)")

    return parser

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
def main(args=None):
    """Main function to run the object calendar plotting"""

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
        



    # DEFAULTS
    elevlim = 20.        # consider only elevations >= elevlim
        
    # no default:
    rao      = float(myargs.RA)%360
    deo      = float(myargs.DEC)
    objlabel = myargs.objlabel

    #- find telescope
    log.info('=====TELESCOPE/INSTRUMENT SETUP================================')
    myTel = cpLib.getTelescope(myargs)
    log.info(myTel)

    #---

    #===========================================================================
    # CONSTELLATIONS
    #===========================================================================
    log.info('=====CONSTELLATIONS==================================')
    CONSTELLATIONS = constLib.findConstellations(myargs.constellations, constFile=myargs.constFile)  # constellationS object

    log.debug(f'Found constellations:\n{CONSTELLATIONS.ToC}')
    log.info(f'Constellations: {CONSTELLATIONS.totSat} satellites in {CONSTELLATIONS.totShells} shells.')


    #===========================================================================
    # CALENDAR
    #===========================================================================
    log.info('=====CALENDAR======================================')

    hstep = 12 # points per hour
    times = np.zeros((365,24*hstep))
    for i in range(365):
        for j in range(24*hstep):
            times[i,j] = 2459580.5 + i + j/24./hstep - (.5 ) # to centre midnight
            

    # Loc ST
    lst = caLib.siderealTimeDeg(times)

    # sun coordinates
    walphas, deltas = caLib.get_sun(times)
    alphas = lst - walphas # alphas is the HA
    elevs = caLib.radec2elev(alphas,deltas,myTel.lat)


    # object coordinates to Az,Elevation
    hao = lst - rao # we work in hour angle
    AzEl = np.array( caLib.radec2azel(hao,deo,myTel.lat))


    densSatAll = np.zeros_like(AzEl[0]) # density of satellites (for all sat, no mag selection)
    densVelAll = np.zeros_like(AzEl[0]) # density of satellite trails (for all sat, no mag selection)
    densSatObs = np.zeros_like(AzEl[0]) # density of observable satellites
    densVelObs = np.zeros_like(AzEl[0]) # density of observable satellite trails
    densSatBloom = np.zeros_like(AzEl[0]) # density of super bright bloomers
    densVelBloom = np.zeros_like(AzEl[0]) # density of super bright bloomer trails
    luminanceUnresolved_tot = np.zeros_like(AzEl[0])  
                                        # total luminance (for all sat)  
    luminanceUnresolved_totNotDet = np.zeros_like(AzEl[0])  
                                        # total luminance 
                                        # (for faint undetected sat)  
    mageffmax = -99
    mageffmin = 99.


    for myShell in CONSTELLATIONS.shells:                 
        # process each shell

        densSi, veli, magi =  myShell.modelOneShell(
            AzEl,myTel.lat, 
            np.reshape(alphas,alphas.shape[0]*alphas.shape[1]),
            np.reshape(deltas,deltas.shape[0]*deltas.shape[1])  )

        # all sat:
        densSatAll += densSi
        densVelAll += densSi * veli


        # effective magnitude
        mageffi = magi  - 2.5*np.log10(myTel.resol/veli/myTel.expt)
        mageffmax = max(mageffmax,np.amax(mageffi))
        mageffmin = min(mageffmin,np.amin(mageffi))

        # only observable ones
        densSobsi = np.copy(densSi)
        densSobsi[ mageffi > myTel.maglim] = 0.
        densSatObs += densSobsi
        densVelObs += densSobsi * veli
        
        # only super bright bloomers
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





    cmap = 'magma'


    # == magnitude selection ==
    # Select effective densities depending on the requested output:
    #   ds: density of satellites;
    #   dv: density of satellite trails (i.e. density of satellites * velocity)
    ds, dv, selectionLab =cpLib.magSelect(myargs.magSelect, 
                    densSatAll, densVelAll, 
                    densSatObs, densVelObs, 
                    densSatBloom, densVelBloom,
                    myTel)


#%#% IN
#%Expected:     dens =  ds* myTel.fovl*myTel.fovw + dv * myTel.fovl * myTel.#%expt    

    # == output selection ==
    # depending on the Output type,
    # set the corresponding logDensity to plot, 
    # the colour bar label, and the title of the plot.
    
    if  myargs.output == "TrailDens":
        log.info("-> Trail Density")
        barLabel = "Trail/deg/sec."
        codeTitle = ["Trail","density"]
        myUnit = "Trail/deg/s"
        logDensity = cpLib.log10Sky( dv )
        barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)

    #--------------------------------------------------------------------------
    elif myargs.output == "satDens":
        log.info("-> Sat Density")
        barLabel = "Number of sat./sq.deg."
        codeTitle = ["Satellite","density"]
        myUnit = "sat./sq.deg."

        logDensity = cpLib.log10Sky( ds )
        barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)

    #--------------------------------------------------------------------------
    elif myargs.output == "skyScattered":
        log.info("-> skyBrightness for scattered light...")
        msg = "skyScattered NOT SUPPORTED"
        log.error(msg)
        return msg

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
            
            logDensity[ luminanceUnresolved_totNotDet == 0.] = logMinValue 
            
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
        logDensity    =  ds* myTel.fovl*myTel.fovw + dv * myTel.fovl * myTel.expt


        # deal with empty sky
        logDensity = cpLib.log10Sky(logDensity)

        # bar info:
        barLabel = "Fraction of FoV lost"
        cmap = gyrd    
        
        logMinValue = -3.5  # log limits for the colour scale
        logMaxValue = 0.  # 2.2

        linTicks = np.array([0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.])
        barTicks = np.array(np.log10(linTicks))
        barTickLabels = [ "{:.2g}".format(x) for x in linTicks]


    #--------------------------------------------------------------------------
    elif myargs.output == "trails":
        log.info("-> number of trails/exp...")
        codeTitle = ["Trails per","exposure"]
    
        # density: number of satellites in the FoV, 
        #          plus number of trails crossing FoV during the exposure.
        #          Magnitude selection already applied in ds and dv.
        logDensity    =  ds* myTel.fovl*myTel.fovw + dv * myTel.fovl * myTel.expt


        # deal with empty sky
        logDensity = cpLib.log10Sky(logDensity)

        # bar info:
        barLabel = "Number of trails per exp."
        barTicks, barTickLabels, logMinValue, logMaxValue = cpLib.setBarLim_standardLog(logDensity)
        
    else:
        msg = ("Invalid output type. Valid options: trails, losses, satDens, TrailDens, skyDiffuse, skyScattered")
        log.error(msg)
        return msg


    #===========================================================================
    # PLOT
    #===========================================================================
    log.info('=====PLOT==========================================')

    # plot setup
    fig = plt.figure(figsize=(12,7))
    plt.rc('font',      size=15) #controls default text size
    plt.rc('axes', titlesize=15) #fontsize of the title
    plt.rc('axes', labelsize=18) #fontsize of the x and y labels
    plt.rc('xtick',labelsize=15) #fontsize of the x tick labels
    plt.rc('ytick',labelsize=15) #fontsize of the y tick labels
    plt.rc('legend',fontsize=15) #fontsize of the legend
    #plt.rcParams['figure.figsize'] = [12, 8]

    # Plot box and ticks and labels
    ax = fig.subplots()
    ax.set_xticks(np.arange(0.,hstep*24, hstep*3))
    ax.set_xticklabels(np.concatenate((np.arange(12,24,3),np.arange(0,12,3))))
    ax.set_xlabel('Local Solar Time')
    ax.set_yticks(np.arange(0.,365,30.5))
    ax.set_yticklabels(['Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'])
    ax.set_ylabel('Calendar')
    ax.set_xlim(4*hstep, 20*hstep) # from 16h to 8h

    # plot content

    # filter logDensity for unobservable areas (elev < elevlim) to avoid colour scale issues; we will fill these areas in grey later
    filteredLogDensity = np.copy(logDensity)

    #- fill in the satellites
    csat = ax.contourf(filteredLogDensity, 
                       levels=np.arange(logMinValue,logMaxValue,.01) , 
                       vmin=logMinValue, vmax=logMaxValue,
                       cmap=cmap, extend='both')

    csat.cmap.set_under('k') # below minimum -> black

    #- unobservable: fill in grey
    osat = np.nan_to_num(logDensity*0) + 1000.*( AzEl[1] < elevlim ) 
    ax.contourf(osat, levels=np.arange(999.,1001.) , cmap='Greys', alpha=1.)

    #- airmass
    wel =   AzEl[1]  *( elevs < 0 )
    cobj = ax.contour(wel, levels=np.arange(elevlim,90.,10.), cmap="summer_r")
    lobj = plt.clabel(cobj, fmt='%.0f$^o$')

    #- daylight
    for myElevMin, myAlpha  in zip([-18., -12., -6., 0.], [.3,.3,.5,1.]):
        ax.contourf(elevs, levels=[myElevMin,90.], 
                    colors='royalblue', alpha=myAlpha)

    #- twilightss
    csun = ax.contour(elevs, levels=[-18.,-12.,-6.,0],
                    linewidths=[1.,1.,1.,5.], 
                    linestyles='solid', 
                    colors='b')
    lsun = plt.clabel(csun, fmt='%.0f$^o$')



    # color bar
    cbar = fig.colorbar(csat)
    cbar.set_ticks(barTicks)
    cbar.set_ticklabels(barTickLabels)

    #cbar.set_ticks(np.log10(np.array([0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,5.,10.,20.,50.,100.,200.])))
    #cbar.set_ticklabels( [ "{:.2g}".format(x) for x in 10.**(cbar.get_ticks())])
    cbar.set_label(barLabel)



    # -labels
    ax.set_title('Object: {} $\\alpha$: {:.2f}$^o$, $\\delta$: {:.2f}$^o$'.format(objlabel,rao,deo))
    #ax.grid()

    ax.text(30.,-32,f'Constellation: {"".join(CONSTELLATIONS.name.split("\n"))} ({CONSTELLATIONS.totSat:.0f} sat.)', 
            fontsize=8, ha='left')

 
    ax.text(30,-40,f'Telescope: {myTel.telescope},    lat= {myTel.lat:.2f}$^o$',fontsize=8, ha='left')

    if myTel.fovl < 1./60.:
        fovll = f'{myTel.fovl_arcsec:.2f}\"'
    elif myTel.fovl < 1./6.:
        fovll = f'{myTel.fovl*60.:.2f}\''
    else:
        fovll = f'{myTel.fovl:.2f}$^o$'

    ax.text(30.,-48,f'Instrument: {myTel.instrument}    FoV= {fovll}   '+
            f'Resolution= {myTel.resol_arcsec:.1g}\"    Exp.T= {myTel.expt:.0f}s',
            fontsize=8, ha='left')

    fig.tight_layout()

    #- save plot
    if objlabel == "": objlabel = "OBJ"
    objlabel = objlabel.replace(" ","_")
    objlabel += f'_{rao:.1f}_{int(deo)}'

    outfileroot = f"{objlabel}_{myargs.code}_{myargs.constellations}_{myargs.output}"
    plt.savefig(outpath+'w.png')

    filename = outpath+outfileroot+myargs.outputformat
    plt.savefig(filename)
    log.info(f'Plot saved as {filename}')
    return filename


if __name__ == "__main__":

    cpLib.init_logger(log)
    log.info('===ObjectCalendar===')


    main()
