#! python3

import argparse, logging

from astropy.time import Time
from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')

import astropy.units as u
import numpy as np

# import ConAn routines
import SatConAnalytic.conan as caLib
import SatConAnalytic.conanplot as cpLib
import SatConAnalytic.satDots as satDots
import SatConAnalytic.skyBrightnessLib as skyLib
import SatConAnalytic.constellations as constLib
import SatConAnalytic.telescopes as telescopeLib
import SatConAnalytic.utils as utLib
from   SatConAnalytic.conanplot import gyrd
from SatConAnalytic.utils import _Dict

log = logging.getLogger('satcon')
log.setLevel(logging.DEBUG)  
for logger_name in ("matplotlib", "matplotlib.pyplot", "matplotlib.font_manager", "matplotlib.backends"):
    logging.getLogger(logger_name).setLevel(logging.WARNING)


#----------------------------------------------------------------------------
def create_argument_parser():
    """Create and return the argument parser"""
    parser = argparse.ArgumentParser(description=
    '''SatCon ''')
    # Constellation
    parser.add_argument('-C','--constellations', default='SLOWGWAK',
                        help="Constellation code (or meta-code); list for a list.")
    parser.add_argument('--constFile', default="constellations.json",
                        help="Constellation: Which constellation file to use (default: constellations.json)")
    parser.add_argument('--mag550', 
                        help="Constellation: overwrite absolute mag for all satellites [Mag at 550km]")


    
    # Observatory and instrument parameters
    parser.add_argument('-D','--date', default="2024-01-01",
                        help="Observatory: Date of the begin of night [YYYY-MM-DD]")
    

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

    parser.add_argument('-l','--lat', 
                        help="Observatory: Latitude of the observatory [deg]")
    parser.add_argument('-r','--resol', 
                        help="Observatory: Resolution of the instrument [deg]")
    parser.add_argument('-t','--expt', 
                        help="Observatory: Exposure time [sec]")
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


    # selection parameters
    parser.add_argument('-M','--magSelect', default="all", 
                        choices=['all', 'detected', 'oversaturated', 'notDetected'],
                        help="selection on magnitude; default=all")

    parser.add_argument('--elev', action='store_true',
                        help="Overplot elevation curves")



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


    # observation catalogue
    parser.add_argument('--cat', default="caldwell.json",
                        help="Path to the observation catalogue JSON file")

    # Plot parameters
    parser.add_argument('--minmax', nargs=2, 
                        help="min, max value for the colorscale")

    parser.add_argument('--noplot', action='store_false',
                        help="Plot: Don't generate the plot (mostly debug)")
    
    parser.add_argument('--noconan', action='store_true',
                        help="Plot: Don't plot the conAn simulation; only the sat dots")

    parser.add_argument('--noshade', action='store_true',
                        help="Plot: Don't shade low elevations")
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


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
class Night():
    '''define the astronomical night,   with attributes:
    - date: date of the night (string)
    - dateT: date of the night (astropy Time)
    - lst_0h_h: local sidereal time at 0h UT (hours)
    - sunAlpha_deg: right ascension of the Sun (degrees)
    - sunDelta_deg: declination of the Sun (degrees)
    - midnightT: time of the astronomical midnight (astropy Time)
    - twilights_eveningT: times of the evening twilights (astropy Time array)
    - twilights_morningT: times of the morning twilights (astropy Time array)
    - twilights_label: labels for the twilights (list of strings)
    '''

    def __init__(self, date, long_deg, lat_deg):

        self.date = date
        self.dateT = Time(date)

        self.long_deg = long_deg
        self.lat_deg = lat_deg

        self.lst_0h_h = caLib.siderealTimeDeg(self.dateT.jd, long_deg)/15.
        self.sunAlpha_deg, self.sunDelta_deg = caLib.get_sun(self.dateT.jd)
        self._set_midnight()
        self._set_twilights()

        log.debug(f'Created Night: {self.date}')

    def _set_midnight(self):
        antisunAlpha_h = ((self.sunAlpha_deg + 180.)/15.)%24.
        self.midnightT = self.dateT + (antisunAlpha_h - self.lst_0h_h)/23.9344444* u.day


    def _set_twilights(self):
        #--- twilight times
        antisunDelta_r = np.radians(-self.sunDelta_deg)

        zdrs = np.radians(90.-np.arange(0.,19.,6.)) 
        # antisun zenith distance at twi
        twilightHAs = caLib.findHAfromZDr(zdrs, 
                                          antisunDelta_r, 
                                          np.radians(self.lat_deg))
        #--- convert to UT
        self.twilights_eveningT  = self.midnightT - twilightHAs /23.9344444* u.day
        self.twilights_morningT  = self.midnightT + twilightHAs /23.9344444* u.day

        self.twilights_label = ['set', 'civil', 'nautical', 'astronomical' ]


    def __repr__(self):
        msg = f'Date: {self.date} \n'
        msg += f'Reference time: {self.dateT.iso} UT\n'
        msg += f'LST at reference time: {self.lst_0h_h:.2f} h\n'
        msg += f'Sun position on {self.date}: alpha={(self.sunAlpha_deg/15)%24:.2f} h, delta={self.sunDelta_deg:.2f} deg \n'
        msg += f'AntiSun=ST at mid-night: {(self.sunAlpha_deg/15 + 12)%24:.2f} h = {(self.sunAlpha_deg+ 180)%360:.2f} deg \n'
        msg += f'Midnight: {self.midnightT.iso} \n'
        msg += f'Twilights: \n'
        for l,et, mt in zip( self.twilights_label, self.twilights_eveningT, self.twilights_morningT):
            msg += f'  {l}: {et.iso}, {mt.iso}\n'
        return msg
    
    
    #-----------------------------------------------------------------------
    def JD2UT(self, jd):
        '''Convert Julian Date to UT hours, relative to the given midnight time'''
        return (jd - int(self.midnightT.jd))*24 +12


#----------------------------------------------------------------------------
def readStars(myFile="caldwell.json"):
    """Read the star catalog from a JSON file and return it as a Table
    
    Args:
        myFile: Path to the JSON file containing the star catalog (default: "caldwell.json")
    
    Returns:
        dict: Dictionary containing the star catalog
    """

    jstars = utLib.read_conanJson(myFile)

    # convert dict of arrays to Table
    starTable = Table(jstars)

    log.debug(f"Read {len(starTable)} stars from {myFile}")

    return starTable
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
def plot_twilights(NIGHT, ax=None):
    if ax is None:  
        fig, ax = plt.subplots(figsize=(12, 8))

    ymin, ymax = ax.get_ylim()
    ymid = (ymin + ymax)/2
    xmin, xmax = ax.get_xlim()

    for i in range(len(NIGHT.twilights_label)):
        et     = NIGHT.JD2UT(NIGHT.twilights_eveningT[i].jd)
        mt     = NIGHT.JD2UT(NIGHT.twilights_morningT[i].jd)
        label = NIGHT.twilights_label[i]

        ax.fill_betweenx([ymin, ymax], xmin, et, color='blue', alpha=0.3) 
        ax.fill_betweenx([ymin, ymax], mt, xmax, color='blue', alpha=0.3)

        # draw vertical lines at the twilight times
        ax.axvline(x=et, color='blue', linewidth=1)
        ax.axvline(x=mt, color='blue', linewidth=1)

        ax.text(et, ymid, f'{NIGHT.twilights_label[i]}', rotation=90, verticalalignment='center', ha='right', color='white')
        ax.text(mt, ymid, f'{NIGHT.twilights_label[i]}', rotation=90, verticalalignment='center', ha='left', color='white')

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

#----------------------------------------------------------------------------
def plot_separators(nightUTs, stars, ax=None):
    for i in range(len(stars)):
        ax.axhline(y=(i+1)/len(stars), color='white', linewidth=0.5, alpha=0.5)
        ax.axhline(y=(i+1)/len(stars), color='k', linestyle=':', linewidth=0.5, alpha=0.5)

#----------------------------------------------------------------------------
def plot_densityBars(elevations, logplotvalues, stars, nightUTs, label=None, ax=None):
    '''
    :param elevations: the elevation values of the stars; shape (nStars, nTimes)
    :param logplotvalues: the values to plot as bars; shape (nStars, nTimes)
    :param stars: Description of the stars; should contain at least 'dec_deg' for labelling
    :param nightUTs: Description of the night times in UT; should be a 1D array of time values corresponding to the columns of logplotvalues
    ---
    :param label: Label for the colorbar
    :param ax: Matplotlib axis object to plot on. If None, a new figure and axis will be created.
    '''


    cmap = gyrd
    cmap = matplotlib.colormaps['magma']

    nStars = len(stars)
    bar_height = 1. / nStars

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    plotmin = np.min(logplotvalues[ logplotvalues > -99. ])  
    plotmax = np.percentile(logplotvalues, 99.9)
    log.info(f"plotValues: min={plotmin:.2e}, 99.9th percentile={plotmax:.2e}")

    for i in range(nStars):
        logvalues = logplotvalues[i, :]
        eleva1 = elevations[i, :]   # elevation in degrees, shape (nTimes,)

        # Build a colored bar segment-by-segment along the time axis
        for j in range(len(nightUTs) - 1):
            color = cmap( 
               np.clip((logvalues[j]-plotmin) / (plotmax-plotmin), 0, 1) ) 

            if eleva1[j] < 0.:   
                color = 'grey'
                #make the color more transparent for low elevations
                #color = (color[0], color[1], color[2], 0.1)  

            ax.barh(
                y      = i * bar_height,
                width  = nightUTs[j+1] - nightUTs[j],
                left   = nightUTs[j],
                height = bar_height,
                color  = color,
                align  = 'edge',
                linewidth = 0,
            )
        


    # Y-axis: one tick per star, labelled with name or NGC id
    ax.set_yticks( (np.arange(nStars) + 0.5) * bar_height )
    #labels = [f'{row["ra_deg"]:.1f}{row["dec_deg"]:+.1f}' for row in stars]
    labels = [row["ngc"] for row in stars]
    ax.set_yticklabels(labels, fontsize=7)

    for i in range(nStars):
        ax.text(nightUTs[0], (i+0.5)*bar_height, 
                f' {stars[i]["ra_deg"]:.1f}{stars[i]["dec_deg"]:+.1f}', 
                verticalalignment='center', horizontalalignment='left', fontsize=7, color='white')


    # Colourbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=plotmin, vmax=plotmax))
    sm.set_array([])

    if label is not None:
        cbar = plt.colorbar(sm, ax=ax, label=label, fraction=0.02, pad=0.01)
    else:
        cbar = plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.01)

    barticks = cbar.get_ticks()
    barticks = barticks[(barticks >= plotmin) & (barticks <= plotmax)]
    cbar.set_ticks(barticks)
    cbar.set_ticklabels([f'{10**x:.1g}' for x in barticks])


#----------------------------------------------------------------------------
def plot_elevationBars(elevations, stars, nightUTs,  ax=None):
    '''
    :param elevations: the values to plot as bars; shape (nStars, nTimes)
    :param stars: Description of the stars; should contain at least 'dec_deg' for labelling
    :param nightUTs: Description of the night times in UT; should be a 1D array of time values corresponding to the columns of elevations
    '''

    cmap = plt.get_cmap('RdYlGn')

    nStars = len(stars)
    bar_height = 1. / nStars

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    for i in range(nStars):
        values = elevations[i, :]   # elevation in degrees, shape (nTimes,)

        ax.plot(nightUTs, (values/180 + i+0.5)*bar_height, color='r', linewidth=0.5   )
        mask = values > 20.
        ax.plot(nightUTs[mask], (values[mask]/180 + i+0.5)*bar_height, color='orange', linewidth=1   )
        mask = values > 30.
        ax.plot(nightUTs[mask], (values[mask]/180 + i+0.5)*bar_height, color='lawngreen', linewidth=1   )



#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
def main(args=None):
    """Main function that can be called with arguments or from command line
    
    Args:
        args: List of arguments (if None, will parse from command line)
    
    Returns:
        dict: Dictionary containing results and output information
    """
    
    print('>>>SatConAnalytic: Planning<<<')
    
    #----- config
    step = 1# 0.1 # hour
    #outpath = "/home/ohainaut/public_html/outsideWorld/"
    outpath = "./"
    
    # Parse arguments
    parser = create_argument_parser()
    if args is None:
        myargs = parser.parse_args()
    else:
        myargs = parser.parse_args(args)


    #===========================================================================
    # TELESCOPE/INSTRUMENT
    #===========================================================================
    log.info('=====TELESCOPE/INSTRUMENT SETUP================================')
    myTel = telescopeLib.getTelescope(myargs)
    log.info(myTel)


    #===========================================================================
    # CONSTELLATIONS
    #===========================================================================
    log.info('=====CONSTELLATIONS==================================')
    CONSTELLATIONS = constLib.findConstellations(myargs.constellations, constFile=myargs.constFile)  # constellationS object

    log.debug(f'Found constellations:\n{CONSTELLATIONS.ToC}')
    log.info(f'Constellations: {CONSTELLATIONS.totSat} satellites in {CONSTELLATIONS.totShells} shells.')


    #===========================================================================
    # NIGHT
    #===========================================================================
    log.info('=====NIGHT=========================================')
    NIGHT = Night(myargs.date, myTel.lon, myTel.lat)

    print(NIGHT.twilights_eveningT[0].jd, NIGHT.twilights_morningT[0].jd)

    nightJDs = np.arange(NIGHT.twilights_eveningT[0].jd -0.05, 
                         NIGHT.twilights_morningT[0].jd +0.05, step/24.)
    nightUTs = NIGHT.JD2UT(nightJDs)
    nightSTs = caLib.siderealTimeDeg(nightJDs, NIGHT.long_deg)

    log.info(NIGHT)

    # sun coordinates
    walphas, nightSunDecs = caLib.get_sun(nightJDs)
    nightSunHAs = nightSTs - walphas # alphas is the HA
    nightSunEls = caLib.radec2elev(nightSunHAs,nightSunDecs,myTel.lat)



    #===========================================================================
    # STARS
    #===========================================================================
    log.info('=====STARS=========================================')
    stars = readStars(myargs.cat)

    # filter stars table to keep only every 10th star
    #stars = stars[::10]

    stars.sort('ra_deg')  # sort by RA for better plotting
    #stars.sort('dec_deg')  # sort by dec for debug

    # only stars that can be above horizon
    log.debug(f"latitude: {NIGHT.lat_deg} deg")
    log.debug(f"limits: {NIGHT.lat_deg - 90.} .. {NIGHT.lat_deg + 90.}")
    log.info(f"Number of stars              : {len(stars)}")
    

    if NIGHT.lat_deg > 0.:
        stars = stars[ stars['dec_deg'] > NIGHT.lat_deg - 90.  ] 
        log.info(f"Number of stars dec > {NIGHT.lat_deg - 90.}  : {len(stars)}")
    else:
        stars = stars[ stars['dec_deg'] < NIGHT.lat_deg + 90.  ] 
        log.info(f"Number of stars dec < {NIGHT.lat_deg + 90.} : {len(stars)}")

    if len(stars) == 0:
        msg = "No stars observable. "
        log.error(msg)
        return msg

    #===========================================================================
    # NIGHT SCAN
    #===========================================================================
    log.info('=====NIGHT SCAN====================================')

    # create a 2D array of HA for each star and each time step
    starHAs = np.zeros((len(stars), len(nightJDs)))
    starDecs = np.zeros((len(stars), len(nightJDs)))
    for i, star in enumerate(stars):
        starHAs[i,:] = (nightSTs - star['ra_deg'] + 360.)%360.
        starDecs[i,:] = star['dec_deg']
    
    AzEl = np.array(caLib.radec2azel(starHAs, starDecs, NIGHT.lat_deg))
    # AzEl shape is (2, nStars, nTimes); [0,:,:] is Az, [1,:,:] is El




    # filter stars that are never above 30 deg elevation
    stars['maxel'] = np.max(AzEl[1,:,:], axis=1)
    filter = stars["maxel"] > 30.
    stars = stars[ filter ]
    AzEl = AzEl[:, filter, :]
    log.info(f"Number of stars above 30 deg elevation: {len(stars)}")
    if len(stars) == 0:
        msg = "No stars observable tonight. "
        log.error(msg)
        return msg


    #===========================================================================
    # CALCULATE CONAN SIMULATION
    #===========================================================================
    log.info('=====CONAN=========================================')

    # initialize output arrays as a SIM_DENSITIES object
    SIM = cpLib.init_output(AzEl)

    # duplicate nightSunHAs, nightSunDecs to match the shape of AzEl
    alphas = np.tile(nightSunHAs, SIM.satAll.shape[0] )
    deltas = np.tile(nightSunDecs,SIM.satAll.shape[0] )

    for myShell in CONSTELLATIONS.shells:                 
        # process each shell
        densSi, veli, magi =  myShell.modelOneShell(
            AzEl,myTel.lat, alphas, deltas  )

        # all sat:
        SIM.satAll += densSi
        SIM.velAll += densSi * veli


        # effective magnitude
        mageffi = magi  - 2.5*np.log10(myTel.resol/veli/myTel.expt)

        # only observable ones
        densSobsi = np.copy(densSi)
        densSobsi[ mageffi > myTel.maglim] = 0.
        SIM.satObs += densSobsi
        SIM.velObs += densSobsi * veli
        
        # only super bright bloomers
        densSbloomi = np.copy(densSi)
        densSbloomi[ mageffi > myTel.magbloom ] = 0.
        SIM.satBloom += densSbloomi
        SIM.velBloom += densSbloomi * veli
        
        #luminance from the satellites in the Element:
        #   f = lumZP*10^(-0.4*mag)   
        #        # flux from one satellite converted in Cd/m^2
        #        # Correct for 1sat/sq.arcsec
        #     * densSi / 3600^2    # n sat per sq.arc 
        lumini = skyLib.luminanceZP *10.**(-0.4*magi) * densSi /3600.**2 
        luminNotDeti = lumini.copy()
        luminNotDeti[ mageffi <= myTel.maglim] = 0.  # only the faint ones
        SIM.luminance       += lumini
        SIM.luminanceNotDet += luminNotDeti


    # == magnitude selection ==
    # Select effective densities depending on the requested output:
    #   ds: density of satellites;
    #   dv: density of satellite trails (i.e. density of satellites * velocity)

    ds, dv, selectionLab = cpLib.magSelect(myargs.magSelect, SIM, myTel)


    # == output selection ==
    # depending on the Output type,
    # set the corresponding logDensity to plot, 
    # the colour bar label, and the title of the plot.
    
    OUT = cpLib.selectOutput(myargs, myTel, ds, dv, SIM)
    logDensity    = OUT.logDensity
    barTicks      = OUT.barTicks
    barTickLabels = OUT.barTickLabels
    barLabel      = OUT.barLabel
    



    #===========================================================================
    # PLOT
    #===========================================================================
    log.info('=====PLOT=========================================')
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlabel('UT (h)')
    ax.set_xlim(nightUTs[0], nightUTs[-1])
    ax.set_ylim(0, 1)

    plot_densityBars( AzEl[1,:,:],  logDensity, stars, nightUTs, 
                     label=barLabel, ax=ax)

    plot_twilights(NIGHT, ax=ax)
    plot_separators(nightUTs, stars, ax=ax) 
    if myargs.elev:
        plot_elevationBars(AzEl[1,:,:], stars, nightUTs, ax=ax)
    # ax.plot(nightUTs, nightSunEls/180.+0.5, color='orange', label='Sun elevation') # sun elevation


    # UT Scale on bottom
    xmin, xmax = ax.get_xlim()
    xticks = np.arange(np.ceil(nightUTs[0]), np.floor(nightUTs[-1])+1, 1.)
    ax.set_xticks(xticks)
    ax.set_xticklabels([f'{int(x%24):02d}' for x in xticks])


    # ST scale on top
    # find 360->0 wrap in nightSTs and add 24 to the values after the wrap 
    nightSTs_h = (nightSTs/15.)
    wrap_indices = np.where(np.diff(nightSTs_h) < -12)[0]
    if len(wrap_indices) > 0:
        nightSTs_h[wrap_indices[0]+1:] += 24.

    ax2 = ax.twiny()
    ax2.set_xlim(xmin, xmax)

    # pick round ST hour values that fall within the night
    _ST0 = np.ceil( nightSTs_h[0])
    _ST1 = np.ceil( nightSTs_h[-1])
    if _ST1 < _ST0:
        _ST1 += 24.
    _STrange = np.arange(_ST0, _ST1+1, 1.)
    ut_for_st = np.interp(_STrange,
                          nightSTs_h,
                          nightUTs )


    mask = (ut_for_st >= nightUTs[0]) & (ut_for_st <= nightUTs[-1])
    ut_for_st = ut_for_st[mask]
    _STrange = _STrange[mask]

    ax2.set_xticks(ut_for_st[:-1])
    ax2.set_xticklabels([f'{int(x%24):02d}' for x in _STrange[:-1]], fontsize=9)
    ax2.set_xlabel('LST (h)')

    ax.text(0.075,0.05,f'Telescope: {myTel.telescope},  long= {myTel.lon:.2f}$^o$ lat= {myTel.lat:.2f}$^o$',fontsize=10, ha='left', 
            transform=plt.gcf().transFigure)

    ax.text(0.075,0.03,f'Instrument: {myTel.instrument}   '+
            f'  Exp.T= {myTel.expt:.0f}s',
            fontsize=10, ha='left',
            transform=plt.gcf().transFigure)

    ax.text(0.925,0.03,f'Constellation: {"".join(CONSTELLATIONS.name.split("\n"))} ({CONSTELLATIONS.totSat:.0f} sat.)', 
            fontsize=10, ha='right',
            transform=plt.gcf().transFigure)


    ax.text(0.925,0.05,r'$\odot$: '+f'{NIGHT.sunAlpha_deg:.2f}$^o$ {NIGHT.sunDelta_deg:+.2f}$^o$', 
            fontsize=10, ha='right',
            transform=plt.gcf().transFigure)


    ax.text(0.475,0.03,NIGHT.date,     
            fontsize=12, ha='center',
            transform=plt.gcf().transFigure)

    ax.grid(True, which='major', axis='x', color='white', linestyle=':', linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    plt.savefig('w.png', dpi=300)
    if plt.get_backend().lower() != 'agg':
        plt.show(block=True)

#----------------------------------------------------------------------------
if __name__ == "__main__":

    utLib.init_logger(log)
    log.info('===Planning===')

    main()