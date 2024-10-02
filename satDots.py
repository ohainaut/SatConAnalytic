#!/usr/bin/env python3
''' SatConAnalytic 

Plot the satellites plot over a skymap, with additional statistics.
This is *not* an analytical representation of the satellite density,
This is a traditional "discrete" representation of the satellites.

Notes
- the orbits are circular
- the motion of the satellites on their orbit is/can be slowed to make
  smoother animations (at natural speed, they are so fast that 
  the time step of the animation must be very small)

- this is a new implementation of vintage satDist
'''

import argparse, logging

import numpy as np
import matplotlib
#matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt

from astropy.table import Table, vstack

# import ConAn routines
import conan as ca
import conanplot as cp
import constellations
import constants as cst

#outpath = "/home/ohainaut/public_html/outsideWorld/"
outpath = "./"

logging.basicConfig(filename='conan.log',  
                    level=logging.INFO,
                    format='[%(levelname)-8s %(name)s/%(funcName)s] %(message)s')

log = logging.getLogger('conan')
log.info('satDots in')
#=====================================================================================
#=====================================================================================


#-----------------------------------------------------------------------------
def propagateTimeAnom(anom0, omega, times):
    '''Rotate the satellites on their orbit
    [rad], [rad/s], [s] -> [rad]'''

    return anom0 + omega * times  

#-----------------------------------------------------------------------------

def propagateTimeNode(node0, times):
    '''rotate the whole constellation (i.e. the Earth is rotatin under it)
    
    [rad], [s] -> [rad] '''

    return node0 + times /86400. *2.*np.pi  
#-----------------------------------------------------------------------------

def elementsToLongLat(noder,incr,anomr):
    '''get long, lat [rad] from propagated elements
     
    all in [rad] 
    '''

    latr   = np.arcsin(   np.sin( anomr ) *np.sin( incr )  )
    longr  = ( noder + np.arctan( np.tan( anomr ) *np.cos( incr )) 
              + np.pi* ( anomr  > np.pi/2. ) 
              + np.pi* ( anomr  > 3*np.pi/2. ) )
    longr  = longr%(2.*np.pi)

    return latr, longr
#-----------------------------------------------------------------------------

def LongLatToGeoXYZ(alt, latr, longr):
    '''geocentric rectangular: rs, xs, ys, zs: geocentric rect coord of the sat

    [same unit as cst.earthRadius, km]'''

    rs = cst.earthRadius + alt
    x = rs* np.cos(latr)* np.cos(longr) 
    y = rs* np.cos(latr)* np.sin(longr)  # note -
    z = rs* np.sin(latr)

    return x,y,z
#-----------------------------------------------------------------------------

def illuminatedSat(xg, yg, zg, sunAlpha, sunDelta):
    '''Which satellites are illuminated ?
    '''
    #== Sun coordinates
    cosalphaSun = np.cos(np.radians(sunAlpha))
    sinalphaSun = np.sin(np.radians(sunAlpha))
    cosdeltaSun = np.cos(np.radians(sunDelta))
    sindeltaSun = np.sin(np.radians(sunDelta))
    #xSun = cosdeltaSun*cosalphaSun
    #ySun = cosdeltaSun*sinalphaSun
    #zSun = sindeltaSun
        

    #- intermediate rotation of the satellites around z to bring X toward sun
    xsw =  xg *cosalphaSun + yg *sinalphaSun
    ysw = -xg *sinalphaSun + yg *cosalphaSun
    zsw =  zg

    #- rotation of the satellites around Y to bring X towards Sun
    xS =  xsw *cosdeltaSun + zsw *sindeltaSun
    yS =  ysw
    zS = -xsw *sindeltaSun + zsw *cosdeltaSun 

    #- distance of the satellite to the radius towards sun (X):
    distS = np.sqrt(yS**2 + zS**2)

    #- satellites in sunlight:
    w =  ((xS>0)*1  #those in front of the Earth
            +(distS > cst.earthRadius)*1) #those behing but out of the shadow
    return w>0 # True for the illuminated satellites
#-----------------------------------------------------------------------------

def geoXYZToTopXYZ(latr, x,y,z ):
    '''Geoctentric xyz to topocentric Az,El'''

    #- rotate to observatory latitude along axis y,
    #       and shift along z to the observatory

    sinLatO = np.sin(latr)
    cosLatO = np.cos(latr)

    xo = x * sinLatO     - z * cosLatO
    yo =              y
    zo = x * cosLatO     + z * sinLatO - cst.earthRadius  # <- shift z to observatory
    return xo, yo, zo
#-----------------------------------------------------------------------------

def topoXYZToDelta(x,y,z):
    '''Distance observer - satellite'''
    return np.sqrt(x**2 + y**2 + z**2)
#-----------------------------------------------------------------------------

def topoXYZToAzrZD(x,y,z, Delta=None):
    '''Az, ZD [rad]'''

    if Delta is None:
        Delta = topoXYZToDelta(x,y,z)

    ZD = np.degrees(np.arccos(z/Delta))             # Zenithal dist in rad
    Azr = np.pi   + np.arctan2(y,x)       # Azimuth # corr apr.18
    return Azr,ZD
#-----------------------------------------------------------------------------

def satMag(x,y,z, alt, mag550=cst.mag550, Delta=None, Sun=None):
    if Delta is None:
        Delta = topoXYZToDelta(x,y,z)

    #== Photometry
    Airm = Delta/alt     # airmass (not 1/cos, for extreme values)
    
    if Sun is None:
        cosalphaSun = 0.
    else:
        cosalphaSun = (x*Sun[0] + y*Sun[1] + z*Sun[2])/Delta

    return  mag550  -13.701 +  5.*np.log10(Delta) + 0.125 * Airm + 2.5* np.log10((1.+cosalphaSun)/2.)
    # -13.701 = -5log(550)
#-----------------------------------------------------------------------------

def magToDotSize(  mag  ):
    '''scale magnitude into a matplotlib dot size'''

    return  0.25*(11.-mag )**3 # size of the dot

#-----------------------------------------------------------------------------
def plot_elevMag(ax, Sv):
    '''plot mag vs zd and histogram of mag vs zd
    
    ax: the axes on which to plot
    Sv: a satellite Table
    '''

    Si = Sv[  Sv["bIlluminated"] ]
    Sb = Si[ Si["mag"] < 7 ]
    ax.set_xlim(0.,90)

    axH = ax.twinx() # second axis
    

    ax.scatter(Si["ZD"], Si["mag"], s=Si["dot"], 
                color="orange", alpha=0.3,
                zorder=20)
    ax.scatter(Sb["ZD"], Sb["mag"], s=Sb["dot"], 
                color="r", alpha=0.8,
                zorder=20)
    ax.set_ylim(12,4)

    sHist = axH.hist( Sv["ZD"],
                    range=(0,90),
                    bins=18,
                    color="k" ,
                    alpha=0.25,
                    zorder=30)
    sHist = axH.hist( Si["ZD"],
                    range=(0,90),
                    bins=18,
                    color="b" ,
                    alpha=0.25,
                    zorder=30)
    #print(sHist)

    #ax.set_ylim(11,4)
    ax.set_xticks( np.arange(0.,91,30.))
    ax.set_xlabel("Zenithal Distance [deg]")
    ax.set_ylabel("Mag")
    axH.set_label("Count")
    #ax.set_yticks( np.arange( 4., 11))

#-----------------------------------------------------------------------------
def plot_legendMag(ax, minMag=5,maxMag=11.1):
    minMag = min(5, minMag)
    maxMag = max(11.1,maxMag)


    def plotit(myMag, color):
        myX = myMag
        myY = np.zeros_like( myMag )
        myDot = magToDotSize( myMag )
        ax.scatter(myY, myX, s=myDot, color=color)
        for i in range(len(myMag)):
            ax.text(  myY[i], myX[i],   f'    {int(myMag[i])}',
                    #size=7.,
                    horizontalalignment='left',
                    verticalalignment= 'center')

    plotit( np.arange(int(minMag),6.1,1.) , 'red')
    plotit( np.arange(7,int(maxMag)) , 'orange')

    #ax.set_title('Mag', size=10)
    ax.set_xlim( -.6,1.5)
    ax.set_ylim(minMag-.5,maxMag+.5)
  

    ax.tick_params( axis='both', which='both', 
                   bottom=False, left=False,
                   labelbottom=False, labelleft=False)

#-----------------------------------------------------------------------------



def makeConstellationStatTable(CONSTELLATIONS, sunAlpha, sunDelta, lat, times):


    # get the orbital element table for all the satellites
    CONSTELLATIONS.Table = vstack(
        [ myShell.orbitalElementsTable() for myShell in CONSTELLATIONS.shells]
    )
    log.info(f'check: {len(CONSTELLATIONS.Table)} lines, totSat = {CONSTELLATIONS.totSat}')
    SAT = CONSTELLATIONS.Table


    # revolution of the sat
    SAT["anomr"] = propagateTimeAnom(SAT["anom0"], SAT["omega"], times/90.)
    #SAT["anomr"] = SAT["anom0"].copy()

    # rotation of the Earth
    SAT["noder"] = propagateTimeNode(SAT["node0"], times/90.)
    #SAT["noder"] = SAT["node0"].copy()

    SAT["latr"],SAT["longr"] = elementsToLongLat(SAT["noder"], SAT["inc"], SAT["anomr"])
    SAT["xg"],SAT["yg"],SAT["zg"]   = LongLatToGeoXYZ(SAT["alt"], SAT['latr'], SAT["longr"])
    SAT["bIlluminated"] = illuminatedSat(SAT["xg"],SAT["yg"],SAT["zg"], sunAlpha, sunDelta)


    # intermediate plot - map of the illuminated satellites etc
    if 0:
        Ill = illuminatedSat(SAT["xg"],SAT["yg"],SAT["zg"], sunAlpha, sunDelta)        
        plt.scatter(np.degrees(SAT["longr"]), np.degrees(SAT["latr"]), c=Ill, s=1, alpha=0.5)
        plt.scatter(sunAlpha,sunDelta,c='r',s=100)
        plt.scatter(0.,lat,c='k',s=100)
        plt.show()

    # intermediate plot - 3d map of the illum satellites
    if 0:
        Ill = illuminatedSat(SAT["xg"],SAT["yg"],SAT["zg"], sunAlpha, sunDelta)        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(SAT["xg"],SAT["yg"],SAT["zg"], s=1, c=Ill )
        plt.show()


    SAT["xt"],SAT["yt"],SAT["zt"]  = geoXYZToTopXYZ(np.radians(lat),
                                                    SAT["xg"],SAT["yg"],SAT["zg"])
    SAT["bVis"] = SAT["zt"] < 0 # observability

    # intermediate plot - 3d map of the observable satellites
    if 0:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(SAT["xt"],SAT["yt"],SAT["zt"], s=1, c=1-SAT["bVis"] )
        plt.show()



    # Select only observable sat
    Sv = SAT[ SAT["zt"] > 0]
    log.info(f'Select Visible: {len(SAT)} -> {len(Sv)}')

    Sv["bIlluminated"] =illuminatedSat(Sv["xg"],Sv["yg"],Sv["zg"], sunAlpha, sunDelta) 
    Sv["Delta"] = topoXYZToDelta( Sv["xt"],Sv["yt"],Sv["zt"] )
    Sv["Azr"],Sv["ZD"] = topoXYZToAzrZD( Sv["xt"],Sv["yt"],Sv["zt"],
                                        Delta=Sv["Delta"]) 

    if 0:
        plt.plot(Sv["Azr"],Sv["ZD"] )
        plt.show()

    Sv["mag"] = satMag(Sv["xt"],Sv["yt"],Sv["zt"] , Sv["alt"],  Delta=Sv["Delta"])
    Sv["dot"] = magToDotSize( satMag(Sv["xt"],Sv["yt"],Sv["zt"] , Sv["alt"], 
                                     Delta=Sv["Delta"]) )

    if 0:
        plt.scatter(Sv["Delta"],  Sv["mag"],       s=1 )
        plt.show()

    return Sv


#=====================================================================================
#=====================================================================================



if __name__ == "__main__":
 
    #--- command line arguments
    parser = argparse.ArgumentParser(
        description='Compute and plot the position and magnitude of each satellite')
    parser.add_argument('-a','--RA',         default=0.,
                        help='''local time [h]''')
    parser.add_argument('-d','--DEC',        default=0.,
                        help='''Declination of the Sun [deg]''')
    parser.add_argument('-n','--objlabel',  default="",
                        help='''Name of the object for label''')
    parser.add_argument('-C','--constellation', default='SLOWGWAK',
                        help="Constellation code (or meta-code)")
    parser.add_argument('-T','--code',  default="FORSimg",
                        help='''Observatory: Telescope/Instrument code''')
    parser.add_argument('-l','--lat',   help='''Observatory: Latitude of the observatory [deg] (OVERWRITE preset)''')
    parser.add_argument('-t','--expt',  help='''Observatory: Individual exposure time [s] (OVERWRITE preset)''')
    parser.add_argument('-r','--resol', help='''Observatory: Resolution element (seeing, pixel) [deg] (OVERWRITE preset)''')
    parser.add_argument('-f','--fovl',  help='''Observatory: Length of the field-of-view [deg] (OVERWRITE preset)''')
    parser.add_argument('-w','--fovw',  help='''Width of the field-of-view [deg] (Default=Fovl; OVERWRITE preset)''')
    parser.add_argument('-k','--trailf',help='''Observatory: Fraction of the exposure destroyed by a trail (1=full) (OVERWRITE preset)''')
    parser.add_argument('-m','--maglim',     help='''Observatory: Limiting magnitude [mag] (detection limit for expTime) (OVERWRITE preset)''')
    parser.add_argument('-M','--magbloom',   help='''Observatory: Saturation magnitude [mag]. Brighter object destroy the full exposure: their trailf=1 (OVERWRITE preset)''')
    parser.add_argument(     '--instrument', help='''Observatory: Name of the instrument for label (OVERWRITE preset)''')
    parser.add_argument(     '--telescope',  help='''Observatory: Name of the telescope for label (OVERWRITE preset)''')
    parser.add_argument(     '--mode',       default="OBS",
                        help='''BRIGHT FAINT ALL OBS EFFECT''')
    myargs = parser.parse_args()


    #- find telescope
    print('TELESCOPE/INSTRUMENT SETUP')
    myTel = cp.getTelescope(myargs)
    print(myTel)

    #- rectangular coord of observatory
    xobs = cst.earthRadius * np.cos(np.radians(myTel.lat))
    yobs = 0.  # by def of xyz
    zobs = cst.earthRadius * np.sin(np.radians(myTel.lat))

    #---

    # satellites
    CONSTELLATIONS = ca.findConstellations(myargs.constellation)
    print('CONSTELLATIONS:')
    print(CONSTELLATIONS.ToC)
    print()


    #- time
    timeh = float(myargs.RA)
    timed = timeh*15. #[deg]
    times = timeh*3600. #[s]

    # sun coordinates
    sunDelta =  float(myargs.DEC) # deg
    sunAlpha = timed -180. # sunAlpha is the HA [deg]
    sunElev = ca.radec2elev(sunAlpha,sunDelta,myTel.lat)  # [deg]
    log.info(f'Sun: {sunAlpha/15.}h, {sunDelta} time={timeh}h elev={sunElev}')


    print('SUN:')
    print(f'Local time: {timeh}h')
    print(f'HA = {sunAlpha}deg  = {sunAlpha/15.:.2f}h, Dec = {sunDelta}')
    print(f'Computed elevation: {sunElev:.2f}')

    # validation
    wxs, wyx, wzs = LongLatToGeoXYZ(150e6, np.radians(sunDelta), np.radians(sunAlpha))
    wxs, wyx, wzs = geoXYZToTopXYZ(np.radians(myTel.lat), wxs, wyx, wzs )
    wAzr, wZD = topoXYZToAzrZD(wxs, wyx, wzs )
    print(f'validation: Az ={np.degrees(wAzr):.2f}, ZD = {wZD:.2f}, el = {90-wZD:.2f} ')

    #==SATELLITES

    Sv = makeConstellationStatTable(CONSTELLATIONS, sunAlpha, sunDelta, 
                                    myTel.lat, times)
    # SV is the table with only the observable satellites

    Sd = Sv[  ~Sv["bIlluminated"] ] # Dark
    Si = Sv[  Sv["bIlluminated"] ] # illuminated
    Sf = Si[ Si["mag"] >= 7 ] # faint
    Sb = Si[ Si["mag"] < 7 ]  # bright



    fig = plt.figure(figsize=(12,8))
    axTitle  = fig.add_subplot(4,2,1)  # just the title and labels
    axW      = fig.add_subplot(4,2,3)  # E-W side view
    axNS     = fig.add_subplot(2,6,10)  # N-S side view
    axLegend = fig.add_subplot(2,6,12)
    axPol    = fig.add_subplot(2,2,3, projection='polar')  # top polar view
    axHist   = fig.add_subplot(2,2,2)  # histogram


    hlimkm = 4000. # km;   limits for the side view
    #- prepare the limb of the Earth for plots.
    wi = np.radians(np.linspace(0,360,360, endpoint=False))
    xearth = cst.earthRadius* np.cos(wi)
    yearth = cst.earthRadius* (np.sin(wi)-1.)

    # NS
    axNS.set( aspect='equal')
    for S, col in zip( [Sd, Sf, Sb], ['k', 'orange', 'red']):
        if len(S) > 0:
            axNS.scatter(S["zt"],S["xt"], 
                     s=S["dot"], c=col, alpha=0.3)
    axNS.set_ylim(hlimkm,-hlimkm) ## corr Apr14## North is negagive /up 
    axNS.set_xlim(-300,1800)
    axNS.plot(yearth,xearth,c='b')
    axNS.scatter(0.,0.,s=30,c='b')

    # EW
    axW.set( aspect='equal')
    for S, col in zip( [Sd, Sf, Sb], ['k', 'orange', 'red']):
        if len(S) > 0:
            axW.scatter(S["yt"],S["zt"], 
                     s=S["dot"], c=col, alpha=0.3)
    axW.set_xlim(-hlimkm,hlimkm) ## 
    axW.set_ylim(-300,1800)
    axW.plot(xearth,yearth,c='b')
    axW.scatter(0.,0.,s=30,c='b')


    if 1:
        #ax =  fig.subplots(1,1,subplot_kw={'projection': 'polar'}) 

        if sunElev > 0:
            pass
        elif sunElev > -18:
            icol = (18+sunElev)/18.
            axPol.set_facecolor( (icol,icol,icol))
        else: 
            axPol.set_facecolor( "k")

        cp.initPolPlot(axPol)

        if myargs.mode == "ALL":
            Sd = Sv[  ~Sv["bIlluminated"] ]
            axPol.scatter(Sd["Azr"],Sd["ZD"], s=Sd["dot"], c="darkblue", alpha=0.5)

        Si = Sv[  Sv["bIlluminated"] ]
        if myargs.mode in ["ALL",  "OBS"] :
            Sb = Si[ Si["mag"] >= 7 ]
            axPol.scatter(Si["Azr"],Si["ZD"], s=Si["dot"], c="orange")

        if myargs.mode in ["ALL", "OBS","BRIGHT"] :
            Sb = Si[ Si["mag"] < 7 ]
            log.info(f"bright: {len(Sb)}")
            axPol.scatter(Sb["Azr"],Sb["ZD"], s=Sb["dot"], c="red")


    # TITLE
    axTitle.set_xlim(0.,1.)
    axTitle.set_ylim(0.,1.)
    axTitle.set_axis_off()
    
    #top left
    x = 0.
    y = 0.9
    dy = 0.2
    axTitle.text(x,y,f'Observatory: {myTel.telescope} Lat.: {myTel.lat:.1f}$^o$')
    y -= dy
    
    loct = (sunAlpha/15.+12.)%24
    loch = int(loct)
    locm = int( (loct-loch)*60.)
    axTitle.text(x,y,f'$\odot$ Sun: Loc.time: {loch:02d}:{locm:02d} '+
             f'$\delta: {sunDelta:.2f}^o$, Elev: {sunElev:.2f}$^o$')
    y -= dy

    axTitle.text(x,y,'Constellation:')
    y -= dy
    axTitle.text(x,y,CONSTELLATIONS.name)
    y -= dy
    axTitle.text(x,y, f'Total {CONSTELLATIONS.totSat:.0f} sat.')
 



    plot_legendMag(axLegend)

    plot_elevMag(axHist, Sv)

    #fig.tight_layout()
    plt.savefig(f'w.png')
    #plt.savefig(f'w{int(times)}.png')
#--
    log.info('============================================================================')