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


option -h for detailed usage
'''

import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt
import numpy as np
import argparse
#sys.path.append(os.path.dirname(__file__)+"/../")

print('conAn ObsSky')

from matplotlib import ticker
from conanplot import gyrd


# import ConAn routines
import conan as ca
import conanplot as cp
import constants as cst
import satDots




#----- config
step = 0.75 #deg >~1. Smaller values take forever
#outpath = "/home/ohainaut/public_html/outsideWorld/"
outpath = "./"


#------Arguments

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
parser.add_argument('-l','--lat', default=-24.6,
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
parser.add_argument('-M','--magSelect', default="OBS", 
                    choices=['ALL', 'OBS', 'BRIGHT', 'FAINT', 'EFFECT'],
                    help="selection on magnitude; default=OBServable")

parser.add_argument('--minmax', nargs=2, 
                    help="min, max value for the colorscale")



parser.add_argument('--noplot', action='store_false',
                    help="Plot: Don't generate the plot (mostly debug)")
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
myargs = parser.parse_args()



# flags
myargs.plotflag      = myargs.noplot
myargs.shadeflag     = myargs.noshade
myargs.scalebarflag  = myargs.noscalebar
myargs.labelplotflag = myargs.nolabel
myargs.almucantar    = myargs.noalmuc
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

#
# OBSERVATORY TELESCOPE INSTRUMENT
#
print('=====TELESCOPE/INSTRUMENT SETUP=======================================')
myTel = cp.getTelescope(myargs)
print(myTel)


#
# SUN
#

sunDelta = float(myargs.deltaSun)
sunElev  = -float(myargs.elevSun)   
if myargs.alphaSun is None:
    sunAlpha = ca.elev2ra(sunElev,sunDelta,myTel.lat) # get sun hourangle for twilight
else:
    sunAlpha = float(myargs.alphaSun)
    sunElev = ca.radec2elev(sunAlpha,sunDelta,myTel.lat)


print('SUN:')
print(f'\tLocal time: {((180+sunAlpha)/15.)%24:.2f}h')
print(f'\tHA = {sunAlpha:.1f}deg  = {(sunAlpha/15.)%24:.2f}h, Dec = {sunDelta:.1f}d')
print(f'\tElevation: {sunElev:.2f}d')
    
sunAz,sunEl = ca.radec2azel(sunAlpha, sunDelta, myTel.lat)
print(f'Validation: az= {sunAz:.2f}, el= {sunEl:.2f}d\n')





# output file
outfileroot  = f'{myargs.code}_{myargs.constellations}_'
outfileroot += f'{myargs.magSelect}_{int(myTel.lat):02d}_{int(sunAlpha):02d}'

#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#
#---  COMPUTE CONSTELLATIONS
#

# expand the constellation Id into a list of real constellations
CONSTELLATIONS = ca.findConstellations(myargs.constellations)
print('=====CONSTELLATIONS===================================================')
print( CONSTELLATIONS.ToC )
print()




# Azimuth-Elevation mesh:
AzEl = ca.fill_AzEl(step)               # mesh of Azimut-Elevation

# Arrays with the various results; same array geometry as AzEl 
densSatAll = np.zeros_like(AzEl[0])    # density of satellites     (all sat)
densVelAll = np.zeros_like(AzEl[0])    # density of satVelocities (all sat)
densSatObs = np.zeros_like(AzEl[0])    # density of Observable satellites (brighter than limiting mag)
densVelObs = np.zeros_like(AzEl[0])    #                       satVel
densSatBloom = np.zeros_like(AzEl[0])  # density of blooming satellites  (brighter than bloom limit)
densVelBloom = np.zeros_like(AzEl[0])  #                     satVel

fluxSatTotal = np.zeros_like(AzEl[0])  # total flux (for all sat)  

magmax = -99.
magmin = 99.
mageffmax = -99
mageffmin = 99.


# Scan the constellation shells
for myShell in CONSTELLATIONS.shells:                 
    # model the shell
    densSi, veli, magi =  ca.modelOneConstMag(AzEl,myTel.lat, sunAlpha,sunDelta,
                            myShell.inc,
                            myShell.alt, 
                            myShell.totSat
                            )

    # extreme magnitudes
    magmax = max(magmax,np.amax(magi))
    magmin = min(magmin,np.amin(magi))
    
    # all sat:
    densSatAll += densSi
    densVelAll += densSi * veli

    # effective magnitude and extremes
    # does the satellite trail more than 1 resolution element during expT:
    trailing =  veli*myTel.expt >= myTel.resol  # bolean for non-zero trailing
    mageffi = magi*1.  # init effective mag for non-trailing
    mageffi[trailing] = magi[trailing]   - 2.5*np.log10(myTel.resol/ (veli[trailing]*myTel.expt)  )
                       # correct for trailing sat
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

    #flux    
    fluxi = 10.**(-0.4*magi) /3600.**2 * densSi
    fluxSatTotal += fluxi

    
# here, we have the following arrays defined:
#   densSatAll: dens in Nsat/sq.deg
#   densVelAll: density if trails in Ntrail/deg/sec
#   densSatObs: same, only for sat with mag < myTel.maglim
#   densVelObs:
#   densSatBloom: same, only for sat with mag < myTel.magbloom
#   densVelBloom
#   fluxSatTotal: flux density for mag/sq.arcsec, for m550 at 550km

if myargs.plotflag:
    print( f'Satellite magnitudes in [{magmax:.2f},{magmin:.2f}]')
    print( f'Satellite eff. mag.  in [{mageffmax:.2f},{mageffmin:.2f}]')
    print( f'\nFolliwing output for zenith: (AzAlt={AzEl[:,-1,-1]})')
    print('- Sat density [n/sq.dg]', round( densSatAll[-1,-1],  4))
    print('- Sat velocity (for last constellation) [deg/s]', round(veli[-1,-1],4))
    print(f'- Diffuse mag [mag/sq/arcsec]: {-2.5*np.log10(fluxSatTotal[-1,-1]):.2f}' )

# sat count for almucantars
elLim = [60.,30.,20., 10.,0.]
elCount = ca.integrateSat(elLim,AzEl,densSatAll)
         # elCount: number of sat higher than elLim


outfile = open(outpath+outfileroot+'.txt','w+')  # store the hist.
outstring =  ('{:6.3f} {:5.1f} '+' {:4.0f}'*16).format(
    sunAlpha,
    sunElev,
    0.,0.,0.,0.,
    elCount[1],    
    elCount[2],    
    elCount[3],    
    elCount[4],
    0.,0.,0.,0.,
    0.,0.,0.,0.)
#will be written to file later

#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#
# PREPARE THE PLOT
#

# select effective densities
if  myargs.magSelect == 'BRIGHT':
    ds = densSatBloom
    dv = densVelBloom

    selectionLab = f'Selection: mag < {myTel.magbloom:.1f}'

elif myargs.magSelect == 'OBS':
    ds = densSatObs
    dv = densVelObs
    selectionLab = f'Selection: mag$_{{eff}}$ < {myTel.maglim:.1f}'

elif myargs.magSelect == 'FAINT':
    ds = densSatAll - densSatBloom
    dv = densVelAll - densVelBloom    
    selectionLab = f'Selection: mag > {myTel.magbloom:.1f}'
    
elif myargs.magSelect == 'ALL':
    ds = densSatAll
    dv = densVelAll

elif myargs.magSelect == 'EFFECT':
    ds = myTel.trailf * densSatObs + (1.-myTel.trailf)* densSatBloom
    dv = myTel.trailf * densVelObs + (1.-myTel.trailf)* densVelBloom
    selectionLab  = 'Selection: all satellites, scaled for effect. '
    selectionLab += f'Detected: V$_{{eff}}$ < {myTel.maglim:.1f} '
    selectionLab += f'Bleeding: V$_{{eff}}$ < {myTel.magbloom:.1f}'


else:
    print(f'invalid mode {myargs.magSelect}')
    exit(1)
    
# colormap
cmap = "magma"

# select what to plot and  limits for the colormap

if myargs.code == "TrailogDensity":
    barLabel = "Number of trails./deg/sec."

    logDensity = np.log10( dv )
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_standardLog(logDensity)
    print("Tdensity")

elif myargs.code == "SatDens":
    barLabel = "Number of sat./sq.deg."

    logDensity = np.log10( ds )
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_standardLog(logDensity)
    print("Sdensity")

elif myargs.code == "skyMag":
    if skyFluxFlag:
        barLabel = r'Surface brightness [$\mu$cd/m$^2$]' # raw string for LaTeX
    
        logDensity =  2.5* np.log10( fluxSatTotal ) # = -1*mag
            # we plot -mag, then we change the scale of the bar
        logMinValue, logMaxValue = cp.getBarLim(logDensity)
        bMin = int(logMinValue*3.)/3. 
        bMax = int(logMaxValue*3. -1)/3.
        barTicks = np.arange(bMin, bMax , .333333)
        muCd = 12e10 * 10**(0.4*barTicks) # conversion to microCandela/m2
        barTickLabels = [ f'{x:.1g}' for x in muCd]
        print("skyFlux")

    else:
        barLabel = "Surface brightness [mag/sq.arcsec]"

        logDensity =  2.5* np.log10( fluxSatTotal ) # = -1*mag
        barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_negMag(logDensity)
        print("skyMag")



elif myargs.code == "skyFrac":  # fraction of the sky surfbrightness
    skymag = 22. 
    barLabel = f'Surface brightness [fraction of sky] ($m_{{sky}} = ${skymag:.2f})'
    
    logDensity =  2.5* np.log10( fluxSatTotal ) + skymag
    #+ 5 # +5 for [%]

    logMinValue, logMaxValue = cp.getBarLim(logDensity)
    bMin = int(logMinValue*3.)/3. 
    bMax = int(logMaxValue*3. -1)/3.
    barTicks = np.arange(bMin, bMax , .333333)
    barTickLabels = [ f'{x:.1g}' for x in 10.**barTicks]
    print("skyFrac")




else: # other specific (including EFFECT)

    dens    =  ds* myTel.fovl*myTel.fovw + dv * myTel.fovl * myTel.expt
             # trailf already accounted for in ds and dv

    # deal with empty sky
    if len( dens[ dens > 0]  ) > 0:
        logDensity   = np.nan_to_num(np.log10(dens ), neginf=-999. )
    else:
        logDensity   = dens*0 - 1000.

    densobs = densSatObs * myTel.fovl*myTel.fovw +  densVelObs* myTel.fovl *myTel.expt
                     # number of observable trails


    # compute average effect:
    EffTot = 0.
    TrailTot = 0.
    icount = 0
    aircut = 20. #30. # start counting at airmass =2
    for i in np.arange(0,len( AzEl[1,:,1]) ):
        if AzEl[1,i,1] >= aircut:
            EffTot   += np.average(dens[i,:]   ) * np.cos(np.radians(AzEl[1,i,1]))
            TrailTot += np.average(densobs[i,:]) * np.cos(np.radians(AzEl[1,i,1]))
            icount += 1
            
    EffTot   = EffTot  /icount
    TrailTot = TrailTot/icount

    print(f'Effect on exposures (at Zenith): Loss fraction: {dens[-1,-1]:.3g}/1.; Trails: {densobs[-1,-1]:.3g}/exp')
    print(f'Effect on exp. (aver above {aircut}): Loss fraction: {EffTot:.3g}/1.; Trails: {TrailTot:.3g}/exp')


    outstring += f' {EffTot} {TrailTot}'

    
    if myargs.magSelect == 'EFFECT':
        #print("Effect")
        barLabel = "Fraction lost"
        cmap = gyrd    
        
        logMinValue = -3.5  # log limits for the colour scale
        logMaxValue = 0.5  # 2.2
        bMin = int(logMinValue*3.)/3. 
        bMax = int(logMaxValue*3. -1)/3.
        barTicks = np.arange(bMin, bMax , .333333)
        barTickLabels = [ f'{x:.1g}' for x in 10.**barTicks]




    else:
        barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_standardLog(logDensity)

        print('HEREl', logMinValue, logMaxValue)
        print('HEREb', barTicks, barTickLabels)
        barLabel = "Number of trails per exp."



#
# PLOT
#

if not myargs.plotflag:
    myargs.labelplotflag  = False
    myargs.scalebarflag = False
    myargs.shadeflag = False
else:
    fig = plt.figure(figsize=(8,8))
    ax =  fig.subplots(1,1,subplot_kw={'projection': 'polar'}) 
    cp.initPolPlot(ax)
    ax.set_facecolor("k")
    
    clab = 'k'
    ccon = 'k'
    tickformat = "{:.1f}".format

    #    logDensity[ logDensity > logMaxValue  ] = logMaxValue
    #    logDensity[ logDensity < logMinValue  ] = logMinValue
 
    cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1], logDensity , 
                      levels=np.linspace(logMinValue,logMaxValue,100),  # NUMBER OF LEVELS
                      vmin=logMinValue, vmax=logMaxValue ,
                      extend='both',
                      cmap=cmap)
    cfd.cmap.set_under('k') # below minimum -> black



#----------------------------------------------------------------------
#Scalebar
if myargs.scalebarflag:
    cbar = fig.colorbar(cfd)
    cbar.set_ticks( barTicks )
    cbar.set_ticklabels( barTickLabels)
    cbar.set_label(barLabel)


#----------------------------------------------------------------------
# airmass shade

if myargs.shadeflag:
    waz = np.arange(0.,2.*np.pi,.01)
    wr1 = 0.*waz +90
    whss = [['xx','x'],['++','+']]
    wcs = ['white','black']

    for whs,wc in zip( whss,wcs ): 
        for wr,wh in zip ([20,30],whs):
            wr2 = wr1 -wr
            if 1 : # hash
                alphashade=0.1
                ax.fill_between(waz,wr1,wr2,
                    facecolor='none',
                            hatch=wh,edgecolor=wc,alpha=alphashade)
                ax.fill_between(waz,wr1,wr2,
                    facecolor='none',
                            hatch=wh,edgecolor=wc,alpha=alphashade)
            else:
                alphashade=0.2
                ax.fill_between(waz,wr1,wr2,
                                color="grey", alpha=alphashade)


#------------------------------------------------------------------------------
# Draw RA,Dec lines
cp.drawHADec(myTel.lat)


#----------------------------------------------------------
#All the labels

if myargs.labelplotflag:

    #Sun symbol on horizon
    plt.text(np.radians(sunAz), 93.,r'$\odot$', va="center", ha='center') # raw string for LaTeX
    

    #top left corner
    x = -1.
    y = 1.2
    dy = 0.08
    
    cp.azlab(ax,x,y,'Observatory: {} Lat.: {:.1f}$^o$'.format(myTel.telescope, myTel.lat))
    y -= dy
    
    
    if myargs.code != "SatDens" and myargs.code != "TrailogDensity" and myargs.code != "skyMag":
        cp.azlab(ax,x,y,'Instrument: {}'.format(myTel.instrument))
        y -= dy

        if myTel.fovl < 1./60.:
            fovll = '{:.2f}\"'.format(myTel.fovl*3600.)
        elif myTel.fovl < 1./6.:
            fovll = '{:.2f}\''.format(myTel.fovl*60.)
        else:
            fovll = '{:.2f}$^o$'.format(myTel.fovl)

        if myTel.fovw < 1./60.:
            fovlw = '{:.2f}\"'.format(myTel.fovw*3600.)
        elif myTel.fovw < 1./6.:
            fovlw = '{:.2f}\''.format(myTel.fovw*60.)
        else:
            fovlw = '{:.2f}$^o$'.format(myTel.fovw)

        cp.azlab(ax,x,y,'Fov: '+fovll+'x'+fovlw)
        y -= dy


        if myTel.expt > 1.:
            cp.azlab(ax,x,y,'Exp.t: {:.0f}s'.format(myTel.expt))
        else:   
            cp.azlab(ax,x,y,'Exp.t: {:.0g}s'.format(myTel.expt))
        y -= dy
        

            
    # bottom left
    x= -1.
    y= -1.08
    cp.azlab(ax,x,y,r'$\odot$ Sun:',14)

    y -= dy
    loct = (sunAlpha/15.+12.)%24
    
    loch = int(loct)
    locm = int( (loct-loch)*60.)
    cp.azlab(ax,x,y,f'Loc.time: {loch:02d}:{locm:02d}')
    y -= dy
    cp.azlab(ax,x,y,r'$\delta: '+f'{sunDelta:.2f}^o$, Elev: {sunElev:.2f}$^o$')
    y -= dy




    # top right
    x=1.
    y=1.2
    cp.azlab(ax,x,y,'Constellation:',14)
    y -= dy
    cp.azlab(ax,x,y,CONSTELLATIONS.name)
    
    y -= dy
    cp.azlab(ax,x,y, f'Total {CONSTELLATIONS.totSat:.0f} sat.')

    #bottom right
    x = 1.
    y = -1.08 -dy
    if 1:
        lab = "Satellite magnitudes: V$_{550km}=$" +\
            f'{cst.mag550:3.1f}' 
#            "{:3.1f}".format(cst.mag550 -5.*np.log10(550./1000.) )
        cp.azlab(ax,x,y,lab)
        y -= dy

    if mageffmin < -1000. or myargs.code == "skyMag":
        lab = "  V$_{sat}$"+" in [{:.1f}, {:.1f}]".format(magmax,magmin)
    else:
        lab = "  V$_{sat}$"+" in [{:.1f}, {:.1f}]".format(magmax,magmin)+\
              "  V$_{eff}$"+" in [{:.1f}, {:.1f}]".format(mageffmax,mageffmin)


    cp.azlab(ax,x,y,lab)
    y -= dy

    if myargs.magSelect == "BRIGHT":
        cp.azlab(ax,x,y,selectionLab)
    elif myargs.magSelect == "OBS":
        wlab = "Selection: mag$_{eff}$ < "+"{:.1f}".format(myTel.maglim)
        cp.azlab(ax,x,y,wlab)
    elif myargs.magSelect == "FAINT":
        cp.azlab(ax,x,y,"Selection: mag > {:.0f}".format(myTel.magbloom))
    elif myargs.magSelect == "EFFECT":
        cp.azlab(ax,x,y,"Selection: all satellites, scaled for effect")
        wlab = "Detected: V$_{eff}$ < "+"{:.1f} ".format(myTel.maglim)
        wlab += "   Bleeding: V$_{eff}$ < "+"{:.1f}".format(myTel.magbloom)
        y -= dy
        cp.azlab(ax,x,y,wlab)


    
# sat count on almucantars
if myargs.almucantar and myargs.plotflag:
    for we, wi in zip(elLim, elCount):
        cp.azlab(ax,-0,(90.-we-5)/90.,'{:.0f} sat.>{:.0f}$^o$:'.format(wi,we),9,0.5*(1.-we/100.))


# Discrete satellites as dots on the plot
if myargs.noDots:
    print("==DOTS==")
    Sv = satDots.makeConstellationStatTable(CONSTELLATIONS, 
                                            sunAlpha, sunDelta, 
                                            myTel.lat, 
                                            (180+sunAlpha)*360.) ## slowed 10x
    if myargs.magSelect == "EFFECT":
        myargs.magSelect = "ALL"

    if myargs.magSelect == "ALL":
        # not illuminated in grey
        Si = Sv[  ~Sv["bIlluminated"] ]
        ax.scatter(Si["Azr"],Si["ZD"], s=Si["dot"], c="grey", alpha=0.2)

    Si = Sv[  Sv["bIlluminated"] ]
    if myargs.magSelect in ["ALL",  "OBS"] :
        # not bright in yellow
        Sb = Si[ Si["mag"] >= 7 ]
        ax.scatter(Si["Azr"],Si["ZD"], s=Si["dot"], c="yellow")


    if myargs.magSelect in ["ALL", "OBS","BRIGHT"] :
        # bright in red
        Sb = Si[ Si["mag"] < 7 ]
        ax.scatter(Sb["Azr"],Sb["ZD"], s=Sb["dot"], c="red")




print('finishing...')

outstring += '\n'
outfile.write(outstring)
outfile.close()
        

# save plot
if myargs.plotflag:
    fig.tight_layout()
    plt.savefig(outpath+'w.png')
    filename = outpath+outfileroot+myargs.outputformat
    print(filename)

    plt.savefig(filename)
#--
print("output in ",outfileroot)


