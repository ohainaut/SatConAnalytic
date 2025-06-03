#!/usr/bin/env python3
'''SatConAnalytic - Satellite Constellation Analytic simulations

Plot debris density over the map of the sky.

option -h for detailed usage

Relies on ESA debris distribution in Debris/spatial_density_1e-*.csv

'''

#import matplotlib
#####matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt
import numpy as np
import argparse

print('conAn ObsSkyDebris')

from conanplot import gyrd


# import ConAn routines
import conan as ca
import conanplot as cp
import constants as cst
import satDots
import debris


#----- config
step = 1.#30. #deg >~1. Smaller values take forever
multiplicationFactor =  1 # magic factor wrt today



outpath = "./"


#------Arguments

parser = argparse.ArgumentParser(description=
'''Debris contribution to the sky background.''')

parser.add_argument('-d','--deltaSun', default=0.,
                    help="Sun: Declination of the Sun [deg]")
parser.add_argument('-a','--alphaSun',
                    help="Sun: Hour Angle of the Sun [deg]. If present, overwrites elevSun")
parser.add_argument('-e','--elevSun', default=24.,
                    help="Sun: Elevation of the Sun BELOW the horizon. Should probably be >0 in most cases [deg]")
parser.add_argument('-M','--mode', default="ratio", choices=['totalMag', 'skyMag', 'ratio', 'debrisMag',  'debrisFlux'],
                    help="Mode: what to plot")
parser.add_argument('-l','--lat', default=-24.6,
                    help="Observatory: Latitude of the observatory [deg]")
parser.add_argument('-s','--telescope',
                    help="Observatory: Name of the telescope")
parser.add_argument('-i','--instrument',
                    help="Observatory: Name of the instrument")
parser.add_argument('--alt', nargs=2, default=[190,1000],
                    help="min  max altitude, in km")
parser.add_argument('--lrad', nargs=2, default=[0,6],
                    help="min max radius in -log (default: 0 6)")


parser.add_argument('-T','--code',
                    help='''Observatory: Predefined telescope/instrument with
                         extptime, FoVl, FoVw, maglim, magbloom, trailf,
                         telescope, instrument, resolution, latitude.
                         Use individual options to overwrite presets.
                         SPECIAL CODES: SatDens for satellite density;
                         TrailDens for trail density;
                         skyMag for sky surface brightness [mag];
                         skyFrac for sky surface brightness as a fraction.
''')
myargs = parser.parse_args()
myargs.outputformat = ".png"

# select debris shell min and max altitude
altmin = float(myargs.alt[0])
altmax = float(myargs.alt[1])

# select debris min and max size (in -log[m], from 0 to 6)
radmin = int(myargs.lrad[0])
radmax = int(myargs.lrad[1])+1

#
# GET TELESCOPE/OBSERVATORY
#
print('TELESCOPE/INSTRUMENT SETUP')
#for backward compatibility (ignore... not used here)
myargs.resol = 1.
myargs.expt = 1.
myargs.fovl = 1.
myargs.fovw = 1.
myargs.maglim = 99.
myargs.magbloom = 99.
myargs.trailf = 1.
myTel = cp.getTelescope(myargs)
print(myTel)

#-------------------------------------------------------------------------------------
#
# SUN position 
#

sunDelta =  float(myargs.deltaSun)
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

wAz, wEl = ca.radec2azel(sunAlpha, sunDelta, myTel.lat)
print(f'Validation: sun az= {wAz:.2f}, el= {wEl:.2f}d\n')


# output file
outfileroot  = f'DEBRIS_{myargs.code}_'
outfileroot += f'{myargs.mode}_{int(myTel.lat):02d}_{int(sunAlpha):02d}'




#
# Read in debris data
#
 
# restricting on requested radius range
print("=Read debris definition")
DEBRIS = debris.make_debrisTab( myrange=np.arange(radmin,radmax))

# Azimut-Elevation mesh
AzEl = ca.fillAzEl(step)
AzElreshape = np.reshape(AzEl,(2,AzEl.shape[1]*AzEl.shape[2]))

#
# Natural Sky brightness
# 

# accounting for the twilight/ sun elevation
magSkyZenith = debris.modelTwilight(sunElev)
# and airmass/elevation 
magSky = np.reshape(debris.modelSkyBrightness(AzElreshape[1,:] , 
                                              mag=magSkyZenith),
                            (AzEl.shape[1],AzEl.shape[2]) )
fluxSky = 10.**(-0.4* magSky)


#
# LOOP THROUGH THE DEBRIS SHELLS
#
fluxDebris = np.zeros(   (AzEl.shape[1],AzEl.shape[2]) )

for i in  np.arange(len(DEBRIS) -1):
    # restrict to requested altitude range
    if DEBRIS['alt'][i] >= altmin and DEBRIS['alt'][i] <= altmax:
        print(i, DEBRIS["alt"][i])
        # model the shell, and summ the contribution
        fluxDebris +=  np.reshape(debris.modelOneConstMag(AzElreshape,myTel.lat, sunAlpha,sunDelta,
                                    DEBRIS['alt'][i],
                                    DEBRIS['surf_a'][i]   ),
                            (AzEl.shape[1],AzEl.shape[2]) )
        # fluxDebris now contains the total flux from all debris for each point of the AzEl mesh

# Multiplication of current values
fluxDebris *= multiplicationFactor

magDebris = -2.5*np.log10(fluxDebris) 
magDebris[ fluxDebris == 0 ] = 9999. # avoid NaN=log(0)

magTotal = -2.5*np.log10(fluxDebris + fluxSky) 

fluxRatio = fluxDebris / fluxSky


#==============================================================================


fig = plt.figure(figsize=(8,8))
ax =  fig.subplots(1,1,subplot_kw={'projection': 'polar'})
cp.initPolPlot(ax)
ax.set_facecolor("k")
nTick = 9

if myargs.mode == "skyMag":
    plotit = -magSky 
    vmin = np.amin( plotit )
    vmax = np.amax( plotit )
    vsign = -1.
    barLabel = "Sky background [Mag/arcsec$^2$]"
    cmap = cp.csunmap
    def barFmt(x):
        return f"{x:.1f}"

elif myargs.mode == "debrisMag":
    plotit = -magDebris
    vmin = np.amin( plotit[plotit > -990] )
    vmax = np.amax( plotit[plotit > -990] )
    vsign = -1.
    barLabel = "Debris [Mag/arcsec$^2$]"
    cmap = cp.csunmap
    def barFmt(x):
        return f"{x:.1f}"
    
elif myargs.mode == "totalMag":
    plotit = -magTotal
    vmin = np.amin( plotit[plotit > -990] )
    vmax = np.amax( plotit[plotit > -990] )
    vsign = -1.
    barLabel = "Debris [Mag/arcsec$^2$]"
    cmap = cp.csunmap
    def barFmt(x):
        return f"{x:.1f}"



elif myargs.mode == "debrisFlux":
    plotit = fluxDebris
    vmin = np.amin( plotit[plotit > 0] )
    vmax = np.amax( plotit[plotit > 0] )
    if vmin == vmax:
        vmin *= .9
        vmax *= 1.1

    vsign = 1.
    barLabel = "Debris [flux/arcsec$^2$]"
    cmap = cp.csunmap
    def barFmt(x):
        return f"{x:.0f}"



elif myargs.mode == "ratio":
    plotit = np.log10(fluxRatio)
    vmin = -5 # np.amin( plotit[fluxRatio > 0] )
    vmax = 1 #np.amax( plotit[fluxRatio > 0] )
    vsign = 1.
    barLabel = "Debris/Sky "
    cmap = cp.gyrd
    nTick = vmax-vmin +1
    def barFmt(x):
        return f"$10^{{{x:.0f}}}$"

print("CONTOURSminmax",  vmin, vmax)

cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1],
                      plotit ,
                      levels=np.linspace(vmin,vmax,nTick*3-2),  # NUMBER OF LEVELS
                      vmin=vmin, vmax=vmax ,
                      extend='both',
                      cmap=cmap)
cfd.cmap.set_under('k') # below minimum -> black


#----------------------------------------------------------------------
#Scalebar
if 1:
    cbar = fig.colorbar(cfd)
    blab = ""
    bnorm = 0. # log of normalization.
    barmag = np.linspace(vmin, vmax, nTick)

    if np.log10(vmin) < 0.:
        print ("LOG")
        bnorm = int( np.log10(vmin) )
        wnorm = int( np.log10( vmax/vmin ) +.5 ) 

        print(vmin, vmax)
        print(bnorm, bnorm+ wnorm)

        barmag = 10.**np.arange(bnorm, bnorm+ wnorm, .5)

        nTick = wnorm+1
        blab = r' $\times 10^'+f'{{{bnorm}}}$'

    print('VMINMAX',vmin, vmax)

    cbar.set_ticks(barmag)
    cbar.set_ticklabels( [ barFmt(x) for x in barmag*10**-bnorm ])
    cbar.set_label(barLabel + blab)


#------------------------------------------------------------------------------
# RA Dec lines
cp.drawHADec(myTel.lat)


#----------------------------------------------------------
#labels



if 1:

    #Sun
    azs,els = ca.radec2azel(sunAlpha, sunDelta, myTel.lat)
    plt.text(np.radians(azs), 93.,r'$\odot$', va="center", ha='center') # raw string for LaTeX


    #top left
    x = -1.
    y = 1.2
    dy = 0.08

    cp.azlab(ax,x,y,'Observatory: {} Lat.: {:.1f}$^o$'.format(myTel.telescope, myTel.lat))
    y -= dy



    #top right
    x = 1.
    y = 1.2
    dy = 0.08

    debris.get_totalShellSurface(DEBRIS)
    totalMass = sum( DEBRIS['surf_m']* DEBRIS['surface'])*multiplicationFactor

    cp.azlab(ax,x,y,f'Albedo: {debris.p_albedo:.2f}')
    y -= dy

    cp.azlab(ax,x,y,f'Debris mass: {totalMass:.2e} kg')
    y -= dy

    cp.azlab(ax,x,y,f'Factor: {multiplicationFactor:.2e} x today')
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
    cp.azlab(ax,x,y,r'Sky brightness at Z: '+f'{magSkyZenith:.1f} mag/sq.arcsec')


    


    # bottom right
    x= 1.
    y= -1.08
    cp.azlab(ax,x,y,f'Altitude: {altmin:.0f} .. {altmax:.0f} km',14)
    y -= dy
    mylabel = f'Size: $10^{{ {-radmin:.0f} }}$ .. $10^{{ {-radmax+1:.0f} }}$ m' # {{: escaped {
    cp.azlab(ax,x,y,mylabel,14)



# save plot
if 1:
    fig.tight_layout()
    plt.savefig(outpath+'w.png')
    filename = outpath+outfileroot+myargs.outputformat
    print(filename)
    plt.savefig(filename)
#--
print("output in ",outfileroot)


