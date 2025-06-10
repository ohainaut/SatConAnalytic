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
import constants
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
parser.add_argument('-l','--lat', default=-24.6,
                    help="Observatory: Latitude of the observatory [deg]")
parser.add_argument('-s','--telescope',
                    help="Observatory: Name of the telescope")
parser.add_argument('-i','--instrument',
                    help="Observatory: Name of the instrument")
parser.add_argument('--alt', nargs=2, default=[0,2000],
                    help="min  max altitude, in km")
parser.add_argument('--lrad', nargs=2, default=[0,6],
                    help="min max radius in -log (default: 0 6)")
parser.add_argument('--density', default=1.,
                    help='''Density threshold [n/sq.dg]: 
                    particle less dense than this don\'t contribute 
                    to the diffuse background''')
parser.add_argument('-T','--code',
                    help='''Observatory: Predefined telescope/instrument with
                         extptime, FoVl, FoVw, maglim, magbloom, trailf,
                         telescope, instrument, resolution, latitude.
                         ''')
parser.add_argument('-M','--mode', default="ratio", 
                    choices=['skyMag', 'skyFlux',
                             'debrisMag',  'debrisFlux', 'debrisCount', 
                             'totalMag', 'ratio', 
                             'surface'],
                    help="Mode: what to plot")

myargs = parser.parse_args()
myargs.outputformat = ".png"


minimumDensity = float(myargs.density)

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
print('Number of altitude shells:',len(DEBRIS))


# Azimut-Elevation mesh
AzEl = ca.fill_AzEl(step)
AzEl_r = np.reshape(AzEl,(2,AzEl.shape[1]*AzEl.shape[2]))
   # all *_r are "reshaped"

#
# Natural Sky brightness
# 

# accounting for the twilight/ sun elevation
magSkyZenith = debris.modelTwilight(sunElev)
# and airmass/elevation 
magSky = np.reshape(debris.modelSkyBrightness(AzEl_r[1,:] , 
                                              mag=magSkyZenith),
                            (AzEl.shape[1],AzEl.shape[2]) )
fluxSky = 10.**(-0.4* magSky)


#
# LOOP THROUGH THE DEBRIS SHELLS
#
fluxDebris  = np.zeros(   (AzEl.shape[1],AzEl.shape[2]) )
countDebris = np.zeros(   (AzEl.shape[1],AzEl.shape[2]) )

countShells = 0

for i in  np.arange(len(DEBRIS) -1):
    # restrict to requested altitude range
    if DEBRIS['alt'][i] >= altmin and DEBRIS['alt'][i] <= altmax:
        countShells += 1

        satAlt = DEBRIS['alt'][i] 

        # compute the geometry for each AzEl pointing
        surface_r, illum_r, sunPhFrac_r, Delta_r = debris.modelOneConst_Geometry(
            AzEl_r, myTel.lat, sunAlpha, sunDelta,  satAlt )

        # model the particle count
        shell_countDebris_r = debris.modelOneConst_Count(
                surface_r, illum_r, sunPhFrac_r, Delta_r,
                DEBRIS['surf_n'][i]   )

        # model the flux
        shell_fluxDebris_r = debris.modelOneConst_Flux(
                surface_r, illum_r, sunPhFrac_r, Delta_r,
                DEBRIS['surf_a'][i]   )
        # account for extinction
        shell_fluxDebris_r *= 10.**( -.4* constants.extinction * (Delta_r/satAlt))


        # contribution to the background is negligible if less than 
        # 1 particle per sq.deg
        shell_fluxDebris_r[ shell_countDebris_r <  minimumDensity ] = 0.


        # summ the contribution
        countDebris +=  np.reshape( shell_countDebris_r,  (AzEl.shape[1],AzEl.shape[2]) )
        fluxDebris +=  np.reshape( shell_fluxDebris_r, (AzEl.shape[1],AzEl.shape[2]) )
        # fluxDebris now contains the total flux from all debris for each point of the AzEl mesh



print(f'Processed {countShells} shells')


# Multiplication of current values
fluxDebris *= multiplicationFactor

# Conversion to sq.arcsec; values are computed for 1sqdeg on sky
fluxDebris /= (3600.*step)**2

magDebris = -2.5*np.log10(fluxDebris) 
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
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_negMag(plotit)
    barLabel = "Sky background [Mag/arcsec$^2$]"
    zenithLabel =''
    cmap = cp.csunmap
    showDebrisLabel = False

if myargs.mode == "skyFlux":
    plotit = -magSky 
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_negMag(plotit)

    # for the unit conversion, we change the labels:
    barTickLabels = [f'{w:.0f}' for w in ca.mag2microcd(-barTicks)] #barTicks is -mag
    barLabel = "Surface brightness [$\mu$cd/m$^2$]"

    zenithLabel =f'Sky flux: {ca.mag2microcd(magSky[-1,0]):.1f} $\mu$cd/m$^2$' 
    cmap = cp.csunmap
    showDebrisLabel = False

elif myargs.mode == "debrisMag":
    plotit = -magDebris
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_negMag(plotit)
    barLabel = "Debris [Mag/arcsec$^2$]"
    cmap = cp.csunmap
    showDebrisLabel = True
    zenithLabel =f'Debris Mag: {-plotit[-1,0]:.1f} mag/arcsec$^2$)' 
        
elif myargs.mode == "totalMag":
    plotit = -magTotal
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_negMag(plotit)
    barLabel = "Sky+Debris [Mag/arcsec$^2$]"
    cmap = cp.csunmap
    showDebrisLabel = True
    zenithLabel =f'Total Mag: {-plotit[-1,0]:.1f} mag/arcsec$^2$)' 


elif myargs.mode == "debrisFlux":
    plotit = -magDebris
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_negMag(plotit)
    
    # for the unit conversion, we change the labels:
    barTickLabels = [f'{w:.2g}' for w in ca.mag2microcd(-barTicks)] #barTicks is -mag
    barLabel = "Surface brightness [$\mu$cd/m$^2$]"
    cmap = cp.csunmap
    showDebrisLabel = True
    zenithLabel =  f'Flux: Debris= {ca.mag2microcd(magDebris[-1,0]):.0f}, ' 
    zenithLabel += f'Sky= {ca.mag2microcd(magSky[-1,0]):.0f} $\mu$cd/m$^2$' 


elif myargs.mode == "debrisCount":
    plotit = np.log10(countDebris)
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_standardLog(plotit)
    barLabel = "Debris [N/deg$^2$]"
    cmap = cp.csunmap
    showDebrisLabel = True
    zenithLabel =f'Debris count: {10.**plotit[-1,0]:.1g} grain/sq.deg' 

elif myargs.mode == "ratio":
    plotit = np.log10(fluxRatio)
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_standardLog(plotit)
    barLabel = "Debris flux/Sky flux "
    cmap = cp.gyrd
    cmap = cp.csunmap

    showDebrisLabel = True
    wpc = 10**(-0.4*(magDebris[-1,0] - magSkyZenith))*100.
    zenithLabel =f'Debris brightness: {magDebris[-1,0]:.1f} = {wpc:.2f}%' 


elif myargs.mode == "surface":
    plotit = np.log10( np.reshape(surface_r,(AzEl.shape[1],AzEl.shape[2]) ) )
    barTicks, barTickLabels, logMinValue, logMaxValue = cp.setBarLim_standardLog(plotit)
    
    barLabel = f'Surface [km2] for a {step**2} sq.deg at alt={satAlt:.0f}km '
    cmap = cp.csunmap

    showDebrisLabel = False
    zenithLabel =f'Surface: {10.**plotit[-1,0]:.1f} km$^2$/sq.deg' 

else:
    print('uh?')

# actual plot

cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1],
                      plotit ,
                      levels=50,
                      vmin=logMinValue, vmax=logMaxValue ,
                      extend='both',
                      cmap=cmap)
cfd.cmap.set_under('k') # below minimum -> black
cfd.cmap.set_over(cmap(1.)) 


#----------------------------------------------------------------------
#Scalebar




if 1:
    cbar = fig.colorbar(cfd)
    cbar.set_ticks( barTicks )
    cbar.set_ticklabels( barTickLabels)
    cbar.set_label(barLabel)


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
    if showDebrisLabel:

        mylabel = f'Debris: Size: $10^{{ {-radmin:.0f} }}$ .. $10^{{ {-radmax+1:.0f} }}$ m' # {{: escaped {
        cp.azlab(ax,x,y,mylabel,14)
        y -= dy

        cp.azlab(ax,x,y,f'Altitude: {min(DEBRIS["alt"]):.0f} .. {max(DEBRIS["alt"]):.0f} km')
        y -= dy

        DEBRIS['surface'] = 4.*np.pi* (constants.earthRadius + DEBRIS['alt'])
        totalMass = sum( DEBRIS['surf_m']* DEBRIS['surface'])*multiplicationFactor
        cp.azlab(ax,x,y,f'Total mass: {totalMass:.2e} kg')
        y -= dy

        cp.azlab(ax,x,y,f'= {multiplicationFactor:.1e} x today')
        y -= dy


        cp.azlab(ax,x,y,f'Albedo: {debris.p_albedo:.2f}')
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


    # bottom right
    x= 1.
    y= -1.08
    cp.azlab(ax,x,y,f'Zenith: ',14)
    
    y -= dy
    cp.azlab(ax,x,y,f'Sky brightness: {magSkyZenith:.1f} mag/sq.arcsec')
    y -= dy
    cp.azlab(ax,x,y,zenithLabel)



# save plot
if 1:
    fig.tight_layout()
    plt.savefig(outpath+'w.png')
    filename = outpath+outfileroot+myargs.outputformat
    print(filename)
    plt.savefig(filename)
#--
print("output in ",outfileroot)


