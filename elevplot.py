#!/usr/bin/env python3
'''SatConAnalytic - Satellite Constellation Analytic simulations

Plot satellite density over the map of the sky.

Alternatively, the plot can show
- velocity density of the satellites
- the number of satellites / detectable satellites / saturating satellites / 
  non-saturating satellites in the field of view of an instrument
- the effect of the satellites on the observations (%loss)
- the sky brighness increase caused by satellites.

obsplot.py -h for usage
'''

import matplotlib
##matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
import os
#sys.path.append(os.path.dirname(__file__)+"/../")


from matplotlib import ticker
from conanplot import gyrd


# import ConAn routines
import conan as ca
import conanplot as cp
from constants import mag550

#----- config
step = 5. #deg >~1. Smaller values take forever
#outpath = "/home/ohainaut/public_html/outsideWorld/"
outpath = "./"


#------Arguments

parser = argparse.ArgumentParser(description=
'''Satellite constellations: sky map of satellite trail density or related information. 
Define the position of the observatory, the position of the Sun, the constellation(s),
the instrument and its characteristics, and the output required.''')
parser.add_argument('-d','--deltaSun', default=0.,
                    help="Sun: Declination of the Sun [deg]")
parser.add_argument('-C','--constellations', default='SLOWGWAK',
                    help="ID of the constellation group; list for a list")
parser.add_argument('-l','--lat', default=-24.6,
                    help="Observatory: Latitude of the observatory [deg]")
parser.add_argument('-r','--resol', 
                    help="Observatory: Resolution of the instrument [deg]")
parser.add_argument('-t','--expt', 
                    help="Observatory: Exposure time [sec]")
parser.add_argument('-f','--fovl',
                    help="Observatory: Field of view of the instrument. Length or diametre [arcsec]")
parser.add_argument('-w','--fovw', 
                    help="Observatory: Field of view of the instrument. Width. Equal to Length if omitted [arcsec]")
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
                         TrailDens for trail density;
                         skyMag for sky surface brightness [mag];
                         skyFrac for sky surface brightness as a fraction.
''')
parser.add_argument('-e','--elevMin', default=20.,
                    help="Elevation limit above which we observe")
parser.add_argument('-M','--mode', default="OBS",
                    help="Plot: ALL (Default), OBS, BRIGHT, FAINT, or EFFECT")
parser.add_argument('--pdf', action='store_true',
                    help="Plot: output file in pdf (default is png)")
myargs = parser.parse_args()


# flags
if myargs.pdf :
    myargs.outputformat = ".pdf"
else:
    myargs.outputformat = ".png"
    




#- find telescope
print('TELESCOPE/INSTRUMENT SETUP')
myTel = cp.getTelescope(myargs)
print(myTel)

# sun
sunDelta = float(myargs.deltaSun)
sunElev  = np.arange(6., -91.,-3.) # sun elevation
sunAlpha = ca.elev2ra(sunElev,sunDelta,myTel.lat) # get sun hourangle 
print('SUNALPHA',sunAlpha)
# get constellation id into a list of real constellations
CONSTELLATIONS = ca.findConstellations(myargs.constellations)
print('CONSTELLATIONS:')
print( CONSTELLATIONS.ToC )
print()



    
######################################################################
#---  COMPUTE CONSTELLATIONS

#fill ElAz:
AzEl = ca.fillAzEl(step)

elLim = [60.,30.,20., 10.,0.] # must be decreasing
densSatAll = np.zeros( (len(elLim), len(sunElev)))
# total density for each elevation, for all shells


elevSat = np.zeros( (len(elLim), len(sunElev))) 


for c in CONSTELLATIONS.list:     
    CONSTELLATIONS.byCode[c].elevSat = np.zeros( (len(elLim), len(sunElev))) 

    for myShell in CONSTELLATIONS.byCode[c].shells:     

        myShell.elevSati = np.zeros( (len(elLim), len(sunElev))) #5 like elLim

        for iSun in range(len(sunElev)): # scan night
        
            # model the shell
            densSi, veli, magi =  ca.modelOneConstMag(AzEl,myTel.lat, 
                                    sunAlpha[iSun],sunDelta,
                                    myShell.inc,
                                    myShell.alt, 
                                    myShell.totSat
                                    )


            # select satellites:
            # default: ALL illuminated (do nothing)
            mageffi = magi  - 2.5*np.log10(myTel.resol/veli/myTel.expt)
            if myargs.mode == "OBS":
                densSi[ mageffi > myTel.maglim] = 0.
            elif myargs.mode == "BRIGHT":
                densSi[ mageffi > myTel.magbloom ] = 0.


            # sat count for almucantars
            myShell.elevSati[:,iSun] = ca.integrateSat(elLim,AzEl,densSi)

        # finished one curve
        iEl = 4
        plt.plot( sunElev, myShell.elevSati[iEl], color='k', alpha=0.1)

        CONSTELLATIONS.byCode[c].elevSat += myShell.elevSati

    plt.plot( sunElev, CONSTELLATIONS.byCode[c].elevSat[iEl],
             label=f'{c} {CONSTELLATIONS.byCode[c].totSat}')
    
    densSatAll += CONSTELLATIONS.byCode[c].elevSat

plt.plot( sunElev, densSatAll[iEl],
             label=f'ALL: {CONSTELLATIONS.totSat}',
             color='k', linewidth=3)

plt.legend()
plt.xlim(6.,-90.)
plt.yscale('log')
plt.ylim(5.,1e4)
plt.ylabel(f'Number of illuinated satellites above {elLim[iEl]}')
plt.show() #block=False)



asdfas
# output file
outfileroot  = f'Elev_{myargs.code}_{myargs.constellations}_'
outfileroot += f'{myargs.mode}_{int(myTel.lat):02d}_{int(-sunElev):02d}'


    #==============================================================================

if myargs.mode == 'BRIGHT':
        ds = densSatBloom
        dv = densVelBloom
        labelmag = True              

elif myargs.mode == 'OBS':
    ds = densSatObs
    dv = densVelObs
    labelmag = True              

elif myargs.mode == 'FAINT':
    ds = densSatAll - densSatBloom
    dv = densVelAll - densVelBloom    
    labelmag = True              

elif myargs.mode == 'ALL':
    ds = densSatAll
    dv = densVelAll

elif myargs.mode == 'EFFECT':
    ds = myTel.trailf * densSatObs + (1.-myTel.trailf)* densSatBloom
    dv = myTel.trailf * densVelObs + (1.-myTel.trailf)* densVelBloom
    labelmag = True
    
else:
    print("valid for -M: BRIGHT OBS FAINT ALL EFFECT")
    exit(1)
    
#==============================================================================
# colormap
cmap = "magma"

# set limits for the colormap

lvmin = -2.5
lvmax = np.log10(30)##< MAXIMUM STANDARD
#lvmax = np.log10(5000)#


print ("telinslabel",myargs.code)
if myargs.code == "TrailDens":
    ldens = np.log10( dv )
    densl = "Number of trails./deg/sec."
    #print("Tdensity")

elif myargs.code == "SatDens":
    ldens = np.log10( ds )
    densl = "Number of sat./sq.deg."
    #print("Sdensity")

elif myargs.code == "skyMag":
    ldens =  2.5* np.log10( fluxSatTotal )
    skybrightmag = False
    lvmin = -30.0
    lvmax = -26.25  ## np.amax(ldens) + 0.5
    lvmax = -24.25  ## np.amax(ldens) + 0.5
    labelmag = True

elif myargs.code == "skyFrac":  # fraction of the sky surfbrightness
    skymag0 = 21.78 + 5 # dark sky, mag/arcsec2, +5 for %
    skymag = skymag0 
    #skymag = 21.00 +5  #  mag/arcsec2    -15deg tw
    #skymag = 19.80 +5  #  mag/arcsec2    -12deg tw
    #skymag = 9.0 +5   #  mag/arcsec2    -9deg tw
    ldens =  2.5* np.log10( fluxSatTotal ) + skymag
    print("pseudomag", np.amax(ldens), np.amin(ldens))
    lvmin = -30.0  + skymag0
    lvmax = -26.25 + skymag0
    labelmag = True


else: # other specific (including EFFECT)
    #print("density: other, specific instrument")
    dens    =  ds* myTel.fovl*myTel.fovw + dv * myTel.fovl * myTel.expt
             # trailf already accounted for in ds and dv
     
    ldens   = np.nan_to_num(np.log10(dens ), neginf=np.log10( np.min( dens[ dens > 0] )))
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


    
    if myargs.mode == 'EFFECT':
        #print("Effect")
        densl = "Fraction lost"
        cmap = gyrd    
        lvmin = -3.5  # log limits for the colour scale
        lvmax = 0.2  # 2.2
        labelmag = True
    else:
        densl = "Number of trails per exp."


#-- plot

if not myargs.plotflag:
    myargs.labelplotflag  = False
    myargs.scalebarflag = False
    myargs.shadeflag = False
else:
    fig = plt.figure(figsize=(8,8))
    ax =  fig.subplots(1,1,subplot_kw={'projection': 'polar'}) 
    _ = cp.initPolPlot(ax)
    
    clab = 'k'
    ccon = 'k'
    tickformat = "{:.1f}".format


    #ldens stat
    lvmax = np.max( ldens ) # np.percentile( ldens, 95.)
    lvmin = np.min(ldens) # np.percentile( ldens, 5.)

    #print(f'log min max {lvmin} {lvmax}')
    
    ldens[ ldens > lvmax  ] = lvmax
    ldens[ ldens < lvmin  ] = lvmin

    cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1], ldens , 
                      levels=np.linspace(lvmin,lvmax,100),  # NUMBER OF LEVELS
                      vmin=lvmin, vmax=lvmax ,
                      extend='both',
                      cmap=cmap)
    cfd.cmap.set_under('k') # below minimum -> black


#----------------------------------------------------------------------
#labels


if myargs.labelplotflag:

    #Sun
    azs,els = ca.radec2azel(sunAlpha, sunDelta, myTel.lat)
    plt.text(np.radians(azs), 93.,"$\odot$", va="center", ha='center')
    

    #top left
    x = -1.
    y = 1.2
    dy = 0.08
    
    cp.azlab(ax,x,y,'Observatory: {} Lat.: {:.1f}$^o$'.format(myTel.telescope, myTel.lat))
    y -= dy
    
    
    if myargs.code != "SatDens" and myargs.code != "TrailDens" and myargs.code != "skyMag":
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


        cp.azlab(ax,x,y,'Exp.t: {:.0f}s'.format(myTel.expt))
        y -= dy
        

            
    # bottom left
    x= -1.
    y= -1.08
    cp.azlab(ax,x,y,'$\odot$ Sun:',14)

    y -= dy
    loct = (sunAlpha/15.+12.)%24
    
    loch = int(loct)
    locm = int( (loct-loch)*60.)
    cp.azlab(ax,x,y,f'Loc.time: {loch:02d}:{locm:02d}')
    y -= dy
    cp.azlab(ax,x,y,f'$\delta: {sunDelta:.2f}^o$, Elev: {sunElev:.2f}$^o$')
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
        lab = "Satellite magnitudes: V$_{1000km}=$" +\
            "{:3.1f}".format(mag550 -5.*np.log10(550./1000.) )
        cp.azlab(ax,x,y,lab)
        y -= dy

    if mageffmin < -1000. or myargs.code == "skyMag":
        lab = "  V$_{sat}$"+" in [{:.1f}, {:.1f}]".format(magmax,magmin)
    else:
        lab = "  V$_{sat}$"+" in [{:.1f}, {:.1f}]".format(magmax,magmin)+\
              "  V$_{eff}$"+" in [{:.1f}, {:.1f}]".format(mageffmax,mageffmin)


    cp.azlab(ax,x,y,lab)
    y -= dy

    if myargs.mode == "BRIGHT":
        cp.azlab(ax,x,y,"Selection: mag < {:.0f}".format(myTel.magbloom))
    elif myargs.mode == "OBS":
        wlab = "Selection: mag$_{eff}$ < "+"{:.1f}".format(myTel.maglim)
        cp.azlab(ax,x,y,wlab)
    elif myargs.mode == "FAINT":
        cp.azlab(ax,x,y,"Selection: mag > {:.0f}".format(myTel.magbloom))
    elif myargs.mode == "EFFECT":
        cp.azlab(ax,x,y,"Selection: all satellites, scaled for effect")
        wlab = "Detected: V$_{eff}$ < "+"{:.1f} ".format(myTel.maglim)
        wlab += "   Bleeding: V$_{eff}$ < "+"{:.1f}".format(myTel.magbloom)
        y -= dy
        cp.azlab(ax,x,y,wlab)


    

if myargs.almucantar and myargs.plotflag:
    # sat count on almucantars
    for we, wi in zip(elLim, elCount):
        cp.azlab(ax,-0,(90.-we-5)/90.,'{:.0f} sat.>{:.0f}$^o$:'.format(wi,we),9,0.5*(1.-we/100.))

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


