#!/usr/bin/env python3
# ConAn
# Constellation Analytic simulations
#
# Generate main plot: density of satellite over the map of the sky

# obsplot.py
#   -h for options

import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt
import numpy as np

import argparse
import sys
import os
sys.path.append(os.path.dirname(__file__)+"/../")


from matplotlib import ticker
from ConAn.conanplot import gyrd


# import ConAn routines
import ConAn.conan as ca
import ConAn.conanplot as cp
from ConAn.constants import mag550

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
parser.add_argument('-e','--elevSun', default=-24.,
                    help="Sun: Elevation of the Sun BELOW the horizon. Should probably be >0 in most cases [deg]")
parser.add_argument('-C','--constellations', default='SLOWGWAK',
                    help="ID of the constellation group; list for a list")
parser.add_argument('-l','--latitude', default=-24.6,
                    help="Observatory: Latitude of the observatory [deg]")
parser.add_argument('-r','--resolution', default=1./3600., 
                    help="Observatory: Resolution of the instrument [deg]")
parser.add_argument('-t','--expTime', default=10.,
                    help="Observatory: Exposure time [sec]")
parser.add_argument('-f','--FoV', default=1.,
                    help="Observatory: Field of view of the instrument. Length or diametre [arcsec]")
parser.add_argument('-w','--FoVw', 
                    help="Observatory: Field of view of the instrument. Width. Equal to Length if omitted [arcsec]")
parser.add_argument('-m','--magLim', default=99.,
                    help="Observatory: Detection limit magnitude of the instrument [5sigma Mag during expTime]")
parser.add_argument('--magBloom', default=-99.,
                    help="Observatory: Magnitude over which the instrument saturates [Mag for expTime]. Default: -99 (no blooming)")
parser.add_argument('-k','--trailFill', default=1.,
                    help="Observatory: Trail filling fraction (width of the trail as fraction of FoV)")
parser.add_argument('-s','--telescope', default=" ",
                    help="Observatory: Name of the telescope")
parser.add_argument('-i','--instrument', default=" ",
                    help="Observatory: Name of the instrument")
parser.add_argument('-T','--system', 
                    help='''Observatory: Predefined telescope/instrument with
                         extptime, FoVl, FoVw, maglim, magbloom, trailf,  
                         telescope, instrument, resolution, latitude. 
                         Use individual options to overwrite presets.
                         SPECIAL CODES: SatDens for satellite density; 
                         TrailDens for trail density;
                         skyMag for sky surface brightness [mag];
                         skyFrac for sky surface brightness as a fraction.
''')
parser.add_argument('-M','--magcut', default="ALL",
                    help="Plot: ALL (Default), OBS, BRIGHT, FAINT, or EFFECT")
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
parser.add_argument('--pdf', action='store_true',
                    help="Plot: output file in pdf (default is png)")
args = parser.parse_args()




if args.system is not None:  #-----presets
    args.telinslabel = args.system
    args.telescope, args.lat, args.instrument, args.expt, args.fovl, args.fovw, args.resol, args.maglim, args.magbloom, args.trailf = cp.findPreset(args.telinslabel)
else:
    args.telinslabel = " "


    # observatory and instrument
    args.lat   = float(args.latitude)
    args.resol = float(args.resolution)

    args.expt = float(args.expTime)
    args.fovl = float(args.FoV)

    if args.FoVw is None:
        args.fovw = args.fovl*1.
    else:
        args.fovw = float(args.FoVw)

    args.maglim = float(args.magLim)
    args.magbloom = float(args.magBloom)
    args.trailf = float(args.trailFill)


# sun
args.deltas = float(args.deltaSun)
args.elevs = float(args.elevSun)   
if args.alphaSun is None:
    args.alphas = ca.elev2ra(args.elevs,args.deltas,args.lat) # get sun hourangle for twilight
else:
    args.alphas = float(args.alphaSun)
    args.elevs = ca.radec2elev(args.alphas,args.deltas,args.lat)



# flags
args.plotflag     = args.noplot
args.shadeflag    = args.noshade
args.scalebarflag = args.noscalebar
args.labelplotflag = args.nolabel
args.almucantar   = args.noalmuc
if args.pdf :
    args.outputformat = ".pdf"
else:
    args.outputformat = ".png"
    



# expand the constellation id into a list of real constellations
constellationsl, constellations = cp.findConstellations(args.constellations)
consLab, consNum, consPla, consNPl, consInc, consAlt = ca.loadConstellations(constellations)
nCons = len(consLab)



if 1:
    print("Telescope:          ", args.telinslabel)
    print("Tel. latitude:      ", args.lat)
    print("Tel. FoV:           ", args.fovl, args.fovw, "deg")
    print("Resolution (pixel):  {:.2g}arcsec".format(args.resol*3600.))
    print("Trail fill fraction: {:.3g}".format(args.trailf))
    print("Sun coord.: HA: {:.1f}, Dec: {:.1f}, Elev: {:.1f} [deg]".format(args.alphas, args.deltas, args.elevs))
    print("Constellations:     ",constellationsl)
    print("Number of sub.constellations:",nCons)
    print("Total num of sat:   ", sum(consNum))
    #print('Number of sat. per const.:',consNum)
    print('Magnitude cut:      ',args.magcut)
    print('plotflag:           ',args.plotflag)
    print('almucantar:         ',args.almucantar)



    
######################################################################
#---  COMPUTE CONSTELLATIONS

#fill ElAz:
AzEl = ca.fillAzEl(step)

densSatAll = np.zeros_like(AzEl[0])
densVelAll = np.zeros_like(AzEl[0])
densSatObs = np.zeros_like(AzEl[0])
densVelObs = np.zeros_like(AzEl[0])
densSatBloom = np.zeros_like(AzEl[0])
densVelBloom = np.zeros_like(AzEl[0])
fluxSatTotal = np.zeros_like(AzEl[0])

magmax = -99.
magmin = 99.
mageffmax = -99
mageffmin = 99.

icons = 0

while icons < len(consInc):
    # constellation
    densSi, veli, magi =  ca.modelOneConstMag(AzEl,args.lat, args.alphas,args.deltas,
                            consInc[icons], consAlt[icons], consNum[icons])

    magmax = max(magmax,np.amax(magi))
    magmin = min(magmin,np.amin(magi))
    
    # all sat:
    densSatAll += densSi
    densVelAll += densSi * veli

    # effective magnitude
    mageffi = magi  - 2.5*np.log10(args.resol/veli/args.expt)
    mageffmax = max(mageffmax,np.amax(mageffi))
    mageffmin = min(mageffmin,np.amin(mageffi))

    # only observable ones
    densSobsi = np.copy(densSi)
    densSobsi[ mageffi > args.maglim] = 0.
    densSatObs += densSobsi
    densVelObs += densSobsi * veli
    
    # only super bright bloomers
    densSbloomi = np.copy(densSi)
    densSbloomi[ mageffi > args.magbloom ] = 0.
    densSatBloom += densSbloomi
    densVelBloom += densSbloomi * veli

    #flux    
    fluxi = 10.**(-0.4*magi) /3600.**2 * densSi
    fluxSatTotal += fluxi

    icons += 1
    
# here, we have the following arrays defined:
#   densSatAll: dens in Nsat/sq.deg
#   densVelAll: density if trails in Ntrail/deg/sec
#   densSatObs: same, only for sat with mag < args.maglim
#   densVelObs:
#   densSatBloom: same, only for sat with mag < args.magbloom
#   densVelBloom
#   fluxSatTotal: flux density for mag/sq.arcsec, for m550 at 550km

if args.plotflag:
    print( "Satellite magnitudes in [{:.2f},{:.2f}]".format(magmax,magmin))
    print( "Satellite eff. mag.  in [{:.2f},{:.2f}]".format(mageffmax,mageffmin))
    print('Values at zenith:', AzEl[:,-1,-1])
    print('Sat density [n/sq.dg]', round(densSatAll[-1,-1],4))
    print('Sat vel (last constellation) [deg/s]', round(veli[-1,-1],4))
    print('Diffuse mag',  -2.5*np.log10(fluxSatTotal[-1,-1]))

# output file
outfileroot = "{}_{}_{}_{:02d}_{:02d}".format(args.telinslabel,args.constellations,args.magcut,int(args.lat),int(-args.elevs))

# sat count for almucantars
elLim = [60.,30.,20., 10.,0.]
elCount = ca.integrateSat(elLim,AzEl,densSatAll)

outfile = open(outpath+outfileroot+'.txt','w+')  # store the hist.
outstring =  ('{:6.3f} {:5.1f} '+' {:4.0f}'*16).format(
    args.alphas,
    args.elevs,
    0.,0.,0.,0.,
    elCount[1],    
    elCount[2],    
    elCount[3],    
    elCount[4],
    0.,0.,0.,0.,
    0.,0.,0.,0.)
#will be written to file later


#==============================================================================

if args.magcut == 'BRIGHT':
    ds = densSatBloom
    dv = densVelBloom
    labelmag = True              

elif args.magcut == 'OBS':
    ds = densSatObs
    dv = densVelObs
    labelmag = True              

elif args.magcut == 'FAINT':
    ds = densSatAll - densSatBloom
    dv = densVelAll - densVelBloom    
    labelmag = True              

elif args.magcut == 'ALL':
    ds = densSatAll
    dv = densVelAll

elif args.magcut == 'EFFECT':
    ds = args.trailf * densSatObs + (1.-args.trailf)* densSatBloom
    dv = args.trailf * densVelObs + (1.-args.trailf)* densVelBloom
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


print ("telinslabel",args.telinslabel)
if args.telinslabel == "TrailDens":
    ldens = np.log10( dv )
    densl = "Number of trails./deg/sec."
    #print("Tdensity")

elif args.telinslabel == "SatDens":
    ldens = np.log10( ds )
    densl = "Number of sat./sq.deg."
    #print("Sdensity")

elif args.telinslabel == "skyMag":
    ldens =  2.5* np.log10( fluxSatTotal )
    skybrightmag = False
    lvmin = -30.0
    lvmax = -26.25  ## np.amax(ldens) + 0.5
    lvmax = -24.25  ## np.amax(ldens) + 0.5
    labelmag = True

elif args.telinslabel == "skyFrac":  # fraction of the sky surfbrightness
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
    dens    =  ds* args.fovl*args.fovw + dv * args.fovl * args.expt
             # trailf already accounted for in ds and dv
     
    ldens   = np.log10(dens)
    densobs = densSatObs * args.fovl*args.fovw +  densVelObs* args.fovl *args.expt
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
    print("Aver. effect on exposures: Loss {:.3g}/1.; Trails: {:.3g}/exp".format( EffTot, TrailTot))
    outstring += " {} {}".format(EffTot,TrailTot)

    
    if args.magcut == 'EFFECT':
        #print("Effect")
        densl = "Fraction lost"
        cmap = gyrd    
        lvmin = -3.5  # log limits for the colour scale
        lvmax = 0.2  # 2.2
        labelmag = True
    else:
        densl = "Number of trails per exp."


#-- plot

if not args.plotflag:
    args.labelplotflag  = False
    args.scalebarflag = False
    args.shadeflag = False
else:
    fig = plt.figure(figsize=(8,8))
    ax =  fig.subplots(1,1,subplot_kw={'projection': 'polar'}) 
    _ = cp.initPolPlot(ax)
    
    clab = 'k'
    ccon = 'k'
    tickformat = "{:.1f}".format

    ldens[ ldens > lvmax  ] = lvmax
    ldens[ ldens < lvmin  ] = lvmin

    cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1], ldens , 
                      levels=np.arange(lvmin,lvmax,0.01),
                      vmin=lvmin, vmax=lvmax ,
                      extend='both',
                      cmap=cmap)
    cfd.cmap.set_under('k') # below minimum -> black


#----------------------------------------------------------------------
#Scalebar
if args.scalebarflag:
    cbar = fig.colorbar(cfd)
    if args.telinslabel == "skyMag":
        # (the plot contains -mag) 
        barmag = np.arange(-30,-23,.5)
        if skybrightmag:
            cbar.set_ticklabels( [ "{:.1f}".format(x) for x in -cbar.get_ticks()])
            densl = "Surface brightness [mag/sq.arcsec]"

        else:
            if 0:
                # mucd/m2
                barlum = 12e10 * 10**(0.4*barmag)   ## CONVERSION mag->mucd
                cbar.set_ticks(barmag)
                cbar.set_ticklabels( [ "{:.1g}".format(x) for x in barlum])
                densl = "Surface brightness [$\mu$cd/m$^2$]"

            else:
                # dark sky:
                barlum = 12e10 * 10**(0.4*barmag)/2.20   ## sky=220; pc= 1/100
                cbar.set_ticks(barmag)
                cbar.set_ticklabels( [ "{:.1g}".format(x) for x in barlum])
                densl = "Surface brightness [ % of sky]"

            
    elif args.telinslabel == "skyFrac":
        barmag = np.arange(-30,-25,.5) + skymag0 # skymag comes with +5mag for %
        barlum = 10**(0.4*barmag) 
        cbar.set_ticks(barmag)
        cbar.set_ticklabels( [ "{:.1g}".format(x) for x in barlum])
        densl = "Surface brightness [ % of sky] ($m_{sky} = $"+"{:.2f}".format(skymag-5)+")]"

    else:
        cbar.set_ticks(np.log10(np.array([0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,5.,10.,20.,50.,100.,200.,500., 1000., 2000.,5000.,])))
#        cbar.set_ticks(np.log10(np.array([0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,5.,10.,20.,50.])))
        cbar.set_ticklabels( [ "{:5.2g}".format(x) for x in 10.**cbar.get_ticks()])


    cbar.set_label(densl)


#----------------------------------------------------------------------
# airmass shade

if args.shadeflag:
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
# RA Dec lines
cp.drawHADec(args.lat)


#----------------------------------------------------------
#labels



if args.labelplotflag:

    #Sun
    azs,els = ca.radec2azel(args.alphas, args.deltas, args.lat)
    plt.text(np.radians(azs), 93.,"$\odot$", va="center", ha='center')
    

    #top left
    x = -1.
    y = 1.2
    dy = 0.08
    
    cp.azlab(ax,x,y,'Observatory: {} Lat.: {:.1f}$^o$'.format(args.telescope, args.lat))
    y -= dy
    
    
    if args.telinslabel != "SatDens" and args.telinslabel != "TrailDens" and args.telinslabel != "skyMag":
        cp.azlab(ax,x,y,'Instrument: {}'.format(args.instrument))
        y -= dy

        if args.fovl < 1./60.:
            fovll = '{:.2f}\"'.format(args.fovl*3600.)
        elif args.fovl < 1./6.:
            fovll = '{:.2f}\''.format(args.fovl*60.)
        else:
            fovll = '{:.2f}$^o$'.format(args.fovl)

        if args.fovw < 1./60.:
            fovlw = '{:.2f}\"'.format(args.fovw*3600.)
        elif args.fovw < 1./6.:
            fovlw = '{:.2f}\''.format(args.fovw*60.)
        else:
            fovlw = '{:.2f}$^o$'.format(args.fovw)

        cp.azlab(ax,x,y,'Fov: '+fovll+'x'+fovlw)
        y -= dy


        cp.azlab(ax,x,y,'Exp.t: {:.0f}s'.format(args.expt))
        y -= dy
        

            
    # bottom left
    x= -1.
    y= -1.08
    cp.azlab(ax,x,y,'$\odot$ Sun:',14)

    y -= dy
    loct = (args.alphas/15.+12.)%24
    
    loch = int(loct)
    locm = int( (loct-loch)*60.)
    cp.azlab(ax,x,y,f'Loc.time: {loch:02d}:{locm:02d}')
    y -= dy
    cp.azlab(ax,x,y,f'$\delta: {args.deltas:.2f}^o$, Elev: {args.elevs:.2f}$^o$')
    y -= dy




    # top right
    x=1.
    y=1.2
    cp.azlab(ax,x,y,'Constellation:',14)
    y -= dy
    cp.azlab(ax,x,y,constellationsl)
    
    y -= dy
    cp.azlab(ax,x,y,"Total {:.0f} sat.".format(sum(consNum)))

    #bottom right
    x = 1.
    y = -1.08 -dy
    if 1:
        lab = "Satellite magnitudes: V$_{1000km}=$" +\
            "{:3.1f}".format(mag550 -5.*np.log10(550./1000.) )
        cp.azlab(ax,x,y,lab)
        y -= dy

    if mageffmin < -1000. or args.telinslabel == "skyMag":
        lab = "  V$_{sat}$"+" in [{:.1f}, {:.1f}]".format(magmax,magmin)
    else:
        lab = "  V$_{sat}$"+" in [{:.1f}, {:.1f}]".format(magmax,magmin)+\
              "  V$_{eff}$"+" in [{:.1f}, {:.1f}]".format(mageffmax,mageffmin)


    cp.azlab(ax,x,y,lab)
    y -= dy

    if args.magcut == "BRIGHT":
        cp.azlab(ax,x,y,"Selection: mag < {:.0f}".format(args.magbloom))
    elif args.magcut == "OBS":
        wlab = "Selection: mag$_{eff}$ < "+"{:.1f}".format(args.maglim)
        cp.azlab(ax,x,y,wlab)
    elif args.magcut == "FAINT":
        cp.azlab(ax,x,y,"Selection: mag > {:.0f}".format(args.magbloom))
    elif args.magcut == "EFFECT":
        cp.azlab(ax,x,y,"Selection: all satellites, scaled for effect")
        wlab = "Detected: V$_{eff}$ < "+"{:.1f} ".format(args.maglim)
        wlab += "   Bleeding: V$_{eff}$ < "+"{:.1f}".format(args.magbloom)
        y -= dy
        cp.azlab(ax,x,y,wlab)


    

if args.almucantar and args.plotflag:
    # sat count on almucantars
    for we, wi in zip(elLim, elCount):
        cp.azlab(ax,-0,(90.-we-5)/90.,'{:.0f} sat.>{:.0f}$^o$:'.format(wi,we),9,0.5*(1.-we/100.))

print('finishing...')

outstring += '\n'
outfile.write(outstring)
outfile.close()
        



# save plot
if args.plotflag:
    fig.tight_layout()
    plt.savefig(outpath+'w.png')
    filename = outpath+outfileroot+args.outputformat
    print(filename)

    plt.savefig(filename)
#--
print("output in ",outfileroot)


