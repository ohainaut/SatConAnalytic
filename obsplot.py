#!/usr/bin/env python3
# SatConAnalytic
# main plot

# obsplot.py
#   -h for options


def helphelp():
    out = """
    obsplot -d deltasun -a alphasun -e elevsun -l latitude*
    -C constellation: YESTURDAY TODAY SL2OW2 SLOW2 OW2r OW2 SLOWr(def)
    -T FlyEye WFI VST EFOSC FORSimg LSST HAWKI MICADO
       FORSspec UVES1h 4MOST ESPRESSO
       VISIR 
       skyMag TrailDens SatDens(def)
          These pre-load default 
                extptime, FoVl, FoVw, maglim, magbloom, trailf, 
                telname, instrName, latitude
    -t exptime* -f Fovl* -m maglim* -k trailf* -r resolution*
              maglim: those with eff.mag fainter are not observable
                      those brigher create a trail
              trailf: fraction of the exposure that is destroyed by trail. 1=full
              magbloom: for those brighter, trailf=1=full
    -s telName* -i instrName*
                (*: overwrites presets)
    -S do not shade low elev
    -L do not label the plot
    -B no scalebar
    -M BRIGHT FAINT ALL OBS EFFECT
         All sats, restrict to bright/faint (wrt magbloom), convert into effect on obs.
    """
    print(out)

import numpy as np
import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt

import sys
import os
sys.path.append(os.path.dirname(__file__)+"/../")

from matplotlib import ticker
from SatConAnalytic.conanplot import gyrd

# import SatConAnalytic routines
import SatConAnalytic.conan as ca
import SatConAnalytic.conanplot as cp

from SatConAnalytic.constants import mag550


#------Arguments

outpath = "/home/ohainaut/public_html/outsideWorld/"
##outpath = "./"

#Default and overwrite
telinslabel= "0"
lat = -24.6
overlat = 0

deltas=0.
elevs=-18
alphas=np.nan

magbloom = -99.
overmagbloom = 0

maglim = 99.
overmaglim = 0

fovl = 1.
fovw = -99.
overfovl = 0
overfovw = 0

resol       = 1./3600. # deg
overresol = 0

exptime = 1
overexptime = 0

trailf = 1
overtrailf = -1

telescope = " "
overtelescope = 0

instrument = " "
overinstrument = 0


constellationsll = "SLOWGWAK"
magcutflag = 'ALL'

almucantar = True # label number of satellites on 0,10,20,30,60deg elevation almucantars
shadeflag = False  # shade airmass above 30 and 20deg
outputformat = '.png'
labelplot= True # write info
scalebar = True # 

# internal
labelmag = False
plotflag = True
step = 0.75 #deg >~1. Smaller values take forever



helpme = False
if len(sys.argv) == 1:
    helpme = True
for i in np.arange(0,len(sys.argv)):
    if sys.argv[i] == '-d': #deltasun
        i += 1
        deltas = float(sys.argv[i])
    if sys.argv[i] == '-a': #alphasun
        i += 1
        alphas = float(sys.argv[i])
    if sys.argv[i] == '-e': #elevation sun
        i += 1
        elevs = float(sys.argv[i])
    if sys.argv[i] == '-l': #latitude
        i += 1
        overlat = float(sys.argv[i])
    if sys.argv[i] == '-r': #resolution
        i += 1
        overresol = float(sys.argv[i])/3600.
    if sys.argv[i] == '-t': #expt
        i += 1
        overexptime = float(sys.argv[i])
    if sys.argv[i] == '-f': #fieldofview
        i += 1
        overfovl = float(sys.argv[i])
    if sys.argv[i] == '-fw': #fieldofview
        i += 1
        overfovw = float(sys.argv[i])
    if sys.argv[i] == '-m': #maglim
        i += 1
        overmaglim = float(sys.argv[i])
    if sys.argv[i] == '-k': #trailf
        i += 1
        overtrailf = float(sys.argv[i])
    if sys.argv[i] == '-s': #telescope
        i += 1
        overtelescope = sys.argv[i]
    if sys.argv[i] == '-i': #instrument
        i += 1
        overinstrument = sys.argv[i]
    if sys.argv[i] == '-C': #constellation ID
        i += 1
        constellationsll = sys.argv[i]
    if sys.argv[i] == '-0': #no plot
        plotflag = False
    if sys.argv[i] == '-S': #do not shade low elevations
        shadeflag = False
    if sys.argv[i] == '-B': #no scalebar
        scalebar = False
    if sys.argv[i] == '-L': #do not plot any label
        labelplot = False
    if sys.argv[i] == '-M': # show only those that are brigher than Mcut
        i += 1
        magcutflag = sys.argv[i]
    if sys.argv[i] == '-T': #telescope/instrument
        i += 1
        telinslabel = sys.argv[i]
    if sys.argv[i] == '-pdf': # output format
        outputformat = ".pdf"
    if sys.argv[i] == '-h' or sys.argv[i] == "--help": #help
        i += 1
        helpme = True         
    i += 1

if helpme:
    helphelp()
    exit(0)

#-----presets
    
# maglim: limiting magnitude
#   maglimoff = 1 = -2.5log(2/5) = offset from 5sig (from table) to 2sig
#   trailf: fraction of the exposure that is destroyed. 1=full
# magbloom: for those brighter, trailf=1, exposure fully destroyed.
#   magbloomoff: -2.5 or -5: offset from saturation to blooming

if telinslabel != "0":
    telescope, lat, instrument, expt, fovl, fovw, resol, maglim, magbloom, trailf =    cp.findPreset(telinslabel)

    
if overexptime:
    expt = overexptime
if overfovl:
    fovl = overfovl
if overfovw:
    fovw = overfovw
if fovw < 0:
    fovw = fovl *1.
    
if overmagbloom:
    magbloom = overmagbloom
if overmaglim:
    maglim = overmaglim
if overresol:
    resol = overresol
if overtrailf >= 0.:
    trailf = overtrailf
if overlat:
    lat = overlat
if overtelescope:
    telescope = overtelescope
if overinstrument:
    instrument = overinstrument

# we have 2 from elevs, alphas, deltas; get the 3rd one.
if np.isnan(alphas):
    alphas = ca.elev2ra(elevs,deltas,lat) # get sun hourangle for twilight
else:
    elevs = ca.radec2elev(alphas,deltas,lat)

constellationsl, constellations = cp.findConstellations(constellationsll)
consLab, consNum, consPla, consNPl, consInc, consAlt = ca.loadConstellations(constellations)
nCons = len(consLab)

if 1:#plotflag:
    print("Telescope:", telinslabel)
    print("Tel. latitude:",lat)
    print("Resolution (pixel): {:.2g}arcsec".format(resol*3600.))
    print("Trail fill fraction: {:.3g}".format(trailf))
    print("Sun coord.: HA: {:.1f}, Dec: {:.1f}, Elev: {:.1f} [deg]".format(alphas, deltas, elevs))
    print("Constellations:",constellationsl)
    print("Number of sub.constellations:",nCons)
    print("Total num of sat:", sum(consNum))
    #print('Number of sat. per const.:',consNum)
    print('Magnitude cut:',magcutflag)

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
    #print("Const: ",icons,    consInc[icons], consAlt[icons], consNum[icons])
                           
    # constellation
    densSi, veli, magi =  ca.modelOneConstMag(AzEl,lat, alphas,deltas,
                            consInc[icons], consAlt[icons], consNum[icons])


    magmax = max(magmax,np.amax(magi))
    magmin = min(magmin,np.amin(magi))
    
    # all sat:
    densSatAll += densSi
    densVelAll += densSi * veli

    # effective magnitude
    mageffi = magi  - 2.5*np.log10(resol/veli/expt)
    mageffmax = max(mageffmax,np.amax(mageffi))
    mageffmin = min(mageffmin,np.amin(mageffi))

    # only observable ones
    densSobsi = np.copy(densSi)
    densSobsi[ mageffi > maglim] = 0.
    densSatObs += densSobsi
    densVelObs += densSobsi * veli
    
    # only super bright bloomers
    densSbloomi = np.copy(densSi)
    densSbloomi[ mageffi > magbloom ] = 0.
    densSatBloom += densSbloomi
    densVelBloom += densSbloomi * veli

    #flux    
    fluxi = 10.**(-0.4*magi) /3600.**2 * densSi
    fluxSatTotal += fluxi

    icons += 1
    
# here, we have the following defined:
# densSatAll: dens in Nsat/sq.deg
# densVelAll: density if trails in Ntrail/deg/sec
# densSatObs: same, only for sat with mag < maglim
# densVelObs:
# densSatBloom: same, only for sat with mag < magbloom
# densVelBloom
# fluxSatTotal: flux density for mag/sq.arcsec, for m550 at 550km


if plotflag:
    print( "Satellite magnitudes in [{:.2f},{:.2f}]".format(magmax,magmin))
    print( "Satellite eff. mag.  in [{:.2f},{:.2f}]".format(mageffmax,mageffmin))

outfileroot = "{}_{}_{}_{:02d}_{:02d}".format(telinslabel,constellationsll,magcutflag,int(deltas),int(-elevs))
outfileroot = "{}_{}_{}_{:02d}_{:02d}".format(telinslabel,constellationsll,magcutflag,int(lat),int(-elevs))

# sat count for almucantars
elLim = [60.,30.,20., 10.,0.]
elCount = ca.integrateSat(elLim,AzEl,densSatAll)

outfile = open(outpath+outfileroot+'.txt','w+')  # store the hist.
outstring =  ('{:6.3f} {:5.1f} '+' {:4.0f}'*16).format(
    alphas,
    elevs,
    0.,0.,0.,0.,
    elCount[1],    
    elCount[2],    
    elCount[3],    
    elCount[4],
    0.,0.,0.,0.,
    0.,0.,0.,0.)
#will be written to file later


#######################################################################

if magcutflag == 'BRIGHT':
    ds = densSatBloom
    dv = densVelBloom
    labelmag = True              

elif magcutflag == 'OBS':
    ds = densSatObs
    dv = densVelObs
    labelmag = True              

elif magcutflag == 'FAINT':
    ds = densSatAll - densSatBloom
    dv = densVelAll - densVelBloom    
    labelmag = True              

elif magcutflag == 'ALL':
    ds = densSatAll
    dv = densVelAll

elif magcutflag == 'EFFECT':
    ds = trailf * densSatObs + (1.-trailf)* densSatBloom
    dv = trailf * densVelObs + (1.-trailf)* densVelBloom
    labelmag = True
    

else:
    print("valid for -M: BRIGHT OBS FAINT ALL EFFECT")
    exit(1)
    

#cmap = copy.copy(mpl.cm.get_cmap("magma"))
cmap = "magma"

lvmin = -2.5
lvmax = np.log10(20)
lvmin = -2.5    #xoxo for satDens
lvmax = np.log10(5) # was 1

print ("telinslabel",telinslabel)
if telinslabel == "TrailDens":
    ldens = np.log10( dv )
    densl = "Number of trails./deg/sec."
    #print("Tdensity")

elif telinslabel == "SatDens":
    ldens = np.log10( ds )
    densl = "Number of sat./sq.deg."
    #print("Sdensity")

elif telinslabel == "skyMag":
    ldens =  2.5* np.log10( fluxSatTotal )
    skybrightmag = False
    lvmin = -30.0
    lvmax = -26.25  ## np.amax(ldens) + 0.5
    labelmag = True

elif telinslabel == "skyFrac":  # fraction of the sky surfbrightness
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
    dens    =  ds* fovl*fovw + dv * fovl * expt
             # trailf already accounted for in ds and dv
     
    ldens   = np.log10(dens)
    densobs = densSatObs * fovl*fovw +  densVelObs* fovl *expt
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

    
    if magcutflag == 'EFFECT':
        #print("Effect")
        densl = "Fraction lost"
        cmap = gyrd    
        lvmin = -3.5  # log limits for the colour scale
        lvmax = 2.2  #
        labelmag = True
    else:
        densl = "Number of trails per exp."


#-- plot

if not plotflag:
    labelplot  = False
    scalebar = False
    shadeflag = False
else:
    #D    print('PLOT ')
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
if scalebar:
    cbar = fig.colorbar(cfd)
    if telinslabel == "skyMag":
        barmag = np.arange(-30,-25,.5)
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
                barlum = 12e10 * 10**(0.4*barmag)/2.20   ## CONVERSION mag->mucd
                cbar.set_ticks(barmag)
                cbar.set_ticklabels( [ "{:.1g}".format(x) for x in barlum])
                densl = "Surface brightness [ % of sky]"

            
    elif telinslabel == "skyFrac":
        barmag = np.arange(-30,-25,.5) + skymag0 # skymag comes with +5mag for %
        barlum = 10**(0.4*barmag) 
        cbar.set_ticks(barmag)
        cbar.set_ticklabels( [ "{:.1g}".format(x) for x in barlum])
        densl = "Surface brightness [ % of sky] ($m_{sky} = $"+"{:.2f}".format(skymag-5)+")]"

    else:
        cbar.set_ticks(np.log10(np.array([1e-5,2e-5,5e-5,1e-4,2e-4,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,5.,10.,20.,50.])))
        cbar.set_ticklabels( [ "{:.2g}".format(x) for x in 10.**cbar.get_ticks()])


    cbar.set_label(densl)


#----------------------------------------------------------------------
# airmass shade

if shadeflag:
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
cp.drawHADec(lat)


#----------------------------------------------------------
#labels



if labelplot:

    #Sun
    azs,els = ca.radec2azel(alphas, deltas, lat)
    plt.text(np.radians(azs), 93.,"$\odot$", va="center", ha='center')
    

    #top left
    x = -1.
    y = 1.2
    dy = 0.08
    
    cp.azlab(ax,x,y,'Observatory: {} Lat.: {:.1f}$^o$'.format(telescope, lat))
    y -= dy
    
    
    if telinslabel != "SatDens" and telinslabel != "TrailDens" and telinslabel != "skyMag":
        cp.azlab(ax,x,y,'Instrument: {}'.format(instrument))
        y -= dy

        if fovl < 1./60.:
            fovll = '{:.2f}\"'.format(fovl*3600.)
        elif fovl < 1./6.:
            fovll = '{:.2f}\''.format(fovl*60.)
        else:
            fovll = '{:.2f}$^o$'.format(fovl)

        if fovw < 1./60.:
            fovlw = '{:.2f}\"'.format(fovw*3600.)
        elif fovw < 1./6.:
            fovlw = '{:.2f}\''.format(fovw*60.)
        else:
            fovlw = '{:.2f}$^o$'.format(fovw)

        cp.azlab(ax,x,y,'Fov: '+fovll+'x'+fovlw)
        y -= dy


        cp.azlab(ax,x,y,'Exp.t: {:.0f}s'.format(expt))
        y -= dy
        

            
    # bottom left
    x= -1.
    y= -1.08
    cp.azlab(ax,x,y,'$\odot$ Sun:',14)

    y -= dy
    loct = (alphas/15.+12.)%24
    
    loch = int(loct)
    locm = int( (loct-loch)*60.)
    cp.azlab(ax,x,y,'Loc.time: {:02d}:{:02d}'.format(loch,locm))
    y -= dy
    #print('$\delta: {:.2f}^o$, Elev: {:.2f}$^o$'.format(deltas, elevs))
    cp.azlab(ax,x,y,'$\delta: {:.2f}^o$, Elev: {:.2f}$^o$'.format(deltas, elevs))
    y -= dy
#    elevs = ca.radec2elev(alphas,deltas,lat)


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

    if mageffmin < -1000. or telinslabel == "skyMag":
        lab = "  V$_{sat}$"+" in [{:.1f}, {:.1f}]".format(magmax,magmin)
    else:
        lab = "  V$_{sat}$"+" in [{:.1f}, {:.1f}]".format(magmax,magmin)+\
              "  V$_{eff}$"+" in [{:.1f}, {:.1f}]".format(mageffmax,mageffmin)


    cp.azlab(ax,x,y,lab)
    y -= dy

    if magcutflag == "BRIGHT":
        cp.azlab(ax,x,y,"Selection: mag < {:.0f}".format(magbloom))
    elif magcutflag == "OBS":
        wlab = "Selection: mag$_{eff}$ < "+"{:.1f}".format(maglim)
        cp.azlab(ax,x,y,wlab)
    elif magcutflag == "FAINT":
        cp.azlab(ax,x,y,"Selection: mag > {:.0f}".format(magbloom))
    elif magcutflag == "EFFECT":
        cp.azlab(ax,x,y,"Selection: all satellites, scaled for effect")
        wlab = "Detected: V$_{eff}$ < "+"{:.1f} ".format(maglim)
        wlab += "   Bleeding: V$_{eff}$ < "+"{:.1f}".format(magbloom)
        y -= dy
        cp.azlab(ax,x,y,wlab)


    

if almucantar and plotflag:
    # sat count on almucantars
    for we, wi in zip(elLim, elCount):
        cp.azlab(ax,-0,(90.-we-5)/90.,'{:.0f} sat.>{:.0f}$^o$:'.format(wi,we),9,0.5*(1.-we/100.))

print('finishing...')

outstring += '\n'
outfile.write(outstring)
outfile.close()
        



# save plot
if plotflag:
    fig.tight_layout()
    plt.savefig(outpath+'w.png')
    filename = outpath+outfileroot+outputformat
    print(filename)

    plt.savefig(filename)
#--
print("output in ",outfileroot)


