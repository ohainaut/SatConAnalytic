#!/usr/bin/env python3
# ConAn
# Constellation Analytic simulations
#
# main plot:  plot the satellite density over a calendar for a given object

import numpy as np
import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote

import matplotlib.pyplot as plt
import sys
sys.path.append('/home/ohainaut/Documents/ESO/E2E/SatelliteConstellations/SIMULATIONS/')

from matplotlib import ticker
from conanplot import gyrd

# import ConAn routines
import conan as ca
import conanplot as cp

#print('ConAn:objplot -in')
def helphelp():
    out = """
    objplot 
    -r RA -d dec -n name
    -T FlyEye WFI VST EFOSC FORSimg LSST HAWKI MICADO
       FORSspec UVES1h 4MOST ESPRESSO
       VISIR 
       skyMag TrailDensity SatDensity(def)
          These pre-load default 
                extptime, FoVl, FoVw, maglim, magbloom, trailf, 
                telname, instrName, latitude
    -l latitude*
    -t exptime* -f Fovl* -m magbloom* -k trailf* 
              maglim: those with eff.mag fainter are not observable
                      those brigher create a trail
              trailf: fraction of the exposure that is destroyed by trail. 1=full
              magbloom: for those brighter, trailf=1=full
    -s telName* -i instrName*
                (*: overwrites presets)
    -M BRIGHT FAINT ALL OBS EFFECT --> HARDCODED
         All sats, restrict to bright/faint (wrt magbloom), convert into effect on obs.
    """
    print(out)

#outpath = "/home/ohainaut/public_html/outsideWorld/"
outpath = "./"

maglimoff = 1.75 # 1.75: 1/5sig ; 1.: 2/5sig
magbloomoff = -2.5
magcutflag = 'OBS'
magcutflag = 'EFFECT' ## 4MOST
fovw = -99.

constellationsll = 'SLOWGWAK'
elevlim = 20.
telinslabel = " "
telescope = " "
overlat = 0
lon = 0.

#default for overwriting values
overtime = 0
overfovl = 0
overfovw = 0
overmagbloom = 0
overmaglim = 0
overresol = 0
overtrailf = -1
overlat = 0
overtelescope = 0
overinstrument = 0
objlabel=" "

helpme = False
if len(sys.argv) == 1:
    helpme = True
for i in np.arange(0,len(sys.argv)):
    if sys.argv[i] == '-r': #ra
        i += 1
        rao = float(sys.argv[i])
    if sys.argv[i] == '-d': #delta
        i += 1
        deo = float(sys.argv[i])
    if sys.argv[i] == '-n': #obj name
        i += 1
        objlabel = sys.argv[i]
    if sys.argv[i] == '-l': #latitude
        i += 1
        overlat = float(sys.argv[i])
    if sys.argv[i] == '-T': #telescope/instrument
        i += 1
        telinslabel = sys.argv[i]
    if sys.argv[i] == '-t': #expt
        i += 1
        overtime = float(sys.argv[i])
    if sys.argv[i] == '-f': #fieldofview
        i += 1
        overfovl = float(sys.argv[i])
    if sys.argv[i] == '-fw': #fieldofview
        i += 1
        overfovw = float(sys.argv[i])
    if sys.argv[i] == '-r': #resolution
        i += 1
        overresol = float(sys.argv[i])/3600.
    if sys.argv[i] == '-m': #maglim
        i += 1
        overmaglim = float(sys.argv[i])


if helpme:
    helphelp()
    exit(0)

    
telescope, lat, instrument, expt, fovl, fovw, resol, maglim, magbloom, trailf =    cp.findPreset(telinslabel)

    
if overtime:
    expt = overtime
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


if fovw < 0:
    fovw = fovl *1.
    


#---


#- fill in the calendar
hstep = 12 # points per hour
times = np.zeros((365,24*hstep))
for i in range(365):
    for j in range(24*hstep):
        times[i,j] = 2459580.5 + i + j/24./hstep - (.5 + lon/360.) # to centre midnight
        

# Loc ST
lst = (280.46061837 + 360.98564736629*(times -2451545.0) + lon)%360

# sun coordinates
walphas, deltas = ca.get_sun(times)
alphas = lst - walphas # alphas is the HA
elevs = ca.radec2elev(alphas,deltas,lat)


# object coordinates
hao = lst - rao
AzEl = np.array( ca.radec2azel(hao,deo,lat))

# satellites

constellationsl, constellations = cp.findConstellations(constellationsll)
consLab, consNum, consPla, consNPl, consInc, consAlt = ca.loadConstellations(constellations)
nCons = len(consLab)

densSatAll = np.zeros_like(AzEl[0])
densVelAll = np.zeros_like(AzEl[0])
densSatObs = np.zeros_like(AzEl[0])
densVelObs = np.zeros_like(AzEl[0])
densSatBloom = np.zeros_like(AzEl[0])
densVelBloom = np.zeros_like(AzEl[0])


consLab, consNum, consPla, consNPl, consInc, consAlt = ca.loadConstellations(constellations)
nCons = len(consLab)

mageffmax = -99
mageffmin = 99.


icons = 0
while icons < len(consInc):                          
    # constellation
    densSi, veli, magi =  ca.modelOneConstMag(AzEl,lat, 
                            np.reshape(alphas,alphas.shape[0]*alphas.shape[1]),
                            np.reshape(deltas,deltas.shape[0]*deltas.shape[1]),
                            consInc[icons], consAlt[icons], consNum[icons])
    
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
    
    icons += 1


cmap = 'magma'
if magcutflag == 'ALL':
    ds = densSatAll
    dv = densVelAll

elif magcutflag == 'OBS':
    ds = densSatObs
    dv = densVelObs
    labelmag = True              
    cblab = "Number of satellite trails per exposure"
    lvmin = -2.5  # log limits for the colour scale
    lvmax = 1.2  #
#    lvmin = 0   # log limits for the colour scale # 4MOST
#    lvmax = 2.7  # 4MOST

elif magcutflag == 'EFFECT':
    ds = trailf * densSatObs + (1.-trailf)* densSatBloom
    dv = trailf * densVelObs + (1.-trailf)* densVelBloom
    labelmag = True
    cblab = "Fraction of exposure lost"
    lvmin = -3.5  # log limits for the colour scale
    lvmax = 2.2  #
    lvmax = .2  #
    lvmin = -3.5  # log limits for the colour scale # 4MOST
#

    cmap = gyrd    
    
else:
    print("valid for -M: BRIGHT OBS FAINT ALL EFFECT")
    exit(1)
    
dens =  ds* fovl*fovw + dv * fovl * expt    # trailf already accounted for


#== plot

fig = plt.figure(figsize=(12,7))
plt.rc('font',      size=15) #controls default text size
plt.rc('axes', titlesize=15) #fontsize of the title
plt.rc('axes', labelsize=18) #fontsize of the x and y labels
plt.rc('xtick',labelsize=15) #fontsize of the x tick labels
plt.rc('ytick',labelsize=15) #fontsize of the y tick labels
plt.rc('legend',fontsize=15) #fontsize of the legend
#plt.rcParams['figure.figsize'] = [12, 8]


ax = fig.subplots()
_ = ax.set_xticks(np.arange(0.,hstep*24, hstep*3))
_ = ax.set_xticklabels(np.concatenate((np.arange(12,24,3),np.arange(0,12,3))))
_ = ax.set_xlabel('Local Solar Time')
_ = ax.set_yticks(np.arange(0.,365,30.5))
_ = ax.set_yticklabels(['Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'])
_ = ax.set_ylabel('Calendar')
_ = ax.set_xlim(4*hstep, 20*hstep) # from 16h to 8h


lsat = np.log10(dens+1.e-10)
lsat = lsat + np.log( elevs < 0 )

#- fill in the satellites
csat = ax.contourf(lsat, levels=np.arange(lvmin,lvmax,.01) , cmap=cmap, extend='both')
csat.cmap.set_under('k') # below minimum -> black

#- unobservable
osat = np.nan_to_num(lsat*0) + 1000.*( AzEl[1] < elevlim ) 
_ = ax.contourf(osat, levels=np.arange(999.,1001.) , cmap='Greys')

#- airmass
wel =   AzEl[1]  *( elevs < 0 )
cobj = ax.contour(wel, levels=np.arange(elevlim,90.,10.), cmap="summer_r")
lobj = plt.clabel(cobj, fmt='%.0f$^o$')

_ = ax.contourf(elevs, levels=[0.,90.], colors='royalblue')

#- daylight
csun = ax.contour(elevs, levels=[-18.,-12.,-6.,0],
                  linewidths=[1.,1.,1.,5.], 
                  linestyles='solid', 
                  colors='b')
lsun = plt.clabel(csun, fmt='%.0f$^o$')



# color bar
cbar = fig.colorbar(csat)
cbar.set_ticks(np.log10(np.array([0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,5.,10.,20.,50.,100.,200.])))
cbar.set_ticklabels( [ "{:.2g}".format(x) for x in 10.**(cbar.get_ticks())])
cbar.set_label(cblab)



# save
_ = ax.set_title('Object: {} $\\alpha$: {:.2f}$^o$, $\\delta$: {:.2f}$^o$'.format(objlabel,rao,deo))
#ax.grid()

#_ = ax.text(0.,-28,'(b)', fontsize=14)
_ = ax.text(30.,-32,'Constellation: {} ({:.0f} sat.)'.format(constellationsl,sum(consNum)), fontsize=8, ha='left')


_ = ax.text(30,-40,'Telescope: {}, lat= {:.2f}$^o$'.format(telescope,lat),fontsize=8, ha='left')

if fovl < 1./60.:
    fovll = '{:.2f}\"'.format(fovl*3600.)
elif fovl < 1./6.:
    fovll = '{:.2f}\''.format(fovl*60.)
else:
    fovll = '{:.2f}$^o$'.format(fovl)

_ = ax.text(30.,-48,'Instrument: {} FoV= {} Resolution= {:.1g}\" Exp.T= {:.0f}s'.format(instrument, fovll, resol*3600., expt), fontsize=8, ha='left')



outfileroot = "{}_{}_{}_{}".format(objlabel,telinslabel,constellationsll,magcutflag)

fig.tight_layout()

plt.savefig(outpath+'w.png')
plt.savefig(outpath+outfileroot+'.png')
