#!/usr/bin/env python3
# SatConAnalytic - Satellite Constellation Analytic simulations
#
# plot the satellite density over a calendar for a given object

import numpy as np
import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt
from conanplot import gyrd
import argparse

# import ConAn routines
import conan as ca
import conanplot as cp


#outpath = "/home/ohainaut/public_html/outsideWorld/"
outpath = "./"



#--- command line arguments
parser = argparse.ArgumentParser(description='Constellation density calendar for an object')
parser.add_argument('-a','--RA',         default=0.,
                    help='''Right Ascension [deg]''')
parser.add_argument('-d','--DEC',        default=0.,
                    help='''Declination [deg]''')
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

# DEFAULTS
maglimoff = 1.75   # 1.75: 1/5sig ; 1.: 2/5sig =offset between limiting magnitude and detection [mag]
plotMode = myargs.mode 
elevlim = 20.        # consider only elevations >= elevlim
      
# no default:
rao      = float(myargs.RA)
deo      = float(myargs.DEC)
objlabel = myargs.objlabel

#- find telescope
print('TELESCOPE/INSTRUMENT SETUP')
myTel = cp.getTelescope(myargs)
print(myTel)

#---

# satellites
CONSTELLATIONS = ca.findConstellations(myargs.constellation)
print('CONSTELLATIONS:')
print(CONSTELLATIONS.ToC)
print()

#- fill in the calendar
hstep = 12 # points per hour
times = np.zeros((365,24*hstep))
for i in range(365):
    for j in range(24*hstep):
        times[i,j] = 2459580.5 + i + j/24./hstep - (.5 ) # to centre midnight
        

# Loc ST
lst = (280.46061837 + 360.98564736629*(times -2451545.0) )%360

# sun coordinates
walphas, deltas = ca.get_sun(times)
alphas = lst - walphas # alphas is the HA
elevs = ca.radec2elev(alphas,deltas,myTel.lat)


# object coordinates to Az,Elevation
hao = lst - rao # we work in hour angle
AzEl = np.array( ca.radec2azel(hao,deo,myTel.lat))


densSatAll = np.zeros_like(AzEl[0])
densVelAll = np.zeros_like(AzEl[0])
densSatObs = np.zeros_like(AzEl[0])
densVelObs = np.zeros_like(AzEl[0])
densSatBloom = np.zeros_like(AzEl[0])
densVelBloom = np.zeros_like(AzEl[0])

mageffmax = -99
mageffmin = 99.


for myShell in CONSTELLATIONS.shells:                 
    # process each shell

    densSi, veli, magi =  ca.modelOneConstMag(AzEl,myTel.lat, 
                            np.reshape(alphas,alphas.shape[0]*alphas.shape[1]),
                            np.reshape(deltas,deltas.shape[0]*deltas.shape[1]),
                            myShell.inc,
                            myShell.alt, 
                            myShell.totSat
                            )
    
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
    


cmap = 'magma'
if plotMode == 'ALL':
    ds = densSatAll
    dv = densVelAll
    lvmin = -2.5  # log limits for the colour scale
    lvmax = 1.2  #
    cblab = "Number of satellite trails per exposure"


elif plotMode == 'OBS':
    ds = densSatObs
    dv = densVelObs
    labelmag = True              
    cblab = "Number of satellite trails per exposure"
    lvmin = -2.5  # log limits for the colour scale
    lvmax = 1.2  #

elif plotMode == 'EFFECT':
    ds = myTel.trailf * densSatObs + (1.)* densSatBloom
    dv = myTel.trailf * densVelObs + (1.)* densVelBloom
    # non-blooming satellites affect myTel.trailf * field,
    # blooming ones affect 100% * field
    labelmag = True
    cblab = "Fraction of exposure lost"
    lvmax = .2  #
    lvmin = -3.5  # log limits for the colour scale # 4MOST
    cmap = gyrd    
    
else:
    raise ValueError("valid Modes: BRIGHT OBS FAINT ALL EFFECT")
    
# effective density:
# ds [spatial density n/deg^2] x (FoVw X FoVl) = instantaneous contribution
#+dv [velocity density n/deg/s] x FoVl [deg] * expt [s] = trail contribution
dens =  ds* myTel.fovl*myTel.fovw + dv * myTel.fovl * myTel.expt    # trailf already accounted for


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
ax.set_xticks(np.arange(0.,hstep*24, hstep*3))
ax.set_xticklabels(np.concatenate((np.arange(12,24,3),np.arange(0,12,3))))
ax.set_xlabel('Local Solar Time')
ax.set_yticks(np.arange(0.,365,30.5))
ax.set_yticklabels(['Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'])
ax.set_ylabel('Calendar')
ax.set_xlim(4*hstep, 20*hstep) # from 16h to 8h


lsat = np.log10(dens+1.e-10)
lsat = lsat + np.log( elevs < 0 )

#- fill in the satellites
csat = ax.contourf(lsat, levels=np.arange(lvmin,lvmax,.01) , cmap=cmap, extend='both')
csat.cmap.set_under('k') # below minimum -> black

#- unobservable
osat = np.nan_to_num(lsat*0) + 1000.*( AzEl[1] < elevlim ) 
ax.contourf(osat, levels=np.arange(999.,1001.) , cmap='Greys', alpha=1.)

#- airmass
wel =   AzEl[1]  *( elevs < 0 )
cobj = ax.contour(wel, levels=np.arange(elevlim,90.,10.), cmap="summer_r")
lobj = plt.clabel(cobj, fmt='%.0f$^o$')

#- daylight
for myElevMin, myAlpha  in zip([-18., -12., -6., 0.], [.3,.3,.5,1.]):
    ax.contourf(elevs, levels=[myElevMin,90.], colors='royalblue', alpha=myAlpha)

#- twilightss
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



# -labels
ax.set_title('Object: {} $\\alpha$: {:.2f}$^o$, $\\delta$: {:.2f}$^o$'.format(objlabel,rao,deo))
#ax.grid()

ax.text(30.,-32,f'Constellation: {CONSTELLATIONS.name} ({CONSTELLATIONS.totSat:.0f} sat.)', 
        fontsize=8, ha='left')


ax.text(30,-40,f'Telescope: {myTel.telescope}, lat= {myTel.lat:.2f}$^o$',fontsize=8, ha='left')

if myTel.fovl < 1./60.:
    fovll = '{myTel.fovl_arcsec:.2f}\"'
elif myTel.fovl < 1./6.:
    fovll = '{myTel.fovl*60.:.2f}\''
else:
    fovll = '{myTel.fovl:.2f}$^o$'

ax.text(30.,-48,f'Instrument: {myTel.instrument} FoV= {fovll} '+
        f'Resolution= {myTel.resol_arcsec:.1g}\" Exp.T= {myTel.expt:.0f}s',
        fontsize=8, ha='left')

fig.tight_layout()

#- save plot
if objlabel == "": objlabel = "OBJ"
outfileroot = f"{objlabel}_{myargs.code}_{myargs.constellation}_{plotMode}"
plt.savefig(outpath+'w.png')
plt.savefig(outpath+outfileroot+'.png')
