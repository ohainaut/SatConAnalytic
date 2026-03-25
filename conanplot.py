#!/usr/bin/env python3
'''SatConAnalytic - Satellite Constellation Analytic simulations

ploting functions supporting conAn 
'''
#==============================================================================

import json
import sys, logging

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import SatConAnalytic.conan as caLib
import SatConAnalytic.telescopes as telescopes

from matplotlib.colors import  LinearSegmentedColormap



log = logging.getLogger('satcon')


# color map for losses
colors = ["black", "lawngreen", "yellow", "orange", "red", "darkred"]
gyrd = LinearSegmentedColormap.from_list("mycmap", colors)

# colors for airmass
colors = ["darkred","red","orange","yellow","green","lawngreen"]
cairmap = LinearSegmentedColormap.from_list("mycmap", colors)

# colors for sun altitude
colors = ['k','darkblue','mediumblue','b','deepskyblue','paleturquoise']
colors = ['#000', '#00f', '#55f', '#aaf', '#eef']
csunmap = LinearSegmentedColormap.from_list("mycmap", colors)


#------------------------------------------------------------------------------
def init_logger(log):
    log.setLevel(logging.DEBUG)
    log.propagate = False
    if log.handlers:
        return

    # Format
    log_format = logging.Formatter('[%(levelname)-8s %(name)s/%(funcName)s] %(message)s')

    # FileHandler for DEBUG and above
    file_handler = logging.FileHandler('satCon.log')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(log_format)

    # create a StreamHandler for INFO and above
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(log_format)

    # Attach handlers to the log
    log.addHandler(file_handler)
    log.addHandler(console_handler)


#------------------------------------------------------------------------------
class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy arrays and scalar types."""
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        return super().default(obj)

#------------------------------------------------------------------------------
def getTelescope(myargs):
    '''return the Telescope object requested by myargs
    The telescope is defined by the 'code' argument,
    and then its parameters are
    overloaded with other myargs specifications.

    myargs is an argparse object.
    '''

    if myargs.code is None: 
        myargs.code = 'DEFAULT'

    myTel = findTelescope(myargs.code)
    debug = False
      

    for what in [
        'telescope',
        'instrument'
    ]:
        if myargs.__dict__[what] is not None:
            myTel.__dict__[what] = myargs.__dict__[what]

    if myargs.expt is not None:
        myTel.expt = float(myargs.expt)
    if myargs.fovl is not None:
        myTel.fovl = float(myargs.fovl)
    if myargs.fovw is not None:
        myTel.fovw = float(myargs.fovw)



    if myargs.fovw is not None:
        myTel.fovw = float(myargs.fovw)
    else:
        myTel.fovw = myTel.fovl *1.
  
    for what in [
        'expt',
        'fovl',
        'magbloom',
        'maglim',
        'resol',
        'trailf',
        'lat'
    ]:
        if what in myargs.__dict__ and myargs.__dict__[what] is not None:
            myTel.__dict__[what] = float(myargs.__dict__[what])

    return myTel

#------------------------------------------------------------------------------
def findTelescope(telinslabel):
    allTel = telescopes.readTelescopes()
    try:
        return allTel.byCode[telinslabel]
    except KeyError:
        if telinslabel != 'list':
            print(f'{telinslabel} not found in telescope list')
        print('Available telescopes are:')
        for x in sorted(allTel.list):
            print(f'  {x}')
        exit(1)




#------------------------------------------------------------------------------
def initPolPlot(ax):
    '''initialize a polar plot'''
    
    plt.rcParams.update({'font.size': 15})
    ax.set_theta_zero_location('N')  
    ax.set_rticks([30,60,70,80,90])
    ax.set_yticklabels([])
    ax.set_xticks(np.radians(range(0,360,90)))
    ax.set_xticklabels(['N','E','S','W'])
    ax.set_rlim(0.,90)
    #ax.set_facecolor("k")

#------------------------------------------------------------------------------
def plotOneDens(ax, AzEl, dens, label):
    '''plot a density function
    IN
    - ax: axis object where to plot
    - AzEl: array with the azimut and elevations
    - dens: the density array to be plotted
    - label: label...
    OUT: updated plot in ax.

    OBSOLETE?
    '''

    cmap = 'magma'
    clab = 'cyan'
    ccon = 'cyan'

    _ = initPolPlot(ax)
            
    vmin = np.amin(dens)
    vmax = np.amax(dens)

    if (vmax-vmin) >0:
        ccd = ax.contour(np.radians(AzEl[0]), 90.-AzEl[1], dens,
             colors=ccon , linewidths=0.5,
             locator=ticker.LogLocator(subs=(.2,.5,1.)))
        ax.clabel(ccd, colors=clab, fmt='%.2g')
    
        cfd = ax.contourf(np.radians(AzEl[0]), 90.-AzEl[1], dens , 
            locator=ticker.LogLocator(subs=10.**np.arange(.1,1.0,.025)),
            cmap=cmap)
    
    ax.set_title(label)
    return

#------------------------------------------------------------------------------
def plotDens(AzEl, nplot, densities, labels, plotLabel):
    '''plot a series of density distributions
    
    OBSOLETE?
    '''
    
    if nplot == 1:
        fig = plt.figure(figsize=(8,8))
        ax =  fig.subplots(1,nplot,subplot_kw={'projection': 'polar'}) 
        dens = densities
        label = labels
        plotOneDens(ax, AzEl, dens, label)
    else:
        fig = plt.figure(figsize=(5*nplot,6))
        axs =  fig.subplots(1,nplot,subplot_kw={'projection': 'polar'}) 
        for ax, dens, label in zip(axs, densities, labels):
            plotOneDens(ax, AzEl, dens, label)
  
    #        cbar0 = fig.colorbar(cfd, ax=axDens)
    #        cbar0.set_label('Trail/exp.')
    
    fig.suptitle(plotLabel,  size=20)
    
    return


#---------------------------------------------------------------------------
def get_skyColor(sunElev):
    """Return sky color based on sun elevation

    Args:
        sunElev (float): Sun elevation in degrees

    Returns:
        str: Color string
    """
    if sunElev >= 0:
        return (0.7,0.8,1.0)  # Day: light blue
    elif sunElev <= -18:
        return (0.0,0.0,0.0)  # Night: black
    else:
        ds = sunElev/18 + 1.
        return (0.5*ds,0.7*ds,ds)    # Twilight: dark blue to black

#---------------------------------------------------------------------------
def cmap_sky(skyColor, cmap_name='magma'):
    '''Create a colormap for the sky based on sun elevation

    Args:
        sunElev (float): Sun elevation in degrees
    '''

    cmap = plt.get_cmap(cmap_name).copy()
    rgb_colors = np.array(cmap.colors)+  np.array(skyColor)
    rgb_colors = np.clip(rgb_colors, 0, 1)
    cmap.colors = list(rgb_colors)
    return cmap


#------------------------------------------------------------------------------
def drawHADec(lat):
    '''draw the HA and Dec lines in an altaz plot'''

    corl = 0.3
    corc = 'r'

    #draw RA
    HAdef = np.arange(-12.,12.01)*15.
    Decdef = np.arange(-90.,91.,1.)
    for i in np.arange(0,len(HAdef)):
        azs,els = caLib.radec2azel(HAdef[i], Decdef, lat)
        az = azs[els>0]
        el = els[els>0]
        plt.plot(np.radians(az), 90.- el, lw=corl, c=corc)
    
    #label RA
    Dec = 0.
    HAdef = np.arange(-12.,12.01,3.)
    az,el = caLib.radec2azel(HAdef*15., Dec, lat)
    for i in np.arange(0,len(HAdef)):
        if el[i] > 0. :
            plt.text(np.radians(az[i]),90. -el[i],'{:.0f}'.format(HAdef[i]), color=corc, fontsize=10)

    #draw Dec
    HAdef = np.arange(-12.,12.01,0.1)*15.
    Decdef = np.arange(-90.,91.,10.)
    for i in np.arange(0,len(Decdef)):
        azs,els = caLib.radec2azel(HAdef, Decdef[i], lat)
        az = azs[els>0]
        el = els[els>0]
        if Decdef[i] == 0.:
            lw = 1.
        else:
            lw = corl
        plt.plot(np.radians(az), 90.- el, lw=lw, c=corc)
    
    #label Dec
    HA = 0.
    Decdef = np.arange(-60.,61.,30.)
    az,el = caLib.radec2azel(HA, Decdef, lat)
    for i in np.arange(0,len(Decdef)):
        if el[i] > 0. :
            plt.text(np.radians(az[i]),90.-el[i],f'{Decdef[i]:.0f}'+r'$^\circ$', color=corc, fontsize=10)

#------------------------------------------------------------------------------
def draw_airmassShade(ax):
    '''Shades below 30 and 20 deg in an altaz plot'''

    azimuth = np.arange(0.,2.*np.pi,.01)
    radius1   = np.full_like(azimuth, 90.)
    
    whss = [['xx','x'],['++','+']]
    wCols = ['white','black']

    for whs,wc in zip( whss,wCols ): 
        for wr,wh in zip ([20,30],whs):
            wr2 = radius1 -wr
            if 1 : # hash
                alphashade=0.1
                ax.fill_between(azimuth,radius1,wr2,
                    facecolor='none',
                            hatch=wh,edgecolor=wc,alpha=alphashade)
                ax.fill_between(azimuth,radius1,wr2,
                    facecolor='none',
                            hatch=wh,edgecolor=wc,alpha=alphashade)
            else:
                alphashade=0.2
                ax.fill_between(azimuth,radius1,wr2,
                                color="grey", alpha=alphashade)




#---
def getBarLim(logDensity):
    '''
    Top and bottom of the bar for logDensity as a standard log

    logMinValue < logMaxValue, ALWAYS
    '''

    if len(logDensity[logDensity > -998] ) == 0:
        #empty sky, set default values; plot will be black
        logMinValue = -4.
        logMaxValue = 0.9
    else:
        logMinValue =  np.amin(logDensity[logDensity > -998])
        logMaxValue =  np.amax(logDensity[logDensity > -998])

    return logMinValue, logMaxValue


def setBarLim_standardLog(logDensity):
    '''prepare Bar limits, ticks and labels for a standard logDensity'''


    logMinValue, logMaxValue = getBarLim(logDensity)
    bMin = np.floor(logMinValue*3.)/3. 
    bMax = np.ceil(logMaxValue*3. )/3.
    barTicks = np.linspace(bMin, bMax, 10)
    barTicks = barTicks[ barTicks >= logMinValue -.35 ]
    barTicks = barTicks[ barTicks <= logMaxValue +.35 ]

    barTickLabels = [""] + [ format_tick(val) for val in 10.**barTicks[1:-1]] + [""] 
    return barTicks, barTickLabels, logMinValue, logMaxValue


def setBarLim_negMag(logDensity):
    '''prepare Bar limits, ticks and labels for a mag plot
    The "logDensity" is -mag'''

    logMinValue, logMaxValue = getBarLim(logDensity)
    bMin = np.floor( logMinValue )
    bMax = np.ceil( logMaxValue )
    barTicks = np.arange(bMin, bMax, .5 )

    barTickLabels = [ f'{x:.1f}' for x in -barTicks]
    
    return barTicks, barTickLabels, logMinValue, logMaxValue

def setBarLim_lin(density):
    '''prepare Bar limits, ticks and labels for a mag plot
    The "density" is linear density (not log)   '''

    logMinValue, logMaxValue = getBarLim(density)
    bMin = np.floor( logMinValue )
    bMax = np.ceil( logMaxValue )
    barTicks = np.arange(bMin, bMax, .5 )

    barTickLabels = [ f'{x:.1f}' for x in -barTicks]
    
    return barTicks, barTickLabels, logMinValue, logMaxValue


#---- format brightness in a nice way
def format_tick(val, rounding=1):
    if val == np.inf:
        return 'inf'

    if 10 < val < 40:
        return f'{val:.1f}'
    if 0.1 < val < 10:
        return f'{val:.2f}'
    else:
        if val < 0.: 
            sign="-"
            val = -val
        else:
            sign=""

        preFormat = f'{val:.1e}' # eg 1.2e+03

        if preFormat[2] == "0":  # eg 1.0e+03
            root = preFormat[0]
        else:
            root = preFormat[:3]
        if preFormat[3:] == "e+00":
            return sign + root
        if preFormat[3:] == "e-01":
            return sign + f'{val:.1f}'
        if preFormat[5] == "0":
            return sign + root + " 10$^{" + preFormat[4] + preFormat[-1] + "}$"
        else:                    
            return sign + root + " 10$^{" + preFormat[-3:] + "}$"
#----



def prepareLabel_instrument(myTel):
    '''prepare a label for the instrument, to be used in plots'''

    label_instrument = [ f'Instrument: {myTel.instrument}']

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
    label_instrument.append(f'FoV: {fovll} x {fovlw}')


    if myTel.expt > 1.:
        label_instrument.append(f'Exp.t: {myTel.expt:.0f}s')
    else:   
        label_instrument.append(f'Exp.t: {myTel.expt:.0g}s')
    
    return label_instrument



#------------------------------------------------------------------------------
def azlab(ax,x,y,lab, **kwargs):
    '''write labels at (x,y) on a polar plot in ax
    
    Additional keyword arguments (e.g. size, alpha, color, fontsize) are
    forwarded to ax.text.  The bbox alpha is taken from the ``alpha`` kwarg
    when provided (default 0).
    '''

    rlab = np.sqrt(x*x + y*y)*90.
    azlab = np.arctan2(-x,y)
       
    if x < 0 :
        halig = 'left'
    elif x == 0:
        halig = 'center'
    else:
        halig = 'right'

    bbox_alpha = kwargs.pop('alpha', 0.)
    t = ax.text(azlab, rlab, lab,
                verticalalignment='top',
                horizontalalignment=halig,
                **kwargs)
    t.set_bbox(dict(alpha=bbox_alpha, facecolor='white', edgecolor='none'))
    return azlab,rlab


def put_labels(ax, x, y, dy, labels, **kwargs):
    '''put labels at (x,y) with a vertical spacing of dy between lines'''

    kwargs.setdefault('fontsize', 12)
    for i, label in enumerate(labels):
        azlab(ax, x, y - i*dy, label, **kwargs)
    return y - len(labels)*dy



#---------------------------------------------------------------------------
def log10Sky( density ):
    '''Logarithm of the density, 
    with a special handling for empty regions.
    '''
    
    if len( density[ density > 0]  ) > 0: # some non-empty elements
        return np.nan_to_num(np.log10(density ), neginf=-999. )
    else: # empty sky
        return np.full_like(density, -1000.)


#------------------------------------------------------------------------------
def magSelect(magSelect, 
            densSatAll,   densVelAll,
            densSatObs  , densVelObs  , 
            densSatBloom, densVelBloom,
              myTel):
    
    '''select the satellites according to the magnitude selection in myTel.magSelect
    Return the corresponding densities of satellites and satellite trails, and a label for the selection.'''

    # ==magSelect==
    # Select effective densities depending on the requested output:
    #   ds: density of satellites;
    #   dv: density of satellite trails (i.e. density of satellites * velocity)


    if  magSelect == 'oversaturated':
        ds = densSatBloom
        dv = densVelBloom
        selectionLab = f'Selection: V$_{{eff}}$  < {myTel.magbloom:.1f}'
        log.debug(f'Mag: selected oversaturated: V_eff < {myTel.magbloom:.1f}')

    elif magSelect == 'detected':
        ds = densSatObs
        dv = densVelObs
        selectionLab = f'Selection: V$_{{eff}}$ < {myTel.maglim:.1f}'
        log.debug(f'Mag: selected detected: V_eff < {myTel.maglim:.1f}')

    elif magSelect == 'notDetected':  
        ds = densSatAll - densSatObs
        dv = densVelAll - densVelObs    
        selectionLab = f'Selection: V$_{{eff}}$ > {myTel.maglim:.1f}'
        log.debug(f'Mag: selected non-detected: V_eff > {myTel.maglim:.1f}')
        
    elif magSelect == 'all':
        ds = densSatAll
        dv = densVelAll
        selectionLab = 'Selection: all satellites'
        log.debug('Mag: selected all satellites')

    else:
        log.error(f'invalid mode {magSelect}')
        raise ValueError(f'invalid mode {magSelect}')
        
    return ds, dv, selectionLab


