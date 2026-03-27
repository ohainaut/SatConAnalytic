#!/usr/bin/env python3
'''SatConAnalytic - Satellite Constellation Analytic simulations

ploting functions supporting conAn 
'''
#==============================================================================

import json
import os, logging
from dataclasses import dataclass, field
from typing import Optional, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import SatConAnalytic.conan as caLib
import SatConAnalytic.skyBrightnessLib as skyLib

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
@dataclass
class SIM_CONVERSION_OUT:
    """Result returned by every output_* handler and by selectOutput.

    Fields
    ------
    logDensity : np.ndarray or None
        Log10 of the quantity to plot (or -999 sentinel for empty sky).
        None when the output type is not supported (e.g. skyScattered).
    barLabel : str
        Human-readable label for the colour bar.
    barTicks : np.ndarray or None
        Tick positions on the colour bar (in log-space unless the handler
        uses a linear/negative-mag scale).
    barTickLabels : list[str] or None
        Formatted labels corresponding to barTicks.
    logMinValue : float or None
        Lower bound for the colour-bar axis.
    logMaxValue : float or None
        Upper bound for the colour-bar axis.
    cmap : str or LinearSegmentedColormap or None
        Colour map for the plot.  May be a Matplotlib named colormap string
        (e.g. ``'magma'``) or a pre-built :class:`LinearSegmentedColormap`
        instance (e.g. ``gyrd``).  ``None`` means "use the caller's default".
    """
    density:       Optional[np.ndarray] = None
    logDensity:    Optional[np.ndarray] = None
    title:         Optional[str]        = None
    unit:          Optional[str]        = None
    barLabel:      str                  = ""
    barTicks:      Optional[np.ndarray] = None
    barTickLabels: Optional[list]       = None
    logMinValue:   Optional[float]      = None
    logMaxValue:   Optional[float]      = None
    cmap:          Union[str, LinearSegmentedColormap, None] = 'magma'


#------------------------------------------------------------------------------
@dataclass
class SIM_DENSITIES:
    """Accumulator for all per-sky-element density arrays produced by the
    constellation simulation loop.

    Initialise with :func:`init_output`; accumulate in-place inside the
    shell loop with ``+=``; then pass directly to :func:`magSelect` and
    :func:`selectOutput`.

    Fields
    ------
    satAll : np.ndarray
        Satellite surface density for *all* satellites [N/sq.deg].
    velAll : np.ndarray
        Trail density for *all* satellites [N/deg/s].
    satObs : np.ndarray
        Satellite surface density for observable satellites
        (V_eff < maglim) [N/sq.deg].
    velObs : np.ndarray
        Trail density for observable satellites [N/deg/s].
    satBloom : np.ndarray
        Satellite surface density for blooming satellites
        (V_eff < magbloom) [N/sq.deg].
    velBloom : np.ndarray
        Trail density for blooming satellites [N/deg/s].
    luminance : np.ndarray
        Total luminance from *all* unresolved satellites [Cd/m²].
    luminanceNotDet : np.ndarray
        Luminance from unresolved satellites that are *not* detected
        (V_eff > maglim) [Cd/m²].
    """
    satAll:         np.ndarray = None
    velAll:         np.ndarray = None
    satObs:         np.ndarray = None
    velObs:         np.ndarray = None
    satBloom:       np.ndarray = None
    velBloom:       np.ndarray = None
    luminance:      np.ndarray = None
    luminanceNotDet: np.ndarray = None


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

#------------------------------------------------------------------------
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


#------------------------------------------------------------------------
def setBarLim_negMag(logDensity):
    '''prepare Bar limits, ticks and labels for a mag plot
    The "logDensity" is -mag'''

    logMinValue, logMaxValue = getBarLim(logDensity)
    bMin = np.floor( logMinValue )
    bMax = np.ceil( logMaxValue )
    barTicks = np.arange(bMin, bMax, .5 )

    barTickLabels = [ f'{x:.1f}' for x in -barTicks]
    
    return barTicks, barTickLabels, logMinValue, logMaxValue

#------------------------------------------------------------------------
def setBarLim_lin(density):
    '''prepare Bar limits, ticks and labels for a mag plot
    The "density" is linear density (not log)   '''

    logMinValue, logMaxValue = getBarLim(density)
    bMin = np.floor( logMinValue )
    bMax = np.ceil( logMaxValue )
    barTicks = np.arange(bMin, bMax, .5 )

    barTickLabels = [ f'{x:.1f}' for x in -barTicks]
    
    return barTicks, barTickLabels, logMinValue, logMaxValue

#------------------------------------------------------------------------
def format_tick(val, rounding=1):
    '''format brightness tick labels in a nice way'''

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



#------------------------------------------------------------------------
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


#------------------------------------------------------------------------
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
    
    _lenDensity = len(density[ density > 0] )  

    if _lenDensity == 0:
        log.debug('All sky elements are empty; flagged as -1000')
        return np.full_like(density, -1000.)
    else:
        _n_zero = np.sum(density == 0.)
        if _n_zero > 0:
            log.debug(f'Some sky elements are empty ({int(_n_zero)} occurrence(s)). That\'s OK ')
        with np.errstate(divide='ignore'):
            return np.nan_to_num(np.log10(density), neginf=-999.)


#------------------------------------------------------------------------------
def magSelect(magSelect, sim, myTel):
    '''Select the satellites according to the magnitude limit in myTel.

    IN

    - magSelect : str
        One of ``'all'``, ``'detected'``, ``'notDetected'``,
        ``'oversaturated'``.
    - sim : SIM_DENSITIES
        Accumulated density arrays from the constellation simulation loop.
    - myTel : telescope object
        Provides ``maglim`` and ``magbloom``.

    OUT
    
    - ds : np.ndarray
        Selected satellite surface density [N/sq.deg].
    - dv : np.ndarray
        Selected trail density [N/deg/s].
    - selectionLab : str
        Human-readable description of the selection applied.
    '''

    if magSelect == 'oversaturated':
        ds = sim.satBloom
        dv = sim.velBloom
        selectionLab = f'Selection: V$_{{eff}}$  < {myTel.magbloom:.1f}'
        log.debug(f'Mag: selected oversaturated: V_eff < {myTel.magbloom:.1f}')

    elif magSelect == 'detected':
        ds = sim.satObs
        dv = sim.velObs
        selectionLab = f'Selection: V$_{{eff}}$ < {myTel.maglim:.1f}'
        log.debug(f'Mag: selected detected: V_eff < {myTel.maglim:.1f}')

    elif magSelect == 'notDetected':
        ds = sim.satAll - sim.satObs
        dv = sim.velAll - sim.velObs
        selectionLab = f'Selection: V$_{{eff}}$ > {myTel.maglim:.1f}'
        log.debug(f'Mag: selected non-detected: V_eff > {myTel.maglim:.1f}')

    elif magSelect == 'all':
        ds = sim.satAll
        dv = sim.velAll
        selectionLab = 'Selection: all satellites'
        log.debug('Mag: selected all satellites')

    else:
        log.error(f'invalid mode {magSelect}')
        raise ValueError(f'invalid mode {magSelect}')

    return ds, dv, selectionLab


#----------------------------------------------------------------------------
def init_output(AzEl):
    """Initialise all density accumulator arrays to zero.

    Returns a :class:`SIM_DENSITIES` instance ready to be accumulated
    into inside the constellation shell loop.
    """
    z = np.zeros_like(AzEl[0])
    return SIM_DENSITIES(
        satAll          = z.copy(),   # N/sq.deg, all satellites
        velAll          = z.copy(),   # N/deg/s,  all satellites
        satObs          = z.copy(),   # N/sq.deg, V_eff < maglim
        velObs          = z.copy(),   # N/deg/s,  V_eff < maglim
        satBloom        = z.copy(),   # N/sq.deg, V_eff < magbloom
        velBloom        = z.copy(),   # N/deg/s,  V_eff < magbloom
        luminance       = z.copy(),   # Cd/m², all unresolved sat
        luminanceNotDet = z.copy(),   # Cd/m², unresolved & undetected sat
    )



#--------------------------------------------------------------------------
def output_TrailDens(dv, **_):
    """Handle output type 'TrailDens': trail density per deg/sec.

    Returns an SIM_CONVERSION_OUT data structure.
    """
    log.info("-> Trail Density")
    logDensity = log10Sky(dv)
    barTicks, barTickLabels, logMinValue, logMaxValue = setBarLim_standardLog(logDensity)
    return SIM_CONVERSION_OUT(
        density       = dv,
        logDensity    = logDensity,
        title         = "Trail Density",
        unit          = "Trail/deg/sec.",
        barLabel      = "Trail/deg/sec.",
        barTicks      = barTicks,
        barTickLabels = barTickLabels,
        logMinValue   = logMinValue,
        logMaxValue   = logMaxValue,
    )


#--------------------------------------------------------------------------
def output_satDens(ds, **_):
    """Handle output type 'satDens': satellite density per sq.deg.

    Returns an SIM_CONVERSION_OUT data structure.
    """
    log.info("-> Sat Density")
    logDensity = log10Sky(ds)
    barTicks, barTickLabels, logMinValue, logMaxValue = setBarLim_standardLog(logDensity)
    return SIM_CONVERSION_OUT(
        density       = ds,
        logDensity    = logDensity,
        title         = "Satellite Density",
        unit          = "sat./sq.deg.",
        barLabel      = "Number of sat./sq.deg.",
        barTicks      = barTicks,
        barTickLabels = barTickLabels,
        logMinValue   = logMinValue,
        logMaxValue   = logMaxValue,
    )


#--------------------------------------------------------------------------
def output_skyScattered(**_):
    """
    Placeholder for output type 'skyScattered': 
    scattered light sky brightness (not supported).

    Scattered light calculation implies integration of the full sky, which is computed only in by obsSky. Not available for other scripts.

    """
    log.info("-> skyBrightness for scattered light...")
    msg = "skyScattered not handled here"
    log.warning(msg)
    return SIM_CONVERSION_OUT()


#--------------------------------------------------------------------------
def output_skyDiffuse(unit, luminanceUnresolved_totNotDet, **_):
    """Handle output type 'skyDiffuse': diffuse sky brightness from unresolved satellites.

    Supported units: 'magarcsec2', 'frac', 'microCandelaPerM2', 'nanoLambert'.

    Returns an SIM_CONVERSION_OUT data structure.
    """
    log.info("-> skyBrightness for unresolved satellites...")

    cmap = 'magma'  # default colormap for sky brightness; may be overridden below

    # luminanceUnresolved_totNotDet in Cd/m^2 — convert to requested unit
    if unit == 'magarcsec2':
        log.info("...SkyMagArcsec2")
        skyBrightness = skyLib.luminance_to_magarc2(luminanceUnresolved_totNotDet)

        logDensity = -skyBrightness
        barTicks, barTickLabels, logMinValue, logMaxValue = setBarLim_negMag(logDensity)
        logDensity[luminanceUnresolved_totNotDet == 0.] = logMinValue
        myUnit = "MpSA"

    elif unit == 'frac':
        log.info("...fraction of sky brightness")
        skyBrightness = skyLib.luminance_to_skyFraction(luminanceUnresolved_totNotDet)
        logDensity = log10Sky(skyBrightness)
        logMinValue = -4.2
        logMaxValue = 0.1
        barMin = int(logMinValue * 3.) / 3.
        barMax = int(logMaxValue * 3. - 1) / 3.
        barTicks = np.arange(barMin, barMax, .333333)
        barTickLabels = [f'{x:.1g}' for x in 10.**barTicks]
        myUnit = 'fraction of dark sky'
        cmap = gyrd  # use the custom colormap for losses

    elif unit == 'microCandelaPerM2':
        log.info("...microCandela/m2")
        skyBrightness = luminanceUnresolved_totNotDet * 1e6   # convert to microCandela/m2
        logDensity = log10Sky(skyBrightness)
        barTicks, barTickLabels, logMinValue, logMaxValue = setBarLim_standardLog(logDensity)
        myUnit = r'$\mu$cd/m$^2$'   # raw string for LaTeX

    elif unit == 'nanoLambert':
        log.info("...nanoLambert")
        skyBrightness = skyLib.luminance_to_lambert(luminanceUnresolved_totNotDet) * 1e9
                                                    # convert to nanoLambert
        logDensity = log10Sky(skyBrightness)
        barTicks, barTickLabels, logMinValue, logMaxValue = setBarLim_standardLog(logDensity)
        myUnit = 'nL'

    else:
        msg = f"Invalid unit '{unit}' for skyDiffuse output."
        log.error(msg)
        return SIM_CONVERSION_OUT(barLabel=msg)

    return SIM_CONVERSION_OUT(
        density       = skyBrightness,
        logDensity    = logDensity,
        title         = "Diffuse Brightness",
        unit          = myUnit,
        barLabel      = 'Surface brightness [' + myUnit + ']',
        barTicks      = barTicks,
        barTickLabels = barTickLabels,
        logMinValue   = logMinValue,
        logMaxValue   = logMaxValue,
        cmap          = cmap
    )


#--------------------------------------------------------------------------
def output_losses(sim, myTel, **_):
    """Handle output type 'losses': fraction of FoV lost to satellite trails.

    LOSSES accounts for the width of the trail (as a fraction of the FoV).
    For blooming satellites trailF = 1 (full FoV affected);
    (1 - trailF) factor avoids double-counting.

    Returns an SIM_CONVERSION_OUT data structure.
    """
    log.info("-> fraction lost...")
    log.info("...overrides magSelect with losses")

    ds = myTel.trailf * sim.satObs  + (1. - myTel.trailf) * sim.satBloom
    dv = myTel.trailf * sim.velObs  + (1. - myTel.trailf) * sim.velBloom

    # density: number of satellites in the FoV,
    #          plus number of trails crossing FoV during the exposure,
    #          accounting for trailF and blooming:
    dens = ds * myTel.fovl * myTel.fovw + dv * myTel.fovl * myTel.expt
    logDensity = log10Sky(dens)

    linTicks = np.array([0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.])

    return SIM_CONVERSION_OUT(
        density       = dens,
        logDensity    = logDensity,
        title         = "Fraction of FoV lost",
        unit          = "Fraction lost",
        barLabel      = "Fraction of FoV lost",
        barTicks      = np.log10(linTicks),
        barTickLabels = ["{:.2g}".format(x) for x in linTicks],
        logMinValue   = -3.5,
        logMaxValue   = 0.,
        cmap          = gyrd,  # use the custom colormap for losses
    )


#--------------------------------------------------------------------------
def output_trails(ds, dv, myTel, **_):
    """Handle output type 'trails': number of trails per exposure.

    Magnitude selection must already have been applied to ds and dv.

    Returns an SIM_CONVERSION_OUT data structure.
    """
    log.info("-> number of trails/exp...")

    # density: number of satellites in the FoV,
    #          plus number of trails crossing FoV during the exposure.
    density = ds * myTel.fovl * myTel.fovw + dv * myTel.fovl * myTel.expt
    logDensity = log10Sky(density)

    barTicks, barTickLabels, logMinValue, logMaxValue = setBarLim_standardLog(logDensity)
    return SIM_CONVERSION_OUT(
        density       = density,
        logDensity    = logDensity,
        title         = "Trails per exposure",
        unit          = "Trails/exp.",
        barLabel      = "Number of trails per exp.",
        barTicks      = barTicks,
        barTickLabels = barTickLabels,
        logMinValue   = logMinValue,
        logMaxValue   = logMaxValue,
    )


#--------------------------------------------------------------------------
# Mapping from output name to the corresponding handler function.
output_HANDLERS = {
    "TrailDens":    output_TrailDens,
    "satDens":      output_satDens,
    "skyScattered": output_skyScattered,
    "skyDiffuse":   output_skyDiffuse,
    "losses":       output_losses,
    "trails":       output_trails,
}


#--------------------------------------------------------------------------
def selectOutput(myargs, myTel, ds, dv, sim):
    """Dispatch to the appropriate output handler based on myargs.output.

    IN:
    - myargs : argparse.Namespace
        Must provide ``output`` and, for ``skyDiffuse``, ``unit``.
    - myTel : telescope object
    - ds, dv : np.ndarray
        Pre-selected satellite and trail densities (from :func:`magSelect`).
    - sim : SIM_DENSITIES
        Full density accumulator; passed to handlers that need fields beyond
        the magnitude-selected ``ds``/``dv`` (e.g. ``losses``).

    OUT:
    - SIM_CONVERSION_OUT
        All fields populated by the matching handler.
        On error, returns an instance whose ``barLabel`` contains the
        error message and all other fields are ``None``.
    """
    handler = output_HANDLERS.get(myargs.output)
    if handler is None:
        msg = ("Invalid output type. "
               "Valid options: trails, losses, satDens, TrailDens, skyDiffuse, skyScattered")
        log.error(msg)
        return SIM_CONVERSION_OUT(barLabel=msg)

    return handler(
        myTel = myTel,
        ds    = ds,
        dv    = dv,
        sim   = sim,
        luminanceUnresolved_totNotDet = sim.luminanceNotDet,
        unit  = getattr(myargs, 'unit', None),
    )

#----------------------------------------------------------------------------
