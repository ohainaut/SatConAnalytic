# Functions for DEBRIS simulations in CONAN

import os
import numpy as np

from astropy.table import join
from astropy.io import ascii


import constants 
import conan

mypath = os.path.dirname(__file__)

p_albedo = .2
rho_density = 0.5* 1000. # km/m3
fluxSun = 10.**(-0.4*constants.magSun)

#----------------------------------------------------------------------
def make_debrisTab(myrange=np.arange(0,7)):
    '''return Table with debris distribution

    Cols: 
    - alt: alt in km
    - {i} with i in "0" to "6": density in 1/km3 for radius = 10^-i
    - {i}_a: cross-section density [m2/km3], with p_albedo
    '''
    for i in np.arange( 0, 7):
        infile = mypath + f'/Debris/spatial_density_1e-{i}.csv'
        T = read_oneDebris(infile)
        T["dens"].name = f'{i}'

        if i == 0:
            Debris = T
        else:
            Debris = join( Debris, T, keys="alt")

        r = 10.**(-i) # radius in m
        a = np.pi * r**2 * p_albedo 
        m = 4./3.* np.pi * r**3 * rho_density 
        print(f'i: {i}, r: {r:.1e}m, a: {a:.2e}m2, m: {m:.2e}kg  ')

        Debris[f'{i}_a'] = a * Debris[f'{i}'] # cross-section in m2
        Debris[f'{i}_m'] = m * Debris[f'{i}'] # mass in kg
        
    # thickness of the shells
    Debris['dAlt'] = 0.
    for i in np.arange(len(Debris)-1):
        Debris['dAlt'][i+1] = (Debris['alt'][i+1] - Debris['alt'][i])
    Debris = Debris[1:]


    # consolidate volume densities over shells into spatial densities
    Debris['surf_n'] = sum(  Debris[f'{i}']   for i in myrange)*Debris['dAlt']
    Debris['surf_a'] = sum(  Debris[f'{i}_a'] for i in myrange)*Debris['dAlt']
    Debris['surf_m'] = sum(  Debris[f'{i}_m'] for i in myrange)*Debris['dAlt']


    return Debris


#----------------------------------------------------------------------
def read_oneDebris(infile):
    '''Read one debris file as Table'''

    T = ascii.read(infile)
    T["Altitude [km]"].name = "alt"
    T["Spatial Density [1/km^3]"].name = "dens"
    return T

#--------------------------------------------------------------------------------------------
def modelSkyBrightness(El, f=0.7, mag=22. ):
    '''return skyBrightness [mag/sqarcsec] 

    IN:
    - El: array of elevations [deg]
    - f: fraction of the sky brightness due to airglow
    - mag: skybrightness at zenith
    
    OUT: mag/sq.arcsec
    
    Based on LaPalma TechReport 115
    '''

    airmass = 1./np.cos(np.radians( 90.- El))
    deltaM = -2.5*np.log10(  f * airmass + (1.-f) )
    return mag+deltaM

def modelTwilight(sunEl):
    if sunEl > 0: return 8.
    elif sunEl < -18: return constants.magSky
    else: return  (8.-constants.magSky)*(sunEl+18.)/18. + constants.magSky


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
def modelOneConst_Geometry(AzEl,obsLatitude, sunAlpha,sunDelta,satAlt,
                           step=1. ):
    '''Geometry of single shell over a set of Az,El pointings
    IN
    - AzEl: lists of [[Azimuths],  [Elevation]]   [deg]
            on which the constellation shall be evaluated. 
    - obsLatitude: latitude of the observer [deg]
    - sunAlpha, sunDelta: HourAngle and Dec. of the Sun [deg]
    - satAlt: altitude [km]

    - step: size [deg] of the "pixel"
    
    OUT
    - surface*illum, [km2], effective illuminated surface 
    - sunPhFrac: solar phase correction [0,1]
    '''


    # geocentric equ. alpha,delta of sat, and   observatory dist, angle 
    alpha, delta, Delta, costheta = conan.AltAz2Delta(obsLatitude,satAlt,AzEl)
    
    # geocentric equatorial rectangular coord of satellite
    xyz = conan.RaDecAlt2xyz(alpha,delta, satAlt)

    # Illumination
    illum = conan.solIllum(xyz, sunAlpha, sunDelta) # 1/0 for illumination of the pixel by the sun


    # area [km2] of the section of shell in the pixel-in-the-sky projected on the shell
    surface = (Delta * np.radians(step))**2/costheta

    #
    # solar phase
    #

    # Sun unit vector
    #   The distance to the sun is >> than the other distances
    #   so, sun_XYZ is also Observer-Sun and Satellite-Sun
    sun_XYZ = conan.RaDecAlt2xyz(-sunAlpha,sunDelta, 1.-constants.earthRadius) #for unit vector...
    #   note -sunAlpha

    # Satellite-Observer vector
    # xyz: geocentric vector to the satellite
    # satellite-observer:
    satObs_XYZ = -xyz
    satObs_XYZ[0,:] += constants.earthRadius
    # |satObs_XYZ| = |topo_XYZ| = Delta, computed above

    # Solar phase correction - Lambertian
    sunPhFrac = (1.+ sun_XYZ@satObs_XYZ / Delta )/2. # |sunXYZ|=1

    if 0:
        for ii in np.arange(len(illum)):
            print(ii, AzEl[:,ii], f' [{alpha[ii]:.2f}, {delta[ii]:.2f}] ', 
              f'[ {xyz[0,ii]-constants.earthRadius:.0f},',
              f'{xyz[1,ii]:.0f}, ',
              f'{xyz[2,ii]:.0f} ],',
               illum[ii])
            
    return surface, illum, sunPhFrac, Delta

#--------------------------------------------------------------------------------------------
def modelOneConst_Flux(surface, illum, sunPhFrac, Delta, 
                       crossSectionDensity,
                       step=1. ):
    '''Model the flux from a single shell over a set of Az,El pointings
    IN
    - surface:  area [km2] 
    - illum: illuminated region [0,1]
    - sunPhFrac: sun phase fraction [0..1]
    - Delta: distance observer-shell [km]
         all these for each AzEl pointing (arrays)
    - cross-section density [m2 / km2]

    - step: size [deg] of the "pixel"

    OUT
    - total flux [(flux unit with ZP=0)] from element, from illuminated particles, corrected for phase
    '''


    # Flux:  
    illumFlux =  surface * illum * sunPhFrac * crossSectionDensity* fluxSun / Delta**2 * 1.E-6
    #  1E-6:  crossSection in m2, Delta2 in km2 


    return    illumFlux


#--------------------------------------------------------------------------------------------
def modelOneConst_Count(surface, illum, sunPhFrac, Delta,
                     numberDensity, 
                     step=1.):
    '''Number of particles of single shell over a set of Az,El pointings
    IN
    - surface:  area [km2] 
    - illum: illuminated region [0,1] - not used
    - sunPhFrac: sun phase fraction [0..1] - not used
    - Delta: distance observer-shell [km] - not used
         all these for each AzEl pointing (arrays)
    - numberDensity:  effective crossSection (in m2) per surface area (km2)
    
    - step: size [deg] of the "pixel"

    OUT
    - flux
    '''

    # Number of particles in the pixel, normalized to sq.deg  
    numberCount =  surface *  numberDensity / step**2

    

    return   numberCount




#------------------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------


if __name__ == "__main__":

    D = make_debrisTab()


