import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps

from astropy.table import Table, join
from astropy.io import ascii


import constants 
import conan

mypath = os.path.dirname(__file__)

p_albedo = .1
rho_density = 0.5* 1000. # km/m3
fluxSun = 10.**(-0.4*constants.magSun)

#----------------------------------------------------------------------
def make_debrisTab(myrange=np.arange(0,7)):
    '''return Table with debris distribution

    Cols: 
    - alt: alt in km
    - i with i in "0" to "6": density in 1/km3 for radius = 10^-i
    - i_a: cross-section density [m2/km3], with p_albedo
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
    Debris['surf_a'] = 0. # consolidated cross section spatial density, at altitude
    Debris['surf_m'] = 0. # consolidated mass          spatial density, at altitude
    ##for ialt in np.arange(len(Debris)): # we drop the last 
    ##    Debris['surf_a'][ialt] = sum(  Debris[f'{i}_a'][ialt] for i in np.arange(0,7))*Debris['dAlt'][ialt]
    ##    Debris['surf_m'][ialt] = sum(  Debris[f'{i}_m'][ialt] for i in np.arange(0,7))*Debris['dAlt'][ialt]

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

#----------------------------------------------------------------------
def surfaceOfSphericalRectangle(   dalpha, delta, ddelta, radius ):
    '''The surface of a section of a shell
    IN
    - dalpha: extent of the section in longitude [deg]
    - delta, ddelta: origin and extent of the section in latitude [deg]
    - radius [unit]

    OUT:
    - surface of the section of the shell [unit^2]


    Notes:
    -  surface of full cap from colat=90-delta to pole
          2pi r2 (1-cos(colat)) = 2pi r2 (1-sin(delta))
    -  surface of full ring from delta1 to delta2:
          2pi r2 (sin(delta1) - sin(delta2))
    - partial ring from  alpha1 to alpha2
          (alpha2-alpha1)/360 *   2pi r2 (sin(delt1)-sin(delta2))


    - to get surface in sq.deg: radius= 180./pi
    -                   sq.arcsec: radius = 3600.*180./pi
    '''

    surface  =  2.*np.pi* radius**2 # 1/2 sphere
    surface *=  (np.sin(np.radians(delta +ddelta)) - np.sin(np.radians(delta))) 
         # difference of callotes; =2 for full range
    surface *=  dalpha/360. # longitude fraction

    return surface

#-------------------------------------------------------------------------
def surfaceOfSection(  dalpha, delta, ddelta, alt):
    return surfaceOfSphericalRectangle(dalpha, delta, ddelta, 
                                       constants.earthRadius + alt )
#-------------------------------------------------------------------------
def sqArcsecOfSection(  dalpha, delta, ddelta):
    return surfaceOfSphericalRectangle(dalpha, delta, ddelta, 
                                       3600.*180./np.pi )



#----------------------------------------------------------------------
def get_totalShellSurface(D):
    # compute surface of the shells
    D['surface']  =  surfaceOfSection(  360., -90., 180., D['alt'])

#----------------------------------------------------------------------
def VolumeOfPixel(  dalpha, delta, ddelta, alt, dradius):
    '''The volume of a section of a shell 
   
    IN
    - dalpha: extent of the section in longitude [deg]
    - delta, ddelta: origin and extent of the section in latitude [deg]
    - alt, dradius: inner altitude and thickness of shell [unit]
    
    OUT
    - volume of the section of the shell [unit^3]
    
    Notes:
    - Volume of the sphere r
          4/3 pi r3
    - shell between r and r+dr
         4/3 pi ( (r+dr)3 - r3 )
         with   (r+dr)3 = r3 + 3r2 dr + 3 r dr2 + dr3
         4/3 pi  (  3r2 dr   +  3r dr2  + dr3 )
         ~   4pi r2 dr = surface * dr
    '''
    surface = surfaceOfSection(  dalpha, delta, ddelta, alt)

    return surface * dradius

#----------------------------------------------------------------------
def get_totalShellVolume(D):
    # compute volume of the shells
    D['volume'] =  VolumeOfPixel(  360., -90., 180., D['alt'], D['dAlt'])

#----------------------------------------------------------------------
def plot_density(D):
    '''plot raw densities'''

    cmap = colormaps['plasma']
    cols = cmap(np.linspace( 0, 7))

    # n
    plt.figure()
    for i in np.arange( 0, 7):
        plt.plot( D['alt'], D[f'{i}'], label=f'10^-{i}', color=cols[i])
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Number of particles/km$^3$')
    plt.legend()
    plt.savefig('debris_n.png')

    # a
    plt.figure()
    for i in np.arange( 0, 7):
        plt.plot( D['alt'], D[f'{i}_a'], label=f'10^-{i}', color=cols[i])
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Effective cross-section [m$^2$]')
    plt.legend()
    plt.savefig('debris_a.png')

    # m
    plt.figure()
    for i in np.arange( 0, 7):
        plt.plot( D['alt'], D[f'{i}_m'], label=f'10^-{i}', color=cols[i])
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Mass [kg]')
    plt.legend()
    plt.savefig('debris_m.png')

#----------------------------------------------------------------------
def plot_integral(D):
    '''convert volume densities to surface densities,
    plot them,
    plot the total (pre-computed)'''

    cmap = colormaps['plasma']
    cols = cmap(np.linspace( 0, 7))

    # compute volume of the shells
    get_totalShellVolume(D)

    # total plots: multiply density by volume
    # ...and compute grand total

    plt.cla()
    nTotTot = 0.
    for i in np.arange( 0, 7):
        D[f'T{i}']   = D[f'{i}'] * D['volume']
        ntot = sum(D[f'T{i}'])
        nTotTot += ntot
        plt.plot( D['alt'], D[f'T{i}'], label=f'10^-{i}: n$_t$={ntot:.2e}', color=cols[i])

    D['T_n'] = D['surf_n'] * D['surface']
    ntot = sum(D[f'T_n'])
    plt.plot(D["alt"], D["T_n"], "r:", linewidth=4, label=f"Total: n$_t$={ntot:.2e}m$^2$")


    plt.title(f'Particle number (total= {nTotTot:.2e})')
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Number of particles in shell')
    plt.legend()
    plt.savefig('debris_Tn.png')
    plt.show()


    plt.cla()
    aTotTot = 0.
    for i in np.arange( 0, 7):
        D[f'T{i}_a'] = D[f'{i}_a'] * D['volume']
        atot = sum(D[f'T{i}_a'])
        aTotTot += atot
        plt.plot( D['alt'], D[f'T{i}_a'], label=f'10^-{i}: a$_t$={atot:.2e}m$^2$', color=cols[i])


    D['T_a'] = D['surf_a'] * D['surface']
    atot = sum(D[f'T_a'])
    plt.plot(D["alt"], D["T_a"], "r:", linewidth=4, label=f"Total: a$_t$={atot:.2e}m$^2$")
    plt.title(f'Cross-section (total= {aTotTot:.2e}m$^2$)')
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Total cross-section in shell [m2]')
    plt.legend()
    plt.savefig('debris_Ta.png')
    plt.show()


    plt.cla()
    mTotTot = 0
    for i in np.arange( 0, 7):
        D[f'T{i}_m'] = D[f'{i}_m'] * D['volume']
        mtot = sum(D[f'T{i}_m'])
        mTotTot += mtot
        plt.plot( D['alt'], D[f'T{i}_m'], label=f'10^-{i}: m$_t$={mtot:.2e}kg', color=cols[i])

    D['T_m'] = D['surf_m'] * D['surface']
    mtot = sum(D[f'T_m'])
    plt.plot(D["alt"], D["T_m"], "r:", linewidth=4, label=f"Total: m$_t$={mtot:.2e}kg")



    plt.title(f'Mass (total= {mTotTot:.2e}kg)')
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Mass in shell [kg]')
    plt.legend()
    plt.savefig('debris_Tm.png')
    plt.show()

#--------------------------------------------------------------------------------------------
def plot_consolidated(D):
    plt.cla()

    get_totalShellSurface(D)

    D['a'] = D['surf_a']* D['surface']
    aTotTot = sum( D['a'])
    plt.plot( D['alt'], D['a'], label='Xsection')
    plt.title(f'Cross-section (total= {aTotTot:.2e}m$^2$)')
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Total cross-section per shell[m2]')
    plt.legend()
    plt.savefig('debris_Ca.png')
    plt.show()

    D['m'] = D['surf_m']* D['surface']
    mTotTot = sum( D['m'])
    plt.plot( D['alt'], D['m'], label='Mass')
    plt.title(f'Mass (total= {mTotTot:.2e}kg$)')
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Total mass per shell[kg]')
    plt.legend()
    plt.savefig('debris_Cm.png')
    plt.show()



#--------------------------------------------------------------------------------------------
def modelSkyBrightness(El, f=0.8, mag=22. ):
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
    elif sunEl < -18: return 22.
    else: return  -14.*(sunEl+18.)/18. + 22.


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

def modelOneConstMag(AzEl,obsLatitude, sunAlpha,sunDelta,
                     satAlt, crossSectionDensity ):
    '''Model one single shell over a set of Az,El pointings
    IN
    - AzEl: lists of [[Azimuths],  [Elevation]]   [deg]
            on which the constellation shall be evaluated. 
    - obsLatitude: latitude of the observer [deg]
    - sunAlpha, sunDelta: HourAngle and Dec. of the Sun [deg]
    - satAlt: altitude [km]
    - crossSectionDensity:  effective crossSection (in m2) per surface area (km2)
    
    OUT
    - flux
    '''

    step=1. # [deg] here we will work in sq.deg


    # geocentric equ. alpha,delta of sat, and   observatory dist, angle 
    alpha, delta, Delta, costheta = conan.AltAz2Delta(obsLatitude,satAlt,AzEl)
    
    # geocentric equatorial rectangular coord of satellite
    xyz = conan.RaDecAlt2xyz(alpha,delta, satAlt)

    # Illumination
    illum = conan.solIllum(xyz, sunAlpha, sunDelta) # 1/0 for illumination of the pixel by the sun


    # area [km2] of the section of sphere in the pixel
    surface =  surfaceOfSection(  step, delta, step, satAlt)

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


    # Flux:  
    illumFlux =  surface * illum * sunPhFrac * crossSectionDensity* fluxSun / Delta**2 * 1.E-6
    #  1E-6:  crossSection in m2, Delta2 in km2 

    #arcsec2 =  180./np.pi*(np.sin( np.radians( AzElreshape[1,:] +1 )) - np.sin( np.radians( AzElreshape[1,:] )))
    arcsec2 = step**2 * 3600.**2
    illumFluxA2 = illumFlux/ arcsec2


    return    illumFluxA2
#------------------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------


if __name__ == "__main__":

    D = make_debrisTab()
    #plot_density(D)
    #plot_integral(D)
    plot_consolidated(D)
    #print( surfaceOfSphericalRectangle(   1., 0., 1., 3600.*180./np.pi ))
