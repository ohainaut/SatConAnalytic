
# Random functions related to CONAN

import numpy as np
import constants 
from matplotlib import colormaps
import matplotlib.pyplot as plt

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
    plt.ylabel('Number of particles [n/km$^3$]')
    plt.legend()
    plt.savefig('debris_n.png')

    # a
    plt.figure()
    for i in np.arange( 0, 7):
        plt.plot( D['alt'], D[f'{i}_a'], label=f'10^-{i}', color=cols[i])
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Effective cross-section [m$^2$/km$^3$]')
    plt.legend()
    plt.savefig('debris_a.png')

    # m
    plt.figure()
    for i in np.arange( 0, 7):
        plt.plot( D['alt'], D[f'{i}_m'], label=f'10^-{i}', color=cols[i])
    plt.yscale('log')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Mass [kg/km$^3$]')
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



