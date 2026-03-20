import numpy as np
from astropy.table import Table
from astropy.io import ascii

import SatConAnalytic.satDots as satDots
import SatConAnalytic.conan as caLib

# CONVERSIONS

#==============================================================================
# Illuminance & Magnitude conversions
#==============================================================================

#illuminanceLuxZeroPoint =  14.18
#REF https://www.vcalc.com/equation/?uuid=47d7099a-ea1c-11e3-b7aa-bc764e2038f2#
# wikipedia
# all from https://web.archive.org/web/20131202231237/http://members.ziggo.nl/jhm.vangastel/Astronomy/Formules.pdf
# too low for sun: 1.07e+05 lux for mag -26.75 star

illuminanceLuxZeroPoint =   13.990065777270708
# based on illuminanceFootCandleZeroPoint = -16.57
# and 1 footcandle = 10.764 lux
# from Krisciunas & Schaefer 1991, Eq.[8]  [checked]
# Validated for Sun: 1.27e+05 lux for mag -26.75 star (expected: 1.2E5 lux)

def magnitude_to_illuminance(mag):    
    """Convert V [magnitude] to illuminance [lux = lm/m²].
    
    Illuminance in lux:  
      10**(-0.4 * (magnitude + illuminanceLuxZeroPoint)) 
    [ok]
    """
    return 10**(-0.4 * (mag + illuminanceLuxZeroPoint))


def magnitude_to_illuminanceFootCandles(magnitude):
    """
    Convert magnitude to illuminance in FootCandles.

    Ref: Krisciunas & Schaefer 1991, Eq.[8]  [checked]

    [ok]
    Validation: 2.128139E-6 lux for mag 0 star. OK with 1fc = 10.764 lux.
    """
    return 10**(-0.4 * (magnitude + 16.57)) #[ok]

def illuminance_to_magnitude(illum):
    """Convert illuminance [lux = lm/m²] to V [magnitude].
    """
    return -2.5 * np.log10(illum) - illuminanceLuxZeroPoint


def illuminance_lux_to_Wm2(lux):
    """Convert illuminance [lux = lm/m²] to flux [W/m²].
    Scaled from Sun: 1.27E+05 lux -> 164.35 W/m²
    1 lux = 1/773 W/m²
    1 W/m² = 773 lux
    [ok]
    """
    return lux / 773.

def Wm2_to_illuminance_lux(flux):
    """Convert flux [W/m²] to illuminance [lux = lm/m²].
    [ok]
    """
    return flux * 773.


def illuminance_lux_to_footcandle(illum):
    """Convert illuminance [lux = lm/m²] to foot-candle [fc = lm/ft²].
    Ref: derived from definition of foot-candle. 
    [checked]"""
    return illum / 10.764

def illuminance_footcandle_to_lux(fc):
    """Convert foot-candle [fc = lm/ft²] to illuminance [lux = lm/m²].
    Ref: derived from definition of foot-candle.  
    [checked]"""
    return fc * 10.764



solarFluxV = 164.346 # W/m² Solar flux through V filter, via integration of solar spectrum * V filter
solarMagV = -26.75  # Solar magnitude in V band
def flux_to_magnitude(flux):
    """Convert flux [W/m^2] to V [magnitude]."""
    return -2.5 * np.log10(flux / solarFluxV) + solarMagV

def magnitude_to_flux(magnitude):
    """Convert V [magnitude] to flux [W/m^2].
    [ok]
    """
    return solarFluxV * 10**(-0.4 * (magnitude - solarMagV))  


solarFluxDensity550nm = 1873.300 # W/m²/μm Solar flux density at 550nm
def fluxDensity550nm_to_magnitude(fluxDensity):
    """Convert flux density at 550nm [W/m²/μm] to V [magnitude]."""
    return -2.5 * np.log10(fluxDensity / solarFluxDensity550nm) + solarMagV

def magnitude_to_fluxDensity550nm(magnitude):
    """Convert V [magnitude] to flux density at 550nm [W/m²/μm].
    """
    return solarFluxDensity550nm * 10**(-0.4 * (magnitude - solarMagV)) 

#==============================================================================
# Luminance and surface brightness conversions
#==============================================================================
darkSky_MpSA=22.0 # mag/arcsec^2 for natural dark sky



#luminanceZP = 10.8E4     # cd/m2 for 0 mag/arcsec2 Allen 1973, p.26
#luminanceZP = 10.8864E4  # cd/m2 Bara et al. 2016 RSocOpenSci 3, Eq.2.3 
luminanceZP = 10.987E4   # cd/m2 Bara 2019 for 5500K blackbody Tab1 Eq10
def magarc2_to_luminance(magnitude_tarcsec2):
    """Convert surface brightness [mag/arcsec^2] to luminance [cd/m^2].
    Reference:
    - Allen, C. W. (1973). Astrophysical Quantities. Allen’s Astrophysical Quantities (3rd ed.) p.26
    - 
    """
    return luminanceZP* 10**(-0.4 * magnitude_tarcsec2)

darkSky_luminance = magarc2_to_luminance(darkSky_MpSA) # cd/m^2 for natural dark sky
print(f'Dark sky luminance: {darkSky_luminance:.2e} cd/m^2')

def luminance_to_magarc2(luminance):
    """Convert luminance [cd/m^2] to surface brightness [mag/arcsec^2]."""
    return -2.5 * np.log10(luminance / luminanceZP)

def luminance_to_skyFraction(luminance):
    """Convert luminance [cd/m^2] to sky fraction."""
    return luminance / darkSky_luminance

def luminance_to_lambert(luminance):
    """Convert luminance [cd/m²] to Lambert [L].
    
    Reference: definition of Lambert: 
      1 L = 1e4/π cd/m² = 3183.1 cd/m²
    """
    return luminance * 1e-4 * np.pi
darkSky_lambert = luminance_to_lambert(darkSky_luminance) # Lambert for natural dark sky
print(f'Dark sky brightness: {darkSky_lambert:.2e} Lambert')


def lambert_to_luminance(lambert):
    """Convert Lambert [L] to luminance [cd/m²].

    Ref: definition of Lambert: 1 L = 1e4/π cd/m²
    """
    return lambert * 1e4 / np.pi # [cd/m²]

def lambert_to_microCandelaPerM2(lambert):
    """Convert Lambert [L] to microCandela per m²."""
    return lambert_to_luminance(lambert) * 1e6 # [μcd/m²]

def magarc2_to_lambert(magnitude_arcsec2):
    """Convert surface brightness [mag/arcsec^2] to Lambert [L].
    """
    return luminance_to_lambert(magarc2_to_luminance(magnitude_arcsec2)) 

def lambert_to_magarc2(lambert):
    """Convert Lambert [L] to surface brightness [mag/arcsec^2]."""
    return luminance_to_magarc2(lambert_to_luminance(lambert))  

def lambert_to_skyFraction(lambert):
    """Convert Lambert [L] to sky fraction."""
    return lambert / darkSky_lambert


def magarc2_to_skyFraction(magnitude_tarcsec2):
    """Convert surface brightness [mag/arcsec^2] to sky fraction."""
    luminance = magarc2_to_luminance(magnitude_tarcsec2)
    darkSky_luminance = magarc2_to_luminance(darkSky_MpSA)
    return luminance / darkSky_luminance

def skyFraction_to_magarc2(sky_fraction):
    """Convert sky fraction to surface brightness [mag/arcsec^2]."""
    luminance = sky_fraction * darkSky_luminance
    return luminance_to_magarc2(luminance)

def magarc2_to_flux(mag):
    """Convert magnitude per arcsec^2 to flux [W/m^2/sr].
    Using solar spectrum
    """
    flux = magnitude_to_flux(mag)  # W/m²/arcsec2
    area = sr_to_squarearcsec(1.)
    return flux * area



def flux_to_magarc2(flux):
    """Convert flux [W/m^2/sr] to magnitude per arcsec^2.
    For V filter using solar spectrum"""

    fluxarcsec2 = flux / sr_to_squarearcsec(1.)
    return flux_to_magnitude(fluxarcsec2)





#==============================================================================
# TRIGONOMETRY / AIRMASS
#==============================================================================

def squarearcsec_to_sr(squarearcsec):
    """
    Convert area from square arcseconds to steradians.
    
    Parameters:
    squarearcsec (float or np.ndarray): Area in square arcseconds.
    
    Returns:
    float or np.ndarray: Area in steradians.
    """
    return squarearcsec * ( (np.pi / (180.0 * 3600.0)) ** 2 )

def sr_to_squarearcsec(sr):
    """
    Convert area from steradians to square arcseconds.
    
    Parameters:
    sr (float or np.ndarray): Area in steradians.
    
    Returns:
    float or np.ndarray: Area in square arcseconds.
    """
    return sr / ( (np.pi / (180.0 * 3600.0)) ** 2 )
def airmassSimple(Z):
    """
    Airmass as a function of Zenithal Distance
    """

    return 1./np.cos(np.radians(Z))

def airmassFull(Z):
    """
    Scattering airmass as a function of Zenithal Distance, accounting for atmospheric refraction and curvature.
    Ref: K&S91 Eq.[3]    
    
    Z: Zenithal Distance in degrees
    """
    return 1./np.sqrt( 1 - 0.96* np.sin(np.radians(Z))**2)

def angular_distance(alt1, az1, alt2, az2):
    """
    Calculate the angular distance between two points on the celestial sphere.
    
    Parameters:
    alt1, az1: Altitude and Azimuth of the first point in degrees.
    alt2, az2: Altitude and Azimuth of the second point in degrees.
    
    Returns:
    float or np.ndarray: Angular distance in degrees.
    """
    # Convert degrees to radians
    alt1_rad = np.radians(alt1)
    az1_rad = np.radians(az1)
    alt2_rad = np.radians(alt2)
    az2_rad = np.radians(az2)
    
    # Haversine formula for angular distance
    delta_az = az2_rad - az1_rad
    delta_alt = alt2_rad - alt1_rad
    
    a = np.sin(delta_alt / 2)**2 + np.cos(alt1_rad) * np.cos(alt2_rad) * np.sin(delta_az / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    return np.degrees(c)



#==============================================================================
# Dark Sky brightness
#==============================================================================

def skyBrightness(Z, Bz=magarc2_to_lambert(22.)*1E9 , k=0.125) :
    """
    Natural Sky brightness at a given zenithal distance Z, scaled from zenithal brightness Bz
    - Z : Zenithal Distance in degrees
    - Bz: Sky brightness at zenith in nanoLamberts
    - k : extinction coefficient in mag/airmass

    IN:
    Z : Zenithal Distance in degrees
    Bz: Sky brightness at zenith in nanoLamberts (default 22 mag/arcsec2)
    k : extinction coefficient in mag/airmass (default 0.125)

    """

    X = airmassFull(Z)

    return Bz*10**(-0.4*k*(X-1))*X


def moonMagnitude(phase_angle):
    """
    Approximate magnitude of the Moon as a function of phase angle (degrees).
    
    Parameters:
    phase_angle (float or np.ndarray): Phase angle in degrees.
    
    Returns:
    float or np.ndarray: Apparent magnitude of the Moon.
    """
    return -12.73 + 0.026 * abs(phase_angle) + 4e-9 * abs(phase_angle)**4

#==============================================================================
# Scattering
#==============================================================================

def scattering_Rayleigh(theta):
    """
    Rayleigh scattering function.
    
    IN:
    theta (float or np.ndarray): Angular distance in degrees.
    
    OUT:
    float or np.ndarray: Rayleigh scattering function value.
    """
    return 10**5.36 * (1.06 + np.cos(np.radians(theta))**2)


def scattering_Mie(theta):
    """
    Mie scattering function (valid for theta > 10 degrees).
    
    Parameters:
    theta (float or np.ndarray): Angular distance in degrees.
    
    Returns:
    float or np.ndarray: Mie scattering function value.
    """
    return 10**(6.15 - theta/40.)


def scattering_total(theta):
    """
    Total scattering function (Rayleigh + Mie).
    
    Parameters:
    theta (float or np.ndarray): Angular distance in degrees.
    
    Returns:
    float or np.ndarray: Total scattering function value.
    """
    return scattering_Rayleigh(theta) + scattering_Mie(theta)

def Scattered_brightness_nanoLambert( magnitude, theta, Z, k=0.125) :
    """
    Scattered sky brightness at angular distance theta from a source of magnitude 'magnitude' at zenithal distance Z

    IN:
    -magnitude : Apparent magnitude of the light source (accounting for atmospheric extinction)
    -theta : Angular distance from the light source in degrees
    -Z : Zenithal Distance in degrees
    -k : extinction coefficient in mag/airmass (default 0.125)

    Out
    -Sky brightness at zenith in nanoLamberts (default 22 mag/arcsec2)
    
    Note: K&S91 Eq.[7] corrects for the extinction. This is not done here, as ConAn  provides magnitudes already corrected for extinction..

    """

    X = airmassFull(Z)

    # original K&S91:
    # IlluminanceOutOfAtmos = magnitude_to_illuminanceFootCandles(magnitude) # Eq.[8]
    # Illuminance = IlluminanceOutOfAtmos * 10**(-0.4*k*X) # Extinction Eq.[7]

    Illuminance = magnitude_to_illuminanceFootCandles(magnitude) # Eq.[8]

    SurfBrightness = Illuminance * scattering_total(theta) 
    SurfBrightness *= ( 1. - 10.**(-0.4*k*X) )  # Eq.[13]->[12] = [15]

    return SurfBrightness

def Scattered_brightness_magarc( magnitude, theta, Z, k=0.125) :
    return luminance_to_magarc2(
        lambert_to_luminance(
            1e-9* Scattered_brightness_nanoLambert(magnitude, theta, Z, k)
            )  )


#=============================================================================

def get_yale_bright_stars():
    from astroquery.vizier import Vizier
    import astropy.units as u

    # Initialize Vizier with no row limit (HR has ~9110 entries)
    # We specifically request J2000 coordinates and the V magnitude.
    v = Vizier(columns=['HR', '_RAJ2000', '_DEJ2000', 'Vmag'], 
               row_limit=-1)

    # Yale Bright Star Catalog (V/50)
    # Catalog V/50/catalog is the main table
    catalogs = v.get_catalogs('V/50/catalog')
    
    if not catalogs:
        print("Catalog not found.")
        return None
    
    hr_table = catalogs[0]
    
    # 3. Clean the data (remove entries without magnitudes or coordinates)
    mask = ~hr_table['Vmag'].mask & ~hr_table['_RAJ2000'].mask
    clean_table = hr_table[mask]
    
    return clean_table



def loadBrightStars(obsLat=-24.,infile=None, sidTime=0.0):
    """
    Load star catalog from file.
    
    Parameters:
    infile (str): Path to the star catalog file.
    sidTime (float): Sidereal time in hours.
    
    Returns:
        astropy.table.Table: Table containing star data with columns for azimuth (radians), zenith distance, magnitude, dot size, and illumination status.
        Table is compatible with satellite data tables.
    """



    if infile == "local":
        import pathlib
        infile = str(pathlib.Path(__file__).parent.resolve() / 'hr.dat')
        starData = ascii.read(infile)

    elif infile is None:
        starData = get_yale_bright_stars()
        starData.rename_column('_RAJ2000', 'RA')
        starData.rename_column('_DEJ2000', 'Dec')
        starData.rename_column('Vmag', 'mag')


    print(f'Loaded {len(starData)} stars from {infile}')
    print(f'faintest star magnitude: {starData["mag"].max():.2f}')

    starData["dot"] = satDots.magToDotSize(starData["mag"]+5.)
    starData["Az"], starData["El"] = caLib.radec2azel(sidTime*15.-starData["RA"], starData["Dec"], obsLat)

    starData["Azr"] = np.radians(starData["Az"])
    starData["ZD"]  = 90. - starData["El"]
    starData["mag"] = starData["mag"]
    starData["dot"] = starData["dot"]
    starData["bIlluminated"] = True  # stars are always illuminated
    return starData[ (starData["ZD"] < 90.0) ]  # only visible stars

#==============================================================================
# LABELS and STUFF
#==============================================================================
def sci_to_label(x):
    """Convert scientific notation number to LaTeX label string."""
    exponent = int(np.floor(np.log10(x)))
    coefficient = x / 10**exponent
    if np.isclose(coefficient, 1.0):
        return f'$10^{{{exponent}}}$'
    else:
        return f'${coefficient:.1f} ~ 10^{{{exponent}}}$' 