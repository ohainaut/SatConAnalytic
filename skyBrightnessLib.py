import numpy as np
from astropy.table import Table
from astropy.io import ascii

import SatConAnalytic.satDots as satDots
import SatConAnalytic.conan as caLib

# CONVERSIONS

def nanoLambert_to_microCandelaperm2(nanoLambert):
    """
    Convert brightness from nanoLamberts to microCandelas per square meter.
    
    Parameters:
    nanoLambert (float or np.ndarray): Brightness in nanoLamberts.
    
    Returns:
    float or np.ndarray: Brightness in microCandelas per square meter.
    """
    return nanoLambert / 3.141592653e-1


def microCandelaperm2_to_nanoLambert(microCandelaperm2):
    """
    Convert brightness from microCandelas per square meter to nanoLamberts.
    Ref: 1 candela/square meter = 	0.3141592654 mL 

    Parameters:
    microCandelaperm2 (float or np.ndarray): Brightness in microCandelas per square meter.
    
    Returns:
    float or np.ndarray: Brightness in nanoLamberts.
    """
    return microCandelaperm2 * 3.141592653e-1


def magarc_to_nanoLambert(mag_arcsec2):
    """
    Convert sky brightness from mag/arcsec^2 to nanoLamberts.

    Parameters:
    mag_arcsec2 (float or np.ndarray): Sky brightness in mag/arcsec^2.

    Returns:
    float or np.ndarray: Sky brightness in nanoLamberts.

    from K&S91, with log10(e)*20.7233 = 9. and log10(3)*0.92104 = 0.4
    """
    return 34.08e9 * 10**(  - 0.4 * mag_arcsec2) 


def nanoLambert_to_magarc(nanoLambert):
    """
    Convert sky brightness from nanoLamberts to mag/arcsec^2.
    Parameters:
    nanoLambert (float or np.ndarray): Sky brightness in nanoLamberts.
    Returns:
    float or np.ndarray: Sky brightness in mag/arcsec^2.
    """
    return -2.5 * np.log10( nanoLambert ) +26.33


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

def magnitude_to_illuminance(magnitude):
    """
    Convert magnitude to illuminance in lux.
    
    Parameters:
    magnitude (float or np.ndarray): Apparent magnitude.
    
    Returns:
    float or np.ndarray: Illuminance in lux:  - 14.18
    Illuminance in foot-candles: + 16·57)
    """
    return 10**(-0.4 * (magnitude + 16.57))


# TRIGONOMETRY / AIRMASS

def airmassSimple(Z):
    """
    Airmass as a function of Zenithal Distance
    """

    return 1./np.cos(np.radians(Z))

def airmassFull(Z):
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



# Dark Sky brightness

def skyBrightness(Z, Bz=magarc_to_nanoLambert(22.) , k=0.125) :
    """
    Natural Sky brightness at a given zenithal distance Z, scaled from zenithal brightness Bz

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

# Scattering

def scattering_Rayleigh(theta):
    """
    Rayleigh scattering function.
    
    Parameters:
    theta (float or np.ndarray): Angular distance in degrees.
    
    Returns:
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
    magnitude : Apparent magnitude of the light source out of atmosphere
    theta : Angular distance from the light source in degrees
    Z : Zenithal Distance in degrees
    k : extinction coefficient in mag/airmass (default 0.125)

    Out
    Sky brightness at zenith in nanoLamberts (default 22 mag/arcsec2)
    
    """

    X = airmassFull(Z)
    IlluminanceOutOfAtmos = magnitude_to_illuminance(magnitude) # Eq.[8]
    Illuminance = IlluminanceOutOfAtmos * 10**(-0.4*k*X) # Extinction Eq.[7]

    SurfBrightness = Illuminance * scattering_total(theta) 
    SurfBrightness *= ( 1. - 10.**(-0.4*k*X) )  # Eq.[13]->[12] = [15]

    return SurfBrightness

def Scattered_brightness_magarc( magnitude, theta, Z, k=0.125) :

    skyBrightness = nanoLambert_to_magarc(
        Scattered_brightness_nanoLambert(magnitude, theta, Z, k))

    return skyBrightness 


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

    starData["dot"] = satDots.magToDotSize(starData["mag"]+5.)
    starData["Az"], starData["El"] = caLib.radec2azel(sidTime*15.-starData["RA"], starData["Dec"], obsLat)

    starData["Azr"] = np.radians(starData["Az"])
    starData["ZD"]  = 90. - starData["El"]
    starData["mag"] = starData["mag"]
    starData["dot"] = starData["dot"]
    starData["bIlluminated"] = True  # stars are always illuminated
    return starData[ (starData["ZD"] < 90.0) ]  # only visible stars


