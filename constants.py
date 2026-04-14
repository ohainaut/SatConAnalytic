# ConAn
# Constellation Analytic simulations
#
# constants.py: a series of useful constants and configurations
#

# from SatAn.constants import magSun, mu, re, we, au, extm

#global magSun, mu, re, we, au, extm

magSun = -26.71  # V magnitude of the Sun
magSky = 21.7    # V magnitude of the sky, mag/sq.arcsec, patat08

mu = 398600.5  # Earth gravitation constant [km^3 / s^2]
earthRad = 6378.1371 # Earth radius [km]
earthOmega = 7.292114992e-5 # Earth rotation rate [rad / s]
au_m = 149597870700. # Astronomical unit [m]
extn = 0.12 # extinction in V  [mag/airmass]

mag550 = 5.7 # magnitude of a satellite at zenith at 550km - 5.7 = 7.0 at 1000km
#mag550 = 4.7 # magnitude of a satellite at zenith at 550km - 5.7 = 6.0 at 1000km
    

ZTFmagAzCut = 40. #  deg, azimuth range in which sun brightens sat - ZTF
ZTFmagAzBright = 0 #  -1.5 # brightening of the sat (i.e must be NEGATIVE) - ZTF
ZTFmagAngPeak = 45. # deg - angle Sat-Sun where the brightening is strongest
