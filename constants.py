'''
 ConAn
 Constellation Analytic simulations

 constants.py: a series of useful constants and configurations
'''


#magSun = -26.75  # V magnitude of the Sun
#magSky = 21.7    # V magnitude of the sky, mag/sq.arcsec, patat08
#au = 149597870.7 # Astronomical unit in km

G             = 6.67430E-11 # Universal Gravity constant
gravityMu     = 398600.5  # Earth gravitation constant [km^3 / s^2]
earthMass     = 5.972E24  # mass of the Earth kg
earthRadius   = 6378.137  # Earth radius [km]
earthRotation = 7.292114992e-5 # Earth rotation rate [rad / s]
extinction    = 0.12 # extinction in V  [mag/airmass]

mag550 = 7. # magnitude of a satellite at zenith at 550km 
#mag550 = 5.7 # magnitude of a satellite at zenith at 550km - 5.7 = 7.0 at 1000km
#mag550 = 4.7 # magnitude of a satellite at zenith at 550km - 5.7 = 6.0 at 1000km
    

ZTFmagAzCut = 40. #  deg, azimuth range in which sun brightens sat - ZTF
ZTFmagAzBright = 0 #  -1.5 # brightening of the sat (i.e must be NEGATIVE) - ZTF
magAngPeak = 45. # deg - angle Sat-Sun where the brightening is strongest
