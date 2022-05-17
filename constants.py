# SatAn: constants
# from SatAn.constants import magSun, mu, re, we, au, extm
global magSun, mu, re, we, au, extm
magSun = -26.75  # V
magSky = 21.7 # V, patat08
mu = 398600.5 # Earth gravitation constant [km^3 / s^2]
re = 6378.1371 # Earth radius [km]
we = 7.292114992e-5 # Earth rotation rate [rad / s]
au = 149597870.7 # Astronomical unit in km
extn = 0.12 # extinction in V  

mag550 = 5.7 # magnitude of a satellite at zenith at 550km
magAzCut = 40. #deg, azimuth range in which sun brightens sat - ZTF
magAzBright = 0 #   -1.5 # brightening of the sat (i.e must be NEGATIVE) - ZTF
magAngPeak = 45. #deg - angle Sat-Sun where the brightening is strongest
