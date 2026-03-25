import numpy as np
import matplotlib.pyplot as plt
from importlib import reload  
import SatConAnalytic.skyBrightnessLib as SBlib
reload(SBlib)

darkSky_mag_arcsec2=22.0
moonSky_mag_arcsec2=18.0


# Test and validate magnitude to illuminance conversions

reload(SBlib)
magnitude = 0.0
luxorig = SBlib.magnitude_to_illuminance(magnitude)
print(f'{luxorig:.3g} lux for mag 0 star ORIG ')

fc = SBlib.magnitude_to_illuminanceFootCandles(magnitude)
print(f'({fc:.3g} footcandle for mag 0 star REF)')

lux = SBlib.illuminance_footcandle_to_lux(fc)
print(f'{lux:.3g} { fc * 10.764:.3g} lux for mag 0 star from footcandle conversion')

print('-----------------------------------------------')


# test and validate flux and lux

reload(SBlib)
#lux to W/m2
def lux_to_Wm2(lux):
    """Convert illuminance [lux = lm/m²] to flux [W/m²].
    Scaled from Sun: 1.27E+05 lux -> 164.35 W/m²
    1 lux = 1/773 W/m²
    1 W/m² = 773 lux
    """
    return lux / 773.

def Wm2_to_lux(flux):
    """Convert flux [W/m²] to illuminance [lux = lm/m²]."""
    return flux * 773.

magnitude = 7#-26.75 # sun
print(f'magnitude:                {magnitude}')

luxorig = SBlib.magnitude_to_illuminance(magnitude)
print(f'magnitude_to_illuminance: {luxorig:.5g} lux' )

fluxorig = SBlib.magnitude_to_flux(magnitude)
print(f'magnitude_to_flux:        {fluxorig:.5g} W/m² ')

flux = SBlib.illuminance_lux_to_Wm2(luxorig)
print(f'illuminance_lux_to_Wm2:   {flux:.5g} W/m² through lux' )

fc = SBlib.magnitude_to_illuminanceFootCandles(magnitude)
print(f'magnitude_to_illuminanceFootCandles: {fc:.5g} footcandle ')

lux = SBlib.illuminance_footcandle_to_lux(fc)
print(f'illuminance_footcandle_to_lux: {lux:.5g} lux from footcandle' )




# validation
#
assert np.isclose(SBlib.magarc2_to_luminance(21.83), 2.04e-4, rtol=1e-2) 
assert np.isclose(SBlib.luminance_to_magarc2(2e-4), 21.83, rtol=1e-2)
#Section 1.3 of Crumey, A. (2014). Human contrast threshold and astronomical visibility. MNRAS 442, 2600–2619

assert np.isclose(SBlib.magarc2_to_skyFraction(17), 100.)
assert np.isclose(SBlib.skyFraction_to_magarc2(100), 17.)


assert np.isclose(SBlib.magarc2_to_skyFraction(22.47), 0.65, rtol=1e-2)
# contribution of airglow to natural dark sky brightness 
assert np.isclose(SBlib.magarc2_to_luminance(22.0), 173.E-6, rtol=1e-2)  # dark sky luminance in cd/m2 - ref: Bara et al. 2016 RSocOpenSci 3
assert np.isclose(SBlib.magarc2_to_luminance(0.0), 1.098E5, rtol=1e-2) , SBlib.magarc2_to_luminance(0.0) # traditional value reported in Barra+20

print(SBlib.illuminance_lux_to_footcandle(1), " footcandle")
print(SBlib.illuminance_footcandle_to_lux(1), " lux")
i = SBlib.magnitude_to_illuminance(0)
f = SBlib.illuminance_lux_to_footcandle(i)
print("mag 0 star ->", i, " lux -> ", f, " footcandle", i/f, " lux/footcandle")
assert np.isclose(i/f, 10.764, rtol=1e-2)



print('dark sky: [mag/arcsec2]:', 22.)
print('Magarc to luminance:  [mucd.m-2]', 1e6*SBlib.magarc2_to_luminance(22))  
print("magarc to L",SBlib.magarc2_to_lambert(22))  
print("          nL",1E9*SBlib.magarc2_to_lambert(22))  

print('Magarc to nLambert: _direct_  [nL]', SBlib.magarc2_to_nanolambert_DIRECT(22))
print(SBlib.magarc2_to_nanolambert_DIRECT(22)/SBlib.magarc2_to_lambert(22))
print('magArc to Flux:  [muW.m-2.sr-1]', SBlib.magarc2_to_flux(22)*1e6)

# more validation

magarc2orig = 0. 
print(f'magarc2:                {magarc2orig} mag/arcsec²')


lux = SBlib.magarc2_to_luminance(magarc2orig)
print(f'magarc2_to_luminance:   {lux:.5g} cd/m² ')

wmarc2  = SBlib.luminance_to_magarc2(lux)
print(f'luminance_to_magarc2:   {wmarc2:.5g} mag/arcsec² was {magarc2orig} mag/arcsec² ')

lambert = SBlib.luminance_to_lambert(lux)
print(f'luminance_to_lambert:   {lambert:.5g} L ')

wlux = SBlib.lambert_to_luminance(lambert)
print(f'lambert_to_luminance:   {wlux:.5g} cd/m² was {lux:.5g} cd/m² ')

lambert2 = SBlib.magarc2_to_lambert(magarc2orig)
print(f'magarc2_to_lambert:     {lambert2:.5g} L was {lambert:.5g} L from lux')



wmarc3 = SBlib.lambert_to_magarc2(lambert2)
print(f'lambert_to_magarc2 (indir): {wmarc3:.3f} mag/arcsec² was {magarc2orig} mag/arcsec² ')

flux = SBlib.magarc2_to_flux(magarc2orig)
print(f'magarc2_to_flux:        {flux:.5g} W/m²/sr ')

wmarc4 = SBlib.flux_to_magarc2(flux)
print(f'flux_to_magarc2:        {wmarc4:.5g} mag/arcsec² was {magarc2orig} mag/arcsec² ')


magarc2 = 22
print(f'Surface brightness:  {magarc2} mag/arcsec²')

flux = SBlib.magnitude_to_flux(magarc2)
print(f'magnitude_to_flux:  {flux:.6g} W/m² / arcsec²')

srarcs = SBlib.sr_to_squarearcsec(1.)
print(f'sr_to_squarearcsec:  {srarcs:.6g} arcsec² / sr')

fluxsr = flux * srarcs
print(f'flux per sr:  {fluxsr:.6g} W/m² / sr')

flux2 = SBlib.magarc2_to_flux(magarc2)
print(f'magarc2_to_flux:  {flux2:.6g} W/m² / sr')

wmarc3 = SBlib.flux_to_magarc2(flux2)
print(f'flux_to_magarc2:  {wmarc3:.6g} mag/arcsec²')



magarc2 = 22
print(f'Surface brightness:  {magarc2} mag/arcsec²')

flux = SBlib.magnitude_to_flux(magarc2)
print(f'magnitude_to_flux:  {flux:.6g} W/m² / arcsec²')

srarcs = SBlib.sr_to_squarearcsec(1.)
print(f'sr_to_squarearcsec:  {srarcs:.6g} arcsec² / sr')

fluxsr = flux * srarcs
print(f'flux per sr:         {fluxsr:.6g} W/m² / sr')

flux2 = SBlib.magarc2_to_flux(magarc2)
print(f'magarc2_to_flux:      {flux2:.6g} W/m² / sr')

wmarc3 = SBlib.flux_to_magarc2(flux2)
print(f'flux_to_magarc2:      {wmarc3:.6g} mag/arcsec²')
print()


sarc2 = SBlib.sr_to_squarearcsec(1.)
print(f'sr_to_squarearcsec:     {sarc2:.6g} arcsec²/sr')
sarc3 =  1/ ( (np.pi / (180.0 * 3600.0)) ** 2 )
print(f'calc sr/arcsec²:        {sarc3:.6g} sr / arcsec²')
print()

flux3a = SBlib.solarFluxV * 10**(-0.4 * (magarc2 - SBlib.solarMagV))  
print(f'calc from solar flux :       {flux3a:.6g} W/m² / arcsec²')
flux3b = SBlib.magnitude_to_flux(magarc2)
print(f'calc from magnitude_to_flux:  {flux3b:.6g} W/m² /arcsec²')
print()


flux3 = SBlib.solarFluxV * 10**(-0.4  * (magarc2 - SBlib.solarMagV)) * (180.0 * 3600.0 / np.pi ) ** 2
print(f'manually :              {flux3:.6g} W/m² / sr')

flux4a = flux3a * sarc2
print(f'calc from steps:        {flux4a:.6g} W/m² / sr')    

flux4 = SBlib.magarc2_to_flux(magarc2)
print(f'magarc2_to_flux:        {flux4:.6g} W/m² / sr')

print(f' flux = {SBlib.solarFluxV} * 10**(-0.4  * (M_V - {SBlib.solarMagV})) * (180.0 * 3600.0 / np.pi ) ** 2')

k = SBlib.solarFluxV * (180.0 * 3600.0 / np.pi ) ** 2
print(f' constant k = {k:.6g} W/m² / sr')

flux5 = k * 10**(-0.4  * (magarc2 - SBlib.solarMagV))
print(f'final check:            {flux5:.6g} W/m² / sr   ')
