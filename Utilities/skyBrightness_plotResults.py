# INPUT EXPECTED IN nanoLambert

import json
import matplotlib.pyplot as plt
import numpy as np

import SatConAnalytic.skyBrightnessLib as skyLib
import SatConAnalytic.conanplot as cpLib

# line arguments for plotting results from batch runs
# define the input json
#input_json = 'skyBrightness_SunLight2.json'
input_json = 'skyBrightness_SLOWGWAK.json'

# Load results from JSON file
with open(input_json, 'r') as f:
    results = json.load(f)

# generic

CONST_name = results[0]["CONST_name"]
CONST_totsat = results[0]["CONST_totsat"]
CONST_range = results[0]["CONST_range"]
CONST_illum = results[0]["CONST_illum"]
unitLabel = results[0]["unitLabel"]

assert unitLabel in ['nanoLambert','nL'], f'Expected unitLabel to be nanoLambert, but got {unitLabel}'



# extract data for plotting
sunAlpha = np.array([ r["sunAlpha"] for r in results ])
sunDelta = np.array([ r["sunDelta"] for r in results ])
sunElev = np.array([ r["sunElev"] for r in results ])
CONST_range = np.array([ r["CONST_range"] for r in results ])
CONST_illum = np.array([ r["CONST_illum"] for r in results ])
skyBrightness_zenith = np.array([ r["skyBrightness_zenith"] for r in results ])
skyBrightness_30 = np.array([ r["skyBrightness_30"] for r in results ])

local_times = sunAlpha/15. + 12. 
local_times[ local_times < 12. ] = local_times[ local_times < 12. ] + 24.



# average brightness for measurements with same sunElev, to get a smoother curve for the fit
unique_sunElev = np.unique(sunElev)
skyBrightness_zenith_avg = np.array([ skyBrightness_zenith[sunElev == elev].mean() for elev in unique_sunElev ])
skyBrightness_30_avg = np.array([ skyBrightness_30[sunElev == elev].mean() for elev in unique_sunElev ])
sunAlpha_avg = np.array([ sunAlpha[sunElev == elev].mean() for elev in unique_sunElev ])
local_times_avg = sunAlpha_avg/15. + 12.
local_times_avg[ local_times_avg < 12. ] = local_times_avg[ local_times_avg < 12. ] + 24.

# stdev for measurements with same sunElev, to get an idea of the scatter
skyBrightness_zenith_std = np.array([ skyBrightness_zenith[sunElev == elev].std() for elev in unique_sunElev ])
skyBrightness_30_std = np.array([ skyBrightness_30[sunElev == elev].std() for elev in unique_sunElev ])


# get Local time of the twilights
twilightCrossings = []
for twilights in [0, -6, -12, -18]:
# interpolate local times where sunElev crosses the twilight levels
    for i in range(len(sunElev) - 1):
        if (sunElev[i] - twilights) * (sunElev[i + 1] - twilights) <= 0:
            # Linear interpolation to find the crossing time
            t1, t2 = local_times[i], local_times[i + 1]
            s1, s2 = sunElev[i], sunElev[i + 1]
            crossing_time = t1 + (twilights - s1) * (t2 - t1) / (s2 - s1)
            twilightCrossings.append((twilights, crossing_time))


# find maximum for plot:

nLmax = max(skyBrightness_zenith.max(), skyBrightness_30.max())
print(f'Max sky brightness: {nLmax:.2e} {unitLabel}')
nLmax = nLmax * 0.9
    


# main plot unit converstion (from nL to myUnit)
myUnit = "microCandelaPerM2"
myUnitLabel = "$\mu$Cd m$^{-2}$"

ymax = np.array([nLmax.copy()])

toBeDone = [skyBrightness_zenith_avg,
            skyBrightness_zenith_std,
            skyBrightness_30_avg,
            skyBrightness_30_std, 
            ymax]



if myUnit == 'nanoLambert':
    darkSkyLevel = 54. # nanoLambert
elif myUnit == 'frac':
    darkSkyLevel = 1.0 # fraction of natural sky brightness
    for what in toBeDone:
        what[:] = skyLib.lambert_to_skyFraction(what*1e-9) # convert from nanoLambert to fraction of natural sky brightness
elif myUnit == 'magarcsec2':
    darkSkyLevel = 22.0 # mag/arcsec^2
    for what in toBeDone:
        what[:] = skyLib.lambert_to_magarcsec2(what*1e-9) # convert from nanoLambert to mag/arcsec^2

elif myUnit == 'microCandelaPerM2':
    darkSkyLevel = skyLib.lambert_to_microCandelaPerM2(54.E-9) # microCandela/m^2
    for what in toBeDone:
        what[:] = skyLib.lambert_to_microCandelaPerM2(what*1e-9) # convert from nanoLambert to microCandela/m^2 
else:
    darkSkyLevel = None





# Plotting
# set default font size
plt.rcParams.update({'font.size': 20})

fig, axLeft = plt.subplots(figsize=(10, 6))

axLeft.set_ylim(0., ymax[0])
axLeft.set_ylabel(f'Sky Brightness [{myUnitLabel}]')



axLeft.plot(local_times_avg, skyBrightness_zenith_avg, '-', color='red')

axLeft.errorbar(local_times_avg, skyBrightness_zenith_avg, 
                yerr=skyBrightness_zenith_std, fmt='-', color='red', alpha=0.5, label='Zenith')


axLeft.plot(local_times_avg+0.051, skyBrightness_30_avg, ':', color='red')
axLeft.errorbar(local_times_avg+0.051, skyBrightness_30_avg, 
                yerr=skyBrightness_30_std, fmt=':', color='red', alpha=0.5, label='30$^o$')


axLeft.legend(loc='upper right')
axLeft.set_xlabel('Local Time (hours)')



# SKY LEVEL
ymin, ymax = axLeft.get_ylim()

# add a tick on the colorbar for the natural dark sky brightness level
for mySky, mySkyCol in zip(np.array([0.01, 0.1, 0.25, 1., 10., 100.]), ['orange', 'r','darkred','k']):
    #print(f'mySky={mySky}x{darkSkyLevel}= {mySky*darkSkyLevel}, ymin={ymin}, ymax={ymax}')
    if darkSkyLevel is not None and mySky*darkSkyLevel >= ymin and mySky*darkSkyLevel <= ymax:  
        axLeft.hlines(mySky*darkSkyLevel, 17, 24, colors=mySkyCol, linestyles='-')
        axLeft.text(23, mySky*darkSkyLevel, f'{mySky*100:.0f}%',  
                    va='bottom', ha='center', fontsize=12)


# TWILIGHTS
ymin, ymax = axLeft.get_ylim()

for twilights, crossing_time in twilightCrossings:

    if crossing_time < 24:
        plt.fill_betweenx([ymin, ymax], 16, crossing_time,  color='blue', alpha=0.2)  
    else:
        plt.fill_betweenx([ymin, ymax], crossing_time, 32, color='blue', alpha=0.2)
    plt.text(crossing_time, ymin + 0.95 * (ymax - ymin), f'{twilights}°', rotation=90, verticalalignment='center')

# tick marks above 24 labeled as 0-12
myTicks = np.arange(18, 32, 1)
plt.xticks(myTicks, [str(int(t) % 24) for t in myTicks])
plt.ylim(ymin, ymax)
plt.xlim(17.5, 24.)

# UNIT CONVERSION

axRight = axLeft.twinx()
# convert left ylimits to right limits: nanoLambert to fraction of natural sky brightness
ymax_right = skyLib.lambert_to_skyFraction(nLmax *1e-9) 
axRight.set_ylim(0., ymax_right)    

tick_values = axRight.get_yticks()
tick_labels =  [ cpLib.format_tick(val) for val in tick_values ] 
axRight.set_yticklabels(tick_labels)

axRight.set_ylabel('[x natural dark sky]')

#plt.title(CONST_name)
plt.tight_layout()
plt.show()
plt.savefig(f'skyBrightness_{CONST_name}_{int(sunDelta[0])}.pdf', dpi=300)  