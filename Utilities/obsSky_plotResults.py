import json
import matplotlib.pyplot as plt
import numpy as np

import SatConAnalytic.conanplot as cpLib
import SatConAnalytic.skyBrightnessLib as skyLib

plt.rcParams.update({'font.size': 16})

# line arguments for plotting results from batch runs
# define the input json
input_json = 'data_skyScattered_SLOWGWAK_0.json'
input_json = 'data_skyScattered_AST_0.json'
input_json = 'data_skyScattered_ASTS_0.json'


# load the skyDiffuse to get the sun elevation values
myFiles = [
    'data_trails_SLOWGWAK_-23.json',
    'data_trails_SLOWGWAK_0.json',
    'data_trails_SLOWGWAK_23.json'
]

coefsList = []
for f in myFiles:
    with open(f, 'r') as f:
        results = json.load(f)
        print(f"Loaded {f} with {len(results)} entries")
        elevs = np.array([ r['sunElev'] for r in results ])
        times = np.array([ r['local_time'] for r in results ])
        times[times < 6] += 24

        plt.plot(times, elevs, '.-', label=f"Sun Delta {results[0]['sunDelta']}°")

        # find a function that gives elevation as time, and plot it
        # from the data, it looks like a linear function is a good fit
        coeffs =    np.polyfit(times, elevs, 5)
        fit_elevs = np.polyval(coeffs, times)
        plt.plot(times, fit_elevs, 'r-')
        coefsList.append(coeffs)

plt.xlabel('Local Time [h]')
plt.ylabel('Sun Elevation [deg]')
plt.title('Sun Elevation vs Local Time for Different Sun Deltas')
plt.grid()
#plt.show()





# main plotting code for skyScattered results

# Load results from JSON file
with open(input_json, 'r') as f:
    results = json.load(f)


# generic
telescope = results[0]['meta']['telescope']
deltaSun = results[0]['sunDelta']
unit = results[0]['value']['unit']

# extract data for plotting
scatter_zenith = [ r['value']['zenith'] for r in results ]
scatter_zenith += [0.]
scatter_30     = [ r['value']['sun30']  for r in results ]
scatter_30     += [0.]
sunElev = [ r['sunElev'] for r in results ]
sunElev += [-89.]  # add a value 


fig = plt.figure(figsize=(14, 6))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 0.04, 0.04, 0.04], hspace=.8)

AX = fig.add_subplot(gs[0])
AX.plot(sunElev, scatter_zenith, 'r:', label='Zenith')
AX.plot(sunElev, scatter_30,     'r-', label='30°')
ymin, ymax = AX.get_ylim()


# TWILIGHTS
for twilight in [0, -6, -12, -18]:
    AX.fill_betweenx([ymin, ymax], 5, twilight, color='blue', alpha=0.2)
    AX.text(twilight, ymax * 0.9, f'{twilight}°', rotation=90, verticalalignment='center')

# TICKS
newUnit = 'microCandelaPerM2'
newUnitLabel = '$\mu$cd/m$^2$'       # raw string for LaTeX
newTicks = np.linspace(0,ymax,5)
assert unit == 'nL', f"Expected input unit to be nL, but got {unit}"
newTickValues = skyLib.lambert_to_microCandelaPerM2(newTicks * 1e-9)
newTickLabels = [ cpLib.format_tick(v) for v in newTickValues]
AX.set_yticks(newTicks)

# make a copy of AX for secondary y-axis
newTickValues2 = skyLib.lambert_to_skyFraction(newTicks * 1e-9)
newTickLabels2 = [ cpLib.format_tick(v) for v in newTickValues2]

AX2 = AX.twinx()
AX2.set_ylim(ymin, ymax)
AX2.set_yticks(newTicks)
AX2.set_yticklabels(newTickLabels2)
AX2.set_ylabel(f'Sky Fraction')




AX.set_xticks( np.arange(6, -90.1, -6))
AX.set_xticklabels( [f'{t:.0f}' for t in np.arange(6, -90.1, -6)] )
AX.set_xlim(5, -90)
AX.text(1.0, 0.1, 'Sun Elevation [deg]',
        verticalalignment='top', horizontalalignment='right',
        transform=AX.transAxes)
AX.set_ylabel(f'Sky Scattered Brightness\n[{newUnitLabel}]')
AX.legend(loc='lower left')


# SECONDARY AXES
time_ticks = np.arange(18, 24.1, 1)
labels = ['Summer', 'Equinox', 'Winter']

for iax, coeffs in enumerate(coefsList):
    ax = fig.add_subplot(gs[iax + 1])

    # hide everything except the bottom spine (= the x axis line)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.set_yticks([])
    ax.patch.set_visible(False)   # transparent background

    # match the x range of the main axis
    ax.set_xlim(5, -90)
    ax.set_ylim(0, 1)  # arbitrary y limits since we won't show y ticks


    # compute the sun-elevation values that correspond to each hour tick
    fit_elevs = np.polyval(coeffs, time_ticks)
    ax.set_xticks(fit_elevs)
    ax.set_xticklabels([f'{t:.0f}' for t in time_ticks])
    ax.tick_params(axis='x', which='major', direction='out', length=8)
    ax.text(1.0, 0.0, labels[iax]+' Local Time [h]',
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes)

    ax.plot([fit_elevs[-1], -90], [.1, .1], 'k-')  
    AX.vlines(fit_elevs[-1], ymin, ymax, colors='gray', linestyles='dashed', alpha=0.5)
    AX.fill_betweenx([ymin, ymax], fit_elevs[-1],-95, color='grey', alpha=0.2)
    AX.text(fit_elevs[-1], ymax, labels[iax]+' midnight ', verticalalignment='top', horizontalalignment='right', rotation=90, color='gray', alpha=0.5, fontsize=12)

AX.set_ylim(ymin, ymax )  # add some vertical space for the twilight labels



#concatenate constellation names into a single string
c_str = ', '.join(results[0]['meta']['constellations']) 
n_sat = f'{results[0]["meta"]["const_stats"]}'
if results[0]['value']['type'] == "skyScattered":
    o_str = 'Scattered light'

AX.set_title(f'{c_str} ({n_sat} satellites): {o_str}', fontsize=16)

plt.tight_layout()
plt.show()

# outfile: replace json by pdf
outfile = input_json.replace('.json', '.pdf')
fig.savefig(outfile)
print(f"Saved plot to {outfile}")