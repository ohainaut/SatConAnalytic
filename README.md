# SatConAnalytic

## Usage

### ObsPlot

Generates a map-of-the-sky

```
obsplot -d deltasun -a alphasun -e elevsun -l latitude*
    -C constellation: YESTURDAY TODAY SL2OW2 SLOW2 OW2r OW2 SLOWr(def)
    -T FlyEye WFI VST EFOSC FORSimg LSST HAWKI MICADO
       FORSspec UVES1h 4MOST ESPRESSO
       VISIR 
       skyMag TrailDens SatDens(def)
          These pre-load default 
                extptime, FoVl, FoVw, maglim, magbloom, trailf, 
                telname, instrName, latitude
    -t exptime* -f Fovl* -m maglim* -k trailf* -r resolution*
              maglim: those with eff.mag fainter are not observable
                      those brigher create a trail
              trailf: fraction of the exposure that is destroyed by trail. 1=full
              magbloom: for those brighter, trailf=1=full
    -s telName* -i instrName*
                (*: overwrites presets)
    -S do not shade low elev
    -L do not label the plot
    -B no scalebar
    -M BRIGHT FAINT ALL OBS EFFECT
         All sats, restrict to bright/faint (wrt magbloom), convert into effect on obs.
```

### ObjPlot

Generates a visibility calendar for an object.

```
objplot 
    -r RA -d dec -n name
    -T FlyEye WFI VST EFOSC FORSimg LSST HAWKI MICADO
       FORSspec UVES1h 4MOST ESPRESSO
       VISIR 
       skyMag TrailDensity SatDensity(def)
          These pre-load default 
                extptime, FoVl, FoVw, maglim, magbloom, trailf, 
                telname, instrName, latitude
    -l latitude*
    -t exptime* -f Fovl* -m magbloom* -k trailf* 
              maglim: those with eff.mag fainter are not observable
                      those brigher create a trail
              trailf: fraction of the exposure that is destroyed by trail. 1=full
              magbloom: for those brighter, trailf=1=full
    -s telName* -i instrName*
                (*: overwrites presets)
    -S do not shade low elev
    -L do not label the plot
    -M BRIGHT FAINT ALL OBS EFFECT
         All sats, restrict to bright/faint (wrt magbloom), convert into effect on obs.
```

## Files

- conan.py	main functions for ConAn
- conanplot.py	plotting functions for ConAn
- constants.py	generic constants for ConAn
- constellations.py	definition of the constellations for ConAn

- obsplot.py	main interface to plot the sky over an observatory
- objplot.py	main interface to plot an object's calendar
- eleccount.py	hack to measure/plot the number of sat as a
		function of elevation
- sequenceplot.py examples of sequential calls
