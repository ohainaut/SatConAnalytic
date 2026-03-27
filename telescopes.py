#!/usr/bin/env python3
#
# SatConAnalytic - Satellite Constellation Analytic simulations

# telescopes.py:
#    Read the parameters of the telescopes and instruments

# telescopes.readTelescope( JSfile ) returns a list of Telescope objects

import logging
import json
import os
log = logging.getLogger('satcon')

from SatConAnalytic.utils import _Dict

def readTelescopeFile(myFile):
     '''read telescope json file'''
     try:
          # Try to open the file as provided (local file or absolute path)
          with open(myFile) as infile:
               tels =  json.load(infile, object_hook=_Dict)
          log.info(f"telescope definitions from {myFile}")
     except FileNotFoundError:
          # If file not found, try to locate it in the same directory as this script
          script_dir = os.path.dirname(os.path.abspath(__file__))
          fallback_path = os.path.join(script_dir, "Data", myFile)
          with open(fallback_path) as infile:
               tels =  json.load(infile, object_hook=_Dict)
          log.info(f"telescope definitions from fallback {fallback_path}")
     return tels



class Telescope():
    '''define a single telescope'''
    def __init__(self, oneTelJS) -> None:
        for x in list(oneTelJS):
            self.__dict__[x] = oneTelJS[x]

        for what in ['resol', 'fovl', 'fovw']:
            if what+'_arcsec' in self.__dict__:
                self.__dict__[what] = self.__dict__[what+'_arcsec']/3600.



        if 'fovl' not in self.__dict__:
             self.fovl = 1./3600 # arcsec
        if 'fovw' not in self.__dict__:
             self.fovw = self.fovl*1.

        if 'resol' not in self.__dict__:
             self.resol = 1./3600. # deg


        if 'trail_arcsec' in self.__dict__:
             self.trailf = self.trail_arcsec/3600./self.fovl
        else:
            if 'trailf' not in self.__dict__:
                self.trailf = 1./3600./self.fovl



        if 'maglim' not in self.__dict__:
             self.maglim = 99. # detect everything

        if 'magbloom' not in self.__dict__:
             self.magbloom = -99. # saturation never a problem

        if 'lat' not in self.__dict__:
             self.lat = -24.62 # VLT...

        if 'expt' not in self.__dict__:
             self.expt = 1. #[s]

        if 'telescope' not in self.__dict__:
             self.telescope = " "

        if 'instrument' not in self.__dict__:
             self.instrument = " "


        for what in ['resol', 'fovl', 'fovw']:
            self.__dict__[what+'_arcsec'] = self.__dict__[what]*3600.
        
        self.trail_arcsec = 3600.*self.fovl*self.trailf
    
        self.ToC = f'{self.code}: {self.telescope} + {self.instrument}'

    def __repr__(self) -> str:
        msg  = f'{self.code}:  '
        msg += f'{self.telescope} {self.instrument}\n'
        msg += f'\tLatitude: \t{self.lat:.1f} deg\n'
        msg += f'\tLongitude: \t{self.lon:.1f} deg\n'
        msg += f'\tExp.time: \t{self.expt} s\n'
        if self.fovl > 0.1 :
            msg += f'\tFoV:    \t{self.fovw:.2f} x {self.fovl:.2f} deg\n'
        else:
            msg += f'\tFoV:    \t{self.fovw_arcsec} x {self.fovl_arcsec} arcsec\n'
        msg += f'\tResolution: \t{self.resol_arcsec:.2f} arcsec\n'
        msg += f'\tTrail width: \t{self.trail_arcsec:.2f} arcsec\n'
        msg += f'\t           = \t{self.trailf*100.:.3f}% of FoV\n'
        msg += f'\tLimiting mag: \t{self.maglim}\n'
        msg += f'\tBloom/sat mag: \t{self.magbloom}\n'
        return msg

class Telescopes():
    '''list of telescope setup'''

    def __init__(self,telJS):
          '''telJS : JS definitions of telescopes
          '''
          self.list = [c.code  for c in telJS]
          self.byCode      = _Dict( { })
          for t in telJS:
                self.byCode[t.code] = Telescope(t)

    def __repr__(self) -> str:
        msg = 'List of telescopes and instruments:\n'
        for t in self.list:
             msg += f'{self.byCode[t]}'
        return msg
         


def readTelescopes( file='telescopes.json'):
     allTel = readTelescopeFile(file)
     return  Telescopes( allTel )


#------------------------------------------------------------------------------
def findTelescope(telinslabel):
     '''Find and return a Telescope object by its code label.'''
     allTel = readTelescopes()
     try:
          return allTel.byCode[telinslabel]
     except KeyError:
          if telinslabel != 'list':
               print(f'{telinslabel} not found in telescope list')
          print('Available telescopes are:')
          for x in sorted(allTel.list):
               print(f'  {x}')
          exit(1)


#------------------------------------------------------------------------------
def getTelescope(myargs):
     '''Return the Telescope object requested by myargs.

     The telescope is defined by the ``code`` argument; individual parameters
     in myargs then overload the presets.  myargs is an argparse Namespace.
     '''

     if myargs.code is None:
          myargs.code = 'DEFAULT'

     myTel = findTelescope(myargs.code)

     for what in ['telescope', 'instrument']:
          if myargs.__dict__[what] is not None:
               myTel.__dict__[what] = myargs.__dict__[what]

     if myargs.expt is not None:
          myTel.expt = float(myargs.expt)
     if myargs.fovl is not None:
          myTel.fovl = float(myargs.fovl)
     if myargs.fovw is not None:
          myTel.fovw = float(myargs.fovw)

     if myargs.fovw is not None:
          myTel.fovw = float(myargs.fovw)
     else:
          myTel.fovw = myTel.fovl * 1.

     for what in ['expt', 'fovl', 'magbloom', 'maglim', 'resol', 'trailf', 'lat']:
          if what in myargs.__dict__ and myargs.__dict__[what] is not None:
               myTel.__dict__[what] = float(myargs.__dict__[what])

     return myTel


if __name__ == "__main__":
     mytels = readTelescopes()
     for t in mytels.list:
          print( mytels.byCode[t] )