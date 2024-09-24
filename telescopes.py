#!/usr/bin/env python3
#
# SatConAnalytic - Satellite Constellation Analytic simulations

# telescopes.py:
#    Read the parameters of the telescopes and instruments

# telescopes.readTelescope( JSfile ) returns a list of Telescope objects

import logging
import json

class _Dict(dict):
    '''convenience: access dict as dict.element.element'''
    __getattr__= dict.__getitem__
    __setattr__= dict.__setitem__
    __delattr__= dict.__delitem__

def readTelescopeFile(myFile):
     '''read constellation json file'''
     with open(myFile) as infile:
          return json.load(infile, object_hook=_Dict)



class Telescope():
    '''define a single shell'''
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
    
        self.ToC = f'{self.code}: \t{self.telescope} {self.instrument}'

    def __repr__(self) -> str:
        msg  = f'{self.code}:  '
        msg += f'{self.telescope} {self.instrument}\n'
        msg += f'\tLatitude: \t{self.lat:.1f} deg\n'
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


if __name__ == "__main__":
     mytel = readTelescopes()
     for t in mytel.list:
          print( mytel.byCode[t].ToC )
