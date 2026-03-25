#!/usr/bin/env python3
#
# SatConAnalytic - Satellite Constellation Analytic simulations

# Constellations.py:
#    Read the parameters of the Constellations
#    Assemble them in a Constellations object,
#    which contains a list of Constellation objects,
#    which each contains a list of Shell objects.
#


import logging, json, os
import numpy as np
from astropy.table import Table
import random

import SatConAnalytic.constants as cCst
import SatConAnalytic.conan as cLib


log = logging.getLogger('satcon')

#----------------------------------------------------------------------------
class _Dict(dict):
    '''convenience: access dict as dict.element.element'''
    __getattr__= dict.__getitem__
    __setattr__= dict.__setitem__
    __delattr__= dict.__delitem__

#---------------------------------------------------------------------------
def findConstellations(constellationsll, constFile="constellations.json"):
    '''Assemble a Constellations object (set of constellations)
    for a list of constellations.

    in: list of constellations ['SL1', 'SL2', 'OWr2']  
        or one of the preset codes defined below

    out: a ConstellationS object
    '''

    log.info(f'looking for constellations with ID={constellationsll} in file {constFile}...')

    metaConstellations = {
        'SL': ['SL1', 'SL2'],
        'OW': ['OW2r'],
        'SLOW': ['SL1', 'SL2','OW2r'],
        'TODAY': ['YESTERDAY', 'TODAYconst'],
        'SLOWGWAK': [#'YESTERDAY',
                    'SL1','SL2',
                    'OW2r',
                    'GW',
                    'AK' ],
        'ALL': ['YESTERDAY',
                'SL1', 'SL2', 
                'OW2r', 
                'GW', 'AK', 'ESP']
    }

    if   constellationsll == 'list' :  
        print(readConstellations(constFile))
        print("=============================================================")
        print("Available preset constellation groups:")
        for c in metaConstellations:
            print(f'  {c}: {metaConstellations[c]}')
            
        exit(0)
    elif constellationsll in metaConstellations:
        constellationsll = metaConstellations[constellationsll]
    else:
        constellationsll = [ constellationsll ]



    return metaConstellation(constellationsll, constFile=constFile)
    

#----------------------------------------------------------------------------
def readConstellationFile(myFile):
     '''read constellation json file
     
     return dict.
     '''

     log.debug(f'Reading in constellation file: {myFile}')

     try:
          # Try to open the file as provided (local file or absolute path)
          with open(myFile) as infile:
               return json.load(infile, object_hook=_Dict)
     except FileNotFoundError:
          # If file not found, try to locate it in the same directory as this script
          script_dir = os.path.dirname(os.path.abspath(__file__))
          fallback_path = os.path.join(script_dir, "Data", myFile)
          with open(fallback_path) as infile:
               return json.load(infile, object_hook=_Dict)

#----------------------------------------------------------------------------
class OneShell():
     '''define a single shell'''
     def __init__(self, oneShellJS) -> None:

          for x in list(oneShellJS):
              self.__dict__[x] = oneShellJS[x]

          self.nSat = int( self.totSat / self.nPlane + 0.5)
          # original number nSat is not reliable.

          if not hasattr(self, 'mag550'):
               self.mag550 = cCst.mag550


     def __repr__(self):
          return f'{self.label}: {self.totSat} = {self.nPlane}*{self.nSat}sat, i={self.inc}'


     def orbitalElementsTable(self):
          '''create the orbital element table for each sat in the shell
          
          'anom0' : anomaly at t0 [rad],
          'node0' : asce.node [rad],
          'inc'   : inclination [rad],
          'alt'   : altitude [km],
          'omega' : angular velocity [rad/s],
          'mag550': abs.mag at 550km
          '''

          # node = ascending nodes of the planes
          nodes =  np.linspace(0.,360., self.nPlane, endpoint=False) # [deg]
          satAnom0 = np.array([0.])  #- initialize some empty vectors
          satNode0 = np.array([0.])
          satInc   = np.array([0.])
          satAlt   = np.array([0.])
          satomega = np.array([0.])


          # create the planes
          for inode, node in enumerate(nodes): # cretes the various planes


               wsatAnom0 = np.linspace(0.,360., self.nSat, endpoint=False) + (20*inode)%(360./self.nSat)
                    # +: so they don't all start at the same place

               wsatInc    = np.full_like(wsatAnom0, self.inc) # [deg]
               wsatNode0  = np.full_like(wsatAnom0, node) # [deg]
               wsatAlt    = np.full_like(wsatAnom0, self.alt) # [deg]
               wsatomega =  np.full_like(wsatAnom0, 
                                         np.sqrt(cCst.G*cCst.earthMass / (1000.*(cCst.earthRadius + self.alt))**3)) #ang.vel, [rad/sec] #r must be in m

               #- and append the results for this plane to the main vectors
               satAnom0 = np.concatenate([satAnom0,   np.radians(wsatAnom0)])      # [RADIANS]
               satNode0  = np.concatenate([satNode0,  np.radians(wsatNode0)])      # [RADIANS]
               satInc    = np.concatenate([satInc,    np.radians(wsatInc)])        # [RADIANS]
               satAlt    = np.concatenate([satAlt,    wsatAlt])                    # [km]
               satomega  = np.concatenate([satomega,  wsatomega])                  # [RADIANS/s]


          self.Table = Table(
               {'anom0' : satAnom0,
                'node0' : satNode0,
                'inc'   : satInc,
                'alt'   : satAlt,
                'omega' : satomega,
                'mag550': np.full(len(satAnom0), self.mag550)    
               }
          )

          log.debug(f'Shell table: {len(self.Table)} -> {self.nSat}sat * {self.nPlane}planes = {self.nSat * self.nPlane} =?= {self.totSat}')
          return self.Table
     

     def modelOneShell(self, AzEl, obsLatitude, sunAlpha,sunDelta):
          '''Model one single shell over a set of Az,El pointings
          IN
          - self: parameters of the satellite constellation shell:
            - satInc: inclination [deg]
            - satAlt: altitude [km]
            - num: number of satellites in the shell
          - AzEl: mesh of [Azimuths,  Elevation]   [deg]
                  on which the constellation shall be evaluated. 
          - obsLatitude: latitude of the observer [deg]
          - sunAlpha, sunDelta: HourAngle and Dec. of the Sun [deg]
          - 

          OUT
          - illuminated satellite number density (same shape as AzEl)
          - illuminated satellite apparent angular velocity (same shape as AzEl)
          - illuminated satellite magnitudes
          '''

          # reshape - some scripts use 2D meshes, some use 3D meshes.
          if len(AzEl.shape) == 3:
               AzElreshape = np.reshape(AzEl,(2,AzEl.shape[1]*AzEl.shape[2]))
               step = AzEl[1,1,0] - AzEl[1,0,0]
          else:
               AzElreshape = AzEl
               step = 1.
               
          # geocentric equ. alpha,delta of sat, and   observatory dist, angle 
          alpha, delta, Delta, costheta = cLib.AltAz2Delta(
               obsLatitude,self.alt,AzElreshape)
          

          # geocentric equ. rect. of satellite
          xyz = cLib.RaDecAlt2xyz(alpha,delta, self.alt)
          

          # Velocities
          wAngularVel = cLib.AzEl2Vel(alpha, delta, Delta,obsLatitude,self.inc,self.alt)  
          
          #Density  at this place on the shell 
          #         adjust angular size
          #         for distance, and
          #         for apparent orientation of the shell on line-of-sight

          numDensity_Shell = cLib.satNumDensity_0(delta,self.inc,self.totSat)
                    # in num/sq.deg on the shell

          areaRatio = (Delta/(cCst.earthRadius+self.alt))**2 / costheta 
               # area on the shell projected on the sky
               # at Zenith: 
               #   areaRatio = alt/(rsat) **2, 
               # = 0.00630 for alt=550km. Verified OK.
               # at horizon: 
               #   sin(theta) = re/rs -> theta = 67dg -> cos(t) = 0.39
               #   Delta2= rs2-re2 -> Delta = 2705 -> D2/rs2 = 0.15
               #   ratio = Delta2/re2 / cos(theta) = 0.3905

          numDensity_Sky = numDensity_Shell * areaRatio


          # Illuminated satellites
          illum = cLib.solIllum(xyz,sunAlpha, sunDelta)
          numDensity_illuminated = numDensity_Sky * illum
          
          #  MAGNITUDE of the satellites:
          
          mag_visual =  self.mag550 + 5.*np.log10(Delta/550.)   # distances
          mag_visual += cCst.extinction*(Delta/self.alt -1.)    # extinction
          

          return \
               np.reshape(numDensity_illuminated,
                           (AzEl.shape[1],AzEl.shape[2]) ) ,\
               np.reshape(wAngularVel, (AzEl.shape[1],AzEl.shape[2]) ),\
               np.reshape(mag_visual,        (AzEl.shape[1],AzEl.shape[2]) )


#----------------------------------------------------------------------------
class Constellation():
     '''define a single constellation

     which is made of one of more shells
     '''
     def __init__(self, oneConstJS):
          for x in list(oneConstJS):
              if x == 'shells':
                   self.shells = [
                        OneShell( oneConstJS.shells[s])
                        for s in   oneConstJS.shells]
              else:
                   self.__dict__[x] = oneConstJS[x]

          self.totSat =  sum( s.totSat for s in self.shells)
          self.totShells = len( self.shells )
          self.ToC = (
               f'\n"{self.name}"\t {self.totSat} sat, {self.totShells} shells'+
               f'\n-----------------------------------------------------------'
               )

          self.vintageTable = [[
                         s.label,
                         s.totSat,
                         s.nPlane,
                         s.nSat,
                         s.inc,
                         s.alt
                         ] for s in self.shells ]


     def __repr__(self):
          msg =  '\n======================================================='
          msg += self.ToC
          for s in self.shells:
               msg += f'\n    {s.label}: \tN= {s.totSat}  '
               msg += f'\ta= {s.alt}km \ti= {s.inc}deg'
          return msg

#----------------------------------------------------------------------------
class Constellations():
     '''define a metaConstellation, list of constellations
     Produced either from
     - an input file or
     - a list of Constellation objects
     '''

     def __init__(self,constJS):
          '''constJS is either a list of
          - JS definitions of Constellations, or
          - Constellation objets'''

          self.list = [c.code  for c in constJS]
          self.name = ",\n".join(self.list)


          self.byCode      = _Dict( { })
          for c in constJS:
               if type(c) == type(_Dict({})):
                    self.byCode[c.code] = Constellation(c)
               else:
                    self.byCode[c.code] = c

          self.totSat = sum(  self.byCode[c].totSat for c in self.list   )
          self.totShells = sum(  self.byCode[c].totShells for c in self.list   )
          self.totConst = len(self.list)
          
          self.ToC =  '\n'.join( [ f'{c} :\t {self.byCode[c].ToC} ' for c in self.list ])
          self.ToC += f'\nTotal N={self.totSat} satellites'
          self.ToC += f'\nover  S={self.totShells} shells'
          self.ToC += f'\nin    C={self.totConst} Constellations'

          #make vintage table
          self.vintageTable = []
          for c in self.list:
               self.vintageTable +=  self.byCode[c].vintageTable

          # consolidate shells
          self.shells = []
          for c in self.list:
               self.shells += self.byCode[c].shells



     def __repr__(self):
          return self.ToC



#----------------------------------------------------------------------------
def readConstellations( constFile='constellations.json'):
     '''Reads a constellation json file into a Constellations object'''

     log.debug(f'file: {constFile}')
     allConst = readConstellationFile(constFile)
     return   Constellations( allConst)

#----------------------------------------------------------------------------
def metaConstellation( cList, constFile='constellations.json'):

     log.debug(f'list: {cList} and file: {constFile}')

     allConstellations = readConstellations(constFile=constFile)

     if cList[0] == 'list':
         return allConstellations


     try:
          return Constellations( [ allConstellations.byCode[c] for c in cList])
     except KeyError:

          cError = [c for c in cList if c not in allConstellations.list]
          log.error(f'{cError} not in constellation list')
          log.debug( allConstellations.list )
          log.error(f'metaConstellation: {cError} not in constellation list: {allConstellations.list}')
          raise ValueError(cError)

if __name__ == "__main__":
     print(readConstellations()  )
