#!/usr/bin/env python3
# SatConAnalytic - Satellite Constellation Analytic simulations
# conan.py: definitions and generic constellation functions
#

import numpy as np
import random

# SatConAn: 
import constants 
import constellations

#---------------------------------------------------------------------------
def get_sun(jd):
    '''Compute coordinates of the Sun

    IN: JD full julian day
    OUT: RA and Dec (both degrees) of the Sun
    '''


    # fast sun
    n = jd - 2451545.0

    eps = 23.439 - 0.0000004*n
    epsr = np.radians(eps)

    # mean longitude
    L = (280.460 + 0.9856474*n)%360.

    # mean anomaly
    g = (357.528 + 0.9856003*n)%360.
    gr = np.radians(g)

    # ecl.long of sun:
    lambdas = L + 1.915*np.sin(gr) + 0.020*np.sin(2.*gr)
    lambdar = np.radians(lambdas)

    alpha = np.degrees(np.arctan2(np.cos(epsr)*np.sin(lambdar), np.cos(lambdar)))
    delta = np.degrees(np.arcsin(np.sin(epsr)*np.sin(lambdar)))

    return alpha, delta


#------------------------------------------------------------------------------
    
def findConstellations(constellationsll):
    '''Assemble a Constellations object (set of constellations)
    for a list of constellations.

    in: list of constellations ['SL1', 'SL2', 'OWr2']  
        or one of the preset codes defined below

    out: a Constellations object
    '''

    if   constellationsll == 'SL' :  constellationsll = ['SL1', 'SL2']
    elif constellationsll == 'OW' :  constellationsll = ['OW2r']
    elif constellationsll == 'TODAY':constellationsll = ['YESTURDAY', 'TODAYconst']
    elif constellationsll == 'SLOWGWAK':constellationsll = ['YESTURDAY',
                                                            'SL1','SL2',
                                                            'OW2r',
                                                            'GW',
                                                            'AK' ]
    elif constellationsll == 'ALL' :  constellationsll = ['YESTURDAY',
                                                          'SL1', 'SL2', 
                                                          'OW2r', 
                                                          'GW', 'AK', 'ESP']
    else:
        constellationsll = [ constellationsll ]

    return constellations.metaConstellation(constellationsll)
    
#---------------------------------------------------------------------------

def velPosAng(delta,satInc):
    '''Compute the velocity position angles for a list of satellites.

    in:
    - delta = lat, latitude of the satellite(s) [deg]
    - satInc: inclination of shell [deg]
    
    out: the two position angles (up and down) [deg]
    
    Note: deals properly with retrograde orbits (with satInc > 90)
    '''

    sintheta = np.cos(np.radians(satInc))/np.cos(np.radians(delta))
    theta1 = np.degrees(np.arcsin(sintheta))
    theta2 = np.degrees(np.arcsin(-sintheta))-180.*np.sign(satInc-90)
    return theta1, theta2
#---------------------------------------------------------------------------

def myarcsin(x):
    '''extended arcsin to any input value

    IN: x, the value from which arcsin must be computed, float, ]-4e4, 4e4[
    OUT: arcsin(x) [radians] 
    '''
    
    myx1 = np.where(x > 1., 1., x)
    myx2 = np.where(myx1 < -1., -1., myx1)
    myarcsin = np.arcsin(myx2)
    
    return myarcsin
#---------------------------------------------------------------------------

def satCount(l1,l2,inc,N):
    '''Number of satellites between two latitudes.

    IN:
    - l1, l2: the two latitudes considered [deg]
    - inc: the inclination of the satellites
    - N: number of satellites in the shell
    OUT:
    - number of satellites with l1<= lat <= l2
    '''
    
    myinc = np.where(inc > 90., 180.-inc, inc) # for retrogr orbits
    return  N/np.pi * (myarcsin(l2/myinc) - myarcsin(l1/myinc))

#---------------------------------------------------------------------------

def satNumDensity(delta1,delta2,satInc,satNum):
    '''Density of satellites in a field

    IN:    
    - delta1,2 = min and max latitude [deg] of the field
    - satNum: total number of sat in the shell
    - satInc: inclination of shell
    
    OUT: the density of satellite at the field [sat/sq.deg]

    Note: accounts for the shrinking sky at higher latitudes.
    '''
    
    satNumDensity = satCount(delta1,delta2,satInc,satNum) \
        / ( 360.*180./np.pi * (np.sin(np.radians(delta2)) - np.sin(np.radians(delta1))) )
        # number of satellites / size of the band
    return satNumDensity
    
#---------------------------------------------------------------------------


def integrateSat(ElLim, AzEl, density ):
    '''count the total number of satellites above an elevation

    IN
    - ElLims: vector of the elevetions above which we want the sat counts [deg]
      The elevations are expected to come sorted by decreasing values (eg [60, 40, 20])
    - AzEl: matrix of Az and El
    - density: n density of satellites over the AzEl matrix
    
    OUT
    - number of satellites above ElLim (vector, same size as ElLim)
    '''
    
    ElCum = np.zeros_like(ElLim)
    Eli = 0
    
    wCum = 0. # integrator

    
    i = len(AzEl[1,:,0]) -1 # we start at zenith
    step = AzEl[1,1,0] - AzEl[1,0,0] # step in elevation
    
    while i >=0 and Eli < len(ElLim): # scan elevation rings
        if AzEl[1,i,0] <= ElLim[Eli]:
            # close one of the requested elevations
            ElCum[Eli] = wCum
            Eli += 1
        
        areaElev = np.degrees(2*np.pi*np.cos(np.radians( AzEl[1,i,0] ))) * step
        averDensity =  np.average(density[i])
                
        wCum += averDensity * areaElev # integrate

        i -= 1 # next elevation ring

    if Eli < len(ElCum):
        ElCum[Eli] = wCum
    return ElCum
#----------------------------------------------------------------------

def Pol2Rec(AzEl,R):
    '''Convert Azimuth,Elevation to rectangular coordinates

    IN:
    - AzEl, an array of [Azimuth, Elevation] (in [deg])
    - R: radius of the points

    OUT:
    - X,Y,Z rectangular coordinates, same unit as R
    '''
    
    Azr = np.radians(AzEl[0])
    Elr = np.radians(AzEl[1])
    cE = np.cos(Elr)
    XYZ = np.array( R*np.array([np.cos(Azr) *cE ,
                    np.sin(Azr) *cE ,
                    np.sin(Elr)
                    ]))
    return XYZ

#---------------------------------------------------------------------------

def Rec2Pol(xyz):
    '''Rectangular to polar coordinate conversion
    IN: on XYZ point as an array
    OUT: 
    - [Az,El], Azimuth and Elevation [deg]
    - R, radius (same unit as XYZ)
    '''
    R = np.linalg.norm(xyz, axis=0)
    Elr  = np.arcsin( xyz[2]/R ) #z
    Az = np.degrees(np.arctan2(xyz[1],xyz[0]))
    return np.array([Az, np.degrees(Elr)]), R

#---------------------------------------------------------------------------

def AltAzEqu(lat,XYZ):
    '''Convert XYZ rectangular coordinates from AltAz to Equatorial 
    (or vice-versa)
    
    IN: 
    - lat: latitude of the observatory [deg]
    - XYZ: array of rectangular coordinates in AltAz (or Eq)

    OUT:
    - XYZ: array of rectangular coordinates in Eq (or AltAz)
    '''

    latr = np.radians(lat)
    sl = np.sin(latr)
    cl = np.cos(latr)
    xyz = np.array([-sl*XYZ[0] + cl*XYZ[2] ,
                    XYZ[1],
                    cl*XYZ[0] + sl*XYZ[2]
                    ])
    return xyz
#---------------------------------------------------------------------------

def AltAz2Delta(lat,alt,AzEl):
    '''Topocentric distance and normal to shell
    in:
    - lat latitude of the site [deg]
    - alt altitude of the shell[km]
    - AzEl: array of topocentric Az, El [deg]
    
    out:
    - alpha: longitude of satellite [deg]
    - delta: latitude of satellite [deg]
    - Delta: topocentric distance [km]
    - costheta: cos of angle between line of sight and normal to shell at satellite.
    '''

    latr = np.radians(lat)
    sl = np.sin(latr)
    cl = np.cos(latr)
    rs = constants.earthRadius+alt
    
    
    # from Az, El to xyz equatorial
    XYZ = Pol2Rec(AzEl,1.)
    xyz = AltAzEqu(lat,XYZ)

    # Delta equation:  Da Delta2 + Db Delta + Dc = 0
    Da = 1.
    Db = 2.*constants.earthRadius * (xyz[0] * np.cos(latr) + xyz[2] * sl)
    Dc = -alt*(alt + 2.* constants.earthRadius)

    # determinant of the equation
    Ddeterm = Db**2 - 4.* Da*Dc

    # solutions
    Delta1 = (np.sqrt(Ddeterm) - Db)/2./Da
    #Delta2 = (-np.sqrt(Ddeterm) - Db)/2./Da
    
    # [7]: extract delta=latitude of satellite
    sindelta = (Delta1* xyz[2] + constants.earthRadius*sl )/ rs
    deltar = (np.arcsin(sindelta))
    delta = np.degrees(deltar)
    cd = np.cos(deltar)
    
    # [5,6]: extract alpha = long
    alphax = (Delta1*xyz[0] + constants.earthRadius*cl)/cd/rs
    alphay = (Delta1*xyz[1]        )/cd/rs
    alpha = np.degrees(np.arctan2(alphay,alphax))
    
    # [10] costheta:
    costheta = (rs**2 + Delta1**2 - constants.earthRadius**2 )/(2.*Delta1*rs)
    
    
    return alpha, delta, Delta1, costheta
#----------------------------------------------------------------------

def fillAzEl(step):
    '''Grid a hemisphere with points in  Elv, Az, with step

    IN: sep, distance beteen points in degrees
    OUT: [ array of Az, array of El]

    First step in elevation is at step/2, so [0, step]
    All Elev rings have the same number of points (so the density 
    at zenith is much higher)
    '''

    El = np.arange(0.+step/2.,90.,step) # so that the 1st one is [0, step]
    Az = np.arange(0,361.,step)
    fillAz, fillEl = np.meshgrid(Az,El)

    return np.array([fillAz, fillEl])
#----------------------------------------------------------------------

def radec2elev(ha,delta,lat):
    '''Elevation from HourAngle, Delta
    
    IN:
    - ha, delta: hour angle (or long), dec (or lat) [deg]
    - lat: latitude of observer [deg]
    OUT
    - elevation [deg]
    '''
    har = np.radians(ha)
    deltar = np.radians(delta)
    latr = np.radians(lat)
    sine = np.sin(latr)*np.sin(deltar) + np.cos(latr)*np.cos(deltar)*np.cos(har)
    el = np.degrees(np.arcsin(sine))

    return el
#----------------------------------------------------------------------

def radec2azel(ha,delta,lat):
    '''Azimut,Elevation from HourAngle,Dec

    IN 
    - ha, delta: hour angle (or long), dec (or lat) [deg]
    - lat: latitude of observer [deg]
    OUT
    - az, elevation [deg], same shape as HA,Delta
    '''
    
    har = np.radians(ha)
    deltar = np.radians(delta)
    latr = np.radians(lat)
    sine = np.sin(latr)*np.sin(deltar) + np.cos(latr)*np.cos(deltar)*np.cos(har)
    elr = np.arcsin(sine)
    cose = np.cos(elr)
    
    azr = np.arctan2(-np.sin(har)*np.cos(deltar)/cose,
                     (np.sin(deltar)-np.sin(latr)*sine)/(np.cos(latr)*cose))

    return np.degrees(azr),np.degrees(elr)
#----------------------------------------------------------------------

def elev2ra(elev,delta,lat):
    # in:
    #   elev:elevation of target
    #   delta: declination of target
    #   lat: latitude of observatory
    #   all in deg
    # out: ra, hourangle. Note that -ra is also a solution
    latr = np.radians(lat)
    deltar = np.radians(delta)
    cosra = (np.sin(np.radians(elev)) - np.sin(latr)*np.sin(deltar))/(np.cos(latr)*np.cos(deltar))
    return np.degrees(np.arccos(cosra))
#----------------------------------------------------------------------


def RaDecAlt2xyz(alpha,delta,alt):
    # input: alpha, delta, altitude of satellites
    # output: xyz equatorial of satellites
    
    rs = constants.earthRadius+alt
    alphar = np.radians(alpha)
    deltar = np.radians(delta)
    xyz = np.array([rs* np.cos(alphar) * np.cos(deltar),
                    rs* np.sin(alphar) * np.cos(deltar),
                    rs* np.sin(deltar) ])

    return xyz
#----------------------------------------------------------------------

def solIllum(xyz,alphas, deltas):
    # input:
    #   xyz: equatorial of satellites,
    #   alphas, deltas [degrees], coordinates of the Sun
    # out: illumination 1/0 for satellites
    
    asr = np.radians(-alphas) ## Sun moves towards West
    dsr = np.radians(deltas)
    cas = np.cos(asr)
    sas = np.sin(asr)
    cds = np.cos(dsr)
    sds = np.sin(dsr)
    re2 = constants.earthRadius*constants.earthRadius
    
    #rotation of alphas along z 
    xyz1 = np.array([xyz[0]* cas + xyz[1]* sas ,
                    -xyz[0]* sas + xyz[1]* cas ,
                    xyz[2] ])

    #rotation of deltas along y2
    xyzs = np.array([ xyz1[0]* cds + xyz1[2]* sds ,
                      xyz1[1] ,
                      -xyz1[0]* sds + xyz1[2]* cds ])

    
    illum = np.zeros_like(xyzs[0]) # init to shadow
    illum[xyzs[0] >= 0] = 1.       # those in front of the Earth are illuminated
    illum[(xyzs[1]**2 + xyzs[2]**2) >= re2 ] = 1.  # those further than constants.earthRadius are illum'd
    
    return illum
#----------------------------------------------------------------------

def satGeoVel(alpha,delta,inc,alt):
    #in:
    #   alpha, delta: longitude and latitude of the satellite, geocentric
    #               equatorial [deg]
    #   inc, alt: orbit inclination [deg] and alt [km]
    #returns:
    #   the two geocentric velocity vectors (xyz geocentric equatorial)
    #   for the two orbits with inc,
    #   alt that cross the alpha delta point.
    
    rs = constants.earthRadius+alt
    
    alphar = np.radians(alpha)
    incr = np.radians(inc)
    si = np.sin(incr)
    ci = np.cos(incr)

    
    # find nodes omega0 and omega1
    longr =  np.arcsin(np.tan(np.radians(delta))/np.tan(np.radians(inc)))
    
    omegar = np.zeros_like([alphar,alphar])
    omegar[0] = alphar - longr
    omegar[1] = alphar + longr + np.pi

    # unit vector normal to orbit
    N = np.array([ np.sin(omegar)*si,
                         -np.cos(omegar)*si,
                         ci + omegar*0.])


    #satellite unit vectors 
    S = Pol2Rec((alpha,delta),1.)

    # satellite velocities [km/s] = Vel * ( N x S ) 
    VS = np.cross(N, S, axis=0) * np.sqrt(constants.gravityMu/rs )

    if 0:
        print("[satGeoVel] -----v")
        print("Nodes:\n", np.degrees(omegar))
        print("Normal:")
        print('[:,0]\n',N[:,0])
        print('[:,1]\n',N[:,1])
        print("|N|\n",np.linalg.norm(N[:,0],axis=0),np.linalg.norm(N[:,1],axis=0))
        print("Sat unit vector:\n", S)
        print("|S|:\n",np.linalg.norm(S,axis=0))
        print("Vsat:\n",VS)
        print("[:,0]:\n",VS[:,0])
        print("[:,1]:\n",VS[:,1])
        print("|V|",np.linalg.norm(VS[:,0],axis=0),np.linalg.norm(VS[:,1],axis=0))
        print("[satGeoVel] -----^")

    return np.nan_to_num(VS)
#----------------------------------------------------------------------


def satTopoVel(VS,lat):
    #In:
    #   VS, geocentric equatorial velocity vectors of the satellites
    #   lat of the observatory
    #OUT
    #   obsvel of topocentric equatorial velocity vector, i.e.
    #     VS corrected for the velocity of the observatory

    
    # observatory velocity
    VO = np.array([0.,constants.earthRotation*constants.earthRadius*np.cos(np.radians(lat)),0.])

    # observed velocity vector
    ObsVel = np.array([VS[0] - VO[0],VS[1] - VO[1],VS[2] - VO[2]])
    
    if 0:
        print("[satTopoVel] --------v")
        print("Vobs:", VO)
        print("|V|: {:.3f} km/s".format(np.linalg.norm(VO)))
        print("apparent Vobs:", ObsVel)
        print("|V|",np.linalg.norm(ObsVel[:,0],axis=0),np.linalg.norm(ObsVel[:,1],axis=0))
        print("[satTopoVel] --------^")
    return ObsVel
#----------------------------------------------------------------------

def AzEl2Vel(alpha, delta, Delta,lat,inc,alt):
    '''Apparent average angular velocity

    IN
    - alpha, delta: geocentric position of the satellite [deg]
    - Delta: distance Observatory-satellite [km]
    - lat latitude of the observatory [deg]
    - inc, alt of the satellites in this shell [deg],[km]

    OUT
    -  AngularVel: apparent (from obs) average (for satellites moving
       up and down) velocity of the satellites. [deg/s]
    '''

    # geocentric coordinates of the sat
    CS = Pol2Rec((alpha,delta), constants.earthRadius+alt) 
    
    # geocentric coordinates of the observatory
    wCO = Pol2Rec((0.,lat),constants.earthRadius)    
    CO = np.array([[wCO[0]], [wCO[1]], [wCO[2]]])
    
    # topocentric coords of sat:
    OS = CS - CO

    #geocentry velocity vector of the sat,
    VS = satGeoVel(alpha, delta, inc, alt)

    #topocentric equ. velocity vector
    ObsVel = satTopoVel(VS,lat)

    #Parallel component of velocity vector:
    Delta2 = Delta*Delta
    NormOV = (OS[0]*ObsVel[0] + OS[1]*ObsVel[1] + OS[2]*ObsVel[2])/Delta2
    ObsVelParallel  = np.array([NormOV[0]*ObsVel[:,0], NormOV[0]*ObsVel[:,1]])
    
    #Perpandicular component of vel vector
    wObsVelPerpan = np.array([
        ObsVel[:,0] - ObsVelParallel[0],
        ObsVel[:,1] - ObsVelParallel[1]
    ])

    #Norm of perp.component of Vel vector
    ObsVelPerpan = np.array([ np.linalg.norm(wObsVelPerpan[0],axis=0), 
                             np.linalg.norm(wObsVelPerpan[1],axis=0) ])

    #apparent angular velocity of satellite  [deg/sec]:
    AngularVel = np.average(np.degrees(np.arctan(ObsVelPerpan/Delta)), axis=0)

    if 0 :
        print("[AzEl2Vel] alpha:",alpha)
        print("[AzEl2Vel] delta:",delta)
        print("[AzEl2Vel] Delta:km",Delta)
        print("CS xyz:", CS)
        print("CO xyz:", CO)
        print("OS xyz:", OS)

        print("Distance: ", np.linalg.norm(OS), Delta)
        print("ObsVel              :", ObsVel)
        print("ObsVel Parallel     :", ObsVelParallel, np.linalg.norm(ObsVelParallel[0]), np.linalg.norm(ObsVelParallel[1]))
        print("ObsVel Perpandicular:", ObsVelPerpan)
        print("AngularVel          :", AngularVel*60, "deg/min")
        
    return AngularVel
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------






def modelOneConstMag(AzEl,lat, alphas,deltas,
                     inc,alt,num ):
    '''Model one single shell over a set of Az,El pointings
    IN
    - AzEl: array [ array of Azimuths, array of Elevations] on which the
      constellation shall be evaluated. Both in [deg]
    - lat: latitude of the observer [deg]
    - alphas, deltas: HourAngle and Dec. of the Sun [deg]
    - inc, alt, num: parameters of the satellite constellation shell:
      - inc: inclination [deg]
      - alt: altitude [km]
      - num: number of satellites in the shell
    OUT
    - illuminated satellite number density (same shape as AzEl)
    - illuminated satellite apparent angular velocity (same shape as AzEl)
    - illuminated satellite magnitudes
    '''


    # sun az, el
    azs,els = radec2azel(alphas,deltas, lat)

    if len(AzEl.shape) == 3:
        AzElreshape = np.reshape(AzEl,(2,AzEl.shape[1]*AzEl.shape[2]))
        step = AzEl[1,1,0] - AzEl[1,0,0]
    else:
        AzElreshape = AzEl
        step = 1.
        
    # geocentric equ. alpha,delta of sat, and   observatory dist, angle 
    alpha, delta, Delta, costheta = AltAz2Delta(lat,alt,AzElreshape)
    
    # geocentric equ. rect. of satellite
    xyz = RaDecAlt2xyz(alpha,delta, alt)

    # Velocities
    wAngularVel = AzEl2Vel(alpha, delta, Delta,lat,inc,alt)  
    
    #Density    
    # get delta of top of field of view
    wAzEl = np.copy(AzElreshape)
    wAzEl[1] += step
    _, deltaTop, _, _ = AltAz2Delta(lat,alt,wAzEl)

    # get delta of bottom of field
    wAzEl = np.copy(AzElreshape)
    wAzEl[1] -= step
    _, deltaBot, _, _ = AltAz2Delta(lat,alt,wAzEl)

    # density at this place
    densitys = satNumDensity(deltaBot, deltaTop,inc,num) \
                   * (Delta/(constants.earthRadius+alt))**2 / costheta 

    # Illuminated satellites
    illum = solIllum(xyz,alphas, deltas)
    wdensityi = densitys * illum

    #  MAGNITUDE of the satellites:

    wmag =  constants.mag550 + 5.*np.log10(Delta/550.)        # distances
    wmag += constants.extinction*(Delta/alt -1.)    # extinction

    ## ZTF brightnening
    # deltaAzs = np.cos(np.radians( AzEl[0] - azs ))
    # mag = 1-np.degrees(np.arccos(DeltaAzs))/constants.ZTFmagAzCut 
    # mag[Dmag < 0] = 0.
    ##    experimental elevation function: peaks when angle(Sat,sun)=45deg
    # el = AzEl[1] - els
    # el[Del>90] = 0.
    # Del = 1-((Del-constants.ZTFmagAngPeak)/constants.ZTFmagAngPeak)**2
    # wmag += Dmag*Del*constants.ZTFmagAzBright     # ZTF brightening
    
    return \
        np.reshape(wdensityi,   (AzEl.shape[1],AzEl.shape[2]) ) ,\
        np.reshape(wAngularVel, (AzEl.shape[1],AzEl.shape[2]) ),\
        np.reshape(wmag,        (AzEl.shape[1],AzEl.shape[2]) )
    


#------------------------------------------------------------------------------
