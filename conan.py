#!/usr/bin/env python3
# ConAn
# Constellation Analytic simulations
#
#
# conan.py: definitions and functions
#==============================================================================
import numpy as np
from astropy.table import Table, join

# ConAn: constants
import constants as Cconst

#---------------------------------------------------------------------------
def get_sun(jd):
    '''
    get coordinates of the Sun

    In: jd, full julian day

    out:
       RA and Dec (both degrees) of the Sun on JD
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


#---------------------------------------------------------------------------
def loadConstellations(constellations):
    '''
    load a collection of constellations into corresponding arrays

    IN: an array of constellations from ConAn/constellations.py

    OUT: set of arrays with the constellation parameters
    - consLab, constellation label
    - consNum, number of satellites in the constellation
    - consPla, number of orbital planse
    - consNPl, number of satellites per plane
    - consInc, orbital inclination
    - consAlt, altitude of the constellation [km]

    Note: consNum = consPla*consNPl, but this is not enforced
    '''

    consLab = np.array([])
    consAlt = np.array([])
    consNum = np.array([])
    consPla = np.array([])
    consNPl = np.array([])
    consInc = np.array([])
    for myCons in constellations:
        consLab = np.append(consLab, myCons[0])
        consNum = np.append(consNum, int(myCons[1]))
        consPla = np.append(consPla, int(myCons[2]))
        consNPl = np.append(consNPl, int(myCons[3]))
        consInc = np.append(consInc, float(myCons[4]))
        consAlt = np.append(consAlt, float(myCons[5]))

    return  consLab, consNum, consPla, consNPl, consInc, consAlt

#---------------------------------------------------------------------------
def velPosAng(delta,satInc):
    '''
    Position Angles (both up and down) of the motion of a satellite at given lat.

    IN
    - delta = lat, latitude of the satellite [deg]
    - satInc: inclination of shell [deg]

    OUT:
    - two position angles (up and down) [deg]
    '''

    # properly deals with retrograde orbits (with satInc > 90)

    sintheta = np.cos(np.radians(satInc))/np.cos(np.radians(delta))
    theta1 = np.degrees(np.arcsin(sintheta))
    theta2 = np.degrees(np.arcsin(-sintheta))-180.*np.sign(satInc-90)
    return theta1, theta2
#---------------------------------------------------------------------------

def myarcsin(x):
    '''
    arcsin function extended beyond [-1,1]
    '''

    myx1 = np.where(x > 1., 1., x)
    myx2 = np.where(myx1 < -1., -1., myx1)
    myarcsin = np.arcsin(myx2)

    return myarcsin
#---------------------------------------------------------------------------

def satCount(l1,l2,inc,N):
    '''
    Number of satellites between two latitudes

    IN:
    - l1, l2: the two latitudes [deg]
    - inc: inclination of the satellite orbit [deg], can be retrograde
    - N: number of satellites in constellation shell

    OUT:
    - number of satellites in [l1, l2]
    '''
    # returns number of satellites between l1 and l2,
    # constellation of inc., N satellites

    myinc = np.where(inc > 90., 180.-inc, inc) # for retrogr orbits
    satCount = N/np.pi * (myarcsin(l2/myinc) - myarcsin(l1/myinc))
    return satCount

#---------------------------------------------------------------------------

def satNumDensity(delta1,delta2,satInc,satNum):
    '''
    Density of satellites in latitude range

    IN:
    - delta1,2 =  latitude of the field [l1,l2] (unit: deg)
    - satInc: orbit inclination of shell [deg]
    - satNum: total number of sat in shell

    OUT:
    - the density of satellite in the field [sat/sq.deg]

    Accounting for the shrinking sky at higher latitudes.

    From surface of spherical cap = 2pi h = 2pi (1-sin(delta)),
    so, surface of the band: 2pi (sin(d2)-sin(d2))

    '''

    satNumDensity = satCount(delta1,delta2,satInc,satNum) \
        / ( 360.*180./np.pi * (np.sin(np.radians(delta2)) - np.sin(np.radians(delta1))) )
                   # size of the band

    return satNumDensity
#---------------------------------------------------------------------------


def integrateSat(ElLim, AzEl, density ):
    '''Number of satellites above elevation

    IN
    - ElLims: vector of the elevetions above which we want the sat counts
    - AzEl: grid of Az and El [deg]
    - density: n density of satellites

    OUT
    -  ElCum: number of satellites above ElLim
    '''

    ElCum = np.zeros_like(ElLim)
    Eli = 0
    wCum = 0. # integrator


    i = len(AzEl[1,:,0]) -1
    step = AzEl[1,1,0] - AzEl[1,0,0]

    while i >=0:
        if AzEl[1,i,0] <= ElLim[Eli]:
            ElCum[Eli] = wCum
            Eli += 1

        areaEl = np.degrees(2*np.pi*np.cos(np.radians( AzEl[1,i,0] ))) * step
        averd =  np.average(density[i])

        wCum += averd * areaEl

        i -= 1

    if Eli < len(ElCum):
        ElCum[Eli] = wCum
    return ElCum
#----------------------------------------------------------------------

def Pol2Rec(AzEl,R):
    '''Polar to rectangular conversion

    IN:
    - AzEl, grid of Az and El [deg]
    - R: radius [unit]

    OUT:
    - XYZ, grid like AzEl of rectangular coord [same unit as R]
    '''

    Azr = np.radians(AzEl[0])
    Elr = np.radians(AzEl[1])
    cE = np.cos(Elr)
    XYZ = np.array(R*np.array([np.cos(Azr) *cE ,
                    np.sin(Azr) *cE ,
                    np.sin(Elr)
                    ]))
    return XYZ
#---------------------------------------------------------------------------


def Rec2Pol(xyz):
    '''Rectangular to spherical conversion

    in:
    - XYZ, array of XYZ

    OUT:
    - AzEl, array of Az, El, same substructure shape as XYZ, [deg]
    - R, array of radii, one per (XYZ)

    Note:
    - El measured from equator
    '''

    R = np.linalg.norm(xyz, axis=0)
    Elr = np.arcsin(xyz[2]/R)
    Az = np.degrees(np.arctan2(xyz[1],xyz[0]))
    return np.array([Az, np.degrees(Elr)]), R
#---------------------------------------------------------------------------

def AltAzEqu(lat,XYZ):
    # AltAz -> Eq
    # OR Equ -> AltAz

    # latitude in deg
    # XYZ in altAz system (or Eq)
    # returns xyz in equatorial system (or AltAz)
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
    '''
    Convert AzEl to Long,Lat and distance

    IN
    - lat: latitude of the observatory [deg]
    - alt: altitude of the constellation [km]
    - AzEl: array topocentric Az, El [deg]

    OUT
    - alpha, delta; longitude, latitude of a satellite
      at Az,El         [deg]
    - Delta: topocentric distance Obs-Satellite [km]
    - costheta: cos of angle between line of sight and normal to shell at satellite.
    '''

    latr = np.radians(lat)
    sl = np.sin(latr)
    cl = np.cos(latr)
    rs = Cconst.earthRad + alt


    # from Az, El to xyz equatorial
    XYZ = Pol2Rec(AzEl,1.)
    xyz = AltAzEqu(lat,XYZ)

    # Delta equation:  Da Delta2 + Db Delta + Dc = 0
    Da = 1.
    Db = 2.*Cconst.earthRad * (xyz[0] * np.cos(latr) + xyz[2] * sl)
    Dc = -alt*(alt + 2.* Cconst.earthRad)

    # determinant of the equation
    Ddeterm = Db**2 - 4.* Da*Dc

    # solutions
    Delta1 = (np.sqrt(Ddeterm) - Db)/2./Da
    #Delta2 = (-np.sqrt(Ddeterm) - Db)/2./Da

    # [7]: extract delta=latitude of satellite
    sindelta = (Delta1* xyz[2] + Cconst.earthRad*sl )/ rs
    deltar = (np.arcsin(sindelta))
    delta = np.degrees(deltar)
    cd = np.cos(deltar)

    # [5,6]: extract alpha = long
    alphax = (Delta1*xyz[0] + Cconst.earthRad*cl)/cd/rs
    alphay = (Delta1*xyz[1]        )/cd/rs
    alpha = np.degrees(np.arctan2(alphay,alphax))

    # [10] costheta:
    costheta = (rs**2 + Delta1**2 - Cconst.earthRad**2 )/(2.*Delta1*rs)


    return alpha, delta, Delta1, costheta
#----------------------------------------------------------------------

def fillAzEl(step):
    '''Create a AzEl grid, [0,360] * [0, 90]

    IN: step [deg] of the array in Az and in El

    OUT: AzEl = [Azimuth, Elevation],
    '''
    # fills in a hemisphere of Elv, Az,
    # in: sep, distance beteen points in degrees
    # out:  fillAz, fillEl

    El = np.arange(0.+step/2.,90.,step) # so that the 1st one is [0, step]
    Az = np.arange(0,360.001,step)
    fillAz, fillEl = np.meshgrid(Az,El)

    return np.array([fillAz, fillEl])
#----------------------------------------------------------------------

def radec2elev(ha,delta,lat):
    #in
    #  ha, delta: hour angle/long, dec/lat [deg]
    #  lat: latitude of observer [deg]
    #out elevation [deg]
    har = np.radians(ha)
    deltar = np.radians(delta)
    latr = np.radians(lat)
    sine = np.sin(latr)*np.sin(deltar) + np.cos(latr)*np.cos(deltar)*np.cos(har)
    el = np.degrees(np.arcsin(sine))
    return el
#----------------------------------------------------------------------

def radec2azel(ha,delta,lat):
    #in
    #  ha, delta: hour angle/long, dec/lat [deg]
    #  lat: latitude of observer [deg]
    #out az, elevation [deg]
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

    rs = Cconst.earthRad+alt
    alphar = np.radians(alpha)
    deltar = np.radians(delta)
    xyz = np.array([rs* np.cos(alphar) * np.cos(deltar),
                    rs* np.sin(alphar) * np.cos(deltar),
                    rs* np.sin(deltar) ])

    return xyz
#----------------------------------------------------------------------

def solIllum(xyz,sunAlpha, sunDelta):
    # input:
    #   xyz: equatorial of satellites,
    #   sunAlpha, sunDelta [degrees], coordinates of the Sun
    # out: illumination 1/0 for satellites

    asr = np.radians(-sunAlpha) ## Sun moves towards West
    dsr = np.radians(sunDelta)
    cas = np.cos(asr)
    sas = np.sin(asr)
    cds = np.cos(dsr)
    sds = np.sin(dsr)
    re2 = Cconst.earthRad * Cconst.earthRad

    #rotation of sunAlpha along z
    xyz1 = np.array([xyz[0]* cas + xyz[1]* sas ,
                    -xyz[0]* sas + xyz[1]* cas ,
                    xyz[2] ])

    #rotation of sunDelta along y2
    xyzs = np.array([ xyz1[0]* cds + xyz1[2]* sds ,
                      xyz1[1] ,
                      -xyz1[0]* sds + xyz1[2]* cds ])


    illum = np.zeros_like(xyzs[0]) # init to shadow
    illum[xyzs[0] >= 0] = 1.       # those in front of the Earth are illuminated
    illum[(xyzs[1]**2 + xyzs[2]**2) >= re2 ] = 1.  # those further than Cconst.earthRad are illum'd

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

    rs = Cconst.earthRad+alt

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
    VS = np.cross(N, S, axis=0) * np.sqrt(Cconst.mu/rs )

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
    VO = np.array([0.,Cconst.earthOmega*Cconst.earthRad*np.cos(np.radians(lat)),0.])

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

def VolumeOfPixel(  dalpha, delta, ddelta, radius, dradius):
    '''The volume of a partial shell
    IN
    - dalpha: extent in longitude [deg]
    - delta, ddelta: origin and extent in latitude [deg]
    - radius, dradius: inner radius and thickness of shell [unit]
    
    OUT
    - volume of the section of the shell [unit^2]
    
    Notes:
    -  surface of full cap from colat=90-delta to pole
          2pi r2 (1-cos(colat)) = 2pi r2 (1-sin(delta))
    -  surface of full ring from delta1 to delta2:
          2pi r2 (sin(delta1) - sin(delta2))
    - partial ring from  alpha1 to alpha2
          (alpha2-alpha1)/360 *   2pi r2 (sin(delt1)-sin(delta2))

    - Volume of the sphere r
          4/3 pi r3
    - shell between r and r+dr
         4/3 pi ( (r+dr)3 - r3 )
         with   (r+dr)3 = r3 + 3r2 dr + 3 r dr2 + dr3
         4/3 pi  (  3r2 dr   +  3r dr2  + dr3 )
         ~   4pi r2 dr = surface * dr
    '''
    
    surface  =  2.*np.pi* radius**2 
    surface *=  (np.sin(np.radians(delta +ddelta)) - np.sin(np.radians(delta))) 
    surface *=  dalpha/360.

    return surface * dradius



#----------------------------------------------------------------------
def AzEl2Vel(alpha, delta, Delta,obsLatitude,inc,alt):
    
    #IN
    #  alpha, delta: geocentric position of the satellite
    #  Delta: distance Observatory-satellite
    #  obsLatitude latitude of the observatory
    #  inc, alt of the satellites
    #OUT
    #  AngularVel: apparent (from obs) average (for satellites moving
    #  up and down) velocity of the satellites. [deg/s]

    # geocentric coordinates of the sat
    CS = Pol2Rec((alpha,delta), Cconst.earthRad+alt)

    # geocentric coordinates of the observatory
    wCO = Pol2Rec((0.,obsLatitude),Cconst.earthRad)
    CO = np.array([[wCO[0]], [wCO[1]], [wCO[2]]])

    # topocentric coords of sat:
    OS = CS - CO

    #geocentry velocity vector of the sat,
    VS = satGeoVel(alpha, delta, inc, alt)

    #topocentric equ. velocity vector
    ObsVel = satTopoVel(VS,obsLatitude)

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

    #apparent angular velocity of satellite (natively: [deg/sec]):
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
def modelOneConstMag(AzEl,obsLatitude, sunAlpha,sunDelta,
                     satInc,satAlt,satNum, mag550=Cconst.mag550 ):
    '''
    IN
    -  AzEl: mesh of Azimuth, Elevation [deg]
    -  obsLatitude: latitude of the observatory 
    -  sunAl, sunDel: coordinates of the Sun [deg]
    -  satInc:  satellite constellation shell Inclination [deg]
    -  satAlt:  satellite constellation shell Altitude [km]
    -  satNum:  number of satellite in constellation shell
    -  mag550:  magnitude of a satellite in the Shell

    OUT
    -  density of illuminated satellites (n/deg2
    -  angular velocity 
    -  magnitude

    '''
    # sun az, el (not used)
    # sunAzimuth,sunElevation = radec2azel(sunAlpha,sunDelta, obsLatitude)

    if len(AzEl.shape) == 3:
        AzElreshape = np.reshape(AzEl,(2,AzEl.shape[1]*AzEl.shape[2]))
        step = AzEl[1,1,0] - AzEl[1,0,0]
    else:
        AzElreshape = AzEl
        step = 1.

    # geocentric equ. alpha,delta of sat, and   observatory dist, angle
    alpha, delta, Delta, costheta = AltAz2Delta(obsLatitude,satAlt,AzElreshape)

    # geocentric equ. rect. of satellite
    xyz = RaDecAlt2xyz(alpha,delta, satAlt)

    # Velocities
    wAngularVel = AzEl2Vel(alpha, delta, Delta,obsLatitude,satInc,satAlt)

    #Density
    # get delta of top of field of view
    wAzEl = np.copy(AzElreshape)
    wAzEl[1] += step
    _, deltaTop, _, _ = AltAz2Delta(obsLatitude,satAlt,wAzEl)

    # get delta of bottom of field
    wAzEl = np.copy(AzElreshape)
    wAzEl[1] -= step
    _, deltaBot, _, _ = AltAz2Delta(obsLatitude,satAlt,wAzEl)

    # density at this place
    # adjust angular size for distance,
    # and for apparent orientation of the shell on line-of-sight
    densitys = satNumDensity(deltaBot, deltaTop,satInc,satNum) \
                   * ( Delta/(Cconst.earthRad+satAlt) )**2 / costheta

    # Illuminated satellites
    illum = solIllum(xyz, sunAlpha, sunDelta)

    wdensityi = densitys * illum




    #  MAGNITUDE of the satellites:

    #   ZTF brightnening
#    DeltaAzs = np.cos(np.radians( AzEl[0] - sunAzimuth ))
#    Dmag = 1-np.degrees(np.arccos(DeltaAzs))/Cconst.ZTFmagAzCut
#    Dmag[Dmag < 0] = 0.
#    #   experimental elevation function: peaks when angle(Sat,sun)=45deg
#    Del = AzEl[1] - sunElevation
#    Del[Del>90] = 0.
#    Del = 1-((Del-Cconst.ZTFmagAngPeak)/magAngPeak)**2


    mag = np.reshape(
        Cconst.mag550
        + 5.*np.log10(Delta/550.) # distances
        + Cconst.extn*(Delta/satAlt -1.)    # extinction
        , (AzEl.shape[1],AzEl.shape[2]))
    # + Dmag*Del*Cconst.ZTFmagAzBright     # ZTF brightening

    return \
        np.reshape(wdensityi,  (AzEl.shape[1],AzEl.shape[2]) ) ,\
        np.reshape(wAngularVel, (AzEl.shape[1],AzEl.shape[2]) ),\
        mag


#-------------------------------------------------------------------------
def modelOneDebrisShell(AzEl,obsLatitude, sunAlpha,sunDelta,
                     debrisDist, iAlt ):
    '''
    Model one shell of debris.

    IN
    -  AzEl: mesh of Azimuth, Elevation [deg]
    -  obsLatitude: latitude of the observatory 
    -  sunAl, sunDel: coordinates of the Sun [deg]

    -  satAlt:  satellite constellation shell Altitude [km]
    '''

    # iterate on the particle sizes
    debris_i   = np.arange(0,7)
    debris_rad = 10.**(-1.*debris_i)
    debris_m550 = absMag_to_mag550( radius_to_absMag( debris_rad ) )

    debris = debrisDist[iAlt]
    satAlt = debris["Altitude [km]"]

    for id in debris_i:
        pass

        



    if len(AzEl.shape) == 3:
        AzElreshape = np.reshape(AzEl,(2,AzEl.shape[1]*AzEl.shape[2]))
        step = AzEl[1,1,0] - AzEl[1,0,0]
    else:
        AzElreshape = AzEl
        step = 1.

    # geocentric equ. alpha,delta of sat, and   observatory dist, angle
    alpha, delta, Delta, costheta = AltAz2Delta(obsLatitude,satAlt,AzElreshape)

    # geocentric equ. rect. of satellite
    xyz = RaDecAlt2xyz(alpha,delta, satAlt)

    
    # get delta of top of pixel
    wAzEl = np.copy(AzElreshape)
    wAzEl[1] += step
    _, deltaTop, _, _ = AltAz2Delta(obsLatitude,satAlt,wAzEl)

    # get delta of bottom of pixel
    wAzEl = np.copy(AzElreshape)
    wAzEl[1] -= step
    _, deltaBot, _, _ = AltAz2Delta(obsLatitude,satAlt,wAzEl)


    # get alpha at left of pix
    wAzEl = np.copy(AzElreshape)
    wAzEl[0] += step
    alphaLeft, _, _, _ = AltAz2Delta(obsLatitude,satAlt,wAzEl)

    # get alpha at left of pix
    wAzEl = np.copy(AzElreshape)
    wAzEl[0] -= step
    alphaRight, _, _, _ = AltAz2Delta(obsLatitude,satAlt,wAzEl)




    # Number of grains in the pixel in the shell
    # angular surface of the grain:
    # full cap from colat=90-delta to pole
    #     2pi r2 (1-cos(colat)) = 2pi r2 (1-sin(delta))
    # full ring from delta1 to delta2:
    #     2pi r2 (sin(delta1) - sin(delta2))
    # partial ring from  alpha1 to alpha2
    #     (alpha2-alpha1)/360 *   2pi r2 (sin(delt1)-sin(delta2))



    # volume of the sphere r
    #     4/3 pi r3
    # shell between r and r+dr
    #     4/3 pi ( (r+dr)3 - r3 )
    #               (r+dr)3 = r3 + 3r2 dr + 3 r dr2 + dr3
    #     4/3 pi  (  3r2 dr   +  3r dr2  + dr3 )
    #     ~   4pi r2 dr = surface * dr

    # volume of the partial ring:
    #      (alpha2-alpha1)/360 *   pi r2 (sin(delt1)-sin(delta2))  * dr









    # density at this place
    # adjust angular size for distance,
    # and for apparent orientation of the shell on line-of-sight
    ##densitys = satNumDensity(deltaBot, deltaTop,satInc,satNum) \
    ##               * ( Delta/(Cconst.earthRad+satAlt) )**2 / costheta










#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def radius_to_absMag( radius_m, albedo=0.2 ):
    '''convert radius and albedo in Absolute Mag m110

    IN:
    - radius_m: radius of the (spherical) particle [m]
    - albedo

    OUT:
    - m110, absolute mag for r=Delta=1au

    (validated)
    '''

    return Cconst.magSun -2.5* np.log10( ( radius_m/Cconst.au_m )**2 * albedo )

def absMag_to_mag550(absMag ):
    '''convert absolute mag m110 in mag at 550km

    IN: absMag, magnitude at r=Delta=1au, alpha=0

    Out: m550, magnitude at 550km at zenith
    '''

    return absMag + 5.*np.log10(  550.e3 / Cconst.au_m)

def absMag_to_mag(absMag, rAU, dAU):
    '''convert an absMag to mag'''

    return absMag + 5.*np.log10( rAU * dAU)

def read_debrisDist():

    logRad = np.arange(0, 7, 1)

    debrisDist = Table()

    for lr in logRad:
        infile = f'spatial_density_1e-{lr}.csv'    
        w = Table.read( infile )
        w["Spatial Density [1/km^3]"].name = f"n_1e-{lr}"
        if lr == 0:
            debrisDist = w
            print('0NEW')
        else:
            print(lr)
            debrisDist = join( debrisDist, w,'Altitude [km]')
        
    return debrisDist




if __name__ == "__main__":

    print(    VolumeOfPixel( 0, 360, -90, 90, 6400.,1) )
