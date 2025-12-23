#!/usr/bin/env python3
# SatConAnalytic - Satellite Constellation Analytic simulations
# conan.py:  generic constellation functions
#
#
# conan.py: definitions and functions
#==============================================================================
import numpy as np
from astropy.table import Table, join
from astropy.io import ascii

# SatConAn: 
import SatConAnalytic.constants as constants
import SatConAnalytic.constellations as constellations

#---------------------------------------------------------------------------
def consolidate_sun(sunAlpha, sunDelta, sunElev, lat):
    ''' using two available parameters of the sun, get the 3rd one.
    
    IN:
    - sunAlpha: HourAngle of the Sun [deg], or None
    - sunDelta: Dec of the sun [deg]
    - sunElev: MINUS elevation of the sun [deg], or None.
    
    OUT:
    - sunAlpha, sunDelta, sunElev'''

    sunDelta = float(sunDelta)

    if sunAlpha is None:
        sunElev = -float(sunElev)
        sunAlpha = elev2ra(sunElev,sunDelta,lat) # get sun hourangle for twilight
    else:
        sunAlpha = float(sunAlpha)
        sunElev = radec2elev(sunAlpha,sunDelta,lat)

    return sunAlpha, sunDelta, sunElev


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
def siderealTimeDeg(timeJD, longitudeDeg=0.0):
    '''Compute the local sidereal time at a given longitude
    IN:
    - timeJD: Julian Day
    - longitudeDeg: longitude of the site [deg], East positive
    OUT:
    - LST: local sidereal time [deg]
    '''

    # Loc ST
    lst = (280.46061837 + 360.98564736629*(timeJD -2451545.0) )%360
    lst += longitudeDeg
    lst = lst %360.
    return lst
#------------------------------------------------------------------------------
def findHAfromZDr(zdr, dr, latr):
    '''Returns HA (symmetric HAmin=-HAmax) based on zenithal Distance.

    In:
    -  dr: declination or a np.array  declinations [RADIANS]
    -  zdr [rad] can be a scalar or a vector like delta [radians]
    -  latr: latitude of the site [rad]

    Out: corresponding HA [hours]; this case is symmetric, so HAmin=-HAmax
    '''

    # generic spherical trigo:
    cosHA = (np.cos(zdr) - np.sin(latr)*np.sin(dr))/(np.cos(latr)*np.cos(dr))

    # case of unreachable dec
    delta0r = latr - zdr # declination that reaches target airmass
                       # at meridian
    cosHA[dr <= delta0r]  = 1.

    delta1r = latr +zdr # declination that reaches target airmass
                      # at meridian on the other side
    cosHA[dr >= delta1r]  = 1.

    # case of dec that are always above airmass
    deltapolr = -np.pi -latr + zdr # declination that is circumpolar for z
    cosHA[dr < deltapolr] = -1.

    HA = np.degrees( np.arccos(cosHA)) /15.
    return HA


#---------------------------------------------------------------------------
def findConstellations(constellationsll):
    '''Assemble a Constellations object (set of constellations)
    for a list of constellations.

    in: list of constellations ['SL1', 'SL2', 'OWr2']  
        or one of the preset codes defined below

    out: a Constellations object
    '''

    metaConstellations = {
        'SL': ['SL1', 'SL2'],
        'OW': ['OW2r'],
        'SLOW': ['SL1', 'SL2','OW2r'],
        'TODAY': ['YESTERDAY', 'TODAYconst'],
        'SLOWGWAK': ['YESTERDAY',
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
        print(constellations.readConstellations())
        print("=============================================================")
        print("Available preset constellation groups:")
        for c in metaConstellations:
            print(f'  {c}: {metaConstellations[c]}')
            
        exit(0)
    elif constellationsll in metaConstellations:
        constellationsll = metaConstellations[constellationsll]
    else:
        constellationsll = [ constellationsll ]

    return constellations.metaConstellation(constellationsll)
    
#---------------------------------------------------------------------------
def velPosAng(delta,satInc):
    '''Compute the velocity position angles for a list of satellites.

    in:
    - delta = latitude of the satellite(s) [deg]
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
    '''
    arcsin function extended beyond [-1,1]

    IN: x, the value from which arcsin must be computed, float, ]-4e4, 4e4[
    OUT: arcsin(x) [radians] 
    '''
    
    myx1 = np.where(x > 1., 1., x)
    myx2 = np.where(myx1 < -1., -1., myx1)
    myarcsin = np.arcsin(myx2)
    
    return myarcsin

#---------------------------------------------------------------------------
def satCount(l1,l2,satInc,N):
    '''Number of satellites between two latitudes.

    IN:
    - l1, l2: the two latitudes considered [deg]
    - satInc: the inclination of the satellites
    - N: number of satellites in the shell
    OUT:
    - number of satellites with l1<= lat <= l2
    '''
    
    myinc = np.where(satInc > 90., 180.-satInc, satInc) # for retrogr orbits
    return  N/np.pi * (myarcsin(l2/myinc) - myarcsin(l1/myinc))

#---------------------------------------------------------------------------
def satNumDensity(delta1,delta2,satInc,satNum):
    '''
    Density of satellites in latitude range

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
    - ElLims: LIST of the elevations above which we want the sat counts [deg]
      The elevations are expected to come sorted by decreasing values (eg [60, 40, 20])
    - AzEl: matrix of Az and El
    - density: n density of satellites over the AzEl matrix
    
    OUT
    - number of satellites above ElLim (vector, same size as ElLim)
    '''
    
    ElCum = np.zeros_like(ElLim)
    Eli = 0    
    wCum = 0. # integrator

    
    i = len(AzEl[1,:,0]) -1          # we start at zenith
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
    Elr  = np.arcsin( xyz[2]/R ) #z
    Az = np.degrees(np.arctan2(xyz[1],xyz[0]))
    return np.array([Az, np.degrees(Elr)]), R

#---------------------------------------------------------------------------
def AltAzEqu(obsLatitude,XYZ):
    '''Convert XYZ rectangular coordinates from AltAz to Equatorial 
    (or vice-versa)
    
    IN: 
    - obsLatitude: latitude of the observatory [deg]
    - XYZ: array of rectangular coordinates in AltAz (or Eq)

    OUT:
    - XYZ: array of rectangular coordinates in Eq (or AltAz)
    '''

    latr = np.radians(obsLatitude)
    sl = np.sin(latr)
    cl = np.cos(latr)
    xyz = np.array([-sl*XYZ[0] + cl*XYZ[2] ,
                    XYZ[1],
                    cl*XYZ[0] + sl*XYZ[2]
                    ])
    return xyz

#---------------------------------------------------------------------------
def AltAz2Delta(obsLatitude,satAlt,AzEl):
    '''Topocentric distance and normal to shell
    in:
    - obsLatitude latitude of the site [deg]
    - satAlt altitude of the shell[km]
    - AzEl: array of topocentric Az, El [deg]
    
    out:
    - alpha: longitude of satellite [deg]
    - delta: latitude of satellite [deg]
    - Delta: topocentric distance [km]
    - costheta: cos of angle between line of sight and normal to shell at satellite.
    '''

    latr = np.radians(obsLatitude)
    sl = np.sin(latr)
    cl = np.cos(latr)
    rs = constants.earthRadius+satAlt
    
    
    # from Az, El to xyz equatorial
    XYZ = Pol2Rec(AzEl,1.)
    xyz = AltAzEqu(obsLatitude,XYZ)

    # Delta equation:  Da Delta2 + Db Delta + Dc = 0
    Da = 1.
    Db = 2.*constants.earthRadius * (xyz[0] * np.cos(latr) + xyz[2] * sl)
    Dc = -satAlt*(satAlt + 2.* constants.earthRadius)

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
def fill_AzEl(step):
    '''Create a AzEl grid, [0,360] * [0, 90]

    IN: step [deg] of the array in Az and in El

    OUT: AzEl = [Azimuth, Elevation],
    '''

    El = np.arange(0.+step/2.,90.,step) # so that the 1st one is [0, step]
    Az = np.arange(0,360.001,step)
    fillAz, fillEl = np.meshgrid(Az,El)

    return np.array([fillAz, fillEl])

#----------------------------------------------------------------------
def surface_AzEl(Az, El, step):
    '''Surface [sq.deg] an element centred on AzEl, widdh=step

    IN: Az, El, step [deg] 

    OUT: surfaceAz [sqDeg]
    '''

    radius = 180./np.pi
    surface = 2.*np.pi* radius**2 # 1/2 sphere
    surface *=  (np.sin(np.radians(El + step/2.)) - np.sin(np.radians(El - step/2.))) 
         # difference of callotes; =2 for full range
    surface *=  step/360. # longitude fraction

    return surface
#----------------------------------------------------------------------
def radec2elev(ha,delta,obsLatitude):
    '''Elevation from HourAngle, Delta
    
    IN:
    - ha, delta: hour angle (or long), dec (or lat) [deg]
    - lat: latitude of observer [deg]
    OUT
    - elevation [deg]
    '''
    har = np.radians(ha)
    deltar = np.radians(delta)
    latr = np.radians(obsLatitude)
    sine = np.sin(latr)*np.sin(deltar) + np.cos(latr)*np.cos(deltar)*np.cos(har)
    el = np.degrees(np.arcsin(sine))

    return el

#----------------------------------------------------------------------
def radec2azel(ha,delta,obsLatitude):
    '''Azimut,Elevation from HourAngle,Dec

    IN 
    - ha, delta: hour angle (or long), dec (or lat) [deg]
    - obsLatitude: latitude of observer [deg]
    OUT
    - az, elevation [deg], same shape as HA,Delta
    '''
    
    har = np.radians(ha)
    deltar = np.radians(delta)
    latr = np.radians(obsLatitude)
    sine = np.sin(latr)*np.sin(deltar) + np.cos(latr)*np.cos(deltar)*np.cos(har)
    elr = np.arcsin(sine)
    cose = np.cos(elr)
    
    azr = np.arctan2(-np.sin(har)*np.cos(deltar)/cose,
                     (np.sin(deltar)-np.sin(latr)*sine)/(np.cos(latr)*cose))

    return np.degrees(azr),np.degrees(elr)

#----------------------------------------------------------------------
def elev2ra(elev,delta,obsLatitude):
    # in:
    #   elev:elevation of target
    #   delta: declination of target
    #   obsLatitude: latitude of observatory
    #   all in deg
    # out: ra, hourangle. Note that -ra is also a solution
    latr = np.radians(obsLatitude)
    deltar = np.radians(delta)
    cosra = (np.sin(np.radians(elev)) - np.sin(latr)*np.sin(deltar))/(np.cos(latr)*np.cos(deltar))
    return np.degrees(np.arccos(cosra))

#----------------------------------------------------------------------
def RaDecAlt2xyz(alpha,delta,satAlt):
    # input: alpha, delta, altitude of satellites
    # output: xyz equatorial of satellites
    
    rs = constants.earthRadius+satAlt
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
    re2 = constants.earthRadius*constants.earthRadius
    
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
def satGeoVel(alpha,delta,satInc,satAlt):
    '''return geocentric velocity vector of satellites
    
    IN:
    -  alpha, delta: longitude and latitude of the satellite,
                      geocentric equatorial [deg]
    - satInc, satAlt: orbit inclination [deg] and alt [km]

    OUT:
    -  the two geocentric velocity vectors (xyz geocentric equatorial)
       for the two orbits with inc, alt that cross the alpha delta point
       (one moving N, the other S)
    '''
    
    rs = constants.earthRadius+satAlt
    
    alphar = np.radians(alpha)
    incr = np.radians(satInc)
    si = np.sin(incr)
    ci = np.cos(incr)

    
    # find nodes omega0 and omega1
    longr =  np.arcsin(np.tan(np.radians(delta))/np.tan(incr))
    
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
def satTopoVel(VS,obsLatitude):
    '''return topocentric velocity vector of satellites
    
    In:
    -   VS, geocentric equatorial velocity vectors of the satellites
    -   obsLatitude of the observatory
    
    OUT:
    -   obsvel of topocentric equatorial velocity vector, i.e.
        VS corrected for the velocity of the observatory
    '''
    
    # observatory velocity
    VO = np.array([0.,constants.earthRotation*constants.earthRadius*np.cos(np.radians(obsLatitude)),0.])

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


#------------------------------------------------------------------------------
def AzEl2Vel(alpha, delta, Delta,obsLatitude,satInc,satAlt):
    '''Apparent average angular velocity

    IN
    - alpha, delta: geocentric position of the satellite [deg]
    - Delta: distance Observatory-satellite [km]
    - obsLatitude latitude of the observatory [deg]
    - satInc, satAlt of the satellites in this shell [deg],[km]

    OUT
    -  AngularVel: apparent (from obs) average (for satellites moving
       up and down) velocity of the satellites. [deg/s]
    '''

    # geocentric coordinates of the sat
    CS = Pol2Rec((alpha,delta), constants.earthRadius+satAlt) 
    
    # geocentric coordinates of the observatory
    wCO = Pol2Rec((0.,obsLatitude), constants.earthRadius)    
    CO = np.array([[wCO[0]], [wCO[1]], [wCO[2]]])


    # topocentric coords of sat:
    OS = CS - CO

    #geocentry velocity vector of the sat,
    VS = satGeoVel(alpha, delta, satInc, satAlt)

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
def mag2microcd(mag):
    '''convert sky brightness in mag/sqarc to microcandela/m2
    
    Based on: Wong+2020 https://doi.org/10.1016/j.jlumin.2020.117256 
    '''

    return 11.E10 * 10.**(-0.4*mag)


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------


