# geodat.py
import numpy as np
import math
import os
from datetime import datetime
import scipy.interpolate as interp
import pyproj


class geodatrxa:

    """ Geodat object - contains information from a geodat file"""

    def __init__(self, file=None, echo=False):
        """ initialize a geodatrxa object, where:
        file\t is optional file name, can be input later with '
        readFile(file=file) echo\t set to true to echo results as they are'
        read in, otherwise no output """
        #
        # set everything to empty values
        self.file = ''
        self.datetime, self.midnight, self.lookdir = [], [], []
        self.t0, self.t1 = None, None
        self.nr, self.na, self.nlr, self.nla = -1, -1, -1, -1
        self.Rn, self.Rf = None, None
        self.ReMajor, self.ReMinor, self.Rc, self.phic = -1, -1, -1, -1
        self.H, self.prf, self.wavelength = -1, -1, -1
        self.corners = []
        self.slpRg, self.slpAz = -1, -1
        self.ascdesc = self.lookdir = []
        self.deltaT, self.deltaR = np.nan, np.nan
        self.rNearSLP = None
        self.stateTime = None
        self.nState, self.tState, self.dTState = -1, -1, -1
        self.position = self.velocity = []
        self.llzToEcef = None

        self.minT, self.maxT = -1, -1
        self.fx, self.fy, self.fz, self.fvx, self.fvy, self.fvz = [None]*6
        # in most cases all or no args would be passe.
        if file is not None:
            self.file = file
            self.readFile(echo=echo)

    # return resolution
    def singleLookResolution(self):
        return self.slpRg, self.slpAz

    def singleLookSize(self):
        return self.nr*self.nlr, self.na*self.nla

    def centerRangem(self):
        return self.Rc*1000.

    def satelliteAltm(self):
        return self.H*1000.

    def TisInImage(self, myTime):
        if myTime >= self.t0 and myTime <= self.t1:
            return True
        return False

    def nearRangem(self):
        if self.Rn is None:
            self.Rn = self.Rc - (self.nr-1)*self.nlr*self.slpRg*0.5*0.001
        return self.Rn*1000.

    def farRangem(self):
        if self.Rf is None:
            self.Rf = self.Rc + (self.nr-1)*self.nlr*self.slpRg*0.5*0.001
        return self.Rf*1000.

    def validRangem(self, R):
        ''' check range in image '''
        if R >= self.nearRangem() and R <= self.farRangem():
            return True
        return False

    def thetaAtTRrad(self, myTime, R, z=0):
        ''' compute theta for given time and range '''
        if not self.TisInImage(myTime) or not self.validRangem(R):
            return -1.0
        # in in image so calculate
        ReH = np.linalg.norm(self.interpPos(myTime))
        # Caution this is using Re from scene center
        Re = self.earthRadm()
        print(ReH)
        #
        print(R, ReH, Re)
        print((R**2 + ReH**2 - (Re+z)**2)/(2*ReH*R))
        print(math.acos((R**2 + ReH**2 - (Re+z)**2)/(2*ReH*R)) * 180./np.pi)
        return math.acos((R**2 + ReH**2 - (Re+z)**2)/(2*ReH*R))

    def thetaCActualrad(self):
        ''' Uses actual ReH rather than nominal H from geodat file '''

        return self.thetaAtTRrad(0.5*(self.t0+self.t1), self.centerRangem())

    def thetaCrad(self):
        alt = self.satelliteAltm()
        r1 = self.centerRangem()
        Re = self.earthRadm()
        theta = r1*r1 + 2.0*alt*Re + alt*alt
        theta = theta/(2.0*r1*(Re+alt))
        theta = math.acos(theta)
        print(Re+alt)
        print('py', r1, alt, Re, self.corners[4, 0])
        return theta

    def earthRadm(self):
        ''' use radius at the center of the image '''
        lat = np.radians(self.corners[4, 0])
        return self.earthRadLatm(lat)

    def earthRadLatm(self, lat):
        ''' use radius at the center of the image '''
        N = (self.ReMajor**2.0) / np.sqrt((self.ReMajor*np.cos(lat))**2.0 +
                                          (self.ReMinor*np.sin(lat))**2.0)
        N = N*1000.
        x = N*np.cos(lat)
        z = (self.ReMinor/self.ReMajor)**2.0 * N * np.sin(lat)
        Re = np.sqrt(x**2 + z**2)
        return Re

    def isSouth(self):
        """ Return True if right looking pass - error if lookdir
        does not exist or has bad value """
        if len(self.corners) == 0:
            print('geodatrxa.isDescending: no corners given')
            exit()
        return True if self.corners[0, 0] < 0 else False

    def isRightLooking(self):
        """ Return True if right looking pass - error if lookdir  does not
        exist or has bad value """
        if (len(self.lookdir) == 0 or not (self.lookdir.lower()
                                           in ['left', 'right'])):
            print('geodatrxa.isDescending: no ascdesc value given')
            exit()
        return True if self.lookdir.lower() == 'right' else False

    def isDescending(self):
        """ Return True if descending pass -
        error if ascdesc does not exist """
        if len(self.ascdesc) == 0 or not (self.ascdesc.lower() in
                                          ['ascending', 'descending']):
            print('geodatrxa.isDescending: no ascdesc value given')
            exit()
        return True if self.ascdesc.lower() == 'descending' else False

    def isAscending(self):
        """ Return True if ascending pass - error if ascdesc does not exist """
        if len(self.ascdesc) == 0 or not (self.ascdesc.lower() in
                                          ['ascending', 'descending']):
            print('geodatrxa.isDescending: no ascdesc value given')
            exit()
        return True if self.ascdesc.lower() == 'ascending' else False

    def readFile(self, file=None, echo=False):
        if file is not None:
            self.file = file
        # check file exists
        if os.path.exists(self.file):
            fp = open(self.file, 'r')
        else:
            print('Attempted to open geodat file that does not exist')
            exit()
        #
        ncorners = 0
        timeSet = False
        nState = 0
        position = True
        self.corners = np.zeros((5, 2))
        #
        # Loop through lines
        #
        for line in fp:
            # image date is the only comment item to extract
            if '; Image date' in line:
                tmp = line.split(':')[-1].strip()
                self.midnight = datetime.strptime(tmp, "%d %b %Y")
            # extract all non-commented data,  assuming its in a set order
            elif 'Skew' in line:
                self.skew, self.squint = [float(x) for x in line.split()[7:9]]
            elif ';' not in line:
                tmp = line.strip().split()
                if self.nr == -1 and len(tmp) == 4:
                    self.nr, self.na = int(tmp[0]), int(tmp[1])
                    self.nlr, self.nla = int(tmp[2]), int(tmp[3])
                    if echo:
                        print('nr,na,nlr,nla ', self.nr, self.na,
                              self.nlr, self.nla)
                elif self.ReMajor == -1 and (len(tmp) == 5 or len(tmp) == 6):
                    self.ReMajor = float(tmp[0])
                    self.ReMinor = float(tmp[1])
                    self.Rc = float(tmp[2])
                    self.phic = float(tmp[3])
                    self.H = float(tmp[4])
                    if len(tmp) == 6:
                        self.deltaR = float(tmp[5])
                    if echo:
                        print('ReMajor, ReMinor, Rc, phic, h, deltaR: ',
                              self.ReMajor, self.ReMinor, self.Rc,
                              self.phic, self.H, self.deltaR)
                elif ncorners < 5 and len(tmp) == 2:
                    self.corners[ncorners, :] = float(tmp[0]), float(tmp[1])
                    ncorners += 1
                    if ncorners == 5 and echo:
                        print(self.corners)
                elif self.slpRg == -1 and len(tmp) == 2:
                    self.slpRg, self.slpAz = float(tmp[0]), float(tmp[1])
                    if echo:
                        print('Single Look Pix Size (r, a)',
                              self.slpRg, self.slpAz)
                elif (tmp[0] in ['descending', 'ascending'] and
                      len(self.ascdesc) == 0):
                    self.ascdesc = tmp[0]
                    if echo:
                        print(self.ascdesc)
                elif tmp[0] in ['right', 'left'] and len(self.lookdir) == 0:
                    self.lookdir = tmp[0]
                    if echo:
                        print('Look direction ', self.lookdir)
                elif tmp[0] in ['state']:
                    if echo:
                        print(tmp[0])
                elif len(tmp) == 3 and not timeSet:
                    hour, minute = int(tmp[0]), int(tmp[1])
                    second = int(float(tmp[2]))
                    microsecond = int((float(tmp[2])-second)*1e6)
                    # this is  a kluge for case where squint time pushes over
                    # 24 hour boundary - not a problem in most casese
                    if hour > 23:
                        hour, minute, second = 23, 59, 59
                    self.datetime = \
                        self.midnight.replace(hour=hour, minute=minute,
                                              second=second,
                                              microsecond=microsecond)
                    self.t0 = (self.datetime - self.midnight).total_seconds()
                    if echo:
                        print('Date/Time ', self.datetime)
                    timeSet = True
                elif len(tmp) == 1 and self.prf == -1:
                    self.prf = float(tmp[0])
                    if echo:
                        print('Prf ', self.prf)
                elif len(tmp) == 1 and self.wavelength == -1:
                    self.wavelength = float(tmp[0])
                    if echo:
                        print('wavelength ', self.wavelength)
                elif len(tmp) == 1 and self.nState == -1:
                    self.nState = int(tmp[0])
                    if echo:
                        print('nState ', self.nState)
                    self.position = np.zeros((self.nState, 3))
                    self.velocity = np.zeros((self.nState, 3))
                elif len(tmp) == 1 and self.tState == -1:
                    self.tState = float(tmp[0])
                    if echo:
                        print('tState ', self.tState)
                elif len(tmp) == 1 and self.dTState == -1:
                    self.dtState = float(tmp[0])
                    if echo:
                        print('dtState ', self.dtState)
                elif nState < self.nState and len(tmp) == 3:
                    if position:
                        self.position[nState, :] = [float(tmp[0]),
                                                    float(tmp[1]),
                                                    float(tmp[2])]
                        position = False
                    else:
                        self.velocity[nState, :] = [float(tmp[0]),
                                                    float(tmp[1]),
                                                    float(tmp[2])]
                        position = True
                        nState += 1
                elif 'deltaT' in line:
                    self.deltaT = float(tmp[1])
                    if echo:
                        print('deltaT ', self.deltaT)
        if echo:
            print('Position (x,y,z): \n', self.position)
            print('Velocity (vx,vy,vz): \n', self.velocity)
        #
        self.stateTime = np.array([self.tState + x*self.dtState
                                   for x in range(0, self.nState)])
        # compute near in slp coordates for geocoding
        self.rNearSLP = self.Rc * 1000.0 - ((self.nr-1)/2) * \
            self.nlr * self.slpRg - (self.nlr-1) * self.slpRg/2.
        self.t1 = self.t0 + (self.na-1) * self.nla/self.prf
    # geocoding methods

    def printState(self):
        print(self.tState, self.dtState, self.position.shape)
        for t, pos, vel in zip(self.stateTime, self.position, self.velocity):
            print(f'{t:10.4f} {pos[0]:8.1f} {pos[1]:8.1f} {pos[2]:8.1f} '
                  f'{vel[0]:8.1f} {vel[1]:8.1f} {vel[2]:8.1f}')

    def setupInterpState(self):
        ''' setup interpolators '''
        kind = 'cubic'
        bError = False
        self.fx = interp.interp1d(self.stateTime, self.position[:, 0],
                                  kind=kind, fill_value=np.nan,
                                  bounds_error=bError)
        self.fy = interp.interp1d(self.stateTime, self.position[:, 1],
                                  kind=kind, fill_value=np.nan,
                                  bounds_error=bError)
        self.fz = interp.interp1d(self.stateTime, self.position[:, 2],
                                  kind=kind, fill_value=np.nan,
                                  bounds_error=bError)
        self.fvx = interp.interp1d(self.stateTime, self.velocity[:, 0],
                                   kind=kind, fill_value=np.nan,
                                   bounds_error=bError)
        self.fvy = interp.interp1d(self.stateTime, self.velocity[:, 1],
                                   kind=kind, fill_value=np.nan,
                                   bounds_error=bError)
        self.fvz = interp.interp1d(self.stateTime, self.velocity[:, 2],
                                   kind=kind, fill_value=np.nan,
                                   bounds_error=bError)

    def interpPos(self, t):
        '''Interp state position vectors; For now only work with scalars '''
        if self.fx is None:
            self.setupInterpState()
        #
        x, y, z = self.fx(t), self.fy(t), self.fz(t)
        if np.isscalar(t):
            return [np.asscalar(x), np.asscalar(y), np.asscalar(z)]
        else:
            return np.array([x, y, z])

    def interpVel(self, t):
        '''Interp state velocity vectors; For now only work with scalars '''
        if self.fvx is None:
            self.setupInterpState()
        #
        vx, vy, vz = self.fvx(t), self.fvy(t), self.fvz(t)
        if np.isscalar(t):
            return [np.asscalar(vx), np.asscalar(vy), np.asscalar(vz)]
        else:
            return np.array([vx, vy, vz])

    def lltoecef(self, lat, lon, zelev):
        ''' convert llz to ecef'''
        #    return pyproj.transform(self.llz, self.ecef, lon, lat, zelev,
        #                            radians=False)
        if self.llzToEcef is None:
            self.ecef = pyproj.Proj(proj='geocent', ellps='WGS84',
                                    datum='WGS84')
            self.llz = pyproj.Proj(proj='latlong', ellps='WGS84',
                                   datum='WGS84')
            self.llzToEcef = pyproj.Transformer.from_proj(self.llz, self.ecef)

        return self.llzToEcef.transform(lat, lon, zelev, radians=False)

    def ReH(self, myTime):
        sPt = np.array(self.interpPos(myTime))
        return np.linalg.norm(sPt)

    def llzPtToRA(self, lat, lon, z, initT=None):
        ''' geocode from lat/lon/z to range/azimuth coordinates '''
        tPt = self.lltoecef(lat, lon, z)
        if initT is None:
            initT = self.t0 + 0.5 * self.na * self.nla / self.prf
        print(initT)
        myTime = initT
        for i in range(0, 50):
            sPt = np.array(self.interpPos(myTime))
            vPt = np.array(self.interpVel(myTime))
            #
            dr = tPt - sPt
            df = np.dot(dr, vPt)
            C1 = -np.dot(vPt, vPt)
            # assume zero dop geom for now
            C2 = 0
            myTime -= df/(C1+C2)
            # print(myTime, np.sqrt(np.dot(der, dr)), i)
            if np.abs(df/C1) < 1e-5:
                break
        # correct to near range for multi-look,  the sub
        r = (np.sqrt(np.dot(dr, dr)) - self.rNearSLP)/self.slpRg
        az = (myTime - self.t0) * self.prf
        print(i)  # print(r,az,i)
        return r, az, myTime
