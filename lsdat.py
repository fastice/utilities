# lsdat.py
import numpy as np
import pyproj
from datetime import datetime, timedelta
from utilities import myerror, mywarning


class lsdat:

    """ lsdat object - contains information from a landsat match dat file"""

    def __init__(self, x0=None, y0=None, xs=None, ys=None, dx=None, dy=None,
                 domain=None, noProjection=False):
        """ Initialize lsdat(x0=None,y0=None,xs=None,ys=None,dx=None,dy=None,
        domain = greenland or antarctica """
        # set everything to zero as default
        self.x0, self.y0, self.dx, self.dy = 0.0, 0.0, 0.0, 0.0
        self.successRate, self.culledRate = -1.0, -1.0
        self.slowFlag = False
        self.xs, self.ys = 0, 0
        self.stepX, self.stepY = 0, 0
        self.dxImage, self.dyImage = 0., 0.
        self.fileEarly, self.fileLate = '', ''
        self.sigmaX, self.sigmaY = 0.0, 0.0
        self.domain = 3413
        # in most cases all or no args would be passed.
        if x0 is not None:
            self.x0 = x0
        if y0 is not None:
            self.y0 = y0
        if xs is not None:
            self.xs = int(xs)
        if ys is not None:
            self.ys = int(ys)
        if dx is not None:
            self.dx = dx
        if dx is not None:
            self.dy = dy
        if domain is not None:
            self.domain = domain.lower()
            domains = [3413, 3031]
            if self.domain not in domains:
                print(f'\n--- setup lsdat invalid domain {self.domain} not '
                      f'in {domains}')
                exit()
        #
        # setup conversions
        # This step is time consuming, so opt out if not needed.
        if not noProjection:
            self.llproj = pyproj.Proj("EPSG:4326")
            if self.domain == 3413:
                self.xyproj = pyproj.Proj("EPSG:3413")
            elif self.domain == 3031:
                self.xyproj = pyproj.Proj("EPSG:3031")
            self.llxyXform = pyproj.Transformer.from_crs("EPSG:4326",
                                                         f'EPSG:{self.domain}')
            self.xyllXform = pyproj.Transformer.from_crs(f'EPSG:{self.domain}',
                                                         "EPSG:4326")
#        print('setting up projections', self.domain)

    def lltoxym(self, lat, lon, scale=None):
        """ lat, lon is an nparray of xy points in units of km, output is x,y
        in meters """
        if scale is None:
            scale = 1
        try:
            sh = lat.shape
        except Exception:
            lat = np.array(lat)
            lon = np.array(lon)
            sh = lat.shape

        lat1 = lat.flatten()
        lon1 = lon.flatten()
        x, y = self.llxyXform.transform(lat1, lon1)
        x = scale * x.reshape(sh)
        y = scale * y.reshape(sh)
        return x, y

    def lltoxykm(self, lat, lon, scale=None):
        """ lat, lon is an nparray of xy points in units of km, output is x,y
        in KILOmeters """
        return self.lltoxym(lat, lon, scale=0.001)

    def xymtoll(self, x, y, scale=None):
        """ x, y  is an nparray of points in PS meters, output lat,lon in
        similar array """
        if scale is None:
            scale = 1
        try:
            sh = x.shape
        except Exception:
            x = np.array(x)
            y = np.array(y)
            sh = x.shape
        xx = x.flatten()
        yy = y.flatten()
        lat, lon = self.xyllXform.transform(xx * scale, yy * scale)
        lat = lat.reshape(sh)
        lon = lon.reshape(sh)
        return lat, lon

    def xykmtoll(self, x, y):
        """ x,y  is an nparray of points in PS KILOmeters, output lat,lon in
        similar array """
        return self.xymtoll(x, y, scale=1000)

    def sizeInPixels(self):
        return self.xs, self.ys

    def sizeInKm(self):
        return self.xs*self.dx*0.001, self.ys*self.dy*0.001

    def sizeInM(self):
        return self.xs*self.dx, self.ys*self.dy

    def originInKm(self):
        return self.x0, self.y0

    def boundsInKm(self):
        return self.x0, self.y0, self.x0 + (self.xs - 1) * self.dx * 0.001, \
            self.y0 + (self.ys - 1) * self.dy * 0.001

    def boundsInM(self):
        return self.x0 * 1000, self.y0 * 1000, \
            (self.x0 + (self.xs - 1) * self.dx * 0.001) * 1000., \
            (self.y0 + (self.ys - 1) * self.dy * 0.001) * 1000.

    def originInM(self):
        return self.x0*1000, self.y0*1000

    def pixSizeInM(self):
        return self.dx, self.dy

    def pixSizeInKm(self):
        return self.dx*0.001, self.dy*0.001

    def xd(self):
        """ return ls information as an 'xd' np matrix """
        xd = np.array([[self.xs, self.ys], [self.dx, self.dy],
                       [self.x0, self.y0]])
        return xd

    def tryOptParam(self, key, myDict, myType):
        ''' check if key in dict before attempting to index.return typed value
        '''
        if key in myDict.keys():
            return myType(myDict[key])
        return None

    def readLSdat(self, LSDatFile, printData=False, warning=False):
        myDict = {}
        try:
            fp = open(LSDatFile, 'r')
        # crude check sum
            for line in fp:
                pieces = line.split('=')
                if len(pieces) == 2:
                    myDict[pieces[0].strip()] = pieces[1].strip()
            fp.close()
            # now stuff values
            self.fileEarly = myDict['fileEarly']
            self.fileLate = myDict['fileLate']
            self.x0 = float(myDict['x0'])/1000.
            self.y0 = float(myDict['y0'])/1000.
            self.dxImage = float(myDict['dx'])
            self.dyImage = float(myDict['dy'])
            self.stepX = int(myDict['stepX'])
            self.stepY = int(myDict['stepY'])
            self.xs = int(myDict['nx'])
            self.ys = int(myDict['ny'])
            self.slowFlag = int(myDict['slowFlag']) > 0
            self.domain = int(myDict['EPSG'])
            self.JD1 = float(myDict['earlyImageJD'])
            self.JD2 = float(myDict['lateImageJD'])
            self.successRate = \
                float(myDict['Success_rate_for_attempted_matches(%)'])
            self.culledRate = \
                self.tryOptParam('Culled_rate_for_attempted_matches(%)',
                                 myDict, float)
            self.sigmaX = self.tryOptParam('Mean_sigmaX', myDict, float)
            self.sigmaY = self.tryOptParam('Mean_sigmaY', myDict, float)
            self.dx = self.dxImage*self.stepX
            self.dy = self.dyImage*self.stepY
            self.date1 = datetime.strptime('2000:01:01', "%Y:%m:%d") + \
                timedelta(days=self.JD1 - 2451544.5)
            self.date2 = datetime.strptime('2000:01:01', "%Y:%m:%d") + \
                timedelta(days=self.JD2 - 2451544.5)
            #
            if printData:
                print(self.fileEarly)
                print(self.fileLate)
                print(self.x0, self.y0)
                print(self.dx, self.dy)
                print(self.stepX, self.stepY)
                print(self.xs, self.ys)
                print(self.domain)
                print(self.date1, self.date2)
                print(self.successRate, self.culledRate)
                print(self.sigmaX, self.sigmaY)
        except Exception:
            if not warning:
                myerror("Problem reading lsdat file "+LSDatFile)
            else:
                mywarning("Problem reading lsdat file "+LSDatFile)
