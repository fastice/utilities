# geodat.py
import numpy as np
import pyproj
from osgeo import gdal
from utilities import myerror
import re


class geodat:

    """ Geodat object - contains information from a geodat file"""

    def __init__(self, x0=None, y0=None, xs=None, ys=None, dx=None, dy=None,
                 domain=None, verbose=True, wkt=None):
        """ Initialize geodat(x0=None,y0=None,xs=None,ys=None,dx=None,dy=None,
        domain=['greenland']) domain = greenland or antarctica """
        # set everything to zero as default
        self.x0, self.y0, self.dx, self.dy = 0.0, 0.0, 0.0, 0.0
        self.xs = self.ys = 0
        self.domain = 'greenland'
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

        self.verbose = verbose

        if domain is not None:
            self.domain = domain.lower()
            domains = ['greenland', 'antarctica']
            if self.domain not in domains:
                myerror(f'\n--- setup geodat invalid domain {self.domain} not'
                        f' in {domains}')
        #
        # setup conversions
        #
        self.llprojEPSG = "EPSG:4326"
        self.llproj = pyproj.Proj(self.llprojEPSG)  # Obsolete
        if wkt is not None:
            #self.domain = re.findall('[a-z,A-Z]+', wkt)[1].lower()
            self.xyprojSRS = wkt  # This could be a wkt or epsg
        else:
            if self.domain == 'greenland':
                self.xyprojSRS = "EPSG:3413"
            elif self.domain == 'antarctica':
                self.xyprojSRS = "EPSG:3031"
        self.xyproj = pyproj.Proj(self.xyprojSRS)
        self.lltoxyXform = pyproj.Transformer.from_crs(self.llprojEPSG,
                                                       self.xyprojSRS)
        self.xytollXform = pyproj.Transformer.from_crs(self.xyprojSRS,
                                                       self.llprojEPSG)
        if self.verbose:
            print('setting up projections', self.domain)

    def lltoxym(self, lat, lon, scale=None):
        """ lat, lon is an nparray of xy points in units of km, output is x, y
        in meters """
        if scale is None:
            scale = 1
        sh = lat.shape
        lat1 = lat.flatten()
        lon1 = lon.flatten()
        x, y = self.lltoxyXform.transform(lat1, lon1)
        x = scale * x.reshape(sh)
        y = scale * y.reshape(sh)
        return x, y

    def lltoImage(self, lat, lon, scale=None):
        x, y = self.lltoxym(lat, lon)
        return self.xymtoImage(x, y)

    def lltoxykm(self, lat, lon):
        """ lat, lon is an nparray of xy points in units of km,  output is x,
        y in KILOmeters """
        return self.lltoxym(lat, lon, scale=0.001)

    def xymtoll(self, x, y, scale=None):
        """ x, y  is an nparray of points in PS meters, output lat,lon in
        similar array """
        if scale is None:
            scale = 1
        sh = x.shape
        xx = x.flatten()
        yy = y.flatten()

        lat, lon = self.xytollXform.transform(xx * scale, yy * scale)
        lat = lat.reshape(sh)
        lon = lon.reshape(sh)
        return lat, lon

    def xymtoImage(self, x, y, scale=None):
        """ x,y  is an nparray of points in PS meters, output xi,yi image
        coordinates """
        if scale is None:
            scale = 1
        xi = (x*scale - self.x0*1000)/self.dx
        yi = (y*scale - self.y0*1000)/self.dy

        return xi, yi

    def xykmtoImage(self, x, y):
        return self.xymtoImage(x, y, scale=1000)

    def xykmtoll(self, x, y):
        """ x, y  is an nparray of points in PS KILOmeters, output lat,lon in
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
        return self.x0, self.y0, self.x0+(self.xs-1)*self.dx*0.001, \
            self.y0+(self.ys-1)*self.dy*0.001

    def boundsInM(self):
        return self.x0 * 1000, self.y0 * 1000, \
            (self.x0 + (self.xs - 1)*self.dx * 0.001) * 1000., \
            (self.y0+(self.ys-1)*self.dy*0.001)*1000.

    def originInM(self):
        return self.x0*1000, self.y0*1000

    def pixSizeInM(self):
        return self.dx, self.dy

    def pixSizeInKm(self):
        return self.dx*0.001, self.dy*0.001

    def xd(self):
        """ return geo information as an 'xd' np matrix """
        xd = np.array([[self.xs, self.ys], [self.dx, self.dy],
                       [self.x0, self.y0]])
        return xd

    def readGeodatFromTiff(self, tiffFile):
        """ Read geoinformation from a tiff file and use it to create geodat
        info - assumes PS coordinates"""
        try:
            gdal.AllRegister()
            ds = gdal.Open(tiffFile)
            self.xs = ds.RasterXSize
            self.ys = ds.RasterYSize
            gt = ds.GetGeoTransform()
            self.dx = abs(gt[1])
            self.dy = abs(gt[5])
            self.x0 = (gt[0]+self.dx/2)*0.001
            if gt[5] < 0:
                self.y0 = (gt[3] - self.ys * self.dy + self.dy/2)*0.001
            else:
                self.y0 = (gt[3] + self.dy/2)*0.001
        except Exception:
            myerror(f"Error trying to readgeodat info from : {tiffFile}")

    def readGeodat(self, geoDatFile):
        """ Read XXX.geodat file - for now skip any projection info """
        try:
            fgeo = open(geoDatFile, 'r')
            skip = [';', '#']
            count = 0
            if self.verbose:
                print('Reading ', geoDatFile)
            for line in fgeo:
                if not any(s in line for s in skip):
                    tmp = line.strip('/n').split()
                    if count == 0:
                        # int(float()) avoids problem with floating point
                        # input - e.g.,  1000.00
                        self.xs = int(float(tmp[0]))
                        self.ys = int(float(tmp[1]))
                        count += 1
                    elif count == 1:
                        self.dx = float(tmp[0])
                        self.dy = float(tmp[1])
                        count += 1
                    elif count == 2:
                        self.x0 = float(tmp[0])
                        self.y0 = float(tmp[1])
                        break
            fgeo.close
        except Exception:
            myerror('Error reading geodat file '+geoDatFile)

    def writeGeodat(self, geoDatFile):
        """ Write a geodat file """
        try:
            fgeo = open(geoDatFile, 'w')
            print('# 2', file=fgeo)
            print(';\n; Image size (pixels) nx ny\n;', file=fgeo)
            print('{:d} {:d}'.format(self.xs, self.ys), file=fgeo)
            print(';\n; Pixel size (m) deltaX deltaY\n;', file=fgeo)
            print('{:.4f} {:.4f}'.format(self.dx, self.dy), file=fgeo)
            print(';\n; Origin, lower left corner (km) Xo  Yo\n;', file=fgeo)
            print('{:.4f} {:.4f}'.format(self.x0, self.y0), file=fgeo)
            print('&', file=fgeo)
            fgeo.close()
        except Exception:
            myerror('Error writing geodat file ', geoDatFile)

    def print(self):
        print('\nSize in pixels:\t', self.xs, self.ys)
        print('Pixel size in m:\t', self.dx, self.dy)
        print('Origin x0, y0:\t\t', self.x0, self.y0, '\n')
