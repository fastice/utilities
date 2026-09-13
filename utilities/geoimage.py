# geoimage.py
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from utilities.readImage import readImage
from utilities.writeImage import writeImage
from utilities.myerror import myerror
from utilities.mywarning import mywarning
from utilities import geodat
import os
from osgeo import gdal, gdal_array, osr
from datetime import datetime

# Per-type configuration: components, file suffixes, geodat suffix, magnitude attribute name.
# 'magAttr' is the attribute that stores the magnitude (None for scalar).
_TYPE_CONFIG = {
    'scalar':     {'components': ['x'],        'suffixes': [''],
                   'geodatSuffix': '.geodat',    'magAttr': None},
    'velocity':   {'components': ['vx', 'vy'], 'suffixes': ['.vx', '.vy'],
                   'geodatSuffix': '.vx.geodat', 'magAttr': 'v'},
    'velocityRA': {'components': ['vr', 'va'], 'suffixes': ['.vr', '.va'],
                   'geodatSuffix': '.vr.geodat', 'magAttr': 'v'},
    'error':      {'components': ['ex', 'ey'], 'suffixes': ['.ex', '.ey'],
                   'geodatSuffix': '.ex.geodat', 'magAttr': 'e'},
    'errorRA':    {'components': ['er', 'ea'], 'suffixes': ['.er', '.ea'],
                   'geodatSuffix': '.er.geodat', 'magAttr': 'e'},
}

# noData values keyed by file suffix (used in writeCloudOptGeo)
_NO_DATA = {
    '.vx': -2.0e9, '.vy': -2.0e9,
    '.vr': -2.0e9, '.va': -2.0e9,
    '.v':  -1.0,
    '.ex': -1.0,   '.ey': -1.0,
    '.er': -1.0,   '.ea': -1.0,
    '.e':  -1.0,
    '':    None,
}


def _tiffName(fileName, suffix):
    ''' build fileName + suffix + '.tif', but avoid a redundant
    '.tif.tif' when fileName is already a .tif and suffix is empty
    (i.e., a scalar geoType with no per-component suffix) '''
    if suffix == '' and fileName.endswith('.tif'):
        return fileName
    return fileName + suffix + '.tif'


class geoimage:
    """ \ngeoimage - object for scalar or velocity PS data + geodat data """

    def __init__(self, x=None, vx=None, vy=None, v=None, ex=None, ey=None,
                 e=None, geoType=None, verbose=True):

        self.x = []
        self.vx, self.vy, self.v = [], [], []
        self.vr, self.va = [], []
        self.ex, self.ey, self.e = [], [], []
        self.geo = []
        self.xx, self.yy = [], []
        self.xGrid, self.yGrid = [], []
        self.geoType = None
        self.velDate = None
        self.fileName = None
        self.fileRoot = None
        if x is not None:
            self.x = x
        if vx is None:
            self.vx = vx
        if vy is None:
            self.vy = vy
        if v is None:
            self.v = v
        if ex is None:
            self.ex = ex
        if ey is None:
            self.ey = ey
        if e is None:
            self.e = e

        self.verbose = verbose
        if geoType is not None:
            self.geoType = geoType
            if self.verbose:
                print('Type ', self.geoType)

    def setGeoType(self, geoType):
        """ setGeoType(geoType) set type to velocity, velocityRA, error, or scalar"""
        if geoType not in _TYPE_CONFIG:
            print(f'\n\n\tgeoImage setType, invalid type: {geoType}\n\n')
            exit()
        self.geoType = geoType

    def xyCoordinates(self):
        """ xyCoordinates - setup xy coordinates in km """
        sx, sy = self.geo.sizeInPixels()
        x0, y0 = self.geo.originInKm()
        dx, dy = self.geo.pixSizeInKm()
        self.xx = np.arange(x0, x0+sx*dx, dx)
        self.yy = np.arange(y0, y0+sy*dy, dy)
        self.xx, self.yy = self.xx[0:sx], self.yy[0:sy]
        self.extent = [min(self.xx), max(self.xx), min(self.yy), max(self.yy)]

    def xyGrid(self):
        if len(self.xx) == 0:
            self.xyCoordinates()
        sx, sy = self.geo.sizeInPixels()
        self.xGrid, self.yGrid = np.zeros((sy, sx)), np.zeros((sy, sx))
        for i in range(0, sy):
            self.xGrid[i, :] = self.xx
        for i in range(0, sx):
            self.yGrid[:, i] = self.yy

    def parseMyMeta(self, metaFile):
        print(metaFile)
        fp = open(metaFile)
        dates = []
        for line in fp:
            if 'MM:DD:YYYY' in line:
                tmp = line.split('=')[-1].strip()
                dates.append(datetime.strptime(tmp, "%b:%d:%Y"))
                if len(dates) == 2:
                    break
        if len(dates) != 2:
            return None
        fp.close()
        return dates[0]+(dates[1]-dates[0])*0.5

    def parseVelCentralDate(self):
        if self.fileName is None:
            metaFile = self.fileName + '.meta'
            if not os.path.exists(metaFile):
                return None
            return self.parseMyMeta(metaFile)
        return None

    def setupInterp(self, method='linear'):
        """ set up interpolation for scalar (xInterp) or two-component types """
        if len(self.xx) < 0:
            myerror('\n\nsetupInterp: x, y limits not set\n\n')
        xy = (self.yy, self.xx)
        cfg = _TYPE_CONFIG[self.geoType]
        for comp in cfg['components']:
            setattr(self, f'{comp}Interp',
                    RegularGridInterpolator(xy, getattr(self, comp), method=method))
        if cfg['magAttr']:
            setattr(self, f'{cfg["magAttr"]}Interp',
                    RegularGridInterpolator(xy, getattr(self, cfg['magAttr']), method=method))

    def interpGeo(self, x, y):
        """ interpolate velocity or x at points x and y, which are in km
        (note x,y is c-r even though data r-c).
        Returns: scalar → single array; others → tuple(comp0, comp1, mag) """
        shapeSave = x.shape
        x1 = x.flatten()
        y1 = y.flatten()
        xgood = np.logical_and(x1 >= self.xx[0], x1 <= self.xx[-1])
        ygood = np.logical_and(y1 >= self.yy[0], y1 <= self.yy[-1])
        igood = np.logical_and(xgood, ygood)
        xy = np.array([y1[igood], x1[igood]]).transpose()

        cfg = _TYPE_CONFIG[self.geoType]

        if self.geoType == 'scalar':
            result = np.full(x1.shape, np.nan)
            result[igood] = self.xInterp(xy)
            return np.reshape(result, shapeSave)

        results = []
        for comp in cfg['components']:
            r = np.full(x1.shape, np.nan)
            r[igood] = getattr(self, f'{comp}Interp')(xy)
            results.append(np.reshape(r, shapeSave))
        if cfg['magAttr']:
            mag = np.full(x1.shape, np.nan)
            mag[igood] = getattr(self, f'{cfg["magAttr"]}Interp')(xy)
            results.append(np.reshape(mag, shapeSave))
        return tuple(results)

    def readGeodat(self, geoFile):
        if self.verbose:
            print(geoFile)
        if self.geo == []:
            self.geo = geodat(verbose=self.verbose)
        if os.path.exists(geoFile):
            self.geo.readGeodat(geoFile)
        else:
            myerror('Missing geodat file '+geoFile)

    def readMyTiff(self, tiffFile, band=1):
        """ read a tiff file and return the array """
        try:
            gdal.AllRegister()
            ds = gdal.Open(tiffFile)
            band = ds.GetRasterBand(band)
            arr = band.ReadAsArray()
            arr = np.flipud(arr)
            ds = None
        except Exception:
            myerror("geoimage.readMyTiff: error reading tiff file "+tiffFile)
        return arr

    def _multibandComponents(self, fileName):
        """ If fileName is itself a single GDAL raster with one band per
        component of the current geoType (matched by band Description, e.g.
        'vx'/'vy' for velocity), return {component: bandIndex}. Otherwise
        return None. Lets a single multi-band file (e.g. a modern VRT
        wrapping a velocity map as vx/vy bands) stand in for the usual
        per-component .tif files -- see readData(). """
        cfg = _TYPE_CONFIG[self.geoType]
        try:
            ds = gdal.Open(fileName)
        except Exception:
            return None
        if ds is None:
            return None
        bandMap = {}
        for b in range(1, ds.RasterCount + 1):
            desc = ds.GetRasterBand(b).GetDescription()
            if desc in cfg['components']:
                bandMap[desc] = b
        ds = None
        if len(bandMap) == len(cfg['components']):
            return bandMap
        return None

    def getWKT_PROJ(self, epsgCode, wktFile):
        ''' get wkt'''
        if epsgCode is None and wktFile is None:
            return None
        if wktFile is not None:
            wkt = self.readWKT(wktFile)
        else:
            sr = osr.SpatialReference()
            sr.ImportFromEPSG(epsgCode)
            wkt = sr.ExportToWkt()
        return wkt

    def readWKT(self, wktFile):
        ''' get wkt from a file '''
        with open(wktFile, 'r') as fp:
            return fp.readline()

    def imageSize(self):
        firstComp = _TYPE_CONFIG[self.geoType]['components'][0]
        ny, nx = getattr(self, firstComp).shape
        return nx, ny

    def computePixEdgeCornersXYM(self):
        nx, ny = self.imageSize()
        x0, y0 = self.geo.originInM()
        dx, dy = self.geo.pixSizeInM()
        xll, yll = x0 - dx/2, y0 - dx/2
        xur, yur = xll + nx * dx, yll + ny * dy
        xul, yul = xll, yur
        xlr, ylr = xur, yll
        corners = {'ll': {'x': xll, 'y': yll}, 'lr': {'x': xlr, 'y': ylr},
                   'ur': {'x': xur, 'y': yur}, 'ul': {'x': xul, 'y': yul}}
        return corners

    def computePixEdgeCornersLL(self):
        corners = self.computePixEdgeCornersXYM()
        llcorners = {}
        for myKey in corners.keys():
            lat, lon = self.geo.xymtoll(np.array([corners[myKey]['x']]),
                                        np.array([corners[myKey]['y']]))
            llcorners[myKey] = {'lat': lat[0], 'lon': lon[0]}
        return llcorners

    def writeMyTiff(self, tiffFile, epsg=None, noDataDefault=None,
                    predictor='YES', noV=False, overviews=None,
                    driverName='COG', wktFile=None, computeStats=True,
                    resampling='AVERAGE', bigTiff=False):
        """ write geotiff(s) for this geoimage.
        tiffFile should not have a '.tif' extension. """
        cfg = _TYPE_CONFIG[self.geoType]
        if wktFile is None:
            epsg = [epsg, 3413][epsg is None]
        firstComp = getattr(self, cfg['components'][0])
        try:
            gdalType = gdal_array.NumericTypeCodeToGDALTypeCode(firstComp.dtype)
        except Exception:
            myerror('writeMyTiff: invalid geoType ' + self.geoType)

        suffixes = list(cfg['suffixes'])
        if cfg['magAttr'] == 'v':
            suffixes.append('.v')
        elif cfg['magAttr'] == 'e':
            suffixes.append('.e')

        try:
            for suffix in suffixes:
                if noV and suffix in ('.v', '.e'):
                    continue
                self.writeCloudOptGeo(tiffFile, suffix, epsg, gdalType,
                                      overviews=overviews, predictor=predictor,
                                      noDataDefault=noDataDefault,
                                      driverName=driverName, wktFile=wktFile,
                                      computeStats=computeStats,
                                      resampling=resampling, bigTiff=bigTiff)
        except Exception:
            myerror(f"geoimage.writeMyTiff: error writing file {tiffFile}")

    def writeCloudOptGeo(self, tiffFile, suffix, epsg, gdalType,
                         overviews=None, predictor='YES', noDataDefault=None,
                         bigTiff=False, driverName='COG', wktFile=None,
                         computeStats=True, resampling='AVERAGE'):
        ''' write a cloud-optimized or plain geotiff.
        Set driverName to GTiff for a plain geotiff '''
        if driverName not in ['COG', 'GTiff']:
            myerror(f'invalid driver for writeCloudOptGeo {driverName}')
        noData = _NO_DATA.get(suffix, noDataDefault)
        driver = gdal.GetDriverByName("MEM")
        nx, ny = self.imageSize()
        dx, dy = self.geo.pixSizeInM()
        dst_ds = driver.Create('', nx, ny, 1, gdalType)
        tiffCorners = self.computePixEdgeCornersXYM()
        dst_ds.SetGeoTransform((tiffCorners['ul']['x'], dx, 0,
                                tiffCorners['ul']['y'], 0, -dy))
        wkt = self.getWKT_PROJ(epsg, wktFile)
        dst_ds.SetProjection(wkt)
        if noData is not None:
            if self.geoType == 'scalar':
                tmp = self.x
            else:
                tmp = getattr(self, suffix.replace('.', ''))
            tmp[np.isnan(tmp)] = noData
            dst_ds.GetRasterBand(1).SetNoDataValue(noData)
        if self.geoType == 'scalar':
            dst_ds.GetRasterBand(1).WriteArray(np.flipud(self.x))
        else:
            dst_ds.GetRasterBand(1).WriteArray(
                np.flipud(getattr(self, suffix.replace('.', ''))))
        if computeStats:
            # An all-noData band (e.g. a velocityStats mean whose only
            # contributors were blank Exclude.pending maps) makes
            # GetStatistics raise -- degrade to a no-stats write rather than
            # letting writeMyTiff's blanket except turn it into a fatal
            # myerror for the whole program.
            try:
                _ = dst_ds.GetRasterBand(1).GetStatistics(0, 1)
            except RuntimeError:
                mywarning(f'writeCloudOptGeo: no valid pixels for statistics '
                          f'in {tiffFile}{suffix} -- writing without stats')
        bigTiffFlag = ["NO", "YES"][bigTiff]
        options = [f'BIGTIFF={bigTiffFlag}', 'COMPRESS=LZW']
        if driverName == 'GTiff':
            if type(predictor) != int and predictor is not None:
                options.append(f'PREDICTOR={1}')
            if overviews is not None:
                if len(overviews) < 2:
                    myerror(f'Overviews {overviews} should be [2, 4, ..])')
                options.append('COPY_SRC_OVERVIEWS=YES')
                dst_ds.BuildOverviews(resampling, overviews)
        else:
            options.append('GEOTIFF_VERSION=1.1')
            options.append(f'RESAMPLING={resampling}')
            if predictor in ['YES', 'NO']:
                options.append(f'PREDICTOR={predictor}')
        dst_ds.FlushCache()
        driver = gdal.GetDriverByName(driverName)
        dst_ds2 = driver.CreateCopy(_tiffName(tiffFile, suffix), dst_ds,
                                    options=options)
        dst_ds2.FlushCache()
        dst_ds, dst_ds2 = None, None

    def writeMyVrt(self, tiffFile, vrtFile=None):
        """ Build a multi-band VRT (tiffFile + '.vrt', or vrtFile if given)
        wrapping the per-component GeoTIFFs already written by writeMyTiff(),
        one band per component (magnitude excluded), with each band's
        Description set to its component name (e.g. 'vr'/'va') -- matches the
        vx+vy/vr+va/ex+ey-pair VRT convention used elsewhere (e.g. mosaic3d's
        write3DFlatVRTs()). tiffFile should not have a '.tif' extension.
        vrtFile lets the VRT be named independently of the tiff basename
        (e.g. a sigma pair sharing froot with the mean pair needs its own
        '<froot>.err.vrt' rather than colliding on '<froot>.vrt'). """
        cfg = _TYPE_CONFIG[self.geoType]
        if vrtFile is None:
            vrtFile = tiffFile + '.vrt'
        srcFiles = [_tiffName(tiffFile, s) for s in cfg['suffixes']]
        # gdal.BuildVRT() already writes a correct relativeToVRT="1" source
        # path when vrtFile and srcFiles share a directory (as they do here,
        # both derived from the same tiffFile base) -- no post-processing needed.
        vrt = gdal.BuildVRT(vrtFile, srcFiles, separate=True)
        for i, comp in enumerate(cfg['components'], start=1):
            vrt.GetRasterBand(i).SetMetadataItem('Description', comp)
        vrt.FlushCache()
        vrt = None

    def getDomain(self, epsg):
        if epsg is None or epsg == 3413:
            domain = 'greenland'
        elif epsg == 3031:
            domain = 'antarctica'
        else:
            myerror('Unexpected epsg code: ' + str(epsg))
        return domain

    def getGeoFile(self, fileName, domain, vxMod=None, geoFile=None,
                   tiff=False, wkt=None, multibandFile=None):
        ''' determine the geodat file name and load geodat info '''
        cfg = _TYPE_CONFIG[self.geoType]
        if not tiff:
            if geoFile is None:
                geoFile = fileName + cfg['geodatSuffix']
        else:
            if geoFile is None:
                if multibandFile is not None:
                    geoFile = multibandFile
                elif vxMod is not None and self.geoType in ('velocity', 'error'):
                    geoFile = fileName + vxMod
                else:
                    geoFile = _tiffName(fileName, cfg['suffixes'][0])
        self.geo = geodat(verbose=self.verbose, domain=domain, wkt=wkt)
        if not tiff:
            self.geo.readGeodat(geoFile)
        else:
            self.geo.readGeodatFromTiff(geoFile)
        return geoFile

    def dataFileNames(self, fileName, tiff=None, vxMod=None):
        ''' compute the file names to be read for all components '''
        cfg = _TYPE_CONFIG[self.geoType]
        if not tiff:
            return [fileName + s for s in cfg['suffixes']]
        # legacy vxMod override for velocity/error tiff paths
        if vxMod is not None and self.geoType in ('velocity', 'error'):
            comps = ['x', 'y'] if self.geoType == 'velocity' else ['x', 'y']
            new = []
            for i, component in enumerate(comps):
                vMod = vxMod.replace('x', component)
                if self.geoType == 'error':
                    vMod = vMod.replace('v', 'e')
                new.append(vMod)
            return [fileName + s for s in new]
        return [_tiffName(fileName, s) for s in cfg['suffixes']]

    def readFiles(self, fileNames, dType, tiff=False, bandMap=None):
        cfg = _TYPE_CONFIG[self.geoType]
        minValue = -2.e9
        sx, sy = self.geo.sizeInPixels()
        for comp, fileName in zip(cfg['components'], fileNames):
            if not tiff:
                myArray = readImage(fileName, sx, sy, dType)
            elif bandMap is not None:
                myArray = self.readMyTiff(fileName, band=bandMap[comp])
            else:
                myArray = self.readMyTiff(fileName)
            if np.sum(np.isnan(myArray)) > 0:
                myArray[np.isnan(myArray)] = np.nan
            elif isinstance(myArray[0, 0], np.floating):
                myArray[myArray <= minValue] = np.nan
            setattr(self, comp, myArray)
        if cfg['magAttr']:
            comps = [getattr(self, c).astype(float) for c in cfg['components']]
            setattr(self, cfg['magAttr'],
                    np.sqrt(sum(c**2 for c in comps)))

    def readData(self, fileName, geoType=None, geoFile=None, dType='>f4',
                 tiff=False, epsg=None, vxMod=None, wktFile=None):
        """ Read geo image data.
        fileName: file basename (suffixes appended automatically)
        geoType: 'velocity', 'velocityRA', 'error', or 'scalar'
        tiff: if True, read GeoTIFF and extract geodat from it.
        If tiff and the usual per-component files (plain suffix.tif, or
        vxMod-substituted names) aren't found on disk, falls back to
        treating fileName itself as a single multi-band raster with one
        band per component, matched by GDAL band Description (e.g. a
        velocity map distributed as one VRT with 'vx'/'vy' bands rather
        than separate .vx.tif/.vy.tif files) -- see _multibandComponents().
        Backwards compatible: this fallback only triggers when the
        existing per-component files are missing, so any setup that
        already works is unaffected. """
        if geoType is not None:
            self.setGeoType(geoType)
        wkt = self.getWKT_PROJ(epsg, wktFile)
        bandMap = None
        if tiff and geoFile is None:
            candidateNames = self.dataFileNames(fileName, tiff=tiff, vxMod=vxMod)
            if not all(os.path.exists(f) for f in candidateNames):
                bandMap = self._multibandComponents(fileName)
        multibandFile = fileName if bandMap is not None else None
        geoFile = self.getGeoFile(fileName, self.getDomain(epsg), wkt=wkt,
                                  geoFile=geoFile, tiff=tiff, vxMod=vxMod,
                                  multibandFile=multibandFile)
        self.xyCoordinates()
        if bandMap is not None:
            fileNames = [fileName] * len(bandMap)
        else:
            fileNames = self.dataFileNames(fileName, tiff=tiff, vxMod=vxMod)
        self.readFiles(fileNames, dType, tiff=tiff, bandMap=bandMap)

    def readDataWindow(self, fileName, xMinKm, xMaxKm, yMinKm, yMaxKm,
                       geoType=None, tiff=True, epsg=None, wktFile=None,
                       padPixels=2):
        """ Read only the part of a per-component GeoTIFF set (fileName +
        suffix + '.tif', as written by writeMyTiff) that covers the km box
        [xMinKm, xMaxKm] x [yMinKm, yMaxKm] (pixel-centre PS coordinates),
        padded by padPixels. Sets self.geo to the WINDOW's geodat and fills the
        components exactly as readData() would for a full read, so interpGeo /
        setupInterp work unchanged on the window. Returns self, or None when
        the box does not overlap the raster or a component file is missing.
        Used by autocleanNISAR against a whole-track velocityStats composite
        that would be far too large to read in full. """
        if geoType is not None:
            self.setGeoType(geoType)
        if not tiff:
            myerror('readDataWindow: only tiff=True is supported')
        wkt = self.getWKT_PROJ(epsg, wktFile)
        cfg = _TYPE_CONFIG[self.geoType]
        fileNames = self.dataFileNames(fileName, tiff=True)
        if not all(os.path.exists(f) for f in fileNames):
            return None
        domain = self.getDomain(epsg)
        full = geodat(verbose=False, domain=domain, wkt=wkt)
        full.readGeodatFromTiff(fileNames[0])
        dxKm, dyKm = full.pixSizeInKm()
        c0 = int(np.floor((xMinKm - full.x0) / dxKm)) - padPixels
        c1 = int(np.ceil((xMaxKm - full.x0) / dxKm)) + padPixels + 1
        r0 = int(np.floor((yMinKm - full.y0) / dyKm)) - padPixels
        r1 = int(np.ceil((yMaxKm - full.y0) / dyKm)) + padPixels + 1
        c0, c1 = max(c0, 0), min(c1, full.xs)
        r0, r1 = max(r0, 0), min(r1, full.ys)
        if c1 <= c0 or r1 <= r0:
            return None
        self.geo = geodat(x0=full.x0 + c0 * dxKm, y0=full.y0 + r0 * dyKm,
                          xs=c1 - c0, ys=r1 - r0, dx=full.dx, dy=full.dy,
                          domain=domain, verbose=False, wkt=wkt)
        self.xyCoordinates()
        minValue = -2.e9
        for comp, f in zip(cfg['components'], fileNames):
            ds = gdal.Open(f)
            # tif rows run top-down; geoimage arrays are bottom-up
            arr = ds.GetRasterBand(1).ReadAsArray(c0, full.ys - r1, c1 - c0, r1 - r0)
            ds = None
            arr = np.flipud(arr)
            if np.sum(np.isnan(arr)) == 0 and isinstance(arr[0, 0], np.floating):
                arr[arr <= minValue] = np.nan
            setattr(self, comp, arr)
        if cfg['magAttr']:
            comps = [getattr(self, c).astype(float) for c in cfg['components']]
            setattr(self, cfg['magAttr'], np.sqrt(sum(c ** 2 for c in comps)))
        return self

    def writeData(self, fileName, geoType=None, geoFile=None, dType='>f4'):
        """ Write binary flat files + geodat sidecars for all components. """
        if geoType is not None:
            self.setGeoType(geoType)

        def writeMyImage(myGeo, NaNVal, x, outName, dType):
            if dType not in ['u1']:
                x[np.isnan(x)] = NaNVal
            writeImage(outName, x, dType)
            myGeo.writeGeodat(f'{outName}.geodat')

        cfg = _TYPE_CONFIG[self.geoType]
        for comp, suffix in zip(cfg['components'], cfg['suffixes']):
            outName = fileName + suffix
            writeMyImage(self.geo, -2.0e9, getattr(self, comp), outName, dType)
