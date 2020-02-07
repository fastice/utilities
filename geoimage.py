# geoimage.py
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from utilities.readImage import readImage
from utilities.writeImage import writeImage
from utilities.myerror import myerror
from utilities import geodat
import os
#import urllib.request
from osgeo.gdalconst import *
from osgeo import gdal, gdal_array, osr
from datetime import datetime

#----------------------------------------------------------------------------------------------------------------------------------------------------
# class defintion for an image object, which covers PS data described by a geodat
#----------------------------------------------------------------------------------------------------------------------------------------------------
class geoimage :
    
    """ \ngeoimage - object for scalar or velocity PS data + geodat data """
    def __init__(self,x=None,vx=None,vy=None,v=None,ex=None,ey=None,e=None,geoType=None,verbose=True) :
        
        self.x=[]
        self.vx,self.vy,self.v,self.ex,self.ey,self.e=[],[],[],[],[],[]
        #self.xs=self.ys=0
        self.geo=[]
        self.xx,self.yy=[],[]
        self.xGrid,self.yGrid=[],[]
        self.geoType=None
        self.velDate=None
        self.fileName=None
        self.fileRoot=None
        if x != None :
            self.x=x
        if vx != None :
            self.vx=vx
        if vy != None :
            self.vy=vy
        if v != None :
            self.v=v
        if ex != None :
            self.ex=ex
        if ey != None :
            self.ey=ey
        if e != None :
            self.e=e   
        
        #if xs != None :
        #    self.xs=xs
        #if ys != None :
        #    self.ys=ys
        self.verbose=verbose
        if geoType != None :
            self.geoType=geoType
            if self.verbose :
                print('Type ',self.geoType)

    #--------------------------------------------------------------------------
    # Set and errorcheck geotype
    #--------------------------------------------------------------------------
    def setGeoType(self,geoType) :
        """ setGeoType(geoType) set type to velocity or scalar"""
        types=['scalar','velocity','error']
        if not geoType in types :
            print('\n\n\tgeoImage setType, invalid type : ',type,'\n\n')
            exit()
        self.geoType=geoType

    #--------------------------------------------------------------------------
    # def setup xy limits
    #--------------------------------------------------------------------------
    def xyCoordinates(self) :
        """ xyCoordinates - setup xy coordinates in km """
        #
        sx,sy=self.geo.sizeInPixels()
        x0,y0=self.geo.originInKm()
        dx,dy=self.geo.pixSizeInKm()
        # remember arange will not generate value for sx*dx (its doing sx-1)
        self.xx=np.arange(x0,x0+sx*dx,dx)
        self.yy=np.arange(y0,y0+sy*dy,dy)
        # force the right length
        self.xx,self.yy=self.xx[0:sx],self.yy[0:sy]   
        

    #--------------------------------------------------------------------------
    # Compute matrix of xy grid points
    #--------------------------------------------------------------------------    
    def xyGrid(self) :
        #
        # if one done grid points not computed, then compute
        if len(self.xx) == 0 :
            self.xyCoordinates()
        sx,sy=self.geo.sizeInPixels()
        #
        # setup array
        self.xGrid,self.yGrid=np.zeros((sy,sx)),np.zeros((sy,sx))
        for i in range(0,sy) :
            self.xGrid[i,:]=self.xx
        for i in range(0,sx) :
            self.yGrid[:,i]=self.yy  
            
    def parseMyMeta(self,metaFile) :
        print(metaFile)
        fp=open(metaFile)
        dates=[]
        for line in fp :
            if 'MM:DD:YYYY' in line :
                tmp=line.split('=')[-1].strip()
                dates.append(datetime.strptime(tmp,"%b:%d:%Y"))
                if len(dates) == 2 :
                    break
        if len(dates) != 2 :
            return None
        fp.close()
        return dates[0]+(dates[1]-dates[0])*0.5
  
    def parseVelCentralDate(self) :

        if(self.fileName != None) :
            metaFile=self.fileName + '.meta'
            if not os.path.exists(metaFile) :
                return None
            return self.parseMyMeta(metaFile)
        return None
        

    #--------------------------------------------------------------------------
    #  setup interpolation functions
    #--------------------------------------------------------------------------
    def setupInterp(self) :
        """ set up interpolation for scalar (xInterp) or velocity/eror (vxInterp,vyInterp,vInterp)  """
    #
        if len(self.xx) < 0 :
            myerror('\n\nsetupInterp: x,y limits not set\n\n')
        #
        # setup interp - flip xy for row colum
        #
        xy=(self.yy,self.xx)
        
        if self.geoType == 'scalar' :
            self.xInterp=RegularGridInterpolator(xy,self.x,method='linear')
        if self.geoType == 'velocity' :
            self.vxInterp=RegularGridInterpolator(xy,self.vx,method='linear')
            self.vyInterp=RegularGridInterpolator(xy,self.vy,method='linear')
            self.vInterp=RegularGridInterpolator(xy,self.v,method='linear')
        if self.geoType == 'error' :
            self.exInterp=RegularGridInterpolator(xy,self.ex,method='linear')
            self.eyInterp=RegularGridInterpolator(xy,self.ey,method='linear')
            self.eInterp=RegularGridInterpolator(xy,self.e,method='linear')            

    #--------------------------------------------------------------------------
    # interpolate geo image
    #--------------------------------------------------------------------------
    def interpGeo(self,x,y) :
        """ interpolate velocity or x at points x and y, which are in km (note x,y is c-r even though data r-c)"""
        # save the original shape
        shapeSave=x.shape
        # flatten and do bounds check
        x1=x.flatten()
        y1=y.flatten()
        xgood=np.logical_and(x1 >= self.xx[0],x1 <= self.xx[-1])
        ygood=np.logical_and(y1 >= self.yy[0],y1 <= self.yy[-1])
        igood=np.logical_and(xgood,ygood)
        #
        # save inbound locations
        xy=np.array([y1[igood],x1[igood]]).transpose()
        #
        if self.geoType == 'scalar' :
            result=np.zeros(x1.transpose().shape)
            result[:]=np.NaN
            result[igood]=self.xInterp(xy)
            result=np.reshape(result,shapeSave)
            return result
        elif self.geoType == 'velocity' :
            # empty arrays the right size
            vxr=np.zeros(x1.transpose().shape)
            vxr[:]=np.NaN
            vyr=vxr.copy()
            vr=vxr.copy()
            # interpolate in bounds
            vxr[igood]=self.vxInterp(xy)
            vyr[igood]=self.vyInterp(xy)
            vr[igood]=self.vInterp(xy)
            vxr=np.reshape(vxr,shapeSave)
            vyr=np.reshape(vyr,shapeSave)
            vr=np.reshape(vr,shapeSave)            
            return vxr,vyr,vr
        elif self.geoType == 'error' :
            # empty arrays the right size
            exr=np.zeros(x1.transpose().shape)
            exr[:]=np.NaN
            eyr=exr.copy()
            er=exr.copy()
            # interpolate in bounds
            exr[igood]=self.exInterp(xy)
            eyr[igood]=self.eyInterp(xy)
            er[igood]=self.eInterp(xy)
            exr=np.reshape(exr,shapeSave)
            eyr=np.reshape(eyr,shapeSave)
            er=np.reshape(er,shapeSave)            
            return exr,eyr,er        
    #
    #
    #
    def readGeodat(self,geoFile) :
        if self.verbose :
            print(geoFile)
        if self.geo == [] :
            self.geo=geodat(verbose=self.verbose)            
        if os.path.exists(geoFile) :
            self.geo.readGeodat(geoFile)            
        else :
            myerror('Missing geodat file '+geoFile)

            
    def readMyTiff(self,tiffFile) :
        """ read a tiff file and return the array """
        try :
            gdal.AllRegister()
            ds = gdal.Open(tiffFile)
            band=ds.GetRasterBand(1)
            arr=band.ReadAsArray()
            arr=np.flipud(arr)
            ds=None
        except :
            myerror("geoimage.readMyTiff: error reading tiff file "+tiffFile)
        return arr
    
    def getWKT_PROJ(self,epsg_code) :
        #url="http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg_code)
        #response=urllib.request.urlopen(url)
        #remove_spaces =str(response.read()).replace(" ","")
        #output=remove_spaces.replace("\\n","").replace("b\'","").replace("\'","")
        # hard code common ones to avoid problems where defs can't be read by too many jobs at once
        wktDict={
            3413 :'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984", \
            SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]], \
            PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],\
            AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],\
            PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],\
            PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]',
            3013 : 'PROJCS["WGS 84 / Antarctic Polar Stereographic",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,\
            AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],\
            UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],\
            PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],\
            PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3031"]]'  }   
        try :
            wkt=wktDict[epsg_code] 
        except  :    
            sr=osr.SpatialReference()
            sr.ImportFromEPSG(epsg_code)
            wkt=sr.ExportToWkt()
        return wkt
    
    def imageSize(self) :
        typeDict={'scalar' : self.x, 'velocity' : self.vx,'error': self.ex}
        ny,nx=typeDict[self.geoType].shape
        return nx,ny
    
    def computePixEdgeCornersXYM(self) :
        nx,ny=self.imageSize()
        x0,y0=self.geo.originInM()
        dx,dy=self.geo.pixSizeInM()
        xll,yll=x0 - dx/2,y0 - dx/2
        xur,yur= xll + nx * dx, yll + ny * dy
        xul,yul=xll,yur
        xlr,ylr=xur,yll
        corners={'ll' : {'x' :xll, 'y' : yll},'lr' : {'x' :xlr, 'y' : ylr}, 'ur' : {'x' :xur, 'y' : yur},'ul' : {'x' :xul, 'y' : yul}}
        return corners
    
    def computePixEdgeCornersLL(self) :
        corners=self.computePixEdgeCornersXYM()
        llcorners={}
        for myKey in corners.keys():
            lat,lon=self.geo.xymtoll(np.array([corners[myKey]['x']]),np.array([corners[myKey]['y']]))
            llcorners[myKey]={'lat' : lat[0],'lon' : lon[0]}
        return llcorners
        
    #--------------------------------------------------------------------------
    # write My Tiff 
    #--------------------------------------------------------------------------
    def writeMyTiff(self,tiffFile,epsg=None,noDataDefault=None,predictor=1,noV=False,overviews=None) :
        """ write a geotiff file  - NEEDS MODIFICATION FOR EPSG AND VX,EX
            Note : tiffFile should not have a ".tif" extension - one will be added.        
        """
        # define various set up stuff
        suffixDict={'scalar' : [''], 'velocity' : ['.vx','.vy','.v'], 'error' : ['.ex','.ey']}
        typeDict={'scalar' : self.x, 'velocity' : self.vx,'error': self.ex}
        predictor=[int(predictor),1][predictor > 3 or predictor < 1]
        epsg=[epsg,3413][epsg == None]
        try :
            suffixes=suffixDict[self.geoType]
            gdalType=gdal_array.NumericTypeCodeToGDALTypeCode(typeDict[self.geoType].dtype) 
        except :
            myerror('writeMyTiff: invalide geoType '+self.geoType)
        #     
        try :      
            # Loop through different components, also write vmag for velocity
            for suffix in suffixes :
                # skip .v if requested
                if noV and suffix == '.v' :
                    continue
                # write the geotiff
                self.writeCloudOptGeo(tiffFile,suffix,epsg,gdalType,overviews=overviews,predictor=predictor,noDataDefault=noDataDefault)
        except :
            myerror("geoimage.writeMyTiff - : error writing tiff file "+tiffFile)
            
            
    def writeCloudOptGeo(self,tiffFile,suffix,epsg,gdalType,overviews=None,predictor=1,noDataDefault=None) :
        ''' write a cloudoptimized geotiff with overviews'''
        # no data info
        noDataDict={'.vx' : -2.0e9 , '.vy' : -2.0e9,'.v' : -1.0, '.ex' : -1.0, '.ey' : -1.0,'' : noDataDefault}
        #
        # use a temp mem driver for CO geo
        driver=gdal.GetDriverByName( "MEM")
        nx,ny=self.imageSize()
        dx,dy=self.geo.pixSizeInM()
        dst_ds=driver.Create('',nx,ny,1,gdalType )
        # set geometry
        tiffCorners=self.computePixEdgeCornersXYM()
        dst_ds.SetGeoTransform((tiffCorners['ul']['x'],dx,0,tiffCorners['ul']['y'],0,-dy))
        #set projection
        wkt=self.getWKT_PROJ(epsg)
        dst_ds.SetProjection(wkt)
        # this seems to force writing full wkt
        #dst_ds.FlushCache()
        # set nodata
        noData=noDataDict[suffix]
        if noData != None :
            if self.geoType=='scalar' :
                tmp=self.x
            else :
                tmp=eval('self'+suffix)
            tmp[np.isnan(tmp)]=noData
            dst_ds.GetRasterBand(1).SetNoDataValue(noData)
        # write data 
        if self.geoType == 'scalar' :
            dst_ds.GetRasterBand(1).WriteArray(np.flipud(self.x))
        else :
            dst_ds.GetRasterBand(1).WriteArray(np.flipud(eval('self'+suffix)))
        #
        if overviews != None :
            dst_ds.BuildOverviews('AVERAGE',overviews)
        # now copy to a geotiff - mem -> geotiff forces correct order for c opt geotiff
        dst_ds.FlushCache()
        driver=gdal.GetDriverByName( "GTiff")
        dst_ds2=driver.CreateCopy(tiffFile+suffix+'.tif',dst_ds,options=['COPY_SRC_OVERVIEWS=YES','COMPRESS=LZW',f'PREDICTOR={predictor}','TILED=YES'] )
        dst_ds=None
        dst_ds2.FlushCache()
        # free memory
        dst_ds,dst_ds2=None,None
        
    def getDomain(self,epsg) :
        if epsg==None or epsg == 3413 :
            domain='greenland'
        elif epsg == 3031 :
            domain='antarctica'
        else :
            myerror('Unexpected epsg code : '+str(epsg))
        return domain
    
    def getGeoFile(self,fileName,domain,vxMod=None,geoFile=None,tiff=False) :
        ''' determine the file name for geodat if not specified and load geodat info'''
        if not tiff :
            suffixes={'scalar' : '.geodat', 'velocity' : '.vx.geodat','error' : '.vx.geodat' } 
        else :
            suffixes={'scalar' : '', 'velocity' : '.vx.tif' ,'error' : '.vx.tif' } 
            if vxMod != None :
                suffixes['error']=vxMod
                suffixes['velocity']=vxMod
        #
        if geoFile == None :
            try :
                geoFile=fileName+suffixes[self.geoType]
            except :
                myerror("geoimage.readData: Invalid type {self.geoType} ")
        #
        self.geo=geodat(verbose=self.verbose,domain=domain)    
        if not tiff :
            self.geo.readGeodat(geoFile)
        else :
            self.geo.readGeodatFromTiff(geoFile)
        return geoFile
    
    
    def dataFileNames(self,fileName,tiff=None,vxMod=None) :
        ''' compute the file names that need to be read'''
        if not tiff :
            suffixes={'scalar' : [''], 'velocity' : ['.vx','.vy'],'error' : ['.ex','.ey'] }[self.geoType]
        else :
            suffixes={'scalar' : [''], 'velocity' : ['.vx.tif','.vy.tif'],'error' : ['.ex.tif','ey.tif'] }[self.geoType] 
            # if vxMod present, use it update the name
            if vxMod != None and self.geoType != 'scalar' :
                for i,component in zip([0,1],['x','y']):
                    vMod=vxMod.replace('x',component)
                    # v -> e for errors
                    if self.geoType == 'error' :
                        vMod.replace('v','e')
                    # replace the standard value with this updated vMod
                    suffixes[i]=suffixes[i].replace(suffixes[i],vMod)
        # now use file root to create file names
        fileNames=[]
        for suffix in suffixes :
            fileNames.append(fileName+suffix)
#        print(fileNames)
        return fileNames
        
    def readFiles(self,fileNames,dType,tiff=False) :
        #  get the values that match the type
        myTypes={'scalar' : ['x'],'velocity' : ['vx','vy'], 'error' : ['ex','ey']}[self.geoType]
        #noData={'scalar' : -2.e9,'velocity' : -2.e9,'error' : -2.e9}[self.geoType]
        minValue={'scalar' : -2.e9,'velocity' : -2.e9,'error' : -2.e9}[self.geoType]
        sx,sy=self.geo.sizeInPixels()
        # loop over arrays and file names to read data
        for myType, fileName in zip(myTypes,fileNames) :
            if not tiff :
                myArray=readImage(fileName,sx,sy,dType)
            else :
                myArray=self.readMyTiff(fileName)
            # handle no data
            myArray[np.isnan(myArray)]=np.nan
            myArray[myArray <=minValue]=np.nan
            exec(f'self.{myType}=myArray')
        #
        # compute mag for velocity and errors
        if self.geoType == 'velocity' :
            self.v=np.sqrt(np.square(self.vx) + np.square(self.vy))
        elif self.geoType == 'error' :
            self.e=np.sqrt(np.square(self.ex) + np.square(self.ey))
            
    #--------------------------------------------------------------------------
    # read geo image data
    #--------------------------------------------------------------------------
    def readData(self,fileName,geoType=None,geoFile=None,dType='>f4',tiff=False,epsg=None,vxMod=None) :
        """ read Data for geo image 
        fileName=filename (or basename if velocity )
        geoType =specify read 'velocity' or 'scalar' data
        geoFile=geodatfile [fileName(.vx).geodat]
        dType=type for scalar ['>f4']   ( 'f4','>f4','>u2','u2','u1','>i2','i2','>u4','u4','>i4','i4')"""    
        #
        # error check type
        if geoType != None :
            self.setGeoType(geoType)
        #
        # read geodat
        domain=self.getDomain(epsg) 
        geoFile=self.getGeoFile(fileName,domain,geoFile=geoFile,tiff=tiff,vxMod=vxMod)
        # compute coordinates for data
        self.xyCoordinates()
        # get size
        sx,sy=self.geo.sizeInPixels()
        #
        # read image (set no data to nan)  
        fileNames=self.dataFileNames(fileName,tiff=tiff,vxMod=vxMod)
        self.readFiles(fileNames,dType,tiff=tiff)
        #exit()
#        if self.geoType == 'scalar' :
#            if tiff :
#                self.x=self.readMyTiff(fileName)
#            else :
#                self.x=readImage(fileName,sx,sy,dType)
#                if 'f' in dType :
#                    missing=self.x <=(-2.0e9+1)
#                    if len(missing) > 0 :
#                        self.x[missing]=np.nan
#        elif self.geoType == 'velocity' :
#            self.fileName=fileName
#            if tiff :
#                self.vx=self.readMyTiff(fileName+'.vx.tif')
#                self.vy=self.readMyTiff(fileName+'.vy.tif')
#            else :
#                self.vx=readImage(fileName+'.vx',sx,sy,'>f4')
#                self.vy=readImage(fileName+'.vy',sx,sy,'>f4')
#            missing=[]
#            self.vx[np.isnan(self.vx)]=-2.0e9
#            missing=self.vx <=(-99998)
#            if len(missing) > 0 :
#                self.vx[missing]=np.nan
#                self.vy[missing]=np.nan
#            self.v=np.sqrt(np.square(self.vx) + np.square(self.vy))
#        elif self.geoType == 'error' :
#            if tiff :
#                self.ex=self.readMyTiff(fileName+'.ex.tif')
#                self.ey=self.readMyTiff(fileName+'.ey.tif')
#            else :
#                self.ex=readImage(fileName+'.ex',sx,sy,'>f4')
#                self.ey=readImage(fileName+'.ey',sx,sy,'>f4')
#            missing=self.ey <=(-99998)
#            if len(missing) > 0 :
#                self.ey[missing]=np.nan
#                self.ex[missing]=np.nan            
#            self.e=np.sqrt(np.square(self.ex) + np.square(self.ey))            


#--------------------------------------------------------------------------
# read geo image data
#--------------------------------------------------------------------------
    def writeData(self,fileName,geoType=None,geoFile=None,dType=None) :
        """ read Data for geo image 
        fileName=filename (or basename if velocity )
        geoType =specify read 'velocity' or 'scalar' data
        geoFile=geodatfile [fileName(.vx).geodat]
        dType=type for scalar ['>f4']   ( 'f4','>f4','>u2','u2','>i2','i2','>u4','u4','>i4','i4')"""
        # default assume byteswapped float
        if dType == None :
            dType='>f4'
        #
        # error check type
        if geoType != None :
            self.setGeoType(geoType)
        #
        # write geodat, with suffix from this dict of possible types
        suffix={'scalar' : '.geodat', 'velocity' : '.vx.geodat', 'error' : '.ex.geodat' }
        if geoFile == None :
            geoFile=fileName+suffix[self.geoType]
        #
        # helper to write image (set no data to nan)
        def writeMyImage(myGeo,NaNVal, x,fileName,dType) :
            x[np.isnan(x)]=-2.e9
            writeImage(fileName,x,dType)
            myGeo.writeGeodat(fileName+'.geodat')
        #   
        # Use above funtion to write image for different types
        if self.geoType == 'scalar' :
            writeMyImage(self.geo,-2.0e9,self.x,fileName,dType)
        elif self.geoType == 'velocity' :
            writeMyImage(self.geo,-2.0e9,self.vx,fileName+'.vx',dType)
            writeMyImage(self.geo,-2.0e9,self.vy,fileName+'.vy',dType)                        
        elif self.geoType == 'error' :
            writeMyImage(self.geo,-2.0e9,self.ex,fileName+'.ex',dType)
            writeMyImage(self.geo,-2.0e9,self.ey,fileName+'.ey',dType)

