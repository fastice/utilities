# lsdat.py
import numpy as np
import pyproj
from datetime import datetime, timedelta
from utilities import myerror

class lsdat :

    """ lsdat object - contains information from a landsat match dat file"""

    def __init__(self,x0=None,y0=None,xs=None,ys=None,dx=None,dy=None,domain=None) :
        """ Initialize lsdat(x0=None,y0=None,xs=None,ys=None,dx=None,dy=None,domain=['greenland'])
        domain = greenland or antarctica """
        # set everything to zero as default
        self.x0,self.y0,self.dx,self.dy=0.0,0.0,0.0,0.0
        self.successRate,self.culledRate=-1.0,-1.0
        self.slowFlag=False
        self.xs,self.ys=0,0
        self.stepX,self.stepY=0,0
        self.dxImage,self.dyImage=0.,0.
        self.fileEarly,self.fileLate='',''
        self.sigmaX,self.sigmaY=0.0,0.0
        self.domain=3413
        # in most cases all or no args would be passed.
        if x0 != None :
            self.x0=x0
        if y0 != None :            
            self.y0=y0
        if xs != None :            
            self.xs=int(xs)
        if ys != None :                        
            self.ys=int(ys)
        if dx != None :            
            self.dx=dx
        if dx != None :            
            self.dy=dy
        if domain != None :
            self.domain=domain.lower()
            domains=[3413,3031] 
            if  not self.domain in domains :
                print('\n--- setup lsdat invalid domain ',self.domain,'not in ',domains)
                exit()
        #
        # setup conversions
        #
        self.llproj=pyproj.Proj("+init=EPSG:4326")
        if self.domain == 3413:
            self.xyproj=pyproj.Proj("+init=EPSG:3413")
        elif self.domain == 3031 :
            self.xyproj=pyproj.Proj("+init=EPSG:3031")
#        print('setting up projections',self.domain)
        
    def lltoxym(self,lat,lon,scale=None) :
        """ lat,lon is an nparray of xy points in units of km, output is x,y in meters """
        if scale == None :
            scale=1
        try :
            sh=lat.shape
        except :
            lat=np.array(lat)
            lon=np.array(lon)
            sh=lat.shape
            
        lat1=lat.flatten()
        lon1=lon.flatten()
        x,y=pyproj.transform(self.llproj,self.xyproj,lon1,lat1)
        x=scale * x.reshape(sh)
        y=scale * y.reshape(sh)
        return x,y

    def lltoxykm(self,lat,lon,scale=None) :
        """ lat,lon is an nparray of xy points in units of km, output is x,y in KILOmeters """
        return self.lltoxym(lat,lon,scale=0.001)
    
    def xymtoll(self,x,y,scale=None) :
        """ x,y  is an nparray of points in PS meters, output lat,lon in similar array """
        if scale == None :
            scale=1
        try :
            sh=x.shape
        except :
            x=np.array(x)
            y=np.array(y)
            sh=x.shape            
        xx=x.flatten()
        yy=y.flatten()                        
        lon,lat=pyproj.transform(self.xyproj,self.llproj,xx*scale,yy*scale)
        #[ self.xytoll.TransformPoint(xx[k]*scale,yy[k]*scale) for k in range(0,len(xx))]
        lat=lat.reshape(sh)
        lon=lon.reshape(sh)
        return lat,lon

    def xykmtoll(self,x,y) :
        """ x,y  is an nparray of points in PS KILOmeters, output lat,lon in similar array """
        return self.xymtoll(x,y,scale=1000)
     
    def sizeInPixels(self):
        return self.xs,self.ys

    def sizeInKm(self):
        return self.xs*self.dx*0.001,self.ys*self.dy*0.001

    def sizeInM(self):
        return self.xs*self.dx,self.ys*self.dy

    def originInKm(self):
        return self.x0,self.y0

    def boundsInKm(self):
        return self.x0,self.y0,self.x0+(self.xs-1)*self.dx*0.001,self.y0+(self.ys-1)*self.dy*0.001

    def boundsInM(self):
        return self.x0*1000,self.y0*1000,(self.x0+(self.xs-1)*self.dx*0.001)*1000.,(self.y0+(self.ys-1)*self.dy*0.001)*1000.
    
    def originInM(self):
        return self.x0*1000,self.y0*1000

    def pixSizeInM(self):
        return self.dx,self.dy
    
    def pixSizeInKm(self):
        return self.dx*0.001,self.dy*0.001

    def xd(self):
        """ return ls information as an 'xd' np matrix """
        xd=np.array([[self.xs,self.ys],[self.dx,self.dy],[self.x0,self.y0]])
        return xd

    def readLSdat(self,LSDatFile,printData=False) :
        try :
            fp=open(LSDatFile,'r')
        # crude check sum
            for line in fp :
                if 'fileEarly' in line :
                    self.fileEarly=line.split('=')[-1].strip()
                elif 'fileLate' in line :
                    self.fileLate=line.split('=')[-1].strip()                
                elif 'x0' in line :
                    self.x0=float(line.split('=')[-1])/1000.
                elif 'y0' in line :
                    self.y0=float(line.split('=')[-1])/1000.
                elif 'dx' in line :
                    self.dxImage=float(line.split('=')[-1])
                elif 'dy' in line :
                    self.dyImage=float(line.split('=')[-1])
                elif 'stepX' in line :
                    self.stepX=int(line.split('=')[-1])
                elif 'step' in line :
                    self.stepY=int(line.split('=')[-1])
                elif 'nx' in line :
                    self.xs=int(line.split('=')[-1])
                elif 'ny' in line :
                    self.ys=int(line.split('=')[-1])
                elif 'slowFlag' in line :
                    self.slowFlag=False
                    if int(line.split('=')[-1]) > 0 :
                        self.slowFlag=True
                elif 'EPSG' in line :
                    self.domain=int(line.split('=')[-1])
                elif 'earlyImageJD' in line :
                    self.JD1=float(line.split('=')[-1])
                elif 'lateImageJD' in line :
                    self.JD2=float(line.split('=')[-1])
                elif 'Success_rate_for_attempted_matches' in line :
                    self.successRate=float(line.split('=')[-1])
                elif 'Culled_rate_for_attempted_matches' in line :
                    self.culledRate=float(line.split('=')[-1])
                elif 'Mean_sigmaX' in line :
                    self.sigmaX=float(line.split('=')[-1])
                elif 'Mean_sigmaY' in line :
                    self.sigmaY=float(line.split('=')[-1])
            fp.close()
            self.dx=self.dxImage*self.stepX
            self.dy=self.dyImage*self.stepY
            self.date1=datetime.strptime('2000:01:01',"%Y:%m:%d") + timedelta(days=self.JD1-2451544.5)
            self.date2=datetime.strptime('2000:01:01',"%Y:%m:%d")+ timedelta(days=self.JD2-2451544.5)                        
            #
            if printData :
                print(self.fileEarly)
                print(self.fileLate)
                print(self.x0,self.y0)
                print(self.dx,self.dy)
                print(self.stepX,self.stepY)
                print(self.xs,self.ys)
                print(self.domain)            
                print(self.date1,self.date2)
                print(self.successRate,self.culledRate)
                print(self.sigmaX,self.sigmaY)

        except :
                myerror("Problem reading lsdat file "+LSDatFile)
        
    
   
