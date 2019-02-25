# lsfit.py
import numpy as np
import pyproj
from osgeo import gdal
from datetime import datetime, timedelta
from utilities import myerror

class lsfit :

    """ lsdat object - contains information from a landsat match dat file"""

    def __init__(self,fitFile=None) :
        # set everything to zero as default
        self.nTiesRejected=np.array([0.0,0.0,0.0])
        self.nTiesUsed=0.0
        self.xFit,self.yFit=np.array([0.0,0.0,0.0]),np.array([0.0,0.0,0.0])
        self.xResidual,self.yResidual=np.array([0.0,0.0,0.0]),np.array([0.0,0.0,0.0])        
        self.Cxx,self.Cyy=np.zeros((3,3)),np.zeros((3,3))
        self.fitFile=None
        if fitFile !=None :
            self.fitFile=fitFile
            self.readLSFit()

#        print('setting up projections',self.domain)
        
    def readLSFit(self,fitFile=None,printData=False) :
        if fitFile != None :
            self.fitFile=fitFile

        if self.fitFile == None :
            myerror('readLSFit: no file specified')

        try :
            fp=open(self.fitFile,'r')
        # crude check sum
            for line in fp :
                if 'N_ties_rejected' in line :
                    self.nTiesRejected =np.array( [int(x) for x in  line.split('=')[-1].strip().split()])
                elif 'nTiesUsed' in line :
                    self.nTiesUsed=int(line.split('=')[-1].strip())
                elif 'Xfit' in line :
                    self.xFit =np.array( [float(x) for x in  line.split('=')[-1].strip().split()])
                elif 'Yfit' in line :
                    self.yFit[:] =np.array( [float(x) for x in  line.split('=')[-1].strip().split()])
                elif 'X_residual' in line :
                    self.xResidual =np.array([float(x) for x in  line.split('=')[-1].strip().split()])
                elif 'Y_residual' in line :
                    self.yResidual[:] =np.array( [float(x) for x in  line.split('=')[-1].strip().split()])                    
                elif 'CX_11_12_13' in line :
                    self.Cxx[0,: ]=np.array([float(x) for x in  line.split('=')[-1].strip().split()])
                elif 'CX_21_22_23' in line :
                    self.Cxx[1,: ]=np.array([float(x) for x in  line.split('=')[-1].strip().split()])
                elif 'CX_31_32_33' in line :
                    self.Cxx[2,:] =np.array([float(x) for x in  line.split('=')[-1].strip().split()])
                elif 'CY_11_12_13' in line :
                    self.Cyy[0,: ]=np.array([float(x) for x in  line.split('=')[-1].strip().split()])
                elif 'CY_21_22_23' in line :
                    self.Cyy[1,: ]=np.array([float(x) for x in  line.split('=')[-1].strip().split()])
                elif 'CY_31_32_33' in line :
                    self.Cyy[2,:] =np.array([float(x) for x in  line.split('=')[-1].strip().split()])                                                                                                  
            #
            fp.close()
            
            if printData :
                print(self.fitFile)
                print(self.nTiesRejected)
                print(self.nTiesUsed)
                print(self.xFit)
                print(self.yFit)
                print(self.xResidual)
                print(self.yResidual)
                print(self.Cxx)
                print(self.Cyy) 

        except :
                myerror("Problem reading lsdat file "+self.fitFile)
        
    
   
