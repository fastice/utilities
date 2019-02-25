# geodat.py
import numpy as np
from utilities.readImage import readImage
from utilities.writeImage import writeImage
from utilities import geodatrxa
from utilities import myerror
import os
import scipy.signal 
from subprocess import call
import ntpath

class offsets :

    """ offsets object - use for offset information """

    def __init__(self,fileRoot=None,rangeFile=None,latlon=None,datFile=None,sigmaAFile=None,sigmaRFile=None,
                 matchTypeFile=None,geodatrxaFile=None,maskFile=None,verbose=True,myPath=None) :
        """ \n\nInitialize offsets\n
        defaults to setting up for azimuth.offsets
        \tfileRoot = specify name to read offsets in from fileRoot (azimuth.offsets, or file.da
        \trangeFile = specify rangeFile Name
        \tlatlon = read lat/lon for latlon.lat/.lon (assumes path from fileRoot
        \tdatFile = specify .dat file 
        \tsigmaAFile/SigmaRFile = specify sigma A/R files, otherwise guess
        \tmatchTypeFile = specify matchTypeFile, otherwise guess\n
        Note  - this will not actually read the offsets - this must be done with readOffsets        """
        #
        # set everything to zero as default
        #
        self.nr=self.na=self.r0=self.a0=self.dr=self.da=self.azErr=0
        self.fileRoot,self.latlon,self.sigmaAFile,self.sigmaRFile,self.datFile,self.latFile,self.lonFile,self.matchTypeFile,self.maskFile=[],[],[],[],[],[],[],[],[]
        self.azOff,self.rgOff,self.sigmaR,self.sigmaA,self.lat,self.lon,self.matchType,self.heading,self.mask=[],[],[],[],[],[],[],[],[]
        self.cc,self.ccFile=[],None
        self.geodatrxa,self.geodatrxaFile=[],[]
        self.slpRg,self.slpAz=-1,-1
        self.rCoord,self.aCoord=[],[]
        self.path='.'
        self.verbose=verbose
        #
        # in most cases all or no args would be passed.
        #
        if fileRoot != None :            
            self.fileRoot=fileRoot
        else :
            self.fileRoot='./azimuth.offsets'
        #
        tmp=fileRoot.split('/')
        if myPath != None :
            self.setMyPath(myPath)
        elif len(tmp) > 1 :
            self.setMyPath('/'.join(tmp[0:-1])+'/')
        else :
            self.setMyPath('./')
        #        
        # these may get overridden,later
        # Note some of these may punt and set to [] if not file exists
        # in which case they need to be explicitly set
        # compute file names
        #
        self.offsetFileNames(self.fileRoot,rangeFile=rangeFile)
        if self.verbose : 
            print('azFile = ',self.azimuthFile)
            print('rgFile = ',self.rangeFile)
        # set or guess datFileName
        self.datFileName(datFile=datFile)
        if self.verbose :         
            print('datFile =',self.datFile)
        # set or guess lat/lon file names
        self.latlonName(latlon=latlon)
        # set or guess sigmaA/R file names
        self.sigmaFileName(sigmaAFile=sigmaAFile,sigmaRFile=sigmaRFile)
        # set or guess 
        self.matchTypeFileName(matchTypeFile=matchTypeFile)
        #
        self.maskFileName(maskFile=maskFile)        
        #
        if geodatrxaFile != None :
            self.geodatrxaFile=self.offFilePath(geodatrxaFile)
            if self.verbose :             
                print(self.geodatrxaFile)
            self.geodatrxa = geodatrxa(file=self.geodatrxaFile,echo=False)
            self.slpRg,self.slpAz=self.geodatrxa.singleLookResolution()
        if self.verbose :
            print('latlon = ',self.latFile,self.lonFile)
            print('sigma files = ',self.sigmaAFile,self.sigmaRFile)
            print('matchType file = ',self.matchTypeFile)
            print('mask file = ',self.maskFile)
            print('path = ',self.path)
            
            
    def checkOffsetFiles(self,fileRoot=None,myLog=None) :
        ''' checks that offsets both exist and are the correct size'''
        if myLog  != None :
            myLog.logEntry('checkOffsetFiles')    
        if fileRoot != None :
            self.fileRoot=fileRoot
            self.datFileName()
        if self.datFile == None or self.fileRoot==None or self.azimuthFile==None :
            myerror('Check offsets called without setting up file names',myLogger=myLog)
        #
        if not os.path.exists(self.datFile) :
            myerror('Missing offsets dat file {0:s}'.format(self.datFile),myLogger=myLog)
        try :
            self.readOffsetsDat() 
        except :
            myerror('Error reading offsets dat file {0:s}'.format(self.datFile),myLog=myLog)
        #
        mySize=self.nr*self.na*4
        myFiles=[self.azimuthFile,self.rangeFile]
        for myFile in myFiles :
            # check exists
            if not os.path.exists(myFile) :
                myerror('Missing offset file {0:s}'.format(myFile),myLogger=myLog)
            # check size
            else :
                statinfo=os.stat(myFile)
                fileSize=statinfo.st_size
                if fileSize != mySize :
                    myerror('Offset file {0:s} should have {1:d} bytes but only has {2:d} bytes'.format(myFile,mySize,fileSize),myLogger=myLog)
        
        if myLog  != None :
           myLog.logReturn('checkOffsetFiles')
           
       
            
    def setMyPath(self,pathName) :
        if pathName == None :
            myerror('sop')
        self.path = pathName
    
    def offFilePath(self,fileName) :
        # check path as not already been appended
        #print(self.path,fileName)
        #print(os.path.join(self.path,ntpath.basename(fileName)))
        #myerror('stop')
        if self.path not in fileName :
            return os.path.join(self.path,ntpath.basename(fileName))
        else :
            return fileName
     
    def areValid(self) :
        return np.logical_and(self.rgOff > -1.99e9, self.azOff > -1.99e9)
        
    def notValid(self) :
        return np.logical_or(self.rgOff < -1.99e9, self.azOff < -1.99e9)    
    # generic minmax  
    def minmax(self,x)   :
        minA,maxA=None,None
        if len(x) > 0 :
            good = self.areValid()
            minA=np.amin(x[good])
            maxA=np.amax(x[good])   
        return minA,maxA
        
        
    def minmaxA(self) :
        return self.minmax(self.azOff)
        
    def minmaxR(self) :
        return self.minmax(self.rgOff)
        
    #----------------------------------------------------------------------------------------
    # read a geodat file
    #----------------------------------------------------------------------------------------
#    def readGeodatxra(self,geodatxraFile) :
#        try :
#            self.geodatxra = geodatrxa(file=self.offFilePath(geodatxraFile))
#        except Exception :
#            myerror('Error loading geodatxra file in offsets')
#        self.slpRg,self.slpAz=self.geodat.singleLookResolution()
       #----------------------------------------------------------------------------------------
    # reads cross corr data (cc) name
    #----------------------------------------------------------------------------------------
    def readCC(self,ccFile=None,datFile=None)  :
        # makes sure there is a name
        if self.ccFile == None or ccFile != None:
            self.ccFileName(ccFile=ccFile)
        # get sizes if needed
        if self.nr < 1 or self.na < 1 :
            self.readOffsetsDat(datFile=datFile)
        self.cc=readImage(self.ccFile,self.nr,self.na,'>f4')
        
    def ccFileName(self,ccFile=None):
        if ccFile !=None :
            self.ccFile=ccFile
            return
        if '.da' in self.fileRoot :
            self.ccFile=self.fileRoot.replace('.da','.cc')
            return
        myerror('do not know how to make cc file name for {0:s}'.format(self.fileRoot))
            
    #----------------------------------------------------------------------------------------
    # process lat/lon name
    #----------------------------------------------------------------------------------------
    def latlonName(self,latlon=None,myPath=None) :
        """" process lat/lon name - if len(latlon)==0, use offsets.lat,offsets.lon """
        if myPath != None :
            self.setMyPath(myPath)
            #self.latFile = self.offFilePath(self.latFile)
            #self.lonFile = self.offFilePath(self.lonFile)
        if  latlon == None:
            self.latFile=self.offFilePath('offsets.lat')
            self.lonFile=self.offFilePath('offsets.lon')
        else :
            self.latFile=self.offFilePath(latlon+'.lat')
            self.lonFile=self.offFilePath(latlon+'.lon')
        # null files names if they don't exist
        if not os.path.exists(self.latFile) or not os.path.exists(self.lonFile) :
            if self.verbose :
                print('Warning: one or of these files do no exist ',self.latFile,self.lonFile)
            self.latFile=[]
            self.lonFile=[]
  
 
            
    #----------------------------------------------------------------------------------------
    # reads lat/lon name
    #----------------------------------------------------------------------------------------
    def readLatlon(self,latlon=None,datFile=None) :
        """ reads lat/lon file - using set filenames, or specifiy with latlon=basefilename"""
        # read datFile if needed. 
        if self.nr < 1 or self.na < 1 :
            try :
                self.readOffsetsDat(datFile=datFile)
            except :
                myerror('Could not read lat/lon dat file {0:s}'.format(datFile))
        #  
        
        if latlon != None :
            self.latlonName(latlon=latlon)
        #
        if os.path.exists(self.latFile) and  os.path.exists(self.lonFile) :
            if self.verbose :
                print('Reading ',self.latFile, ' and ' ,self.lonFile,' with size ',self.nr,self.na)
            self.lat=readImage(self.latFile,self.nr,self.na,'>f8')
            self.lon=readImage(self.lonFile,self.nr,self.na,'>f8')
        else :
            myerror( 'Offsets tried to read an invalid lat or lon file - {0:s} {1:s}'.format(self.latFile,self.lonFile))
    #----------------------------------------------------------------------------------------
    # return lat/lon dat  - force read if need.
    #----------------------------------------------------------------------------------------
    def getLatLon(self,latlon=None) :
      
        if len(self.lat) <= 0 :
            self.readLatlon(latlon=latlon)
       
        # check again in case it faile 
        if len(self.lat) <= 0 :
            myerror(' Could not read lat/lon - check files exist ')

        # force conversion to 64 bit for coordinate transforms
        return self.lat.astype(float),self.lon.astype(float)
    

    #----------------------------------------------------------------------------------------
    # create matchType file names
    #----------------------------------------------------------------------------------------
    def matchTypeFileName(self,matchTypeFile=None,myPath=None) :
        """ compute matchType file names or set with matchTypeFile=filename """
        #
        # compute defaults first
        if '.da' in self.azimuthFile :
            self.matchTypeFile=self.azimuthFile.replace('.da','.mt')
        # now over ride
        if matchTypeFile != None :
            self.matchTypeFile=matchTypeFile
        # now check path
        # return if couldn't form name
        if len(self.matchType) ==0 :
            return
        # have name, so check, cancel if doesn't exit
        if not os.path.exists(self.matchTypeFile):
            if self.verbose :
                print('*** warning no matchType files exist ****')               
                self.matchTypeFile=[]
                return
        # ammend path if necessary
        if myPath != None :
            self.setMyPath(myPath)
        self.matchTypeFile = self.offFilePath(self.matchTypeFile) 
        return

    #----------------------------------------------------------------------------------------
    # create matchType file names
    #----------------------------------------------------------------------------------------
    def maskFileName(self,maskFile=None,myPath=None) :
        """ compute mask file names or set with matchTypeFile=filename """
        #
        # compute defaults first
        if '.da' in self.azimuthFile :
            self.maskFile=self.azimuthFile.replace('.da','.mask')
        # now over ride
        if maskFile != None :
            self.maskFile=maskFile
        # now check path
        # return if couldn't form name
        if len(self.maskFile) ==0 :
            return
        # have name, so check, cancel if doesn't exit
        #if not os.path.exists(self.maskFile):
        #    if self.verbose :
        #        print('*** warning no mask type files exist ****')
        #    self.maskFile=[]
        # Override path
        # ammend path if necessary
        if myPath != None :
            self.setMyPath(myPath)
        self.maskFile = self.offFilePath(self.maskFile)               
        return    

 #----------------------------------------------------------------------------------------
    # read matchType data
    #----------------------------------------------------------------------------------------
    def readMask(self,maskFile=None) :
        """ reads mask files - assumes filenames set, or specify here with maskFile=filename """
        if maskFile != None :
            self.maskFileName(maskFile=maskFile)
        if len(self.maskFile) == 0:
            myerror('\n\nError : Offsets tried to read an blank mask file - '+self.maskFile)
        if os.path.exists(self.maskFile) :
            if self.verbose :            
                print('Reading ',self.maskFile,'with size ',self.nr,self.na)
            self.mask=readImage(self.maskFile,self.nr,self.na,'u1')
        else :
            myerror('\n\nError : Offsets tried to read an invalid mask file - '+self.maskFile)
        return self.mask
#----------------------------------------------------------------------------------------
    # read matchType data
    #----------------------------------------------------------------------------------------
    def getMask(self,maskFile=None) :
        """ return mask - force read if not already ready """
        if len(self.mask) <= 0 :
            self.readMask(maskFile)
        # check again in case it faile 
        if len(self.mask) <= 0 :
            self.mask=np.zeros((self.na,self.nr),dtype=np.byte)
            myerror(' Could not read mask - check files exist - returning 0s')
        # 
        return self.mask
    
    #----------------------------------------------------------------------------------------
    # read matchType data
    #----------------------------------------------------------------------------------------
    def readMatchType(self,matchTypeFile=None) :
        """ reads matchType files - assumes filenames set, or specify here with matchTypeFile=filename """
        if matchTypeFile != None :
            self.matchTypeFileName(matchTypeFile=matchTypeFile)
        if len(self.matchTypeFile) == 0:
            myerror('\n\nError : Offsets tried to read an blank matchType file - '+self.matchTypeFile)
        if os.path.exists(self.matchTypeFile) :
            if self.verbose :
                print('Reading ',self.matchTypeFile,'with size ',self.nr,self.na)
            self.matchType=readImage(self.matchTypeFile,self.nr,self.na,'u1')
        else :
            myerror('\n\nError : Offsets tried to read an invalid matchType file - '+self.matchTypeFile)
        return self.matchType
 #----------------------------------------------------------------------------------------
    # get matchType data
    #----------------------------------------------------------------------------------------
    def getMatchType(self,matchTypeFile=None) :
        """ return matchType - force read if not already ready """
        if len(self.matchType) <= 0 :
            self.readMatchType(matchTypeFile)
        # check again in case it faile 
        if len(self.matchType) <= 0 :
            print(' Could not read lat/lon - check files exist ')
            exit()
        return self.matchType

    #----------------------------------------------------------------------------------------
    # create sigma file names
    #----------------------------------------------------------------------------------------
    def sigmaFileName(self,sigmaAFile=None,sigmaRFile=None,myPath=None) :
        """ compute sigma file names or set values sigmaAFile,sigmaRFile """
        #
        # compute defaults first
        if 'azimuth' in self.azimuthFile :
            self.sigmaAFile=self.azimuthFile+'.sa'
            self.sigmaRFile=self.rangeFile+'.sr'
        elif '.da' in self.azimuthFile :
            self.sigmaAFile=self.azimuthFile.replace('.da','.sa')
            self.sigmaRFile=self.rangeFile.replace('.dr','.sr')
        # now over ride
        if sigmaAFile != None :
            self.sigmaAFile=sigmaAFile
        if sigmaRFile != None :
            self.sigmaRFile=sigmaRFile
        # ammend path if necessary
        if myPath != None :
            self.setMyPath(myPath)
        self.sigmaRFile = self.offFilePath(self.sigmaRFile) 
        self.sigmaAFile = self.offFilePath(self.sigmaAFile) 
        # now check path
        if not os.path.exists(self.sigmaAFile) or not os.path.exists(self.sigmaRFile) :
            if self.verbose :
                print('*** warning no sigma files exist ****')
        return

    #----------------------------------------------------------------------------------------
    # read sigma data
    #----------------------------------------------------------------------------------------
    def readSigma(self,sigmaAFile=None,sigmaRFile=None) :
        """ reads sigmaA/R files - assumes filenames set - but can override with sigmaAFile,sigmaRFile """
        #
        # set file name only if specified, assumes defaults tried in init
        if sigmaAFile != None or sigmaRFile != None :
            self.sigmaFileName(sigmaAFile=sigmaAFile,sigmaRFile=sigmaRFile)
        # read fresults after checking existence
        if os.path.exists(self.sigmaAFile) and  os.path.exists(self.sigmaRFile) :
            if self.verbose :
                print('Reading ',self.sigmaAFile, ' and ' ,self.sigmaRFile,' with size ',self.nr,self.na)
            self.sigmaA=readImage(self.sigmaAFile,self.nr,self.na,'>f4')
            self.sigmaR=readImage(self.sigmaRFile,self.nr,self.na,'>f4')
        else :
            print('\n\nError : Offsets tried to read an invalid sigmaA/R file - ',self.sigmaAFile,self.sigmaRFile)
            exit()
        return self.sigmaR,self.sigmaA
    
    #----------------------------------------------------------------------------------------
    # create data file names
    #----------------------------------------------------------------------------------------
    def datFileName(self,datFile=None,myPath=None) :
        """ compute dat file name - follows rules to set up based on standard names, but can specify overide name """
        if len(self.fileRoot) == 0 :
            print('\n\ndatFile Error: offsets no file Root selected')
            exit()
        # force to use specified value
        elif datFile != None :
            self.datFile=datFile
        # else use default
        elif 'azimuth' in self.fileRoot :
            self.datFile=self.fileRoot+'.dat'
        elif '.da'  in self.fileRoot and '.dat' not in self.fileRoot: 
            self.datFile=self.fileRoot.replace('.da','.da.dat')   
        # Override path
        if myPath != None :
            self.setMyPath(myPath)
        self.datFile = self.offFilePath(self.datFile) 
            
    def setFileRoot(self,fileRoot) :
        self.fileRoot=fileRoot
            
    #----------------------------------------------------------------------------------------
    # Offset file names
    #----------------------------------------------------------------------------------------
    def offsetFileNames(self,fileRoot,rangeFile=None,myPath=None,updateDatFileName=False) :
        """ Assumes fileRoot is the base name - estimates range name based on rules, or specify directly with rangeFile"""
    #
    # default azimuth.range
        self.setFileRoot(fileRoot)
        if len(fileRoot)== 0 :
            self.fileRoot='./azimuth.offsets'
            self.azimuthFile='./azimuth.offsets'
            self.rangeFile='./range.offsets'
        #
        # range file specified
        elif rangeFile != None :
            self.rangeFile=rangeFile
            self.azimuthFile=fileRoot
        #
        # assume azimuth/range offsets
        #
        elif 'azimuth' in fileRoot :
            self.azimuthFile=fileRoot
            self.rangeFile=fileRoot.replace('azimuth','range')
        # root is based on something ".da"
        #
        elif '.da' in fileRoot :
            self.azimuthFile=fileRoot
            self.rangeFile=fileRoot.replace('.da','.dr')
        #
        # update path if requested
        print('--',self.azimuthFile,self.rangeFile,myPath,self.path)
        if myPath != None :
            self.setMyPath(myPath)
        # join filename with path, this is the only place this should happen 
        self.azimuthFile = self.offFilePath(self.azimuthFile) 
        self.rangeFile = self.offFilePath(self.rangeFile)  
        #
        if updateDatFileName :
            self.datFileName()
            self.maskFileName()
#----------------------------------------------------------------------------------------
# Read the offset dat file
#----------------------------------------------------------------------------------------
    def readOffsetsDat(self,datFile=None) :
        """ read Offset dat file with either names setup earlier, or file set here """
        #
        # read offsets first
        if len(self.datFile) == 0 or datFile !=None :
            self.datFileName(datFile=datFile)
        # read the dat file
        if not os.path.exists(self.datFile) :
            print('\nError readOffsetsDat : tried to open an non-existent data file ',self.datFile)
            exit()
        fdat=open(self.datFile,'r')
        for line in fdat :
            if len(line) > 15 and (not ';' in line) :
                a=line.split()
                self.r0=int(a[0])
                self.a0=int(a[1])
                self.nr=int(a[2])
                self.na=int(a[3])
                self.dr=int(a[4])
                self.da=int(a[5])
                if len(a) == 7 :
                    self.azErr=a[6]
                break

#----------------------------------------------------------------------------------------
# Read the actual range and offset data
#----------------------------------------------------------------------------------------        
    def readOffsets(self,fileRoot=None,rangeFile=None,datFile=None) :
        """ Read da/dr offset file at previously set up offset names, or specify directly with values as defined for init"""
        # over ride names if needed.

        if len(self.azimuthFile) == 0 or rangeFile != None or fileRoot !=None :
            self.offsetFileNames(fileRoot,rangeFile=None)

        # read datFile if needed. 
        if self.nr < 1 or self.na < 1 :
            self.readOffsetsDat(datFile=datFile)
        #
        if os.path.exists(self.azimuthFile)  and os.path.exists(self.rangeFile) :
            if self.verbose :
                print('Reading ',self.azimuthFile,' with size ',self.nr,self.na)
            self.azOff=readImage(self.azimuthFile,self.nr,self.na,'>f4')
            self.azOff[np.isnan(self.azOff)]=-2.e9
            if self.verbose :
                print('Reading ',self.rangeFile,' with size ',self.nr,self.na)
            self.rgOff=readImage(self.rangeFile,self.nr,self.na,'>f4')
            self.rgOff[np.isnan(self.rgOff)]=-2.e9
        else :
            print('\n\nError : Offsets tried to read offset files - ',self.azimuthFile,self.rangeFile)
            exit()
        return self.rgOff.astype(float),self.azOff.astype(float)            

#----------------------------------------------------------------------------------------
# Read the actual range and offset data
#----------------------------------------------------------------------------------------
    def getOffsets(self,fileRoot=None,rangeFile=None,datFile=None) :
        """ Return offsets values as float """
        if len(self.rgOff) <= 0 :
            self.readOffsets(fileRoot,rangeFile,datFile) 
            # check again in case it faile
        if len(self.rgOff) <= 0 :
            print(' Could not read lat/lon - check files exist ')
            exit()
            # force conversion to 64 bit for coordinate transforms
        return self.rgOff.astype(float),self.azOff.astype(float)


#----------------------------------------------------------------------------------------
# remove offsets as specified by an logical array in a list (ie., everythign flattened)
#----------------------------------------------------------------------------------------
    def removeList(self,toRemove) :
    #
    #
        shapeSave=self.rgOff.shape
        if len(toRemove) <= 0 :
           return
        rg=self.rgOff.flatten()
        az=self.azOff.flatten()

        inRange= np.logical_and( toRemove >= 0, toRemove < len(rg) )
        use=toRemove[inRange]
        if len(use) > 0 :
           rg[use]=-2.e9
           az[use]=-2.e9
        self.rgOff=rg.reshape(shapeSave)
        self.azOff=az.reshape(shapeSave)
        
        if(len(self.matchType) > 0) :
           m=self.matchType.flatten()
           m[use]=0
           self.matchType=m.reshape(shapeSave)
#----------------------------------------------------------------------------------------
# remove offsets as specified by an logical array of same size
#----------------------------------------------------------------------------------------
    def remove(self,toRemove) :
    #
    #
        if(len(self.rgOff) > 0) :
            if toRemove.shape != self.rgOff.shape :
                print('removeOffsets : toRemove shape doesn not equal offset set shape')
                exit()
            self.rgOff[toRemove]=-2.e9
            self.azOff[toRemove]=-2.e9
        else :
            print('warning - no offsets remove because offsets not specified')
            return
        # do other types 
        if(len(self.sigmaR) > 0) :
            self.sigmaR[toRemove]=-2.e9
            self.sigmaA[toRemove]=-2.e9
        if(len(self.matchType) > 0) :
            self.matchType[toRemove]=0
#
#
#----------------------------------------------------------------------------------------
# write Offsets
#----------------------------------------------------------------------------------------
    def writeOffsets(self,fileRoot=None,rangeFile=None,datFile=None) :
        # check names defined
        if len(self.azimuthFile) < 1 or fileRoot != None :
            if fileRoot != None :
                self.offsetFileNames(fileRoot,rangeFile=rangeFile)
            else :
                self.offsetFileNames(self.fileRoot,rangeFile=rangeFile)
        azFile=self.azimuthFile
        rgFile=self.rangeFile
        if len(self.datFile) < 1 or datFile != None : 
            self.datFileName(datFile=datFile)
        # write dat file
        datFile=self.datFile
        fpDat=open(datFile,'w')
        if self.verbose :
            print(self.datFile) 
        if 'azimuth.offsets' not in self.datFile :
            # no azimuth, so output 6 entries
            print(self.r0, self.a0,self.nr,self.na,self.dr,self.da,file=fpDat)
        else :
            # azimuth so output azErr (even if zero)
            print(self.r0, self.a0,self.nr,self.na,self.dr,self.da,self.azErr,file=fpDat) 
        fpDat.close()
        # check if az/range pair so azError treated properly
        if 'azimuth.offsets' in self.datFile :
            # if azimuth, write the range with no azErr
            fpDat=open(datFile.replace('azimuth','range'),'w')
            print(self.r0, self.a0,self.nr,self.na,self.dr,self.da,file=fpDat)
            fpDat.close()
        # otherwise just make a copy
        else :
            call('cp '+datFile+' '+datFile.replace('.dat','.da.dat'),shell=True)
            call('cp '+datFile+' '+datFile.replace('.dat','.dr.dat'),shell=True)        
        # write offsets
        writeImage(rgFile,self.rgOff,'>f4')        
        writeImage(azFile,self.azOff,'>f4')
        if len(self.mask) > 0 and len(self.maskFile) > 0 :
            writeImage(self.maskFile,self.mask,'u1')
        if len(self.sigmaAFile) > 0 and len(self.sigmaRFile) > 0 and len(self.sigmaA) > 0 and len(self.sigmaR) > 0:
            if self.verbose :
                print('writing sigmaA ',self.sigmaAFile,' SigmaR ',   self.sigmaRFile)
            writeImage(self.sigmaAFile,self.sigmaA,'>f4')
            writeImage(self.sigmaRFile,self.sigmaR,'>f4')
#----------------------------------------------------------------------------------------
# compute satellite heading
# this program returns the satellite heading relative to north (in NH) or south (in SH)
#----------------------------------------------------------------------------------------
    def computeHeading(self,fileRoot=None) :
        if len(self.lat) ==0  or len(self.lon) == 0:
            print('Error in offsets.computeHeading : called without latitude loaded')
            exit()
        if self.slpAz < 0 :
            print('Error in offsets.computeHeading : no azimuth size, check geodatrxa specified')
        #
        if not self.geodatrxa.isRightLooking :
            print('offsets.computeHeading : Left looking heading not implemented yet')
            exit()
        #
        lattom=110947.
        skip=1
        azsp=self.slpAz * self.da*skip
        dlat=np.zeros(np.shape(self.lat))
        dlat[skip: ,:]=(self.lat[skip:,:]-self.lat[:-skip,:])*lattom 
        # slowly varying so replicate bottom line
        dlat[0:1,:]=dlat[skip,:]
        # make relative to south
        if self.geodatrxa.isSouth() :
            #if not self.geodatrxa.isDescending() :
                # for ascending, this will make dlat negative since moving away from South
                # for descending, this will make dlat positive since moving toward south
                dlat*=-1
        else :
            # for NH, this will make descending +
            if self.geodatrxa.isDescending() :
               dlat *= -1.0   
        #median filter
        self.heading=scipy.signal.medfilt2d((np.arccos(np.clip(dlat/azsp,-1.,1.))),kernel_size=5)
        # compute median and sigma
        med=np.median(self.heading)
        sig=np.std(self.heading)
        # clip to avoid extrem values
        self.heading=np.clip(self.heading,med-sig,med+sig)
        # rotate back to north (this seems to fix cos issue)
        if not self.geodatrxa.isSouth() and self.geodatrxa.isDescending()  :
            self.heading -= np.pi
        #print(3,dlat[0,0],azsp,self.heading[10,10])
        return self.heading

 #----------------------------------------------------------------------------------------
 # return whether south looking
 #----------------------------------------------------------------------------------------
    def isSouth(self) :
        return self.geodatrxa.isSouth()

#----------------------------------------------------------------------------------------
# compute single look pixel coordinates
#----------------------------------------------------------------------------------------
    def getRACoords(self) :
        if self.na < 1 or self.nr < 1 :
            print('offsets.slpCoords : tried to define coordinate with no size input - make sure .dat file read')
            exit()
        if len(self.rCoord) > 1 :
            return self.rCoord,self.aCoord
        
        rc=np.arange(0,self.nr*self.dr,self.dr)+self.r0
        ac=np.arange(0,self.na*self.da,self.da)+self.a0
        self.rCoord=np.zeros((self.na,self.nr))
        self.aCoord=np.zeros((self.na,self.nr))
        # range coord
        for i in range(0,self.nr):
            self.aCoord[:,i]=ac.copy()
        # azimuth coord
        for i in range(0,self.na) :
            self.rCoord[i,:]=rc.copy()
        # return
        return self.rCoord,self.aCoord
    

#----------------------------------------------------------------------------------------
# convert slp coords to pixels in offset image
#----------------------------------------------------------------------------------------

    def slpRAtoOffsetCoords(self,r,a) :
        """convert slp coords to pixels in offset image"""        
        ro = np.round((r-self.r0)/self.dr).astype(np.int)
        ao = np.round((a-self.a0)/self.da).astype(np.int)
        # check bounds
        rcheck=np.logical_and(ro >= 0,ro < self.nr)
        acheck=np.logical_and(ao >= 0,ao < self.na)
        good=np.logical_and(rcheck,acheck)
        # return        
        if len(good) > 0 :
            ro=ro[good]
            ao=ao[good]
            index=ao*self.nr + ro
            return ro,ao,index
        else :
            return [],[],[]


