# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 15:10:07 2018

@author: ian
"""
from datetime import datetime 
from utilities import myerror
class logger :

    """ Start a log file and allow continued writes """

    def __init__(self,fileRoot=None,echo=False) :
        #
        if fileRoot == None :
            fileRoot='log'
        self.logFile='log.{0:s}.{1:s}'.format(fileRoot,datetime.now().strftime('%Y-%m-%d.%H-%M-%S'))
        try :
            print('Opening log file: '+self.logFile)
            self.fp=open(self.logFile,'w')
        except :
            myerror('Could not open logfile {0:s}'.format(self.logFile))
            
    def logLine(self,myMsg,logTime=True) :
        # print message
        myTime=''
        if logTime :
            myTime=datetime.now().strftime('%Y-%m-%d.%H:%M:%S')
            print('{0:s} -- {1:s}'.format(myTime,myMsg),file=self.fp)
            self.fp.flush()
        else :
            print(myMsg,file=self.fp)
            self.fp.flush()
            
    def logfile(self) :
        return self.logFile
        
    def logError(self,myMsg) :
        self.logLine('*** Terminating With Error : {0:s}'.format(myMsg),logTime=True)
        self.closeLog()
            
    def closeLog(self) :
        self.logLine('Log closed',logTime=True)
        self.fp.close()
        
    def logEntry(self,progName,logTime=True) :
        self.logLine('Entering {0:s} '.format(progName),logTime=logTime)
        
    def logRunProgram(self,progName,logTime=True) :
        self.logLine('Calling {0:s} '.format(progName),logTime=logTime)
        
    def logReturn(self,progName,logTime=True) :
        self.logLine('Returning from {0:s} '.format(progName),logTime=logTime)
        
    def logArgs(self,args) :
        ''' record arguments from arg parser '''
        argDict=vars(args)
        self.logLine('------------- Arguments---------------')
        #print(argDict)
        for myArg in argDict :
            if type(argDict[myArg]) == list :
                argList=' '.join([str(x) for x in argDict[myArg]])
            else :
                argList=str(argDict[myArg])
            print(myArg,argList)
            self.logLine('{0:s} = {1:s}'.format(myArg,argList))
        self.logLine('-------------End Arguments---------------')