# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:06:53 2018

@author: ian
"""
import subprocess
import os

def callMyProg(myProg,myInputs=None,myArgs=None,screen=False,logger=None,env=None) :
    ''' call a program with optional 'myInputs' for user input, 'myArgs' for 
    command line args, and 'screen' to output io to the screen. Returns the programs stdout'''
    myCommand=[myProg]
    if myArgs != None :
        myCommand += myArgs
    #
    if logger != None :
        logger.logRunProgram(myProg)
        logger.logLine('\t'+' '.join(myCommand))
    #print(myCommand)
    if env == None :
        env=os.environ.copy()
    if screen :
        proc=subprocess.Popen(myCommand,stdin=subprocess.PIPE,env=env)
    else :
        proc=subprocess.Popen(myCommand,stdin=subprocess.PIPE,stdout=subprocess.PIPE,env=env)
    
    if myInputs != None :
        fullInput=''
        for myInput in myInputs :
            fullInput+=' '.join([str(x) for x in myInput])+'\n'
        myOutput=proc.communicate(fullInput.encode('utf-8'))[0].decode('utf-8')
        #print(myOutput)
    else :
        myOutput=proc.communicate()
        #print(type(myOutput))
        #print(myOutput)
        if myOutput[0] != None :
            myOutput=myOutput[0].decode('utf-8')
    #print(myOutput,type(myOutput))
    if logger != None :
        logger.logReturn(myProg)     
    return myOutput