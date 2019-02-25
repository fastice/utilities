#!/usr/bin/python3
import utilities as u
import sys
import os
from subprocess import call

def runvel(dem,flags,inputFile,outputFile) :
    command='mosaic3d  -center ' + flags +  ' ' + inputFile+' ' + dem+' '+ outputFile
    print(command)
    command=command.split()

    try :
        fout=open(outputFile+'.stdout','w')
        ferr=open(outputFile+'.stderr','w')
        call(command,stdout=fout,stderr=ferr)
        fout.close()
        ferr.close()
    except Exception:
        print('error running ',command);
        exit()
        
    
