# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 14:33:23 2018

@author: ian
"""
from utilities import myerror
import numpy as np
class argusMesh :
    # nodes
    def __init__(self,expFile=None,meshType='Triangle') :
        self.nNodes,self.nElements=0,0
        self.nodes,self.elements=[],[]
        self.nodeValues,self.elementValues=[],[]
        self.nNVal,self.nEVal=0,0
        self.meshpy=None
        self.meshtype=meshType
        
        if expFile !=None :
            self.readExp(expFile)
        
    #
    #   Read an argus mesh file       
    def readExp(self,expFile) :
        try :
            fpExp=open(expFile,'r')
        except :
            myerror('argusMesh.readExp: could not open meshfile for read')
        nodesRead,elementsRead=0,0
        for line in fpExp : 
            # process first line
            pieces=line.split()
            nPieces=len(pieces)
            if nPieces == 0 :
                continue
            #
            # parse header
            if self.nNodes == 0 and nPieces == 4 :
                self.nNodes=int(pieces[1])
                self.nElements=int(pieces[0])
                self.nNVal=int(pieces[3])
                self.nEVal=int(pieces[2])
                self.nodeValues=[[] for i in range(self.nNVal)]
                self.elementValues=[[] for i in range(self.nEVal)]
                print(self.nodeValues)
     
                print('Number of Nodes : ',self.nNodes,'Number of Elements : ',self.nElements,'Number of node/element values',self.nNVal,self.nEVal)  
            #
            # parse nodes 
            elif pieces[0] == 'N' :
                # append coordinates for a point 
                self.nodes.append((float(pieces[2]),float(pieces[3])))
                for i in range(self.nNVal) :
                    self.nodeValues[i].append(float(pieces[4+i]))
                #print(self.nodeValues)
                nodesRead+=1
            #
            # parse elements
            elif  pieces[0] == 'E' :
                # parse elements, substract one to zero index
                self.elements.append( [ int(pieces[2]),int(pieces[3]),int(pieces[4])] )
                for i in range(self.nEVal) :
                    self.elementValues[i].append(float(pieces[5+i]))
                elementsRead+=1
        print(nodesRead,self.nNodes,elementsRead,self.nElements)
        #print(self.nodeValues)
        return
        
    #
    # convert Meshpy
    def toGmsh(self,gmshfile) :
        #
        # init
        try :
            fp=open(gmshfile,'w')
        except :
            myerror('error opening gmsh output file')
        #
        # Format
        print('$MeshFormat',file=fp)
        print('2.2 0 8',file=fp)
        print('$EndMeshFormat',file=fp)
        #
        # Nodes
        print('$Nodes',file=fp)
        print(self.nNodes,file=fp)
        for i in range(1,self.nNodes+1) :
            print(i,'{} {}'.format(*self.nodes[i-1]),0.0,file=fp)
            #print(i,'{} {}'.format(*self.nodes[i-1]),0.0)
        print('$EndNodes',file=fp)
        #
        # Elements
        print('$Elements',file=fp)
        print(self.nElements,file=fp)
        if self.meshtype == 'Triangle' :
            for i in range(1,self.nElements+1) :
                print(i,2,2,99,2,'{} {} {}'.format(*self.elements[i-1]),file=fp)         
        else :
            myerror('Invalid mesh type')
        print('$EndElements',file=fp)
        #
        # Node data
        print('$NodeData',file=fp)
        # string tag
        print(1,file=fp)
        print('"Argus Node data"',file=fp)
        # 1 dummy real tag
        print(1,file=fp)
        print(0.0,file=fp)
        # integer tags
        print(3,file=fp)
        print(0,file=fp) # dummy time step
        print(self.nNVal,file=fp) # number of values per node
        print(self.nNodes,file=fp) # number of nodes; assume all
        data=np.array(self.nodeValues)
        for i in range(1,self.nNodes+1) :
            print(i,*data[:,i-1] ,file=fp)
        print('$EndNodeData',file=fp)
        #
        # Element data
        print('$ElementData',file=fp)
        # string tag
        print(1,file=fp)
        print('"Argus Element data"',file=fp)
        # 1 dummy real tag
        print(1,file=fp)
        print(0.0,file=fp)
        # integer tags
        print(3,file=fp)
        print(0,file=fp) # dummy time step
        print(self.nEVal,file=fp) # number of values per elements
        print(self.nElements,file=fp) # number of elements; assume all
        dataE=np.array(self.elementValues)
        for i in range(1,self.nElements+1) :
            print(i,*dataE[:,i-1] ,file=fp)
        print('$EndNodeData',file=fp)    
        
   
            
            