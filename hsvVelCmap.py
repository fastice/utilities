# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 13:50:06 2017

@author: ian
"""
import matplotlib.pylab as plt
import matplotlib.colors as colors
import numpy as np

def hsvVelCmap(transparency=None):
    # default transparency
    if transparency == None :
        transparency=0.9
    # get hsv to start with
    cmap=plt.get_cmap('hsv',360)
    # pull the part we want
    cmapPart=cmap(np.arange(23,360,1))
    # stretch the map
    cmap2=colors.LinearSegmentedColormap.from_list('velcol',cmapPart)
    # set the transparency
    cmap2._init()
    cmap2._lut[:,-1]=0.9
    cmap2._lut[0,-1]=0 
    #
    return cmap2