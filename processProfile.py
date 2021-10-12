# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 07:31:31 2019

@author: ian
"""
import numpy as np


def processProfile(xl, yl, deltaX):
    ''' For a series lines segments forming a polyline (xl,yl), uniformly
    sample with spacing deltaX'''
    # compute number of segments
    nSegs = len(xl) - 1
    xdom, ydom, dist = [], [], []
    # approximate uniform sampling of points
    for i in range(0, nSegs):
        dist = np.sqrt((xl[i+1] - xl[i])**2 + (yl[i+1] - yl[i])**2)
        nPts = int(round(dist/deltaX))
        if nPts > 0:
            dxN = (xl[i+1]-xl[i])/nPts
            dyN = (yl[i+1]-yl[i])/nPts
            for j in range(0, nPts):
                xdom.append(xl[i]+j*dxN)
                ydom.append(yl[i]+j*dyN)
    xdom.append(xl[-1])
    ydom.append(yl[-1])
    #
    # Compute distances along profiles
    d = [0]
    print(xdom[0], ydom[0], d[0])
    for i in range(1, len(xdom)):
        delta = np.sqrt((xdom[i] - xdom[i-1])**2 + (ydom[i] - ydom[i-1])**2)
        d.append(d[i-1] + delta)
        # print(xdom[i],ydom[i],d[i])

    return np.array(xdom), np.array(ydom), np.array(d)
