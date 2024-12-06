# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 07:31:31 2019

@author: ian
"""
import numpy as np


def processProfile(xls, yls, deltaX):
    ''' For a series lines segments forming a polyline (xl,yl),
    sample with nearly uniform spacing deltaX'''
    # Setup return values
    xdom, ydom, dist = np.array([]), np.array([]), np.array([])
    # Compute distances for each segment
    dxs, dys = np.diff(xls), np.diff(yls)
    distances = np.cumsum(np.sqrt(dxs**2 + dys**2))
    pointCounts = np.rint(distances/deltaX)
    # Process each segment
    for dist, dx, dy, xl, yl, nPts in zip(distances, dxs, dys, xls, yls,
                                          pointCounts):
        if nPts > 0:
            # Step size
            dxN = dx/nPts
            dyN = dy/nPts
            #
            xdom = np.concatenate((xdom, xl + dxN*np.arange(0, nPts)))
            ydom = np.concatenate((ydom, yl + dyN*np.arange(0, nPts)))
    # append last point
    xdom = np.concatenate((xdom, xls[-1:]))
    ydom = np.concatenate((ydom, yls[-1:]))
    # Compute distances along profiles
    d = np.cumsum(np.concatenate(([0], np.sqrt(np.diff(xdom)**2 +
                                               np.diff(ydom)**2))))
    return xdom, ydom, d