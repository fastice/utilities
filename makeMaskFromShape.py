# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 08:04:19 2017

@author: ian
"""
import shapefile
import numpy as np
try:
    from PIL import Image, ImageDraw
except Exception:
    print('couldnot open PIL - may cause problems')
from utilities import myerror
import os


def makeMaskFromShape(geo, shpfile):
    '''
    For a geodat specified by geo, produce a corresponding mask from a shpfile
    (inlucd .shp)
    mask is 0 for exclude, 1 for keep
    If mask file is zero length, create mask of ones
    '''
    #
    # loop through features
    nx, ny = geo.sizeInPixels()
    # create PIL image, which is widthxheight, L indicates 8-bit
    mask = Image.new('L', (nx, ny), 1)
    # e.g., < .shp
    if len(shpfile) < 4:
        return []
    if not os.path.exists(shpfile):
        myerror('Shape file : ' + shpfile+' does not exist')
    #
    # Open shape file
    shape = shapefile.Reader(shpfile)
    #
    # Loop over each feature in shape file
    for feature in shape.shapeRecords():
        # extract polygon from feature
        poly = feature.shape.__geo_interface__
        xyPoly = np.array(poly['coordinates'])[0]
        # check it 2D
        if len(xyPoly.shape) == 2:
            # get image coords for polygon
            xi, yi = np.rint(geo.xymtoImage(xyPoly[:, 0], xyPoly[:, 1]))
            # create array of tuple: poly=[ (x1,y1), (x2,y2) ...]
            poly = list(zip(xi, yi))
            # burn mask in polygon
            ImageDraw.Draw(mask).polygon(poly, outline=0, fill=0)
    # create np array of mask
    # goes from [x][y] to [y][x]
    mask = np.array(mask, dtype='u1')
    return mask
