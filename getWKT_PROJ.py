# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 10:50:56 2017

@author: ian
"""
#  urllib.request
import pyproj


def getWKT_PROJ(epsg_code):
    '''
        Get wkt for projection in gdal format, which is needed to use
        in shapefile for QGIS
    '''
    crs = pyproj.CRS.from_epsg(f'{epsg_code}')
    return crs.to_wkt(pyproj.enums.WktVersion.WKT1_GDAL)
    # Obsolete
    # output.replace("\\n","").replace("b\'","").replace("\'","")
    # return  output
    # url="http://spatialreference.org/ref/epsg/{0}/prettywkt/".
    # format(epsg_code)
    # print(url)
    # response=urllib.request.urlopen(url)
    # remove_spaces =str(response.read()).replace(" ","")
    # output=remove_spaces
    # return output
