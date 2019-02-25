# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 10:50:56 2017

@author: ian
"""
import urllib.request

def getWKT_PROJ(epsg_code) :
    url="http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg_code)
    print(url)
    response=urllib.request.urlopen(url)
    remove_spaces =str(response.read()).replace(" ","")
    output=remove_spaces.replace("\\n","").replace("b\'","").replace("\'","")
    return output