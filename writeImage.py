#writeImage.py
import numpy as np

def writeImage(fileName,x,dataType) :
    """ write a binary image of size nx by ny with dataType = to one ['f4','>f4','>u2','u2','>i2','i2','>u4','u4','>i4','i4','u1'] """
#
# reads several types of binary images and creates a numpy matrix
#
    types=['f4','>f4','f8','>f8','>u2','u2','>i2','i2','>u4','u4','>i4','i4','u1']
    if not dataType in types :
        print('\nError readImage: Specified data type, \033[1m',dataType,'\033[0m not in accepted types :\n\n\t',types,'\n')
        exit()
    x1=x
    if 'f4' in dataType :
        x1=x.astype('float32')
    if 'f8' in dataType :
        x1=x.astype('float64')        
    elif 'i2' in dataType :
        x1=x.astype('int16')
    elif 'i4' in dataType :
        x1=x.astype('int32')
    elif 'u2' in dataType :
        x1=x.astype('uint16')
    elif 'u4' in dataType :
        x1=x.astype('uint32')        
    
    if '>' in dataType :
        x1.byteswap(True)

    sh=x1.shape
  #  x1=x1.transpose()
    dt=np.dtype(dataType)
    fOut=open(fileName,'w')
    x1.tofile(fOut)
    fOut.close()
    return
