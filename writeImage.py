from utilities import myerror


def writeImage(fileName, x, dataType):
    """ write a binary image of size nx by ny with dataType = to one
    ['f4','>f4','>u2','u2','>i2','i2','>u4','u4','>i4','i4','u1'] """
#
# reads several types of binary images and creates a numpy matrix
#
    types = ['f4', '>f4', 'f8', '>f8', '>u2', 'u2', '>i2', 'i2', '>u4', 'u4',
             '>i4', 'i4', 'u1']
    if dataType not in types:
        myerror(f'writeImage: invalid data type - {dataType}')

    x1 = x
    types = {'f4': 'float32', 'f8': 'float64', 'i2': 'int16', 'i4': 'int32',
             'u1': 'uint8', 'u2': 'uint16', 'u4': 'uint32'}
    x1 = x1.astype(types[dataType.replace('>', '')])
    if '>' in dataType:
        x1.byteswap(True)

    fOut = open(fileName, 'w')
    x1.tofile(fOut)
    fOut.close()
    return
