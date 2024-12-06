import numpy as np


def readImage(fileName, nx, ny, dataType):
    """ read a binary image of size nx by ny with dataType = to one
    ['f4','>f4','>u2','u2','>i2','i2','>u4','u4','>i4','i4','u1','>f8','f8']
    if f4 variant, conver to float for internal use   """
#
# reads several types of binary images and creates a numpy matrix
#
    types = ['f8', '>f8', 'f4', '>f4', '>u2', 'u2', '>i2', 'i2', '>u4',
             'u4', '>i4', 'i4', 'u1']
    if dataType not in types:
        print(f'\nError readImage: Specified data type,\033[1m'
              f'{dataType}\033[0m not in accepted types :\n\n\t{types}\n')
        exit()

    dt = np.dtype(dataType)
    x = np.fromfile(fileName, dtype=dt)[0:nx*ny]
    x = np.reshape(x, [ny, nx])
    # print('Data Type ',dataType)
    # swap data so its in native format
    if '>' in dataType:
        x = x.astype(np.dtype(dataType.replace('>', '')))
    # print(x.dtype)
    return x
