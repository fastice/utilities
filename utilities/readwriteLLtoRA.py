import numpy as np
import os


def writeLLtoRAformat(lat, lon, tempfile=None):
    """ Write lat,lon array in binary format for lltora in32 nr ncol,
    lat,lon pairs as >fp4"""
    if tempfile is None:
        tempfile = 'temp.ll'
    # flatten arrays to columns
    lat1 = lat.flatten()
    lon1 = lon.flatten()
    # make nrowx2col array
    ll = np.zeros((len(lat1), 2), dtype='float64')
    ll[:, 0] = np.copy(lat1)
    ll[:, 1] = np.copy(lon1)
    ll = ll.byteswap(True)
    # open file
    fp = open(tempfile, 'w')
    # write header
    s = ll.shape
    s = np.int32(s)
    s = s.byteswap(True)
    s.tofile(fp)
    # write data
    ll.tofile(fp)
    # clean up
    fp.close()


def readLLtoRA(tempfile=None, returnZ=False):
    """ read r,a from the output of lltora """
    if tempfile is None:
        tempfile = 'temp.ra'
    # check file exits
    if not os.path.exists(tempfile):
        print('readLLtoRA : error ra file does not exist ')
        exit()
    # open file
    with open(tempfile, 'r') as fp:
        # size
        sz = np.fromfile(fp, count=2, dtype='>i4')
        # data
        data = np.fromfile(fp, dtype='>f8')
        # reshap8
        data = np.reshape(data, (sz[0], sz[1]))
        r = np.copy(data[:, 0])
        a = np.copy(data[:, 1])
        if returnZ:
            z = np.copy(data[:, 2])
            return r, a, z
    # return result
    return r, a
