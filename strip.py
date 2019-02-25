# strip.py
# This function will read columns of data. WIth the format
#  # n  (n is number of columns)
#   ; comment (; delimited comment)
#   x y z (whitespace separated numbers)
#   & (& indicates end of data)
#
import numpy as np


def strip(file):
    # define data
    data = []
    # Open file
    f = open(file, 'r')
    i = 0
    ncol = 0
    # cycle through lines
    for line in iter(f):
        # break if end of data
        if "&" in line:
            break
        # process lines that are not comments
        if ";" not in line:
            tmp = line.split()
            if '#' in line:
                print(tmp[0], tmp[1])
                ncol = int(tmp[1])
            elif ncol > 0:
                tmp = tmp[0:ncol]
                data.append(tmp)
                i += 1

    # make a numpy array out of the result
    data1 = np.array(data)
    data1=data1.astype(np.float)
    return data1
