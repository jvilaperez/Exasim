import numpy as np
from pathlib import Path
from mkmaster import mkmaster 

def readbin(filename):
    """
    Reads a binary file of doubles, equivalent to MATLAB readbin().
    """
    with open(filename, "rb") as f:
        data = np.fromfile(f, dtype=np.float64)
    return data

def getmaster(filename):


    tm = readbin(filename)

    # Step 1: read header sizes
    sz = int(tm[0])

    k1 = 1
    k2 = k1 + sz
    nsize = tm[k1:k2].astype(int)

    k1 = k2
    k2 = k1 + nsize[0]
    ndims = tm[k1:k2].astype(int)

    nd = int(ndims[0])
    elemtype = int(ndims[1])
    nodetype = int(ndims[2])
    porder = int(ndims[3])
    pgauss = int(ndims[4])

    master = mkmaster(nd,porder,pgauss,elemtype,nodetype)

    return master