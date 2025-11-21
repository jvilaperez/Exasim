import numpy as np
from pathlib import Path

def readbin(filename):
    """
    Reads a binary file of doubles, equivalent to MATLAB readbin().
    """
    with open(filename, "rb") as f:
        data = np.fromfile(f, dtype=np.float64)
    return data


def readxdgsol(filesol,ne):
    """
    Python equivalent of MATLAB readxdg.m

    Parameters
    ----------
    filesol : str
        Path to the binary solution file.

    Returns
    -------
    xdg : ndarray
        Array of shape (npe, ncx, ne)
    """

    tm = readbin(filesol)

    # Step 1: read header sizes
    sz = int(tm[0])

    k1 = 1
    k2 = k1 + sz
    nsize = tm[k1:k2].astype(int)

    k1 = k2
    k2 = k1 + nsize[0]
    ndims = tm[k1:k2].astype(int)

    k1 = k2
    k2 = k1 + nsize[1]
    xdg = tm[k1:k2]

    k1 = k2
    k2 = k1 + nsize[2]
    udg = tm[k1:k2]

    k1 = k2
    k2 = k1 + nsize[3]
    vdg = tm[k1:k2]


    # Step 2: reshape according to dimensions
    ncx = ndims[0]
    npe = int(len(xdg) / (ncx * ne))
    xdg = np.reshape(xdg, (npe, ncx, ne), order='F')

    ncu = ndims[1]
    udg = np.reshape(udg, (npe, ncu, ne), order='F')

    ncv = ndims[2]
    vdg = np.reshape(vdg, (npe, ncv, ne), order='F')


    
    return xdg, udg, vdg

def readxdgsolmpi(base, elemparts):
    """
    Python equivalent of MATLAB readxdgmpi.m

    Parameters
    ----------
    base : str
        Base filename (e.g., 'sol' for files like 'sol1.bin', 'sol2.bin', ...).
    elemparts : list
        List of element partitions.

    Returns
    -------
    xdg : ndarray
        Concatenated array of shape (npe, ncx, total_ne)
    """

    nprocs = len(elemparts)
    total_ne = len(np.unique(np.concatenate(elemparts[:])))
    
    if nprocs == 1:
        ne = len(elemparts[0])
        filesol = f"{base}.bin"
        if not Path(filesol).exists():
            raise FileNotFoundError(f"Cannot open file: {filesol}")
        xdg = readxdgsol(filesol)
        return xdg

    xdg0 = None
    udg0 = None
    vdg0 = None
    for i in range(1, nprocs + 1):
        ne = len(elemparts[i-1])
        filesol = f"{base}{i}.bin"
        if not Path(filesol).exists():
            raise FileNotFoundError(f"Cannot open file: {filesol}")
        xdgi, udgi, vdgi = readxdgsol(filesol,ne)
        xdgi = xdgi[:, :, :ne]
        udgi = udgi[:, :, :ne]
        vdgi = vdgi[:, :, :ne]

        if xdg0 is None:
            xdg0 = xdgi
        else:
            xdg0 = np.concatenate((xdg0, xdgi), axis=2)

        if udg0 is None:
            udg0 = udgi
        else:
            udg0 = np.concatenate((udg0, udgi), axis=2)

        if vdg0 is None:
            vdg0 = vdgi
        else:
            vdg0 = np.concatenate((vdg0, vdgi), axis=2)

    npe, ncx, _ = xdg0.shape
    _, ncu, _ = udg0.shape
    _, ncv, _ = vdg0.shape

    xdg = None if ncx == 0 else np.zeros((npe, ncx, total_ne))
    udg = None if ncu == 0 else np.zeros((npe, ncu, total_ne))
    vdg = None if ncv == 0 else np.zeros((npe, ncv, total_ne))

    ind0 = 0
    for i in range(nprocs):
        ne = len(elemparts[i])
        if xdg is not None:
            xdg[:, :, elemparts[i]] = xdg0[:, :, ind0:ind0+ne]
        if udg is not None:
            udg[:, :, elemparts[i]] = udg0[:, :, ind0:ind0+ne]
        if vdg is not None:
            vdg[:, :, elemparts[i]] = vdg0[:, :, ind0:ind0+ne] 
        ind0 += ne


    return xdg, udg, vdg