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

    # Step 2: reshape according to dimensions
    ncx = ndims[0]
    npe = int(len(xdg) / (ncx * ne))

    xdg = np.reshape(xdg, (npe, ncx, ne), order='F')
    return xdg

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
    for i in range(1, nprocs + 1):
        ne = len(elemparts[i-1])
        filesol = f"{base}{i}.bin"
        if not Path(filesol).exists():
            raise FileNotFoundError(f"Cannot open file: {filesol}")
        xdgi = readxdgsol(filesol,ne)
        xdgi = xdgi[:, :, :ne]

        if xdg0 is None:
            xdg0 = xdgi
        else:
            xdg0 = np.concatenate((xdg0, xdgi), axis=2)

    npe, ncx, _ = xdg0.shape
    xdg = np.zeros((npe, ncx, total_ne))
    ind0 = 0
    for i in range(nprocs):
        ne = len(elemparts[i])
        xdg[:, :, elemparts[i]] = xdg0[:, :, ind0:ind0+ne]
        ind0 += ne

    return xdg