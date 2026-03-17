import numpy as np


def gradu(shapen, dgnodes, udg):
    """
    Compute gradient of solution.
    
    Args:
        shapen:   (npe, npe, nd)   - nodal shape functions by dimension
        dgnodes:  (npe, nd, ne)    - DG node coordinates
        udg:      (npe, ne)   - solution values
    
    Returns:
        qdg:      (npe, ncu*nd, ne)
    """
    ncu = 1
    npe, ne = udg.shape
    nd = shapen.shape[2]

    Xx = volgeom(shapen, dgnodes)  # (npe*ne, nd, nd)

    tmp = np.zeros((npe, ne, nd))
    for i in range(nd):
        tmp[:, :, i] = shapen[:, :, i] @ udg

    tmp = tmp.reshape(npe * ne, 1, nd)

    qdg = np.zeros((npe * ne, 1, nd))
    for k in range(nd):
        for i in range(nd):
            qdg[:, 0, k] += tmp[:, 0, i] * Xx[:, k, i]

    qdg = qdg.reshape(npe, ne, nd)
    # permute(qdg,[0 2 1]) -> (npe, nd, ne)
    qdg = np.transpose(qdg, (0, 2, 1))

    return qdg


def volgeom(shapen, dgnodes):
    """
    Compute inverse Jacobian (scaled by determinant) at quadrature points.
    
    Args:
        shapen:   (npe, npe, nd)
        dgnodes:  (npe, nd, ne)
    
    Returns:
        Xx:       (npe*ne, nd, nd)  - cofactor matrix / Jacobian determinant
    """
    npe, nd, ne = dgnodes.shape

    # MATLAB: permute(dgnodes,[1 3 2]) -> (npe, ne, nd)
    dgnodes_p = np.transpose(dgnodes, (0, 2, 1))

    # MATLAB: reshape(permute(shapen,[1 3 2]),[npe*nd, npe])
    # permute(shapen,[1 3 2]): (npe, nd, npe) -> reshape: (npe*nd, npe)
    dshapvt = np.transpose(shapen, (0, 2, 1)).reshape(npe * nd, npe)

    # Jg computation: (npe*nd, npe) @ (npe, ne*nd) -> (npe*nd, ne*nd)
    Jg = dshapvt @ dgnodes_p.reshape(npe, ne * nd)

    # MATLAB: permute(reshape(Jg,[npe nd ne nd]),[1 3 2 4])
    Jg = Jg.reshape(npe, nd, ne, nd)
    Jg = np.transpose(Jg, (0, 2, 1, 3))   # (npe, ne, nd, nd)
    Jg = Jg.reshape(npe * ne, nd, nd)

    if nd == 1:
        jac = Jg[:, 0, 0]
        Xx = np.ones((npe * ne, 1, 1))

    elif nd == 2:
        jac = Jg[:, 0, 0] * Jg[:, 1, 1] - Jg[:, 0, 1] * Jg[:, 1, 0]
        Xx = np.zeros((npe * ne, 2, 2))
        Xx[:, 0, 0] =  Jg[:, 1, 1]   # dxi/dx
        Xx[:, 1, 0] = -Jg[:, 1, 0]   # dxi/dy
        Xx[:, 0, 1] = -Jg[:, 0, 1]   # deta/dx
        Xx[:, 1, 1] =  Jg[:, 0, 0]   # deta/dy

    elif nd == 3:
        jac = (
            Jg[:, 0, 0] * Jg[:, 1, 1] * Jg[:, 2, 2] - Jg[:, 0, 0] * Jg[:, 2, 1] * Jg[:, 1, 2]
          + Jg[:, 1, 0] * Jg[:, 2, 1] * Jg[:, 0, 2] - Jg[:, 1, 0] * Jg[:, 0, 1] * Jg[:, 2, 2]
          + Jg[:, 2, 0] * Jg[:, 0, 1] * Jg[:, 1, 2] - Jg[:, 2, 0] * Jg[:, 1, 1] * Jg[:, 0, 2]
        )
        Xx = np.zeros((npe * ne, 3, 3))
        Xx[:, 0, 0] = Jg[:, 1, 1] * Jg[:, 2, 2] - Jg[:, 1, 2] * Jg[:, 2, 1]
        Xx[:, 1, 0] = Jg[:, 1, 2] * Jg[:, 2, 0] - Jg[:, 1, 0] * Jg[:, 2, 2]
        Xx[:, 2, 0] = Jg[:, 1, 0] * Jg[:, 2, 1] - Jg[:, 1, 1] * Jg[:, 2, 0]
        Xx[:, 0, 1] = Jg[:, 0, 2] * Jg[:, 2, 1] - Jg[:, 0, 1] * Jg[:, 2, 2]
        Xx[:, 1, 1] = Jg[:, 0, 0] * Jg[:, 2, 2] - Jg[:, 0, 2] * Jg[:, 2, 0]
        Xx[:, 2, 1] = Jg[:, 0, 1] * Jg[:, 2, 0] - Jg[:, 0, 0] * Jg[:, 2, 1]
        Xx[:, 0, 2] = Jg[:, 0, 1] * Jg[:, 1, 2] - Jg[:, 0, 2] * Jg[:, 1, 1]
        Xx[:, 1, 2] = Jg[:, 0, 2] * Jg[:, 1, 0] - Jg[:, 0, 0] * Jg[:, 1, 2]
        Xx[:, 2, 2] = Jg[:, 0, 0] * Jg[:, 1, 1] - Jg[:, 0, 1] * Jg[:, 1, 0]

    else:
        raise ValueError("Dimension not implemented")

    Xx = Xx / jac[:, np.newaxis, np.newaxis]

    return Xx