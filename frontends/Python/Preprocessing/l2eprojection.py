import numpy as np


def l2eprojection(mesh, master, f, nc):
    """
    L2 element projection.

    Args:
        mesh:   dict with keys 'ne', 'dgnodes' (npv, nd+1, ne)
        master: dict with keys 'npv', 'ngv', 'nd', 'shapvt', 'shapvgdotshapvl', 'shapvg'
        f:      field at dg nodes (ngv, nc, ne)
        nc:     number of solution components

    Returns:
        UDG:    (npv, nc, ne)
    """
    ngv,nc,ne  = f.shape
    npv = master['npe']
    nd  = master['nd']

    # shapvt(:,:,1) -> index 0 in last axis, shape (ngv, npv)
    shapvt = np.squeeze(master['shapvt'][:, :, 0])              # (ngv, npv)

    # MATLAB: reshape(permute(shapvt(:,:,2:nd+1), [1 3 2]), [ngv*nd, npv])
    # permute [1 3 2] on (ngv, npv, nd) -> (ngv, nd, npv)
    # then reshape to (ngv*nd, npv)
    dshapvt = np.transpose(
        master['shapvt'][:, :, 1:nd+1], (0, 2, 1)
    ).reshape(ngv * nd, npv, order='F')                         # (ngv*nd, npv)

    UDG = np.zeros((npv, nc, ne), order='F')

    for i in range(ne):
        dg = mesh['dgnodes'][:, :, i]                           # (npv, nd+1)

        # Jacobian at Gauss points: (ngv*nd, npv) @ (npv, nd) -> (ngv*nd, nd)
        Jg = (dshapvt @ dg[:, :nd]).reshape(ngv, nd, nd, order='F')   # (ngv, nd, nd)
        jac = volgeom(Jg)                                              # (ngv,)

        pg = shapvt @ dg                                        # (ngv, nd+1)
        fg = f[:, :, i]                                         # (ngv, nc)

        # (npv, ngv) * (ngv,) -> (npv, ngv), then reshape to (npv, npv)
        M = (master['shapvgdotshapvl'][:, :, 0] * jac).reshape(npv, npv, order='F')

        for j in range(nc):
            F = master['shapvg'][:, :, 0] @ (jac * fg[:, j])   # (npv,)
            UDG[:, j, i] = np.linalg.solve(M, F)

    return UDG


def volgeom(Jg):
    """
    Compute Jacobian determinant at Gauss points.

    Args:
        Jg: (ngv, nd, nd)

    Returns:
        jac: (ngv,)
    """
    nd = Jg.shape[1]

    if nd == 1:
        jac = Jg[:, 0, 0]

    elif nd == 2:
        jac = Jg[:, 0, 0] * Jg[:, 1, 1] - Jg[:, 0, 1] * Jg[:, 1, 0]

    elif nd == 3:
        jac = (
            Jg[:, 0, 0] * Jg[:, 1, 1] * Jg[:, 2, 2] - Jg[:, 0, 0] * Jg[:, 2, 1] * Jg[:, 1, 2]
          + Jg[:, 1, 0] * Jg[:, 2, 1] * Jg[:, 0, 2] - Jg[:, 1, 0] * Jg[:, 0, 1] * Jg[:, 2, 2]
          + Jg[:, 2, 0] * Jg[:, 0, 1] * Jg[:, 1, 2] - Jg[:, 2, 0] * Jg[:, 1, 1] * Jg[:, 0, 2]
        )

    else:
        raise ValueError("Dimension not implemented")

    return jac