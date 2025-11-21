import numpy as np

def getnodes(xdg,elemtype,porder):

    npe, nsd, ne = xdg.shape

    if nsd==2:
        if elemtype==0:
            # error: not defined
            raise ValueError(f"elemtype {elemtype} not defined for nsd=2")
        elif elemtype==1:
            nodestruct = [0,porder,(porder+1)**2-1,(porder+1)*porder]        
    elif nsd==3:
        if elemtype==0:
            # error: not defined
            raise ValueError(f"elemtype {elemtype} not defined for nsd=3")
        elif elemtype==1:
            nodestruct = [0,porder,(porder+1)**2-1,(porder+1)*porder, porder*(porder+1)**2, porder*(porder+1)**2 + porder, (porder+1)**3-1, porder*(porder+1)*(porder+2)]
        
    p = xdg[(nodestruct),:,:].transpose(1, 0, 2).reshape(nsd, -1)
    tol = 1e-6  # tolerance for uniqueness
    p1 = np.round(p / tol) * tol
    p, t = np.unique(p1.T, axis=0, return_inverse=True)
    p = p.T
    t = t.reshape(len(nodestruct), -1)

    return p, t