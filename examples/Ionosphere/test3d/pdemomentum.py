from numpy import array, reshape, hstack, zeros
from sympy import sqrt, exp, log

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0])
    return m

def flux(u, q, w, v, x, t, mu, eta):
    gam = mu[0]
    
    u1 = u[0]
    u2 = u[1]
    u3 = u[2]
    
    T = v[10]
    p = T/gam

    fi = array([u1*u1+p, u1*u2, u1*u3,
                u1*u2, u2*u2+p, u2*u3,
                u1*u3, u2*u3, u3*u3+p])
    
    f = fi
    f = reshape(f,(3,3),'F')
    return f

def source(u, q, w, v, x, t, mu, eta):
    gam = mu[0]
    Ro = mu[1]

    R0 = mu[8]
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    r  = sqrt(x1**2 + x2**2 + x3**2)
    g0 = 1/gam
    g = g0*(R0/r)**2

    u1 = u[0]
    u2 = u[1]
    u3 = u[2]
    T = v[10]
    p = T/gam

    u1x = -q[0]
    u2x = -q[1]
    u3x = -q[2]

    u1y = -q[3]
    u2y = -q[4]
    u3y = -q[5]

    u1z = -q[6]
    u2z = -q[7]
    u3z = -q[8]

    # acceleration
    ax = -g*x1/r + x1/Ro**2 + 2*u2/Ro
    ay = -g*x2/r + x2/Ro**2 - 2*u1/Ro
    az = -g*x3/r

    # divergence of velocity
    divV = u1x + u2y + u3z

    # Lorentz force: E + vxB
    piE = mu[2]
    piB = mu[3]
    Bx = v[0]
    By = v[1]
    Bz = v[2]
    Ex = v[3]
    Ey = v[4]
    Ez = v[5]
    fLx = piE*Ex + piB*(u2*Bz - u3*By)
    fLy = piE*Ey + piB*(u3*Bx - u1*Bz)
    fLz = piE*Ez + piB*(u1*By - u2*Bx)

    # Density gradient
    lnx = v[7]
    lny = v[8]
    lnz = v[9]

    # Ion drag
    piNu = mu[4]
    nu = v[11]
    un1 = v[12]
    un2 = v[13]
    un3 = v[14]
    id1 = piNu*nu*(un1-u1)
    id2 = piNu*nu*(un2-u2)
    id3 = piNu*nu*(un3-u3)


    # source terms
    su1 = ax + divV*u1 - p*lnx + fLx + id1
    su2 = ay + divV*u2 - p*lny + fLy + id2
    su3 = az + divV*u3 - p*lnz + fLz + id3

    s = array([su1, su2, su3])
    return s

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = zeros((3,2))
    return ub

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fb = zeros((3,2))
    return fb

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    # Lower boundary: vi = un
    un1 = v[12]
    un2 = v[13]
    un3 = v[14]

    fl = 0*u
    fl[0] = un1 - uhat[0]
    fl[1] = un2 - uhat[1]
    fl[2] = un3 - uhat[2]

    # Upper boundary:
    # vperp = ExB/B^2 & d(vpar)/db = 0
    Bx = v[0]
    By = v[1]
    Bz = v[2]
    Ex = v[3]
    Ey = v[4]
    Ez = v[5]

    B2 = Bx*Bx + By*By + Bz*Bz
    B = sqrt(B2)
    b1 = Bx/B
    b2 = By/B
    b3 = Bz/B

    drift1 = (Ey*Bz - Ez*By)/B2
    drift2 = (Ez*Bx - Ex*Bz)/B2
    drift3 = (Ex*By - Ey*Bx)/B2

    vhatb = uhat[0]*b1 + uhat[1]*b2 + uhat[2]*b3
    Fperp1 = uhat[0] - vhatb*b1 - drift1
    Fperp2 = uhat[1] - vhatb*b2 - drift2
    Fperp3 = uhat[2] - vhatb*b3 - drift3


    u1x = -q[0]
    u2x = -q[1]
    u3x = -q[2]

    u1y = -q[3]
    u2y = -q[4]
    u3y = -q[5]

    u1z = -q[6]
    u2z = -q[7]
    u3z = -q[8]

    bGb = b1*(u1x*b1 + u2x*b2 + u3x*b3) + b2*(u1y*b1 + u2y*b2 + u3y*b3) + b3*(u1z*b1 + u2z*b2 + u3z*b3)
    vb = u[0]*b1 + u[1]*b2 + u[2]*b3

    Fpar = bGb + tau[0]*(vb-vhatb)
    Fpar1 = Fpar*b1
    Fpar2 = Fpar*b2
    Fpar3 = Fpar*b3

    fu = 0*u
    fu[0] = Fperp1 + Fpar1
    fu[1] = Fperp2 + Fpar2
    fu[2] = Fperp3 + Fpar3

    fb = hstack((fl, fu))
    fb = reshape(fb,(3,2),'F')
    return fb

def initu(x, mu, eta):
    u0 = array([0.0, 0.0, 0.0])
    return u0
