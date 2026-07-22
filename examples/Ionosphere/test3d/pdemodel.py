from numpy import ones, reshape, array, zeros, hstack
from sympy import sqrt, exp

def mass(u, q, w, v, x, t, mu, eta):
    m = ones(shape=10)
    return m

def flux(u, q, w, v, x, t, mu, eta):
    m_op = mu[3] # molar mass of O+, 16
    m_o2p = mu[4] # molar mass of O2+, 32
    m_np = mu[5] # molar mass of N+, 28
    m_n2p = mu[6] # molar mass of N2+, 28
    m_nop = mu[7] # molar mass of NO+, 40
    gamma_i = mu[8] # ion adiabatic index, 5/3
    gamma_e = mu[9] # electron adiabatic index, 5/3
    lamda = mu[10] # kB*T0*t0**2/(mp*H0**2)

    B1 = v[5]
    B2 = v[6]
    B3 = v[7]
    kappa_i = v[26] # ion thermal conduction coefficient
    kappa_e = v[27] # electron thermal conduction coefficient

    logop = u[0]
    logo2p = u[1]
    lognp = u[2]
    logn2p = u[3]
    lognop = u[4]
    v1 = u[5]
    v2 = u[6]
    v3 = u[7]
    ti = u[8]
    te = u[9]

    tix = -q[8]
    tex = -q[9]
    tiy = -q[18]
    tey = -q[19]
    tiz = -q[28]
    tez = -q[29]

    Bmag = sqrt(B1**2 + B2**2 + B3**2)
    b1 = B1/Bmag
    b2 = B2/Bmag
    b3 = B3/Bmag
    op = exp(logop)
    o2p = exp(logo2p)
    np = exp(lognp)
    n2p = exp(logn2p)
    nop = exp(lognop)
    ne = op + o2p + np + n2p + nop
    mbar = (op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop) / ne # mean molar mass
    p = lamda*(ti + te)/mbar
    qi = -(gamma_i-1)*kappa_i*ti**2*sqrt(ti)*(b1*tix + b2*tiy + b3*tiz) # ion heat flux
    qe = -(gamma_e-1)*kappa_e*te**2*sqrt(te)*(b1*tex + b2*tey + b3*tez) # electron heat flux

    # logn*v, v*v+p.I, ne*ti*v+qi*b, qe*b
    f = reshape([logop*v1, logo2p*v1, lognp*v1, logn2p*v1, lognop*v1, v1**2+p, v1*v2, v1*v3, ne*ti*v1+qi*b1, qe*b1,
                logop*v2, logo2p*v2, lognp*v2, logn2p*v2, lognop*v2, v1*v2, v2**2+p, v2*v3, ne*ti*v2+qi*b2, qe*b2,
                logop*v3, logo2p*v3, lognp*v3, logn2p*v3, lognop*v3, v1*v3, v2*v3, v3**2+p, ne*ti*v3+qi*b3, qe*b3],
                shape=(10, 3), order='F')
    return f

def source(u, q, w, v, x, t, mu, eta):
    r0 = mu[0] # earth radius (non-dimensionalized): Re/H0
    g0 = mu[1] # gravity at earth surface (non-dimensionalized): g0*t0**2/H0
    wrot = mu[2] # earth rotation (non-dimensionalized): wrot*t0
    m_op = mu[3] # molar mass of O+, 16
    m_o2p = mu[4] # molar mass of O2+, 32
    m_np = mu[5] # molar mass of N+, 28
    m_n2p = mu[6] # molar mass of N2+, 28
    m_nop = mu[7] # molar mass of NO+, 40
    gamma_i = mu[8] # ion adiabatic index, 5/3
    gamma_e = mu[9] # electron adiabatic index, 5/3
    lamda = mu[10] # kB*T0*t0**2/(mp*H0**2)
    omega = mu[11] # e*B0*t0/mp

    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    bgb = v[0] # b.grad(B)
    tn = v[1]
    u1 = v[2]
    u2 = v[3]
    u3 = v[4]
    B1 = v[5]
    B2 = v[6]
    B3 = v[7]
    E1 = v[8]
    E2 = v[9]
    E3 = v[10]
    prod_op = v[11]
    prod_o2p = v[12]
    prod_np = v[13]
    prod_n2p = v[14]
    prod_nop = v[15]
    loss_op = v[16]
    loss_o2p = v[17]
    loss_np = v[18]
    loss_n2p = v[19]
    loss_nop = v[20]
    nu_op = v[21]
    nu_o2p = v[22]
    nu_np = v[23]
    nu_n2p = v[24]
    nu_nop = v[25]
    kappa_i = v[26]
    kappa_e = v[27]
    heat_i = v[28]
    heat_e = v[29]
    qei = v[30]
    qin = v[31]
    qen = v[32]

    logop = u[0]
    logo2p = u[1]
    lognp = u[2]
    logn2p = u[3]
    lognop = u[4]
    v1 = u[5]
    v2 = u[6]
    v3 = u[7]
    ti = u[8]
    te = u[9]

    logopx = -q[0]
    logo2px = -q[1]
    lognpx = -q[2]
    logn2px = -q[3]
    lognopx = -q[4]
    v1x = -q[5]
    tix = -q[8]
    tex = -q[9]
    logopy = -q[10]
    logo2py = -q[11]
    lognpy = -q[12]
    logn2py = -q[13]
    lognopy = -q[14]
    v2y = -q[16]
    tiy = -q[18]
    tey = -q[19]
    logopz = -q[20]
    logo2pz = -q[21]
    lognpz = -q[22]
    logn2pz = -q[23]
    lognopz = -q[24]
    v3z = -q[27]
    tiz = -q[28]
    tez = -q[29]

    r = sqrt(x1**2 + x2**2 + x3**2)
    Bmag = sqrt(B1**2 + B2**2 + B3**2)
    b1 = B1/Bmag
    b2 = B2/Bmag
    b3 = B3/Bmag
    op = exp(logop)
    o2p = exp(logo2p)
    np = exp(lognp)
    n2p = exp(logn2p)
    nop = exp(lognop)
    ne = op + o2p + np + n2p + nop
    mbar = (op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop) / ne

    # divergence of velocity
    divV = v1x + v2y + v3z

    source_op = prod_op/op - loss_op + (logop-1)*divV
    source_o2p = prod_o2p/o2p - loss_o2p + (logo2p-1)*divV
    source_np = prod_np/np - loss_np + (lognp-1)*divV
    source_n2p = prod_n2p/n2p - loss_n2p + (logn2p-1)*divV
    source_nop = prod_nop/nop - loss_nop + (lognop-1)*divV

    # gravity + centrifugal + Coriolis
    g = g0*(r0/r)**2
    a1 = -g*x1/r + wrot**2*x1 + 2*wrot*v2
    a2 = -g*x2/r + wrot**2*x2 - 2*wrot*v1
    a3 = -g*x3/r

    # pressure gradient
    pg = lamda*(ti + te)/(ne*mbar)
    pg1 = -pg*(op*logopx + o2p*logo2px + np*lognpx + n2p*logn2px + nop*lognopx)
    pg2 = -pg*(op*logopy + o2p*logo2py + np*lognpy + n2p*logn2py + nop*lognopy)
    pg3 = -pg*(op*logopz + o2p*logo2pz + np*lognpz + n2p*logn2pz + nop*lognopz)

    # Lorentz force: E + vxB
    em1 = omega*(E1 + v2*B3 - v3*B2)/mbar
    em2 = omega*(E2 + v3*B1 - v1*B3)/mbar
    em3 = omega*(E3 + v1*B2 - v2*B1)/mbar

    # collision
    col = (op*m_op*nu_op
            + o2p*m_o2p*nu_o2p
            + np*m_np*nu_np
            + n2p*m_n2p*nu_n2p
            + nop*m_nop*nu_nop) / (ne*mbar)

    source_v1 = a1 + v1*divV + pg1 + em1 + col*(u1-v1)
    source_v2 = a2 + v2*divV + pg2 + em2 + col*(u2-v2)
    source_v3 = a3 + v3*divV + pg3 + em3 + col*(u3-v3)

    # chemical heating
    hchem = ti*(prod_op-loss_op*op
                + prod_o2p-loss_o2p*o2p
                + prod_np-loss_np*np
                + prod_n2p-loss_n2p*n2p
                + prod_nop-loss_nop*nop)

    # heat conduction
    hicond = kappa_i*bgb*ti**2*sqrt(ti)*(b1*tix + b2*tiy + b3*tiz)
    hecond = kappa_e*bgb*te**2*sqrt(te)*(b1*tex + b2*tey + b3*tez)

    # heat transfer
    hitrans = heat_i + qei*(te-ti) - qin*(ti-tn)
    hetrans = heat_e - qei*(te-ti) - qen*(te-tn)

    # mechanical work
    hiwork = -ne*ti*divV

    source_Ti = hchem + (gamma_i-1)*(hicond + hitrans + hiwork)
    source_Te = (gamma_e-1)*(hecond + hetrans)

    s = array([source_op, source_o2p, source_np, source_n2p, source_nop,
                source_v1, source_v2, source_v3, source_Ti, source_Te])
    return s

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = zeros(shape=(10, 2))
    return ub

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fb = zeros(shape=(10, 2))
    return fb

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    tn = v[1]
    u1 = v[2]
    u2 = v[3]
    u3 = v[4]
    prod_op = v[11]
    prod_o2p = v[12]
    prod_np = v[13]
    prod_n2p = v[14]
    prod_nop = v[15]
    loss_op = v[16]
    loss_o2p = v[17]
    loss_np = v[18]
    loss_n2p = v[19]
    loss_nop = v[20]
    heat_i = v[28]
    heat_e = v[29]
    qei = v[30]
    qin = v[31]
    qen = v[32]

    B1 = v[5]
    B2 = v[6]
    B3 = v[7]
    E1 = v[8]
    E2 = v[9]
    E3 = v[10]

    v1x = -q[0]
    v2x = -q[1]
    v3x = -q[2]
    v1y = -q[3]
    v2y = -q[4]
    v3y = -q[5]
    v1z = -q[6]
    v2z = -q[7]
    v3z = -q[8]

    Bmag2 = B1**2 + B2**2 + B3**2
    Bmag = sqrt(Bmag2)
    b1 = B1/Bmag
    b2 = B2/Bmag
    b3 = B3/Bmag

    ExB1 = (E2*B3 - E3*B2)/Bmag2
    ExB2 = (E3*B1 - E1*B3)/Bmag2
    ExB3 = (E1*B2 - E2*B1)/Bmag2

    # thermal equilibrium temperature
    ti_eq = (heat_i + qei*te + qin*tn)/(qei + qin)
    te_eq = (heat_e + qei*ti + qen*tn)/(qei + qen)

    # Lower boundary:
    # chemical equilibrium n = P/L
    # v = u
    # thermal equilibrium T = t_eq
    fl = array([prod_op/loss_op - uhat[0],
                prod_o2p/loss_o2p - uhat[1],
                prod_np/loss_np - uhat[2],
                prod_n2p/loss_n2p - uhat[3],
                prod_nop/loss_nop - uhat[4],
                u1 - uhat[5], u2 - uhat[6], u3 - uhat[7],
                ti_eq - uhat[8], te_eq - uhat[9]])

    # Upper boundary:
    # chemical equilibrium n = P/L
    # vperp = ExB/B^2 & d(vpar)/db = 0
    # thermal equilibrium T = t_eq
    vhatb = uhat[5]*b1 + uhat[6]*b2 + uhat[7]*b3
    fperp1 = uhat[5] - vhatb*b1 - ExB1
    fperp2 = uhat[6] - vhatb*b2 - ExB2
    fperp3 = uhat[7] - vhatb*b3 - ExB3

    bGb = (v1x*b1**2 + v2x*b1*b2 + v3x*b1*b3
            + v1y*b1*b2 + v2y*b2**2 + v3y*b2*b3
            + v1z*b1*b3 + v2z*b2*b3 + v3z*b3**2)
    vb = u[5]*b1 + u[6]*b2 + u[7]*b3
    fpar = bGb + tau[0]*(vb-vhatb)
    fpar1 = fpar*b1
    fpar2 = fpar*b2
    fpar3 = fpar*b3

    fu = array([prod_op/loss_op - uhat[0],
                prod_o2p/loss_o2p - uhat[1],
                prod_np/loss_np - uhat[2],
                prod_n2p/loss_n2p - uhat[3],
                prod_nop/loss_nop - uhat[4],
                fperp1 + fpar1, fperp2 + fpar2, fperp3 + fpar3,
                ti_eq - uhat[8], te_eq - uhat[9]])

    fb = reshape(hstack(tup=(fl, fu)), shape=(10, 2), order='F')
    return fb
