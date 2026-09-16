from numpy import ones, reshape, array, zeros, hstack
from sympy import sqrt, exp

def mass(u, q, w, v, x, t, mu, eta):
    m = ones(shape=10)
    return m

def flux(u, q, w, v, x, t, mu, eta):
    m_op = mu[3] # molar mass of O+, 16
    m_o2p = mu[4] # molar mass of O2+, 32
    m_np = mu[5] # molar mass of N+, 14
    m_n2p = mu[6] # molar mass of N2+, 28
    m_nop = mu[7] # molar mass of NO+, 30
    gamma_i = mu[8] # ion adiabatic index, 5/3
    gamma_e = mu[9] # electron adiabatic index, 5/3
    lamda = mu[10] # kB*T0*t0**2/(mp*H0**2)
    omega = mu[11] # e*B0*t0/mp

    E1 = v[15]
    E2 = v[16]
    E3 = v[17]
    kappa_i = v[23]
    kappa_e = v[24]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]
    Bmag = v[32]

    op = u[0]
    o2p = u[1]
    np = u[2]
    n2p = u[3]
    nop = u[4]
    phi1 = u[5]
    phi2 = u[6]
    phi3 = u[7]
    pi = u[8]
    pe = u[9]

    opx = -q[0]
    o2px = -q[1]
    npx = -q[2]
    n2px = -q[3]
    nopx = -q[4]
    pix = -q[8]
    pex = -q[9]
    opy = -q[10]
    o2py = -q[11]
    npy = -q[12]
    n2py = -q[13]
    nopy = -q[14]
    piy = -q[18]
    pey = -q[19]
    opz = -q[20]
    o2pz = -q[21]
    npz = -q[22]
    n2pz = -q[23]
    nopz = -q[24]
    piz = -q[28]
    pez = -q[29]

    ne = op + o2p + np + n2p + nop
    nm = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    nx = opx + o2px + npx + n2px + nopx
    ny = opy + o2py + npy + n2py + nopy
    nz = opz + o2pz + npz + n2pz + nopz
    ti = pi/ne
    te = pe/ne
    tix = (pix - nx*ti)/ne
    tiy = (piy - ny*ti)/ne
    tiz = (piz - nz*ti)/ne
    tex = (pex - nx*te)/ne
    tey = (pey - ny*te)/ne
    tez = (pez - nz*te)/ne

    fphi11 = phi1**2/nm + lamda*(pi + pe*b1**2)
    fphi12 = phi1*phi2/nm + lamda*pe*b1*b2
    fphi13 = phi1*phi3/nm + lamda*pe*b1*b3
    fphi21 = phi2*phi1/nm + lamda*pe*b2*b1
    fphi22 = phi2**2/nm + lamda*(pi + pe*b2**2)
    fphi23 = phi2*phi3/nm + lamda*pe*b2*b3
    fphi31 = phi3*phi1/nm + lamda*pe*b3*b1
    fphi32 = phi3*phi2/nm + lamda*pe*b3*b2
    fphi33 = phi3**2/nm + lamda*(pi + pe*b3**2)

    bgti = b1*tix + b2*tiy + b3*tiz
    ci = (gamma_i-1)*kappa_i*ti**2.5*bgti
    fpi1 = pi*phi1/nm - ci*b1
    fpi2 = pi*phi2/nm - ci*b2
    fpi3 = pi*phi3/nm - ci*b3

    Ep1 = E1 + lamda*pex/(omega*ne)
    Ep2 = E2 + lamda*pey/(omega*ne)
    Ep3 = E3 + lamda*pez/(omega*ne)
    ve1 = (Ep2*b3 - Ep3*b2)/Bmag
    ve2 = (Ep3*b1 - Ep1*b3)/Bmag
    ve3 = (Ep1*b2 - Ep2*b1)/Bmag
    bgte = b1*tex + b2*tey + b3*tez
    ce = (gamma_e-1)*kappa_e*te**2.5*bgte
    fpe1 = pe*ve1 - ce*b1
    fpe2 = pe*ve2 - ce*b2
    fpe3 = pe*ve3 - ce*b3

    f = reshape([op*phi1/nm, o2p*phi1/nm, np*phi1/nm, n2p*phi1/nm, nop*phi1/nm, fphi11, fphi12, fphi13, fpi1, fpe1,
                 op*phi2/nm, o2p*phi2/nm, np*phi2/nm, n2p*phi2/nm, nop*phi2/nm, fphi21, fphi22, fphi23, fpi2, fpe2,
                 op*phi3/nm, o2p*phi3/nm, np*phi3/nm, n2p*phi3/nm, nop*phi3/nm, fphi31, fphi32, fphi33, fpi3, fpe3],
                shape=(10, 3), order='F')
    return f

def source(u, q, w, v, x, t, mu, eta):
    r0 = mu[0] # earth radius (non-dimensionalized): Re/H0
    g0 = mu[1] # gravity at earth surface (non-dimensionalized): g0*t0**2/H0
    wrot = mu[2] # earth rotation (non-dimensionalized): wrot*t0
    m_op = mu[3] # molar mass of O+, 16
    m_o2p = mu[4] # molar mass of O2+, 32
    m_np = mu[5] # molar mass of N+, 14
    m_n2p = mu[6] # molar mass of N2+, 28
    m_nop = mu[7] # molar mass of NO+, 30
    gamma_i = mu[8] # ion adiabatic index, 5/3
    gamma_e = mu[9] # electron adiabatic index, 5/3
    lamda = mu[10] # kB*T0*t0**2/(mp*H0**2)
    omega = mu[11] # e*B0*t0/mp
    dampfac = mu[12] # damping coefficient (non-dimensionalized): dampfac*t0

    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    prod_op = v[0]
    prod_o2p = v[1]
    prod_np = v[2]
    prod_n2p = v[3]
    prod_nop = v[4]
    loss_op = v[5]
    loss_o2p = v[6]
    loss_np = v[7]
    loss_n2p = v[8]
    loss_nop = v[9]
    nu_op = v[10]
    nu_o2p = v[11]
    nu_np = v[12]
    nu_n2p = v[13]
    nu_nop = v[14]
    E1 = v[15]
    E2 = v[16]
    E3 = v[17]
    heat_i = v[18]
    heat_e = v[19]
    qei = v[20]
    qin = v[21]
    qen = v[22]
    kappa_i = v[23]
    kappa_e = v[24]
    tn = v[25]
    u1 = v[26]
    u2 = v[27]
    u3 = v[28]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]
    Bmag = v[32]
    lnBx = v[33]
    lnBy = v[34]
    lnBz = v[35]
    bglnB = v[36] # b.grad(lnB)

    op = u[0]
    o2p = u[1]
    np = u[2]
    n2p = u[3]
    nop = u[4]
    phi1 = u[5]
    phi2 = u[6]
    phi3 = u[7]
    pi = u[8]
    pe = u[9]

    opx = -q[0]
    o2px = -q[1]
    npx = -q[2]
    n2px = -q[3]
    nopx = -q[4]
    phi1x = -q[5]
    pix = -q[8]
    pex = -q[9]
    opy = -q[10]
    o2py = -q[11]
    npy = -q[12]
    n2py = -q[13]
    nopy = -q[14]
    phi2y = -q[16]
    piy = -q[18]
    pey = -q[19]
    opz = -q[20]
    o2pz = -q[21]
    npz = -q[22]
    n2pz = -q[23]
    nopz = -q[24]
    phi3z = -q[27]
    piz = -q[28]
    pez = -q[29]

    r = sqrt(x1**2 + x2**2 + x3**2)
    ne = op + o2p + np + n2p + nop
    nm = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    nx = opx + o2px + npx + n2px + nopx
    ny = opy + o2py + npy + n2py + nopy
    nz = opz + o2pz + npz + n2pz + nopz
    vi1 = phi1/nm
    vi2 = phi2/nm
    vi3 = phi3/nm
    ti = pi/ne
    te = pe/ne
    tix = (pix - nx*ti)/ne
    tiy = (piy - ny*ti)/ne
    tiz = (piz - nz*ti)/ne
    tex = (pex - nx*te)/ne
    tey = (pey - ny*te)/ne
    tez = (pez - nz*te)/ne

    net_op = prod_op - loss_op*op
    net_o2p = prod_o2p - loss_o2p*o2p
    net_np = prod_np - loss_np*np
    net_n2p = prod_n2p - loss_n2p*n2p
    net_nop = prod_nop - loss_nop*nop

    # gravity + centrifugal
    g = g0*(r0/r)**2
    a1 = -g*x1/r + wrot**2*x1
    a2 = -g*x2/r + wrot**2*x2
    a3 = -g*x3/r

    # along magnetic fields: damping + pressure gradient (not included for now)
    alongb = -dampfac*(phi1*b1 + phi2*b2 + phi3*b3) #- 2*lamda*pe*bglnB

    # chemical
    chem = m_op*net_op + m_o2p*net_o2p + m_np*net_np + m_n2p*net_n2p + m_nop*net_nop

    # Lorentz
    em1 = nm*E1 + Bmag*(phi2*b3 - phi3*b2)
    em2 = nm*E2 + Bmag*(phi3*b1 - phi1*b3)
    em3 = nm*E3 + Bmag*(phi1*b2 - phi2*b1)

    # collision
    col = op*m_op*nu_op + o2p*m_o2p*nu_o2p + np*m_np*nu_np + n2p*m_n2p*nu_n2p + nop*m_nop*nu_nop

    source_phi1 = nm*a1 + alongb*b1 + (chem*phi1 + omega*em1 + col*(nm*u1-phi1))/nm + 2*wrot*phi2
    source_phi2 = nm*a2 + alongb*b2 + (chem*phi2 + omega*em2 + col*(nm*u2-phi2))/nm - 2*wrot*phi1
    source_phi3 = nm*a3 + alongb*b3 + (chem*phi3 + omega*em3 + col*(nm*u3-phi3))/nm

    # chemical heating
    hchem = net_op + net_o2p + net_np + net_n2p + net_nop

    # mechanical work
    nmx = opx*m_op + o2px*m_o2p + npx*m_np + n2px*m_n2p + nopx*m_nop
    nmy = opy*m_op + o2py*m_o2p + npy*m_np + n2py*m_n2p + nopy*m_nop
    nmz = opz*m_op + o2pz*m_o2p + npz*m_np + n2pz*m_n2p + nopz*m_nop
    divvi = (phi1x + phi2y + phi3z - nmx*vi1 - nmy*vi2 - nmz*vi3)/nm
    hiwork = -pi*divvi
    Ep1 = E1 + lamda*pex/(omega*ne)
    Ep2 = E2 + lamda*pey/(omega*ne)
    Ep3 = E3 + lamda*pez/(omega*ne)
    ve1 = (Ep2*b3 - Ep3*b2)/Bmag
    ve2 = (Ep3*b1 - Ep1*b3)/Bmag
    ve3 = (Ep1*b2 - Ep2*b1)/Bmag
    det = nx*pey*b3 + ny*pez*b1 + nz*pex*b2 - nx*pez*b2 - ny*pex*b3 - nz*pey*b1
    divve = -lamda*det/(omega*ne**2*Bmag) - 2*(ve1*lnBx + ve2*lnBy + ve3*lnBz)
    hework = -pe*divve

    # heat conduction (not included for now)
    bgti = b1*tix + b2*tiy + b3*tiz
    hicond = kappa_i*ti**2.5*bglnB*bgti
    bgte = b1*tex + b2*tey + b3*tez
    hecond = kappa_e*te**2.5*bglnB*bgte

    # heat transfer
    hitrans = (ne*heat_i + qei*(pe-pi) - qin*(pi-ne*tn))/ne
    hetrans = (ne*heat_e - qei*(pe-pi) - qen*(pe-ne*tn))/ne

    source_pi = pi*hchem/ne + (gamma_i-1)*(hiwork + hitrans)
    source_pe = pe*hchem/ne + (gamma_e-1)*(hework + hetrans)

    s = array([net_op, net_o2p, net_np, net_n2p, net_nop,
               source_phi1, source_phi2, source_phi3, source_pi, source_pe])
    return s

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = zeros(shape=(10, 2))
    return ub

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fb = zeros(shape=(10, 2))
    return fb

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    m_op = mu[3] # molar mass of O+, 16
    m_o2p = mu[4] # molar mass of O2+, 32
    m_np = mu[5] # molar mass of N+, 14
    m_n2p = mu[6] # molar mass of N2+, 28
    m_nop = mu[7] # molar mass of NO+, 30
    psi_op = mu[13] # O+ number flux at upper boundary
    psi_o2p = mu[14] # O2+ number flux at upper boundary
    psi_np = mu[15] # N+ number flux at upper boundary
    psi_n2p = mu[16] # N2+ number flux at upper boundary
    psi_nop = mu[17] # NO+ number flux at upper boundary
    qi = mu[18] # ion heat flux at upper boundary
    qe = mu[19] # electron heat flux at upper boundary

    prod_op = v[0]
    prod_o2p = v[1]
    prod_np = v[2]
    prod_n2p = v[3]
    prod_nop = v[4]
    loss_op = v[5]
    loss_o2p = v[6]
    loss_np = v[7]
    loss_n2p = v[8]
    loss_nop = v[9]
    E1 = v[15]
    E2 = v[16]
    E3 = v[17]
    heat_i = v[18]
    heat_e = v[19]
    qei = v[20]
    qin = v[21]
    qen = v[22]
    kappa_i = v[23]
    kappa_e = v[24]
    tn = v[25]
    u1 = v[26]
    u2 = v[27]
    u3 = v[28]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]
    Bmag = v[32]
    b1x = v[37]
    b1y = v[38]
    b1z = v[39]
    b2x = v[40]
    b2y = v[41]
    b2z = v[42]
    b3x = v[43]
    b3y = v[44]
    b3z = v[45]

    op = u[0]
    o2p = u[1]
    np = u[2]
    n2p = u[3]
    nop = u[4]
    phi1 = u[5]
    phi2 = u[6]
    phi3 = u[7]
    pi = u[8]
    pe = u[9]

    opx = -q[0]
    o2px = -q[1]
    npx = -q[2]
    n2px = -q[3]
    nopx = -q[4]
    phi1x = -q[5]
    phi2x = -q[6]
    phi3x = -q[7]
    pix = -q[8]
    pex = -q[9]
    opy = -q[10]
    o2py = -q[11]
    npy = -q[12]
    n2py = -q[13]
    nopy = -q[14]
    phi1y = -q[15]
    phi2y = -q[16]
    phi3y = -q[17]
    piy = -q[18]
    pey = -q[19]
    opz = -q[20]
    o2pz = -q[21]
    npz = -q[22]
    n2pz = -q[23]
    nopz = -q[24]
    phi1z = -q[25]
    phi2z = -q[26]
    phi3z = -q[27]
    piz = -q[28]
    pez = -q[29]

    ophat = uhat[0]
    o2phat = uhat[1]
    nphat = uhat[2]
    n2phat = uhat[3]
    nophat = uhat[4]
    phihat1 = uhat[5]
    phihat2 = uhat[6]
    phihat3 = uhat[7]
    pihat = uhat[8]
    pehat = uhat[9]

    ne = op + o2p + np + n2p + nop
    nm = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    nx = opx + o2px + npx + n2px + nopx
    ny = opy + o2py + npy + n2py + nopy
    nz = opz + o2pz + npz + n2pz + nopz
    nmx = opx*m_op + o2px*m_o2p + npx*m_np + n2px*m_n2p + nopx*m_nop
    nmy = opy*m_op + o2py*m_o2p + npy*m_np + n2py*m_n2p + nopy*m_nop
    nmz = opz*m_op + o2pz*m_o2p + npz*m_np + n2pz*m_n2p + nopz*m_nop

    # chemical equilibrium
    op_eq = prod_op/loss_op
    o2p_eq = prod_o2p/loss_o2p
    np_eq = prod_np/loss_np
    n2p_eq = prod_n2p/loss_n2p
    nop_eq = prod_nop/loss_nop

    nm_eq = op_eq*m_op + o2p_eq*m_o2p + np_eq*m_np + n2p_eq*m_n2p + nop_eq*m_nop

    # thermal equilibrium temperature
    ne_eq = op_eq + o2p_eq + np_eq + n2p_eq + nop_eq
    denom = qei*qen + qei*qin + qen*qin
    pi_eq = ne_eq*((qei+qen)*(heat_i + qin*tn) + qei*(heat_e + qen*tn))/denom
    pe_eq = ne_eq*((qei+qin)*(heat_e + qen*tn) + qei*(heat_i + qin*tn))/denom

    # Lower boundary: chemical equilibrium, v = u, thermal equilibrium
    fl = array([ophat - op_eq,
                o2phat - o2p_eq,
                nphat - np_eq,
                n2phat - n2p_eq,
                nophat - nop_eq,
                phihat1 - nm_eq*u1,
                phihat2 - nm_eq*u2,
                phihat3 - nm_eq*u3,
                pihat - pi_eq,
                pehat - pe_eq])

    # Upper boundary:
    # number flux nv = psi0*b
    # vperp = ExB/B^2 & d(vpar)/ds = 0
    # heat flux gradT = q0*b
    bn = b1*n[0] + b2*n[1] + b3*n[2]
    phin = phi1*n[0] + phi2*n[1] + phi3*n[2]
    phihatn = phihat1*n[0] + phihat2*n[1] + phihat3*n[2]

    vd1 = (E2*b3 - E3*b2)/Bmag
    vd2 = (E3*b1 - E1*b3)/Bmag
    vd3 = (E1*b2 - E2*b1)/Bmag
    phihatb = phihat1*b1 + phihat2*b2 + phihat3*b3
    fperp1 = phihat1 - phihatb*b1 - nm*vd1
    fperp2 = phihat2 - phihatb*b2 - nm*vd2
    fperp3 = phihat3 - phihatb*b3 - nm*vd3

    bgphib = (phi1x*b1**2 + phi2x*b1*b2 + phi3x*b1*b3
            + phi1y*b2*b1 + phi2y*b2**2 + phi3y*b2*b3
            + phi1z*b3*b1 + phi2z*b3*b2 + phi3z*b3**2)
    bgbphi = (b1x*b1*phi1 + b2x*b1*phi2 + b3x*b1*phi3
            + b1y*b2*phi1 + b2y*b2*phi2 + b3y*b2*phi3
            + b1z*b3*phi1 + b2z*b3*phi2 + b3z*b3*phi3) # not included for now
    bgn = nmx*b1 + nmy*b2 + nmz*b3
    phib = phi1*b1 + phi2*b2 + phi3*b3
    fpar = bgphib - bgn*phib/nm + tau[5]*(phib-phihatb)

    ti = pi/ne
    te = pe/ne
    gradpin = pix*n[0] + piy*n[1] + piz*n[2]
    gradpen = pex*n[0] + pey*n[1] + pez*n[2]
    gradnn = nx*n[0] + ny*n[1] + nz*n[2]

    fu = array([op*phin + tau[0]*(op*phin - ophat*phihatn) - nm*psi_op*bn,
                o2p*phin + tau[1]*(o2p*phin - o2phat*phihatn) - nm*psi_o2p*bn,
                np*phin + tau[2]*(np*phin - nphat*phihatn) - nm*psi_np*bn,
                n2p*phin + tau[3]*(n2p*phin - n2phat*phihatn) - nm*psi_n2p*bn,
                nop*phin + tau[4]*(nop*phin - nophat*phihatn) - nm*psi_nop*bn,
                fperp1 + fpar*b1, fperp2 + fpar*b2, fperp3 + fpar*b3,
                gradpin-gradnn*ti + tau[6]*(pi-pihat) + ne*qi*bn/(kappa_i*ti**2.5),
                gradpen-gradnn*te + tau[7]*(pe-pehat) + ne*qe*bn/(kappa_e*te**2.5)])

    fb = reshape(hstack(tup=(fl, fu)), shape=(10, 2), order='F')
    return fb

def initu(x, mu, eta):
    u0 = zeros(shape=10)
    return u0
