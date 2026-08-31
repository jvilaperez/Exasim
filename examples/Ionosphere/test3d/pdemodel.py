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
    omega = mu[11] # e*B0*t0/mp

    B1 = v[4]
    B2 = v[5]
    B3 = v[6]
    Bmag = v[7]
    E1 = v[12]
    E2 = v[13]
    E3 = v[14]
    kappa_i = v[30] # ion thermal conduction coefficient
    kappa_e = v[31] # electron thermal conduction coefficient

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

    Bmag = sqrt(B1**2 + B2**2 + B3**2)
    b1 = B1/Bmag
    b2 = B2/Bmag
    b3 = B3/Bmag
    ne = op + o2p + np + n2p + nop
    nm = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    nx = opx + o2px + npx + n2px + nopx
    ny = opy + o2py + npy + n2py + nopy
    nz = opz + o2pz + npz + n2pz + nopz

    fphi11 = phi1**2/nm + lamda*(pi + pe*b1**2)
    fphi12 = phi1*phi2/nm + lamda*pe*b1*b2
    fphi13 = phi1*phi3/nm + lamda*pe*b1*b3
    fphi21 = phi2*phi1/nm + lamda*pe*b2*b1
    fphi22 = phi2**2/nm + lamda*(pi + pe*b2**2)
    fphi23 = phi2*phi3/nm + lamda*pe*b2*b3
    fphi31 = phi3*phi1/nm + lamda*pe*b3*b1
    fphi32 = phi3*phi2/nm + lamda*pe*b3*b2
    fphi33 = phi3**2/nm + lamda*(pi + pe*b3**2)

    qix = ne*pix - pi*nx
    qiy = ne*piy - pi*ny
    qiz = ne*piz - pi*nz
    bqi = b1*qix + b2*qiy + b3*qiz
    fpi_1 = pi*phi1/nm - (gamma_i-1)*kappa_i*b1*bqi/ne**2
    fpi_2 = pi*phi2/nm - (gamma_i-1)*kappa_i*b2*bqi/ne**2
    fpi_3 = pi*phi3/nm - (gamma_i-1)*kappa_i*b3*bqi/ne**2

    vex = (E2*b3 - E3*b2)/Bmag + lamda*(pey*b3 - pez*b2)/(omega*ne*Bmag)
    vey = (E3*b1 - E1*b3)/Bmag + lamda*(pez*b1 - pex*b3)/(omega*ne*Bmag)
    vez = (E1*b2 - E2*b1)/Bmag + lamda*(pex*b2 - pey*b1)/(omega*ne*Bmag)
    qex = ne*pex - pe*nx
    qey = ne*pey - pe*ny
    qez = ne*pez - pe*nz
    bqe = b1*qex + b2*qey + b3*qez
    fpe_1 = pe*vex - (gamma_e-1)*kappa_e*b1*bqe/ne**2
    fpe_2 = pe*vey - (gamma_e-1)*kappa_e*b2*bqe/ne**2
    fpe_3 = pe*vez - (gamma_e-1)*kappa_e*b3*bqe/ne**2

    f = reshape([op*phi1/nm, o2p*phi1/nm, np*phi1/nm, n2p*phi1/nm, nop*phi1/nm, fphi11, fphi12, fphi13, fpi_1, fpe_1,
                 op*phi2/nm, o2p*phi2/nm, np*phi2/nm, n2p*phi2/nm, nop*phi2/nm, fphi21, fphi22, fphi23, fpi_2, fpe_2,
                 op*phi3/nm, o2p*phi3/nm, np*phi3/nm, n2p*phi3/nm, nop*phi3/nm, fphi31, fphi32, fphi33, fpi_3, fpe_3],
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
    dampfac = mu[12] # damping coefficient (non-dimensionalized): dampfac*t0

    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    tn = v[0]
    u1 = v[1]
    u2 = v[2]
    u3 = v[3]
    B1 = v[4]
    B2 = v[5]
    B3 = v[6]
    Bmag = v[7]
    Bx = v[8]
    By = v[9]
    Bz = v[10]
    bgradB = v[11] # b.grad(Bmag)
    E1 = v[12]
    E2 = v[13]
    E3 = v[14]
    prod_op = v[15]
    prod_o2p = v[16]
    prod_np = v[17]
    prod_n2p = v[18]
    prod_nop = v[19]
    loss_op = v[20]
    loss_o2p = v[21]
    loss_np = v[22]
    loss_n2p = v[23]
    loss_nop = v[24]
    nu_op = v[25]
    nu_o2p = v[26]
    nu_np = v[27]
    nu_n2p = v[28]
    nu_nop = v[29]
    kappa_i = v[30]
    kappa_e = v[31]
    heat_i = v[32]
    heat_e = v[33]
    qei = v[34]
    qin = v[35]
    qen = v[36]

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
    b1 = B1/Bmag
    b2 = B2/Bmag
    b3 = B3/Bmag
    ne = op + o2p + np + n2p + nop
    nm = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    nx = opx + o2px + npx + n2px + nopx
    ny = opy + o2py + npy + n2py + nopy
    nz = opz + o2pz + npz + n2pz + nopz

    net_op = prod_op - loss_op*op
    net_o2p = prod_o2p - loss_o2p*o2p
    net_np = prod_np - loss_np*np
    net_n2p = prod_n2p - loss_n2p*n2p
    net_nop = prod_nop - loss_nop*nop

    # gravity + centrifugal
    g = g0*(r0/r)**2
    a1 = nm*(-g*x1/r + wrot**2*x1)
    a2 = nm*(-g*x2/r + wrot**2*x2)
    a3 = nm*(-g*x3/r)

    # Coriolis
    cor1 = 2*wrot*phi2
    cor2 = -2*wrot*phi1

    # damping
    phib = -dampfac*(phi1*b1 + phi2*b2 + phi3*b3)

    # pressure gradient
    pg = -2*lamda*pe*bgradB

    # chemical
    chem = m_op*net_op + m_o2p*net_o2p + m_np*net_np + m_n2p*net_n2p + m_nop*net_nop

    # Lorentz
    em1 = omega*(nm*E1 + phi2*B3 - phi3*B2)
    em2 = omega*(nm*E2 + phi3*B1 - phi1*B3)
    em3 = omega*(nm*E3 + phi1*B2 - phi2*B1)

    # collision
    col = (op*m_op*nu_op
            + o2p*m_o2p*nu_o2p
            + np*m_np*nu_np
            + n2p*m_n2p*nu_n2p
            + nop*m_nop*nu_nop)

    source_phi1 = a1 + cor1 + phib*b1 + pg*b1 + (chem*phi1 + em1 + col*(nm*u1-phi1))/nm
    source_phi2 = a2 + cor2 + phib*b1 + pg*b2 + (chem*phi2 + em2 + col*(nm*u2-phi2))/nm
    source_phi3 = a3 + phib*b1 + pg*b3 + (chem*phi3 + em3 + col*(nm*u3-phi3))/nm

    # chemical heating
    hchem = net_op + net_o2p + net_np + net_n2p + net_nop

    # heat conduction
    qix = ne*pix - pi*nx
    qiy = ne*piy - pi*ny
    qiz = ne*piz - pi*nz
    bqi = b1*qix + b2*qiy + b3*qiz
    hicond = kappa_i*bgradB*bqi/ne**2
    qex = ne*pex - pe*nx
    qey = ne*pey - pe*ny
    qez = ne*pez - pe*nz
    bqe = b1*qex + b2*qey + b3*qez
    hecond = kappa_e*bgradB*bqe/ne**2

    # heat transfer
    hitrans = (ne*heat_i + qei*(pe-pi) - qin*(pi-ne*tn))/ne
    hetrans = (ne*heat_e - qei*(pe-pi) - qen*(pe-ne*tn))/ne

    # mechanical work
    nmx = opx*m_op + o2px*m_o2p + npx*m_np + n2px*m_n2p + nopx*m_nop
    nmy = opy*m_op + o2py*m_o2p + npy*m_np + n2py*m_n2p + nopy*m_nop
    nmz = opz*m_op + o2pz*m_o2p + npz*m_np + n2pz*m_n2p + nopz*m_nop
    hiwork = -pi*(nm*(phi1x + phi2y + phi3z) - (phi1*nmx + phi2*nmy + phi3*nmz))/nm**2
    vex = (E2*b3 - E3*b2)/Bmag + lamda*(pey*b3 - pez*b2)/(omega*ne*Bmag)
    vey = (E3*b1 - E1*b3)/Bmag + lamda*(pez*b1 - pex*b3)/(omega*ne*Bmag)
    vez = (E1*b2 - E2*b1)/Bmag + lamda*(pex*b2 - pey*b1)/(omega*ne*Bmag)
    det = nx*pey*b3 + ny*pez*b1 + nz*pex*b2 - nx*pez*b2 - ny*pex*b3 - nz*pey*b1
    hework = pe*(lamda*det/(omega*ne**2*Bmag) + 2*(vex*Bx + vey*By + vez*Bz))

    source_pi = pi*hchem/ne + (gamma_i-1)*(hicond + hitrans + hiwork)
    source_pe = pe*hchem/ne + (gamma_e-1)*(hecond + hetrans + hework)

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
    m_np = mu[5] # molar mass of N+, 28
    m_n2p = mu[6] # molar mass of N2+, 28
    m_nop = mu[7] # molar mass of NO+, 40

    tn = v[0]
    u1 = v[1]
    u2 = v[2]
    u3 = v[3]
    B1 = v[4]
    B2 = v[5]
    B3 = v[6]
    Bmag = v[7]
    E1 = v[12]
    E2 = v[13]
    E3 = v[14]
    prod_op = v[15]
    prod_o2p = v[16]
    prod_np = v[17]
    prod_n2p = v[18]
    prod_nop = v[19]
    loss_op = v[20]
    loss_o2p = v[21]
    loss_np = v[22]
    loss_n2p = v[23]
    loss_nop = v[24]
    heat_i = v[32]
    heat_e = v[33]
    qei = v[34]
    qin = v[35]
    qen = v[36]

    op = u[0]
    o2p = u[1]
    np = u[2]
    n2p = u[3]
    nop = u[4]
    phi1 = u[5]
    phi2 = u[6]
    phi3 = u[7]

    phi1x = -q[5]
    phi2x = -q[6]
    phi3x = -q[7]
    phi1y = -q[15]
    phi2y = -q[16]
    phi3y = -q[17]
    phi1z = -q[25]
    phi2z = -q[26]
    phi3z = -q[27]

    b1 = B1/Bmag
    b2 = B2/Bmag
    b3 = B3/Bmag
    ne = op + o2p + np + n2p + nop
    nm = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop

    vd1 = (E2*b3 - E3*b2)/Bmag
    vd2 = (E3*b1 - E1*b3)/Bmag
    vd3 = (E1*b2 - E2*b1)/Bmag

    v1 = phi1/nm
    v2 = phi2/nm
    v3 = phi3/nm
    v1x = phi1x/nm
    v2x = phi2x/nm
    v3x = phi3x/nm
    v1y = phi1y/nm
    v2y = phi2y/nm
    v3y = phi3y/nm
    v1z = phi1z/nm
    v2z = phi2z/nm
    v3z = phi3z/nm

    # thermal equilibrium temperature
    ti_eq = ((qei+qen)*(heat_i + qin*tn) + qei*(heat_e + qen*tn))/(qei*qen + qei*qin + qen*qin)
    te_eq = ((qei+qin)*(heat_e + qen*tn) + qei*(heat_i + qin*tn))/(qei*qen + qei*qin + qen*qin)

    # Lower boundary:
    # chemical equilibrium n = P/L
    # v = u
    # thermal equilibrium T = t_eq
    fl = array([prod_op/loss_op - uhat[0],
                prod_o2p/loss_o2p - uhat[1],
                prod_np/loss_np - uhat[2],
                prod_n2p/loss_n2p - uhat[3],
                prod_nop/loss_nop - uhat[4],
                nm*u1 - uhat[5], nm*u2 - uhat[6], nm*u3 - uhat[7],
                ne*ti_eq - uhat[8], ne*te_eq - uhat[9]])

    # Upper boundary:
    # chemical equilibrium n = P/L
    # vperp = ExB/B^2 & d(vpar)/ds = 0
    # thermal equilibrium T = t_eq
    vhatb = uhat[5]*b1 + uhat[6]*b2 + uhat[7]*b3
    fperp1 = uhat[5] - vhatb*b1 - vd1
    fperp2 = uhat[6] - vhatb*b2 - vd2
    fperp3 = uhat[7] - vhatb*b3 - vd3

    bGb = (v1x*b1**2 + v2x*b1*b2 + v3x*b1*b3
            + v1y*b1*b2 + v2y*b2**2 + v3y*b2*b3
            + v1z*b1*b3 + v2z*b2*b3 + v3z*b3**2)
    vb = v1*b1 + v2*b2 + v3*b3
    fpar = bGb + tau[0]*(vb-vhatb)
    fpar1 = fpar*b1
    fpar2 = fpar*b2
    fpar3 = fpar*b3

    fu = array([prod_op/loss_op - uhat[0],
                prod_o2p/loss_o2p - uhat[1],
                prod_np/loss_np - uhat[2],
                prod_n2p/loss_n2p - uhat[3],
                prod_nop/loss_nop - uhat[4],
                nm*(fperp1 + fpar1), nm*(fperp2 + fpar2), nm*(fperp3 + fpar3),
                ne*ti_eq - uhat[8], ne*te_eq - uhat[9]])

    fb = reshape(hstack(tup=(fl, fu)), shape=(10, 2), order='F')
    return fb

def initu(x, mu, eta):
    u0 = zeros(shape=10)
    return u0
