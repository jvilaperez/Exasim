from numpy import ones, reshape, array, zeros, hstack
from sympy import sqrt, exp

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    return m

def flux(u, q, w, v, x, t, mu, eta):
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

    Pe = mu[4] # Peclet number (inverse)
    m_op  = mu[5] # molar mass of O+, 16
    m_o2p = mu[6] # molar mass of O2+, 32
    m_np  = mu[7] # molar mass of N+, 14
    m_n2p = mu[8] # molar mass of N2+, 28
    m_nop = mu[9] # molar mass of NO+, 30
    gamma_i = mu[0] # ion adiabatic index, 5/3
    gamma_e = mu[0] # electron adiabatic index, 5/3

    ne = op + o2p + np + n2p + nop
    rho = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    Ti = pi/ne
    Te = pe/ne
    ne  = smoothpos(ne, 1e-6)
    rho = smoothpos(rho, 1e-4)
    Ti  = smoothpos(Ti, 1e-3)
    Te  = smoothpos(Te, 1e-3)

    vi1 = phi1/rho
    vi2 = phi2/rho
    vi3 = phi3/rho

    nx = opx + o2px + npx + n2px + nopx
    ny = opy + o2py + npy + n2py + nopy
    nz = opz + o2pz + npz + n2pz + nopz
    Tix = (pix - nx*Ti)/ne
    Tiy = (piy - ny*Ti)/ne
    Tiz = (piz - nz*Ti)/ne
    Tex = (pex - nx*Te)/ne
    Tey = (pey - ny*Te)/ne
    Tez = (pez - nz*Te)/ne

    kappa_i0 = v[23]
    kappa_e0 = v[24]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]

    kappa_i = kappa_i0*Ti**2.5
    kappa_e = kappa_e0*Te**2.5
    bgti = b1*Tix + b2*Tiy + b3*Tiz
    bgte = b1*Tex + b2*Tey + b3*Tez

    hfi1 = (gamma_i-1)*Pe*kappa_i*bgti*b1
    hfi2 = (gamma_i-1)*Pe*kappa_i*bgti*b2
    hfi3 = (gamma_i-1)*Pe*kappa_i*bgti*b3
    hfe1 = (gamma_e-1)*Pe*kappa_e*bgte*b1
    hfe2 = (gamma_e-1)*Pe*kappa_e*bgte*b2
    hfe3 = (gamma_e-1)*Pe*kappa_e*bgte*b3

    f = array([op*vi1, o2p*vi1, np*vi1, n2p*vi1, nop*vi1, phi1*vi1 + (pi+pe), phi2*vi1, phi3*vi1, pi*vi1 - hfi1, pe*vi1 - hfe1,
                 op*vi2, o2p*vi2, np*vi2, n2p*vi2, nop*vi2, phi1*vi2, phi2*vi2 + (pi+pe), phi3*vi2, pi*vi2 - hfi2, pe*vi2 - hfe2,
                 op*vi3, o2p*vi3, np*vi3, n2p*vi3, nop*vi3, phi1*vi3, phi2*vi3, phi3*vi3 + (pi+pe), pi*vi3 - hfi3, pe*vi3 - hfe3])
    f = reshape(f,(10,3),'F')
    return f

def source(u, q, w, v, x, t, mu, eta):
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

    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    Ro = mu[1] # Rossby number (inverse)
    Kn = mu[2] # Knudsen number (inverse)
    Omega = mu[3] # ratio between gyro frequency and reference frequency
    Pe = mu[4] # Peclet number (inverse)
    m_op  = mu[5] # molar mass of O+, 16
    m_o2p = mu[6] # molar mass of O2+, 32
    m_np  = mu[7] # molar mass of N+, 14
    m_n2p = mu[8] # molar mass of N2+, 28
    m_nop = mu[9] # molar mass of NO+, 30
    gamma_i = mu[0] # ion adiabatic index, 5/3
    gamma_e = mu[0] # electron adiabatic index, 5/3

    ne = op + o2p + np + n2p + nop
    rho = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    Ti = pi/ne
    Te = pe/ne
    ne  = smoothpos(ne, 1e-6)
    rho = smoothpos(rho, 1e-4)
    Ti  = smoothpos(Ti, 1e-3)
    Te  = smoothpos(Te, 1e-3)

    vi1 = phi1/rho
    vi2 = phi2/rho
    vi3 = phi3/rho
    rhox = opx*m_op + o2px*m_o2p + npx*m_np + n2px*m_n2p + nopx*m_nop
    rhoy = opy*m_op + o2py*m_o2p + npy*m_np + n2py*m_n2p + nopy*m_nop
    rhoz = opz*m_op + o2pz*m_o2p + npz*m_np + n2pz*m_n2p + nopz*m_nop
    v1x = (phi1x-vi1*rhox)/rho
    v2y = (phi2y-vi2*rhoy)/rho
    v3z = (phi3z-vi3*rhoz)/rho
    divV = v1x + v2y + v3z

    # Input fields
    nu_op  = v[10]
    nu_o2p = v[11]
    nu_np  = v[12]
    nu_n2p = v[13]
    nu_nop = v[14]
    E1 = v[15]
    E2 = v[16]
    E3 = v[17]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]
    Bmag = v[32]
    u1 = v[26]
    u2 = v[27]
    u3 = v[28]

    # Accelerations
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    r = sqrt(x1**2 + x2**2 + x3**2)
    R0 = mu[18]
    g = (R0/r)**2
    ax = rho*(-g*x1/r + x1*Ro**2 + 2*u2*Ro)
    ay = rho*(-g*x2/r + x2*Ro**2 - 2*u1*Ro)
    az = rho*(-g*x3/r)

    # Collisions
    col = op*m_op*nu_op + o2p*m_o2p*nu_o2p + np*m_np*nu_np + n2p*m_n2p*nu_n2p + nop*m_nop*nu_nop
    id1 = Kn*col*(u1-vi1)
    id2 = Kn*col*(u2-vi2)
    id3 = Kn*col*(u3-vi3)

    # Lorentz force: E + vxB
    fLx = Omega*ne*(E1 + u2*Bmag*b3 - u3*Bmag*b2)
    fLy = Omega*ne*(E2 + u3*Bmag*b1 - u1*Bmag*b3)
    fLz = Omega*ne*(E3 + u1*Bmag*b2 - u2*Bmag*b1)

    # collision
    col = op*m_op*nu_op + o2p*m_o2p*nu_o2p + np*m_np*nu_np + n2p*m_n2p*nu_n2p + nop*m_nop*nu_nop

    s_op = 0
    s_o2p = 0
    s_np = 0
    s_n2p = 0
    s_nop = 0
    source_phi1 = ax + fLx + id1
    source_phi2 = ay + fLy + id2
    source_phi3 = az + fLz + id3
    source_pi = -(gamma_i-1)*pi*divV
    source_pe = -(gamma_e-1)*pe*divV

    s = array([s_op, s_o2p, s_np, s_n2p, s_nop,
               source_phi1, source_phi2, source_phi3, source_pi, source_pe])
    return s

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = zeros(shape=(10, 2))
    return ub

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fb = zeros(shape=(10, 2))
    return fb

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):

    ophat   = uhat[0]
    o2phat  = uhat[1]
    nphat   = uhat[2]
    n2phat  = uhat[3]
    nophat  = uhat[4]
    phihat1 = uhat[5]
    phihat2 = uhat[6]
    phihat3 = uhat[7]
    pihat   = uhat[8]
    pehat   = uhat[9]

    op   = u[0]
    o2p  = u[1]
    np   = u[2]
    n2p  = u[3]
    nop  = u[4]
    phi1 = u[5]
    phi2 = u[6]
    phi3 = u[7]
    pi   = u[8]
    pe   = u[9]

    opx   = -q[0]
    o2px  = -q[1]
    npx   = -q[2]
    n2px  = -q[3]
    nopx  = -q[4]
    phi1x = -q[5]
    phi2x = -q[6]
    phi3x = -q[7]
    pix   = -q[8]
    pex   = -q[9]
    opy   = -q[10]
    o2py  = -q[11]
    npy   = -q[12]
    n2py  = -q[13]
    nopy  = -q[14]
    phi1y = -q[15]
    phi2y = -q[16]
    phi3y = -q[17]
    piy   = -q[18]
    pey   = -q[19]
    opz   = -q[20]
    o2pz  = -q[21]
    npz   = -q[22]
    n2pz  = -q[23]
    nopz  = -q[24]
    phi1z = -q[25]
    phi2z = -q[26]
    phi3z = -q[27]
    piz   = -q[28]
    pez   = -q[29]

    Pe = mu[4] # Peclet number (inverse)
    gamma_i = mu[0] # ion adiabatic index, 5/3
    gamma_e = mu[0] # electron adiabatic index, 5/3
    m_op    = mu[5] # molar mass of O+, 16
    m_o2p   = mu[6] # molar mass of O2+, 32
    m_np    = mu[7] # molar mass of N+, 14
    m_n2p   = mu[8] # molar mass of N2+, 28
    m_nop   = mu[9] # molar mass of NO+, 30
    psi_op  = mu[10] # O+ number flux at upper boundary
    psi_o2p = mu[11] # O2+ number flux at upper boundary
    psi_np  = mu[12] # N+ number flux at upper boundary
    psi_n2p = mu[13] # N2+ number flux at upper boundary
    psi_nop = mu[14] # NO+ number flux at upper boundary
    qi      = mu[15] # ion heat flux at upper boundary
    qe      = mu[16] # electron heat flux at upper boundary

    ne = op + o2p + np + n2p + nop
    rho = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    Ti = pi/ne
    Te = pe/ne
    ne  = smoothpos(ne, 1e-6)
    rho = smoothpos(rho, 1e-4)
    Ti  = smoothpos(Ti, 1e-3)
    Te  = smoothpos(Te, 1e-3)

    nx = opx + o2px + npx + n2px + nopx
    ny = opy + o2py + npy + n2py + nopy
    nz = opz + o2pz + npz + n2pz + nopz
    rhox = opx*m_op + o2px*m_o2p + npx*m_np + n2px*m_n2p + nopx*m_nop
    rhoy = opy*m_op + o2py*m_o2p + npy*m_np + n2py*m_n2p + nopy*m_nop
    rhoz = opz*m_op + o2pz*m_o2p + npz*m_np + n2pz*m_n2p + nopz*m_nop

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
    kappa_i0 = v[23]
    kappa_e0 = v[24]
    Tn = v[25]
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

    un = u1*n[0] + u2*n[1] + u3*n[2]



    # Lower boundary: nhat=n, v = (I-nn)*u, T=Tn
    fl = array([ophat - op,
                o2phat - o2p,
                nphat - np,
                n2phat - n2p,
                nophat - nop,
                phihat1 - rho*(u1-un*n[0]),
                phihat2 - rho*(u2-un*n[1]),
                phihat3 - rho*(u3-un*n[2]),
                pihat - ne*Tn,
                pehat - ne*Tn])

    # Upper boundary:
    # number flux nhat=n
    # vperp = ExB/B^2 & d(vpar)/ds = 0
    # heat flux gradT = q
    vd1 = (E2*b3 - E3*b2)/Bmag
    vd2 = (E3*b1 - E1*b3)/Bmag
    vd3 = (E1*b2 - E2*b1)/Bmag
    phihatb = phihat1*b1 + phihat2*b2 + phihat3*b3
    fperp1 = phihat1 - phihatb*b1 - rho*vd1
    fperp2 = phihat2 - phihatb*b2 - rho*vd2
    fperp3 = phihat3 - phihatb*b3 - rho*vd3

    bgphib = (phi1x*b1**2 + phi2x*b1*b2 + phi3x*b1*b3
            + phi1y*b2*b1 + phi2y*b2**2 + phi3y*b2*b3
            + phi1z*b3*b1 + phi2z*b3*b2 + phi3z*b3**2)
    bgbphi = (b1x*b1*phi1 + b2x*b1*phi2 + b3x*b1*phi3
            + b1y*b2*phi1 + b2y*b2*phi2 + b3y*b2*phi3
            + b1z*b3*phi1 + b2z*b3*phi2 + b3z*b3*phi3)
    bgn = rhox*b1 + rhoy*b2 + rhoz*b3
    phib = phi1*b1 + phi2*b2 + phi3*b3
    fpar = bgphib + bgbphi - bgn*phib/rho +  tau[0]*(phib-phihatb)

    gradpin = pix*n[0] + piy*n[1] + piz*n[2]
    gradpen = pex*n[0] + pey*n[1] + pez*n[2]
    gradnn = nx*n[0] + ny*n[1] + nz*n[2]

    gradTin = (gradpin-gradnn*Ti)/ne
    gradTen = (gradpen-gradnn*Te)/ne
    kappa_i = kappa_i0*Ti**2.5
    kappa_e = kappa_e0*Te**2.5

    fu = array([op - ophat,
                o2p - o2phat,
                np - nphat,
                n2p - n2phat,
                nop - nophat,
                fperp1 + fpar*b1,
                fperp2 + fpar*b2,
                fperp3 + fpar*b3,
                (gamma_i-1)*Pe*kappa_i*gradTin - qi + tau[0]*(pi-pihat),
                (gamma_e-1)*Pe*kappa_e*gradTen - qe + tau[0]*(pe-pehat)])

    fb = hstack((fl, fu))
    fb = reshape(fb,(10,2),'F')
    return fb

def initu(x, mu, eta):
    u0 = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return u0

def smoothpos(a, eps):
    # smooth positive part: ~a for a >> eps, ~eps/2 for a <= 0, always > 0
    return 0.5*(a + sqrt(a*a + eps*eps))
