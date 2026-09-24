from numpy import ones, reshape, array, zeros, hstack
from sympy import sqrt, exp

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    return m

def flux(u, q, w, v, x, t, mu, eta):
    ln_op = u[0]
    ln_o2p = u[1]
    ln_np = u[2]
    ln_n2p = u[3]
    ln_nop = u[4]
    vi1 = u[5]
    vi2 = u[6]
    vi3 = u[7]
    Ti = u[8]
    Te = u[9]

    ln_opx = -q[0]
    ln_o2px = -q[1]
    ln_npx = -q[2]
    ln_n2px = -q[3]
    ln_nopx = -q[4]
    vi1x = -q[5]
    vi2x = -q[6]
    vi3x = -q[7]
    Tix = -q[8]
    Tex = -q[9]
    ln_opy = -q[10]
    ln_o2py = -q[11]
    ln_npy = -q[12]
    ln_n2py = -q[13]
    ln_nopy = -q[14]
    vi1y = -q[15]
    vi2y = -q[16]
    vi3y = -q[17]
    Tiy = -q[18]
    Tey = -q[19]
    ln_opz = -q[20]
    ln_o2pz = -q[21]
    ln_npz = -q[22]
    ln_n2pz = -q[23]
    ln_nopz = -q[24]
    vi1z = -q[25]
    vi2z = -q[26]
    vi3z = -q[27]
    Tiz = -q[28]
    Tez = -q[29]

    Pe = mu[4] # Peclet number (inverse)
    m_op  = mu[5] # molar mass of O+, 16
    m_o2p = mu[6] # molar mass of O2+, 32
    m_np  = mu[7] # molar mass of N+, 14
    m_n2p = mu[8] # molar mass of N2+, 28
    m_nop = mu[9] # molar mass of NO+, 30
    gamma_i = mu[0] # ion adiabatic index, 5/3
    gamma_e = mu[0] # electron adiabatic index, 5/3

    op = exp(ln_op)
    o2p = exp(ln_o2p)
    np = exp(ln_np)
    n2p = exp(ln_n2p)
    nop = exp(ln_nop)

    ne = op + o2p + np + n2p + nop
    rho = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    mean_m = rho/ne

    ve1 = vi1
    ve2 = vi2
    ve3 = vi3

    kappa_i0 = v[23]
    kappa_e0 = v[24]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]

    kappa_i = kappa_i0*Ti**2.5/ne
    kappa_e = kappa_e0*Te**2.5/ne
    bgti = b1*Tix + b2*Tiy + b3*Tiz
    bgte = b1*Tex + b2*Tey + b3*Tez

    # hfi1 = (gamma_i-1)*Pe*kappa_i*bgti*b1
    # hfi2 = (gamma_i-1)*Pe*kappa_i*bgti*b2
    # hfi3 = (gamma_i-1)*Pe*kappa_i*bgti*b3
    # hfe1 = (gamma_e-1)*Pe*kappa_e*bgte*b1
    # hfe2 = (gamma_e-1)*Pe*kappa_e*bgte*b2
    # hfe3 = (gamma_e-1)*Pe*kappa_e*bgte*b3

    hfi1 = (gamma_i-1)*Pe*kappa_i*Tix
    hfi2 = (gamma_i-1)*Pe*kappa_i*Tiy
    hfi3 = (gamma_i-1)*Pe*kappa_i*Tiz
    hfe1 = (gamma_e-1)*Pe*kappa_e*Tex
    hfe2 = (gamma_e-1)*Pe*kappa_e*Tey
    hfe3 = (gamma_e-1)*Pe*kappa_e*Tez

    f = array([ ln_op*vi1, ln_o2p*vi1, ln_np*vi1, ln_n2p*vi1, ln_nop*vi1, vi1*vi1 + (Ti+Te)/mean_m, vi2*vi1, vi3*vi1, Ti*vi1 - hfi1, Te*ve1 - hfe1,
                ln_op*vi2, ln_o2p*vi2, ln_np*vi2, ln_n2p*vi2, ln_nop*vi2, vi1*vi2, vi2*vi2 + (Ti+Te)/mean_m, vi3*vi2, Ti*vi2 - hfi2, Te*ve2 - hfe2,
                ln_op*vi3, ln_o2p*vi3, ln_np*vi3, ln_n2p*vi3, ln_nop*vi3, vi1*vi3, vi2*vi3, vi3*vi3 + (Ti+Te)/mean_m, Ti*vi3 - hfi3, Te*ve3 - hfe3])

    f = reshape(f,(10,3),'F')
    return f

def source(u, q, w, v, x, t, mu, eta):
    ln_op = u[0]
    ln_o2p = u[1]
    ln_np = u[2]
    ln_n2p = u[3]
    ln_nop = u[4]
    vi1 = u[5]
    vi2 = u[6]
    vi3 = u[7]
    Ti = u[8]
    Te = u[9]

    ln_opx = -q[0]
    ln_o2px = -q[1]
    ln_npx = -q[2]
    ln_n2px = -q[3]
    ln_nopx = -q[4]
    vi1x = -q[5]
    vi2x = -q[6]
    vi3x = -q[7]
    Tix = -q[8]
    Tex = -q[9]
    ln_opy = -q[10]
    ln_o2py = -q[11]
    ln_npy = -q[12]
    ln_n2py = -q[13]
    ln_nopy = -q[14]
    vi1y = -q[15]
    vi2y = -q[16]
    vi3y = -q[17]
    Tiy = -q[18]
    Tey = -q[19]
    ln_opz = -q[20]
    ln_o2pz = -q[21]
    ln_npz = -q[22]
    ln_n2pz = -q[23]
    ln_nopz = -q[24]
    vi1z = -q[25]
    vi2z = -q[26]
    vi3z = -q[27]
    Tiz = -q[28]
    Tez = -q[29]

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

    # Input fields
    nu_op  = v[10]
    nu_o2p = v[11]
    nu_np  = v[12]
    nu_n2p = v[13]
    nu_nop = v[14]
    E1 = v[15]
    E2 = v[16]
    E3 = v[17]
    kappa_i0 = v[23]
    kappa_e0 = v[24]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]
    Bmag = v[32]
    u1 = v[26]
    u2 = v[27]
    u3 = v[28]

    op = exp(ln_op)
    o2p = exp(ln_o2p)
    np = exp(ln_np)
    n2p = exp(ln_n2p)
    nop = exp(ln_nop)

    ne = op + o2p + np + n2p + nop
    rho = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    mean_m = rho/ne

    lnx = (op*m_op*ln_opx + o2p*m_o2p*ln_o2px + np*m_np*ln_npx + n2p*m_n2p*ln_n2px + nop*m_nop*ln_nopx)/rho
    lny = (op*m_op*ln_opy + o2p*m_o2p*ln_o2py + np*m_np*ln_npy + n2p*m_n2p*ln_n2py + nop*m_nop*ln_nopy)/rho
    lnz = (op*m_op*ln_opz + o2p*m_o2p*ln_o2pz + np*m_np*ln_npz + n2p*m_n2p*ln_n2pz + nop*m_nop*ln_nopz)/rho

    lnex = (op*ln_opx + o2p*ln_o2px + np*ln_npx + n2p*ln_n2px + nop*ln_nopx)/ne
    lney = (op*ln_opy + o2p*ln_o2py + np*ln_npy + n2p*ln_n2py + nop*ln_nopy)/ne
    lnez = (op*ln_opz + o2p*ln_o2pz + np*ln_npz + n2p*ln_n2pz + nop *ln_nopz)/ne

    ve1 = vi1
    ve2 = vi2
    ve3 = vi3

    divVi = vi1x + vi2y + vi3z
    divVe = divVi

    kappa_i = kappa_i0*Ti**2.5/ne
    kappa_e = kappa_e0*Te**2.5/ne

    # Accelerations
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    r = sqrt(x1**2 + x2**2 + x3**2)
    R0 = mu[18]
    g = (R0/r)**2
    ax = (-g*x1/r + x1*Ro**2 + 2*u2*Ro)
    ay = (-g*x2/r + x2*Ro**2 - 2*u1*Ro)
    az = (-g*x3/r)

    # Collisions
    col = op*m_op*nu_op + o2p*m_o2p*nu_o2p + np*m_np*nu_np + n2p*m_n2p*nu_n2p + nop*m_nop*nu_nop
    mean_nu = col/rho
    id1 = Kn*mean_nu*(u1-vi1)
    id2 = Kn*mean_nu*(u2-vi2)
    id3 = Kn*mean_nu*(u3-vi3)

    # Lorentz force: E + vxB
    fLx = Omega*(E1 + vi2*Bmag*b3 - vi3*Bmag*b2)/mean_m
    fLy = Omega*(E2 + vi3*Bmag*b1 - vi1*Bmag*b3)/mean_m
    fLz = Omega*(E3 + vi1*Bmag*b2 - vi2*Bmag*b1)/mean_m

    # collision
    lneTi = lnex*Tix + lney*Tiy + lnez*Tiz
    lneTe = lnex*Tex + lney*Tey + lnez*Tez

    s_op = (ln_op - 1)*divVi
    s_o2p = (ln_o2p - 1)*divVi
    s_np = (ln_np - 1)*divVi
    s_n2p = (ln_n2p - 1)*divVi
    s_nop = (ln_nop - 1)*divVi
    source_vi1 = ax + divVi*vi1 - (Ti+Te)*lnx/mean_m + fLx + id1
    source_vi2 = ay + divVi*vi2 - (Ti+Te)*lny/mean_m + fLy + id2
    source_vi3 = az + divVi*vi3 - (Ti+Te)*lnz/mean_m + fLz + id3
    source_Ti = (2-gamma_i)*Ti*divVi + (gamma_i-1)*Pe*kappa_i*lneTi
    source_Te = (2-gamma_e)*Te*divVe + (gamma_e-1)*Pe*kappa_e*lneTe

    s = array([s_op, s_o2p, s_np, s_n2p, s_nop,
               source_vi1, source_vi2, source_vi3, source_Ti, source_Te])
    return s

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = zeros(shape=(10, 2))
    return ub

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fb = zeros(shape=(10, 2))
    return fb

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):

    ln_ophat   = uhat[0]
    ln_o2phat  = uhat[1]
    ln_nphat   = uhat[2]
    ln_n2phat  = uhat[3]
    ln_nophat  = uhat[4]
    vihat1 = uhat[5]
    vihat2 = uhat[6]
    vihat3 = uhat[7]
    Tihat   = uhat[8]
    Tehat   = uhat[9]

    ln_op   = u[0]
    ln_o2p  = u[1]
    ln_np   = u[2]
    ln_n2p  = u[3]
    ln_nop  = u[4]
    vi1 = u[5]
    vi2 = u[6]
    vi3 = u[7]
    Ti   = u[8]
    Te   = u[9]

    lnopx   = -q[0]
    lno2px  = -q[1]
    lnnpx   = -q[2]
    lnn2px  = -q[3]
    lnnopx  = -q[4]
    vi1x = -q[5]
    vi2x = -q[6]
    vi3x = -q[7]
    Tix   = -q[8]
    Tex   = -q[9]
    lnopy   = -q[10]
    lno2py  = -q[11]
    lnnpy   = -q[12]
    lnn2py  = -q[13]
    lnnopy  = -q[14]
    vi1y = -q[15]
    vi2y = -q[16]
    vi3y = -q[17]
    Tiy   = -q[18]
    Tey   = -q[19]
    lnopz   = -q[20]
    lno2pz  = -q[21]
    lnnpz   = -q[22]
    lnn2pz  = -q[23]
    lnnopz  = -q[24]
    vi1z = -q[25]
    vi2z = -q[26]
    vi3z = -q[27]
    Tiz   = -q[28]
    Tez   = -q[29]

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

    op = exp(ln_op)
    o2p = exp(ln_o2p)
    np = exp(ln_np)
    n2p = exp(ln_n2p)
    nop = exp(ln_nop)

    ne = op + o2p + np + n2p + nop
    rho = op*m_op + o2p*m_o2p + np*m_np + n2p*m_n2p + nop*m_nop
    mean_m = rho/ne

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
    fl = array([ln_ophat - ln_op,
                ln_o2phat - ln_o2p,
                ln_nphat - ln_np,
                ln_n2phat - ln_n2p,
                ln_nophat - ln_nop,
                vihat1 - (u1-un*n[0]),
                vihat2 - (u2-un*n[1]),
                vihat3 - (u3-un*n[2]),
                Tihat - Tn,
                Tehat - Tn])

    # Upper boundary:
    # number flux nhat=n
    # vperp = ExB/B^2 & d(vpar)/ds = 0
    # heat flux gradT = q
    vd1 = (E2*b3 - E3*b2)/Bmag
    vd2 = (E3*b1 - E1*b3)/Bmag
    vd3 = (E1*b2 - E2*b1)/Bmag
    vihatb = vihat1*b1 + vihat2*b2 + vihat3*b3
    fperp1 = vihat1 - vihatb*b1 - vd1
    fperp2 = vihat2 - vihatb*b2 - vd2
    fperp3 = vihat3 - vihatb*b3 - vd3

    bgvib = (vi1x*b1**2 + vi2x*b1*b2 + vi3x*b1*b3
            + vi1y*b2*b1 + vi2y*b2**2 + vi3y*b2*b3
            + vi1z*b3*b1 + vi2z*b3*b2 + vi3z*b3**2)
    bgbvi = (b1x*b1*vi1 + b2x*b1*vi2 + b3x*b1*vi3
            + b1y*b2*vi1 + b2y*b2*vi2 + b3y*b2*vi3
            + b1z*b3*vi1 + b2z*b3*vi2 + b3z*b3*vi3)
    vib = vi1*b1 + vi2*b2 + vi3*b3
    fpar = bgvib + bgbvi +  tau[0]*(vib-vihatb)

    gradTin = Tix*n[0] + Tiy*n[1] + Tiz*n[2]
    gradTen = Tex*n[0] + Tey*n[1] + Tez*n[2]
    kappa_i = kappa_i0*Ti**2.5/ne
    kappa_e = kappa_e0*Te**2.5/ne

    fu = array([ln_op - ln_ophat,
                ln_o2p - ln_o2phat,
                ln_np - ln_nphat,
                ln_n2p - ln_n2phat,
                ln_nop - ln_nophat,
                fperp1 + fpar*b1,
                fperp2 + fpar*b2,
                fperp3 + fpar*b3,
                (gamma_i-1)*Pe*kappa_i*gradTin - qi + tau[0]*(Ti-Tihat),
                (gamma_e-1)*Pe*kappa_e*gradTen - qe + tau[0]*(Te-Tehat)])

    fb = hstack((fl, fu))
    fb = reshape(fb,(10,2),'F')
    return fb

def initu(x, mu, eta):
    u0 = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return u0

def smoothpos(a, amin):
    # smooth max(a, amin): >= amin always, == a for a >> amin
    return 0.5*(a + amin + sqrt((a - amin)**2 + amin**2))
