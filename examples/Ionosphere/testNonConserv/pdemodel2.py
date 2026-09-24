from numpy import ones, reshape, array, zeros, hstack
from sympy import sqrt, exp

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0, 1.0, 1.0])
    return m

def flux(u, q, w, v, x, t, mu, eta):
    ln_op = u[0]
    vi1 = u[1]
    vi2 = u[2]
    vi3 = u[3]
    Ti = u[4]

    Tix = -q[4]
    Tiy = -q[9]
    Tiz = -q[14]

    Pe = mu[4] # Peclet number (inverse)
    m_op  = mu[5] # molar mass of O+, 16
    gamma_i = mu[0] # ion adiabatic index, 5/3

    kappa_i0 = v[23]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]

    nop = exp(ln_op)

    # kappa_i = kappa_i0*Ti**2.5
    kappa_i = 1

    hfi1 = (gamma_i-1)*Pe*kappa_i*Tix/nop
    hfi2 = (gamma_i-1)*Pe*kappa_i*Tiy/nop
    hfi3 = (gamma_i-1)*Pe*kappa_i*Tiz/nop


    f = array([ ln_op*vi1, vi1*vi1 + Ti/m_op, vi2*vi1, vi3*vi1, Ti*vi1 - hfi1,
                ln_op*vi2, vi1*vi2, vi2*vi2 + Ti/m_op, vi3*vi2, Ti*vi2 - hfi2,
                ln_op*vi3, vi1*vi3, vi2*vi3, vi3*vi3 + Ti/m_op, Ti*vi3 - hfi3])

    f = reshape(f,(5,3),'F')
    return f

def source(u, q, w, v, x, t, mu, eta):
    ln_op = u[0]
    vi1 = u[1]
    vi2 = u[2]
    vi3 = u[3]
    Ti = u[4]

    ln_opx = -q[0]
    vi1x = -q[1]
    vi2x = -q[2]
    vi3x = -q[3]
    Tix = -q[4]
    ln_opy = -q[5]
    vi1y = -q[6]
    vi2y = -q[7]
    vi3y = -q[8]
    Tiy = -q[9]
    ln_opz = -q[10]
    vi1z = -q[11]
    vi2z = -q[12]
    vi3z = -q[13]
    Tiz = -q[14]

    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    Ro = mu[1] # Rossby number (inverse)
    Kn = mu[2] # Knudsen number (inverse)
    Omega = mu[3] # ratio between gyro frequency and reference frequency
    Pe = mu[4] # Peclet number (inverse)
    m_op  = mu[5] # molar mass of O+, 16
    gamma_i = mu[0] # ion adiabatic index, 5/3

    # Input fields
    nu_op  = v[10]
    E1 = v[15]
    E2 = v[16]
    E3 = v[17]
    kappa_i0 = v[23]
    b1 = v[29]
    b2 = v[30]
    b3 = v[31]
    Bmag = v[32]
    u1 = v[26]
    u2 = v[27]
    u3 = v[28]

    nop = exp(ln_op)
    divVi = vi1x + vi2y + vi3z

    # kappa_i = kappa_i0*Ti**2.5
    kappa_i = 1

    # Accelerations
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    r = sqrt(x1**2 + x2**2 + x3**2)
    R0 = mu[18]
    g = (R0/r)**2
    # ax = (-g*x1/r + x1*Ro**2 + 2*u2*Ro)
    # ay = (-g*x2/r + x2*Ro**2 - 2*u1*Ro)
    # az = (-g*x3/r)

    ax = -g*x1/r
    ay = -g*x2/r
    az = -g*x3/r

    # Collisions
    id1 = Kn*nu_op*(u1-vi1)
    id2 = Kn*nu_op*(u2-vi2)
    id3 = Kn*nu_op*(u3-vi3)

    # Lorentz force: E + vxB
    fLx = Omega*(E1 + vi2*Bmag*b3 - vi3*Bmag*b2)/m_op
    fLy = Omega*(E2 + vi3*Bmag*b1 - vi1*Bmag*b3)/m_op
    fLz = Omega*(E3 + vi1*Bmag*b2 - vi2*Bmag*b1)/m_op

    # diff source
    lneTi = ln_opx*Tix + ln_opy*Tiy + ln_opz*Tiz
    # continuity source
    s_op = (ln_op - 1)*divVi

    source_vi1 = ax + divVi*vi1 - Ti*ln_opx/m_op + fLx + id1
    source_vi2 = ay + divVi*vi2 - Ti*ln_opy/m_op + fLy + id2
    source_vi3 = az + divVi*vi3 - Ti*ln_opz/m_op + fLz + id3
    source_Ti = (2-gamma_i)*Ti*divVi + (gamma_i-1)*Pe*kappa_i*lneTi/nop

    s = array([s_op, source_vi1, source_vi2, source_vi3, source_Ti])
    return s

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = zeros(shape=(5, 2))
    return ub

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fb = zeros(shape=(5, 2))
    return fb

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):

    ln_ophat   = uhat[0]
    vihat1 = uhat[1]
    vihat2 = uhat[2]
    vihat3 = uhat[3]
    Tihat   = uhat[4]

    ln_op   = u[0]
    vi1 = u[1]
    vi2 = u[2]
    vi3 = u[3]
    Ti   = u[4]

    lnopx   = -q[0]
    vi1x = -q[1]
    vi2x = -q[2]
    vi3x = -q[3]
    Tix   = -q[4]
    lnopy   = -q[5]
    vi1y = -q[6]
    vi2y = -q[7]
    vi3y = -q[8]
    Tiy   = -q[9]
    lnopz   = -q[10]
    vi1z = -q[11]
    vi2z = -q[12]
    vi3z = -q[13]
    Tiz   = -q[14]

    Pe = mu[4] # Peclet number (inverse)
    gamma_i = mu[0] # ion adiabatic index, 5/3
    m_op    = mu[5] # molar mass of O+, 16
    psi_op  = mu[10] # O+ number flux at upper boundary
    qi      = mu[15] # ion heat flux at upper boundary

    nop = exp(ln_op)

    E1 = v[15]
    E2 = v[16]
    E3 = v[17]
    kappa_i0 = v[23]
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
                vihat1 - (u1-un*n[0]),
                vihat2 - (u2-un*n[1]),
                vihat3 - (u3-un*n[2]),
                Tihat - Tn])

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
    # kappa_i = kappa_i0*Ti**2.5
    kappa_i = 1

    fu = array([ln_op - ln_ophat,
                fperp1 + fpar*b1,
                fperp2 + fpar*b2,
                fperp3 + fpar*b3,
                (gamma_i-1)*Pe*kappa_i*gradTin/nop - qi + tau[0]*(Ti-Tihat)])

    fb = hstack((fl, fu))
    fb = reshape(fb,(5,2),'F')
    return fb

def initu(x, mu, eta):
    u0 = array([0.0, 0.0, 0.0, 0.0, 0.0])
    return u0

def smoothpos(a, amin):
    # smooth max(a, amin): >= amin always, == a for a >> amin
    return 0.5*(a + amin + sqrt((a - amin)**2 + amin**2))
