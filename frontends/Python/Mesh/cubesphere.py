from numpy import *
import Preprocessing, Mesh

def cubesphere(order, R0, R1, n, m):
    pc, tc = Mesh.cubemesh(n, n, m, 1)
    pc[2, :] = loginc(pc[2, :], 3)

    pc[0, :] = -pi / 4 + pi * pc[0, :] / 2
    pc[1, :] = -pi / 4 + pi * pc[1, :] / 2
    pc[2, :] = R0 + (R1 - R0) * pc[2, :]

    dgnod = Preprocessing.createdgnodes(pc, tc, 0, [], [], order)
    npe, nd, ne = dgnod.shape
    dgnod = dgnod.transpose(1, 0, 2)
    dgnod = dgnod.reshape((nd, npe * ne), order='F') 

    offe = tc.shape[1]
    offp = pc.shape[1]
    offd = npe * ne

    ph = zeros((3, 6 * offp))
    th = zeros((8, 6 * offe), dtype=int)
    dg = zeros((3, 6 * npe * ne))

    p = transform(pc)
    g = transform(dgnod)

    ph[:, 0:offp] = p
    dg[:, 0:offd] = g
    th[:, 0:offe] = tc

    R = rotmat(array([0, 0, 1]))
    for i in range(1, 4):
        p = R @ p
        g = R @ g
        ph[:, i * offp:(i + 1) * offp] = p
        dg[:, i * offd:(i + 1) * offd] = g
        th[:, i * offe:(i + 1) * offe] = tc + i * offp

    R = rotmat(array([0, 1, 0]))
    p5 = R.T @ ph[:, :offp]
    g5 = R.T @ dg[:, :offd]
    ph[:, 4 * offp:5 * offp] = p5
    dg[:, 4 * offd:5 * offd] = g5
    th[:, 4 * offe:5 * offe] = tc + 4 * offp

    p6 = R @ ph[:, :offp]
    g6 = R @ dg[:, :offd]
    ph[:, 5 * offp:6 * offp] = p6
    dg[:, 5 * offd:6 * offd] = g6
    th[:, 5 * offe:6 * offe] = tc + 5 * offp

    ph, th = fixmesh(ph.T, th.T)
    p = ph.T
    t = th.T

    dg = dg.reshape((nd, npe, 6 * ne), order='F') 
    dgnodes = dg.transpose(1, 0, 2)

    # dgnodes = dg.reshape(3, npe, 6 * ne).transpose(1, 0, 2)
    return p, t, dgnodes

def cubesphere_scaleheight(order, n, nDiv, hL, hT, H0):
    
    z = scale_height_vertical_grid(z_bottom=hL/1e3, z_top=hT/1e3, nDiv=nDiv)
    m = len(z)-1
    
    pc, tc = Mesh.cubemesh(n, n, m, 1)
    pc[2, :m*(n+1)**2] = repeat(z[:m], (n+1)**2)

    Re = 6378*1e3
    R0 = (Re + hL)/H0
    R1 = (Re + hT)/H0

    pc[0, :] = -pi / 4 + pi * pc[0, :] / 2
    pc[1, :] = -pi / 4 + pi * pc[1, :] / 2
    pc[2, :] = R0 + (R1 - R0) * pc[2, :]

    dgnod = Preprocessing.createdgnodes(pc, tc, 0, [], [], order)
    npe, nd, ne = dgnod.shape
    dgnod = dgnod.transpose(1, 0, 2)
    dgnod = dgnod.reshape((nd, npe * ne), order='F') 

    offe = tc.shape[1]
    offp = pc.shape[1]
    offd = npe * ne

    ph = zeros((3, 6 * offp))
    th = zeros((8, 6 * offe), dtype=int)
    dg = zeros((3, 6 * npe * ne))

    p = transform(pc)
    g = transform(dgnod)

    ph[:, 0:offp] = p
    dg[:, 0:offd] = g
    th[:, 0:offe] = tc

    R = rotmat(array([0, 0, 1]))
    for i in range(1, 4):
        p = R @ p
        g = R @ g
        ph[:, i * offp:(i + 1) * offp] = p
        dg[:, i * offd:(i + 1) * offd] = g
        th[:, i * offe:(i + 1) * offe] = tc + i * offp

    R = rotmat(array([0, 1, 0]))
    p5 = R.T @ ph[:, :offp]
    g5 = R.T @ dg[:, :offd]
    ph[:, 4 * offp:5 * offp] = p5
    dg[:, 4 * offd:5 * offd] = g5
    th[:, 4 * offe:5 * offe] = tc + 4 * offp

    p6 = R @ ph[:, :offp]
    g6 = R @ dg[:, :offd]
    ph[:, 5 * offp:6 * offp] = p6
    dg[:, 5 * offd:6 * offd] = g6
    th[:, 5 * offe:6 * offe] = tc + 5 * offp

    ph, th = fixmesh(ph.T, th.T)
    p = ph.T
    t = th.T

    dg = dg.reshape((nd, npe, 6 * ne), order='F') 
    dgnodes = dg.transpose(1, 0, 2)

    # dgnodes = dg.reshape(3, npe, 6 * ne).transpose(1, 0, 2)
    return p, t, dgnodes


def rotmat(u):
    W = array([[0, -u[2], u[1]],
               [u[2], 0, -u[0]],
               [-u[1], u[0], 0]])
    return eye(3) + W + W @ W


def fixmesh(p, t):
    snap = 100 * max(max(p, axis=0) - min(p, axis=0)) * 1024 * finfo(float).eps
    rounded = round(p / snap) * snap
    _, ix, jx = unique(rounded, axis=0, return_index=True, return_inverse=True)
    p = p[ix]
    t = jx[t]
    if t.ndim == 1:
        t = t[newaxis, :]

    pix, ix, jx = unique(t, return_index=True, return_inverse=True)
    t = jx.reshape(t.shape)
    p = p[pix]
    
    return p, t


def transform(p):
    tanth = tan(p[0, :])
    tanph = tan(p[1, :])
    x = p[2, :] / sqrt(1 + tanth**2 + tanph**2)
    return vstack((x, x * tanth, x * tanph))


def loginc(x, alpha):
    a = min(x)
    b = max(x)
    return a + (b - a) * (exp(alpha * (x - a) / (b - a)) - 1) / (exp(alpha) - 1)


def scale_height(z, h_0=100):
    """Scale height in km. Assumes Bates temperature profile and mass profile transitioning from N2/O2 to O."""
    T_0 = 200
    T_inf = 1200
    kB = 1.38e-23
    g0 = 9.81
    Re= 6378
    m_low=28.9 # average mass of N2 and O2
    m_high=16.0 # mass of O
    
    # Gravity
    g = g0 * (Re / (Re + z))**2

    # Mass profile
    z_turbo=120; width=40
    m_amu = 1.67e-27*(m_high + (m_low - m_high) * 0.5 * (1 - tanh((z - z_turbo) / width)))

    # Bates temperature profile
    s=0.02
    T = T_inf - (T_inf - T_0) * exp(-s * (z - h_0))

    # Scale height
    H = kB * T / (m_amu * g)
    return H / 1e3  # in km


def scale_height_vertical_grid(z_bottom=100, z_top=600, nDiv=2):
    """Build altitude grid in [0,1] with nDiv divisions per scale height"""
    z = [z_bottom]
    while z[-1] < z_top:
        H = scale_height(z[-1], h_0=z_bottom)
        dz = H/nDiv
        z.append(z[-1] + dz)
    # z[-1] = z_top  # snap last point to top
    z = array(z)
    # z = (z_top - z_bottom) * (z - z_bottom) / (z[-1] - z_bottom) + z_bottom
    z = (z - z_bottom) / (z[-1] - z_bottom) # normalized to [0, 1]
    return z

