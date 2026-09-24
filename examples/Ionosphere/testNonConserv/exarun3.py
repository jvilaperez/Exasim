# ---------------------------------------------------------------------------
# exarun3.py -- MINIMAL VERIFICATION CASE for the non-conservative O+ model
#
# Isothermal hydrostatic atmosphere at rest.  With Ro = Kn = Omega = Pe = 0
# this is an EXACT steady solution of the 5-equation system in pdemodel2.py:
#
#     psi(r) = psi0 - (m_op*R0/T)*(1 - R0/r)     (log number density)
#     v      = 0
#     Ti     = T                                  (constant)
#
# Check:  (1/rho) grad(p) = (T/m_op) grad(psi) = -(R0/r)^2 rhat = gravity.
# The code should hold this state.  Any drift is discretisation error, so the
# steady residual is a direct measure of whether the scheme is sound.
#
# Self-contained: no netCDF input.  Same mesh and same pdemodel2.py as exarun2.
# ---------------------------------------------------------------------------
from operator import ne
import numpy, os

cdir = os.getcwd(); ii = cdir.find("Exasim")
exec(open(cdir[0:(ii+6)] + "/install/setpath.py").read())

import Preprocessing, Postprocessing, Gencode, Mesh

pde,mesh = Preprocessing.initializeexasim()

pde['model'] = "ModelD"
pde['modelfile'] = "pdemodel2"   # unchanged -- this tests the real model
pde['hybrid'] = 1

pde['platform'] = "gpu"
pde['mpiprocs'] = 8
pde['cpucompiler']  = "CC"
pde['mpicompiler'] = "CC"
pde['gpucompiler'] = cdir[0:(ii+6)] + "/kokkos/bin/nvcc_wrapper"

pde['nd'] = 3
pde['porder'] = 3
pde['pgauss'] = 2*pde['porder']
pde['tau'] = numpy.array([5.0])

pde['GMRESrestart'] = 249
pde['GMRESortho'] = 1
pde['linearsolvertol'] = 1e-8
pde['linearsolveriter'] = 250
pde['preconditioner'] = 1
pde['ppdegree'] = 0
pde['RBdim'] = 10
pde['NLtol'] = 1e-10       # absolute; the equilibrium residual is ~1e-3 so
pde['NLiter'] = 5         # Newton should stop on the iteration count
pde['matvectol'] = 1e-10

flagGencode = False

# ------------------------------- constants ---------------------------------
dt0 = 1.0          # s
tf  = 60.0         # s  (short: we only want to see whether it stays put)
nsave = 1
nt = int(tf/dt0)

Re = 6378e3
hL = 100e3; hT = 600e3

gam = 5/3
mp = 1.67e-27
kB = 1.38e-23
e  = 1.602e-19
g0 = 9.81*Re*Re/(Re + hL)**2

n0 = 1e12
T0 = 200.0
H0 = kB*T0/(mp*g0)
v0 = numpy.sqrt(gam*kB*T0/mp)
t0 = H0/v0

R0 = (Re + hL)/H0
R1 = (Re + hT)/H0
m_op = 16

# THE knob: uniform ion temperature, in units of T0 = 200 K.
# T_hat = 5 -> 1000 K -> psi drops by 8.55 across the domain (8.5 scale heights)
T_hat = 5.0

# every forcing term switched off -- this is what makes the state exact
Ro = 0.0; Kn = 0.0; Omega = 0.0; Pe = 0.0
psi_op = 0.0; qi = 0.0; qe = 0.0; dampfac = 0.0
m_o2p = 32; m_np = 14; m_n2p = 28; m_nop = 30

pde['physicsparam'] = numpy.array([gam, Ro, Kn, Omega, Pe, m_op, m_o2p, m_np, m_n2p, m_nop,
                                   psi_op, 0.0, 0.0, 0.0, 0.0, qi, qe, dampfac,
                                   R0, R1, Re, H0])

pde['torder'] = 2
pde['nstage'] = 2
pde['dt'] = (dt0/t0)*numpy.ones(nt)
pde['soltime'] = numpy.arange(1,pde['dt'].size+1,nsave)
pde['visdt'] = 1.0

# --------------------------------- mesh ------------------------------------
hDiv = 18
nDiv = 1
mesh['p'], mesh['t'], mesh['dgnodes'] = Mesh.cubesphere_scaleheight(pde['porder'], hDiv, nDiv, hL, hT, H0)
mesh['boundaryexpr'] = [
    lambda p: abs(p[0,:]**2 + p[1,:]**2 + p[2,:]**2 - R0**2) < 1e-6,
    lambda p: abs(p[0,:]**2 + p[1,:]**2 + p[2,:]**2 - R1**2) < 1e-6]
mesh['boundarycondition'] = numpy.array([1, 2])
npe, nd, ne = mesh['dgnodes'].shape

# ------------------------- analytic initial state --------------------------
x1 = mesh['dgnodes'][:,0,:]; x2 = mesh['dgnodes'][:,1,:]; x3 = mesh['dgnodes'][:,2,:]
r  = numpy.sqrt(x1**2 + x2**2 + x3**2)

psi = -(m_op*R0/T_hat)*(1.0 - R0/r)      # psi0 = 0 at the lower boundary

mesh['udg'] = numpy.zeros((npe,5,ne))
mesh['udg'][:,0,:] = psi                 # ln(n_O+)
mesh['udg'][:,1,:] = 0.0                 # vi1
mesh['udg'][:,2,:] = 0.0                 # vi2
mesh['udg'][:,3,:] = 0.0                 # vi3
mesh['udg'][:,4,:] = T_hat               # Ti

# -------------------------- external fields (vdg) --------------------------
# Only b, Bmag, grad(b) and Tn are read by pdemodel2 in this configuration;
# everything else multiplies a zeroed coefficient.  b is taken radial.
mesh['vdg'] = numpy.zeros((npe,46,ne))
b1 = x1/r; b2 = x2/r; b3 = x3/r
mesh['vdg'][:,25,:] = T_hat              # Tn  -> lower BC Tihat = Tn is satisfied
mesh['vdg'][:,26,:] = 0.0                # u1  -> lower BC vihat = 0 is satisfied
mesh['vdg'][:,27,:] = 0.0                # u2
mesh['vdg'][:,28,:] = 0.0                # u3
mesh['vdg'][:,29,:] = b1
mesh['vdg'][:,30,:] = b2
mesh['vdg'][:,31,:] = b3
mesh['vdg'][:,32,:] = 1.0                # Bmag (only ever a divisor; E = 0 so vd = 0)

xpe,telem,xpf,tface,perm = Preprocessing.masternodes(pde['porder'],pde['nd'],1)
gpe, gwe = Preprocessing.gaussnodes(pde['pgauss'],pde['nd'],pde['elemtype'])
gpe = numpy.array(gpe,float)
shapeg = Preprocessing.mkshape(pde['porder'],xpe,gpe,1)
for d in range(0,pde['nd']+1):
    shapeg[:,:,d] = shapeg[:,:,d].transpose()

mesh['vdg'][:,37:40,:] = Preprocessing.gradu(shapeg[:,:,1:4], mesh['dgnodes'], b1)
mesh['vdg'][:,40:43,:] = Preprocessing.gradu(shapeg[:,:,1:4], mesh['dgnodes'], b2)
mesh['vdg'][:,43:46,:] = Preprocessing.gradu(shapeg[:,:,1:4], mesh['dgnodes'], b3)

print("hydrostatic IC:  T_hat = %g (%.0f K),  psi in [%.3f, %.3f]  (%.1f scale heights)"
      % (T_hat, T_hat*T0, psi.min(), psi.max(), psi.max()-psi.min()))
print("                 n in [%.3e, %.3e] n0,   r in [%.3f, %.3f]"
      % (numpy.exp(psi).min(), numpy.exp(psi).max(), r.min(), r.max()))

# ------------------------------- run ---------------------------------------
pde = Gencode.setcompilers(pde)
pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh)
if flagGencode:
    Gencode.gencode(pde)
compilerstr = Gencode.cmakecompile(pde)
pde['mpirun'] = "mpirun"
runstr = Gencode.runcode(pde, 1)
print("Done!")
