# import external modules
from operator import ne

import numpy, os

# Add Exasim to Python search path
cdir = os.getcwd(); ii = cdir.find("Exasim")
exec(open(cdir[0:(ii+6)] + "/install/setpath.py").read())

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh
import netCDF4

# Create pde object and mesh object
pde,mesh = Preprocessing.initializeexasim()

# Define a PDE model: governing equations and boundary conditions
pde['model'] = "ModelD";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model
pde['hybrid'] = 1;          # 0 -> LDG, 1 -> HDG

# Choose computing platform and set number of processors
pde['platform'] = "gpu"
pde['mpiprocs'] = 4
pde['cpucompiler']  = "CC"
pde['mpicompiler'] = "CC"
pde['gpucompiler'] = cdir[0:(ii+6)] + "/kokkos/bin/nvcc_wrapper" # e.g. /path_to_kokkos/bin/nvcc_wrapper

# Set discretization parameters, physical parameters, and solver parameters
pde['nd'] = 3
pde['porder'] = 2;         # polynomial degree
pde['pgauss'] = 2*pde['porder'] # number of Gauss points for quadrature
pde['tau'] = numpy.array([2.0]);            # DG stabilization parameter
pde['GMRESrestart'] = 250; # GMRES restart parameter
pde['GMRESortho'] = 1
pde['linearsolvertol'] = 1e-6; # linear solver tolerance
pde['linearsolveriter'] = 500
pde['preconditioner'] = 1
pde['ppdegree'] = 0; # degree of polynomial precontiditoning
pde['RBdim'] = 0; # reduced basis dimension for preconditioner
pde['NLtol'] = 1e-7; # nonlinear solver tolerance
pde['NLiter'] = 10; # number of nonlinear iterations
pde['matvectol'] = 1e-6; # matrix-vector product tolerance

# Physical parameters
dt0 = 60.0  # time step size in seconds

Re = 6378e3
hL = 100e3
hT = 600e3

t0 = 10 # temporal scale (second)
H0 = 10e3 # spatial scale (meter)
n0 = 1e12 # typical electron density (per cubic meter)
T0 = 1e3 # typical temperature (Kelvin)
B0 = 1e-5 # average magnetic field (Tesla)

g0 = 9.81 # gravity at earth surface
wrot = 2*numpy.pi/86400 # earth rotation
mp = 1.67e-27 # proton mass
kB = 1.38e-23 # Boltzmann constant
e = 1.602e-19 # elementary charge

# Nondimensional parameters
R0 = (Re + hL)/H0
R1 = (Re + hT)/H0
G0 = g0*t0**2/H0
Wrot = wrot*t0
lamda = kB*T0*t0**2/(mp*H0**2) # ratio between thermal energy and kinetic energy
omega = e*B0*t0/mp # time scale of gyro motion
dampfac = 1e-6*t0 # damping factor for numerical stability

# earth radius, gravity at earth surface, earth rotation,
# molar mass of O+, O2+, N+, N2+, NO+,
# ion adiabatic index, electron adiabatic index,
# ratio between thermal energy and kinetic energy, time scale of gyro motion
pde['physicsparam'] = numpy.array([Re/H0, G0, Wrot, 16, 32, 14, 28, 30, 5/3, 5/3, lamda, omega, dampfac])

# Time-stepping parameters
pde['torder'] = 1;          # time-stepping order of accuracy
pde['nstage'] = 1;          # time-stepping number of stages
pde['dt'] = (dt0/t0)*numpy.ones(10);   # time step sizes
pde['soltime'] = numpy.arange(1,pde['dt'].size+1,1); # steps at which solution are collected
pde['visdt'] = 1.0; # visualization timestep size

# Mesh generation
hDiv = 18 #horizontal divisions
nDiv = 1 # number of grid elements per scale height
mesh['p'], mesh['t'], mesh['dgnodes'] = Mesh.cubesphere_scaleheight(pde['porder'], hDiv, nDiv, hL, hT, H0)
# mesh['p'], mesh['t'], mesh['dgnodes'] = Mesh.cubesphere(pde['porder'],R0,R1,15,10)
mesh['boundaryexpr'] = [
    lambda p: abs(p[0,:]**2 + p[1,:]**2 + p[2,:]**2 - R0**2) < 1e-6,
    lambda p: abs(p[0,:]**2 + p[1,:]**2 + p[2,:]**2 - R1**2) < 1e-6]
mesh['boundarycondition'] = numpy.array([1, 2]); # Set boundary condition for each boundary
npe, nd, ne = mesh['dgnodes'].shape

# # Defining dipole magnetic field (test)
# x1 = mesh['dgnodes'][:,0,:]; x2 = mesh['dgnodes'][:,1,:]; x3 = mesh['dgnodes'][:,2,:]
# r = numpy.sqrt(x1**2 + x2**2 + x3**2)
# Bx = 3*x1*x3/r**2; By = 3*x2*x3/r**2; Bz = (3*x3**2 - r**2)/r**2
# Bx = Bx*(R0/r)**3; By = By*(R0/r)**3; Bz = Bz*(R0/r)**3


data = netCDF4.Dataset(filename='/glade/work/haonan/data.nc')
nalt = data.dimensions['alt'].size
nnode = data.dimensions['node'].size
nodeidx = data['nodeidx_out'][:].filled() - 1
altidx = data['altidx_out'][:].filled() - 1
tn = data['TN'][:].filled()/T0
u1 = data['UN'][:].filled()*t0/H0
u2 = data['VN'][:].filled()*t0/H0
u3 = data['WN'][:].filled()*t0/H0
B1 = data['BX'][:].filled()
B2 = data['BY'][:].filled()
B3 = data['BZ'][:].filled()
Bmag = numpy.sqrt(B1**2 + B2**2 + B3**2)
b1 = B1/Bmag
b2 = B2/Bmag
b3 = B3/Bmag
Bmag = Bmag/B0
E1 = data['EX'][:].filled()*t0/(B0*H0)
E2 = data['EY'][:].filled()*t0/(B0*H0)
E3 = data['EZ'][:].filled()*t0/(B0*H0)
prod_op = data['PROD_OP'][:].filled()*t0/n0
prod_o2p = data['PROD_O2P'][:].filled()*t0/n0
prod_np = data['PROD_NP'][:].filled()*t0/n0
prod_n2p = data['PROD_N2P'][:].filled()*t0/n0
prod_nop = data['PROD_NOP'][:].filled()*t0/n0
loss_op = data['LOSS_OP'][:].filled()*t0
loss_o2p = data['LOSS_O2P'][:].filled()*t0
loss_np = data['LOSS_NP'][:].filled()*t0
loss_n2p = data['LOSS_N2P'][:].filled()*t0
loss_nop = data['LOSS_NOP'][:].filled()*t0
nu_op = data['NU_OP'][:].filled()*t0
nu_o2p = data['NU_O2P'][:].filled()*t0
nu_np = data['NU_NP'][:].filled()*t0
nu_n2p = data['NU_N2P'][:].filled()*t0
nu_nop = data['NU_NOP'][:].filled()*t0
kappa_i = data['KAPPA_I'][:].filled()*data['TI'][:].filled()**2.5*t0/(n0*kB*H0**2)
kappa_e = data['KAPPA_E'][:].filled()*data['TE'][:].filled()**2.5*t0/(n0*kB*H0**2)
heat_i = data['HEAT_I'][:].filled()*t0/(n0*kB*T0)
heat_e = data['HEAT_E'][:].filled()*t0/(n0*kB*T0)
qei = data['QEI'][:].filled()*t0/(n0*kB)
qin = data['QIN'][:].filled()*t0/(n0*kB)
qen = data['QEN'][:].filled()*t0/(n0*kB)
op = data['OP'][:].filled()/n0
o2p = data['O2P'][:].filled()/n0
np = data['NP'][:].filled()/n0
n2p = data['N2P'][:].filled()/n0
nop = data['NOP'][:].filled()/n0
eden = op + o2p + np + n2p + nop
nm = op*16 + o2p*32 + np*14 + n2p*28 + nop*30
phi1 = data['UI'][:].filled()*nm*t0/(n0*H0)
phi2 = data['VI'][:].filled()*nm*t0/(n0*H0)
phi3 = data['WI'][:].filled()*nm*t0/(n0*H0)
pi = data['TI'][:].filled()*eden/(n0*T0)
pe = data['TE'][:].filled()*eden/(n0*T0)
data.close()

mesh['vdg'] = numpy.zeros((npe,37,ne))
mesh['udg'] = numpy.zeros((npe,10,ne))
for icell in range(ne):
    for inode in range(npe):
        # 0: tn, 1: u1, 2: u2, 3: u3, 4: b1, 5: b2, 6: b3, 7: Bmag
        # 8: lnBx, 9: lnBy, 10: lnBz, 11: b.grad(lnB), 12: E1, 13: E2, 14: E3
        # 15: prod_op, 16: prod_o2p, 17: prod_np, 18: prod_n2p, 19: prod_nop
        # 20: loss_op, 21: loss_o2p, 22: loss_np, 23: loss_n2p, 24: loss_nop
        # 25: nu_op, 26: nu_o2p, 27: nu_np, 28: nu_n2p, 29: nu_nop
        # 30: kappa_i, 31: kappa_e, 32: heat_i, 33: heat_e, 34: qei, 35: qin, 36: qen
        mesh['vdg'][inode, 0, icell] = tn[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 1, icell] = u1[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 2, icell] = u2[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 3, icell] = u3[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 4, icell] = b1[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 5, icell] = b2[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 6, icell] = b3[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 7, icell] = Bmag[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 12, icell] = E1[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 13, icell] = E2[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 14, icell] = E3[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 15, icell] = prod_op[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 16, icell] = prod_o2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 17, icell] = prod_np[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 18, icell] = prod_n2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 19, icell] = prod_nop[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 20, icell] = loss_op[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 21, icell] = loss_o2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 22, icell] = loss_np[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 23, icell] = loss_n2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 24, icell] = loss_nop[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 25, icell] = nu_op[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 26, icell] = nu_o2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 27, icell] = nu_np[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 28, icell] = nu_n2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 29, icell] = nu_nop[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 30, icell] = kappa_i[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 31, icell] = kappa_e[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 32, icell] = heat_i[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 33, icell] = heat_e[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 34, icell] = qei[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 35, icell] = qin[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 36, icell] = qen[altidx[inode, icell], nodeidx[inode, icell]]

        mesh['udg'][inode, 0, icell] = op[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 1, icell] = o2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 2, icell] = np[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 3, icell] = n2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 4, icell] = nop[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 5, icell] = phi1[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 6, icell] = phi2[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 7, icell] = phi3[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 8, icell] = pi[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 9, icell] = pe[altidx[inode, icell], nodeidx[inode, icell]]

xpe,telem,xpf,tface,perm = Preprocessing.masternodes(pde['porder'],pde['nd'],1)

gpe, gwe = Preprocessing.gaussnodes(pde['pgauss'],pde['nd'],pde['elemtype'])
gpe = numpy.array(gpe,float)
shapeg = Preprocessing.mkshape(pde['porder'],xpe,gpe,1)

for d in range(0,pde['nd']+1):
    shapeg[:,:,d] = shapeg[:,:,d].transpose()

# calculate b.grad(lnB)
lnB = numpy.log(mesh['vdg'][:,7,:])
gradlnB = Preprocessing.gradu(shapeg[:,:,1:4], mesh['dgnodes'], lnB)
mesh['vdg'][:,8,:] = gradlnB[:,0,:]
mesh['vdg'][:,9,:] = gradlnB[:,1,:]
mesh['vdg'][:,10,:] = gradlnB[:,2,:]
mesh['vdg'][:,11,:] = b1*gradlnB[:,0,:] + b2*gradlnB[:,1,:] + b3*gradlnB[:,2,:]

# search compilers and set options
pde = Gencode.setcompilers(pde)

# generate input files and store them in datain folder
pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh)

# generate source codes and store them in app folder
Gencode.gencode(pde)

# compile source codes to build an executable file and store it in build folder
compilerstr = Gencode.cmakecompile(pde)

# Run code
pde['mpirun'] = "mpirun"; # command to run MPI programs
runstr = Gencode.runcode(pde, 1)

# # # get solution from output files in dataout folder
# pde['vistime'] = [];
# sol = Postprocessing.fetchsolution(pde,master,dmd, pde['buildpath'] + "/dataout");
# x = mesh['dgnodes'][:,0,:]; y = mesh['dgnodes'][:,1,:]; z = mesh['dgnodes'][:,2,:]
# # uexact = numpy.sin(numpy.pi*x)*numpy.sin(numpy.pi*y)*numpy.sin(numpy.pi*z); # exact solution
# # uh = sol[:,0,:];  # numerical solution
# # print("Maximum absolute error: %g\n" % max(abs(uh.flatten()-uexact.flatten())));

print("Done!")
