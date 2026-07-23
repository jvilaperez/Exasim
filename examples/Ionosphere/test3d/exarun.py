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

# earth radius, gravity at earth surface, earth rotation,
# molar mass of O+, O2+, N+, N2+, NO+,
# ion adiabatic index, electron adiabatic index,
# ratio between thermal energy and kinetic energy, time scale of gyro motion
pde['physicsparam'] = numpy.array([Re/H0, G0, Wrot, 16, 32, 14, 28, 30, 5/3, 5/3, lamda, omega])

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
ux = data['UN'][:].filled()*t0/H0
uy = data['VN'][:].filled()*t0/H0
uz = data['WN'][:].filled()*t0/H0
bx = data['BX'][:].filled()/B0
by = data['BY'][:].filled()/B0
bz = data['BZ'][:].filled()/B0
ex = data['EX'][:].filled()*t0/(B0*H0)
ey = data['EY'][:].filled()*t0/(B0*H0)
ez = data['EZ'][:].filled()*t0/(B0*H0)
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
kappa_i = data['KAPPA_I'][:].filled()*T0**2.5*t0/(n0*kB*H0**2)
kappa_e = data['KAPPA_E'][:].filled()*T0**2.5*t0/(n0*kB*H0**2)
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
vx = data['UI'][:].filled()*t0/H0
vy = data['VI'][:].filled()*t0/H0
vz = data['WI'][:].filled()*t0/H0
ti = data['TI'][:].filled()/T0
te = data['TE'][:].filled()/T0
data.close()

mesh['vdg'] = numpy.zeros((npe,33,ne))
mesh['udg'] = numpy.zeros((npe,10,ne))
for icell in range(ne):
    for inode in range(npe):
        # 0: b.grad(B), 1: tn, 2: ux, 3: uy, 4: uz, 5: bx, 6: by, 7: bz, 8: ex, 9: ey, 10: ez
        # 11: prod_op, 12: prod_o2p, 13: prod_np, 14: prod_n2p, 15: prod_nop
        # 16: loss_op, 17: loss_o2p, 18: loss_np, 19: loss_n2p, 20: loss_nop
        # 21: nu_op, 22: nu_o2p, 23: nu_np, 24: nu_n2p, 25: nu_nop
        # 26: kappa_i, 27: kappa_e, 28: heat_i, 29: heat_e, 30: qei, 31: qin, 32: qen
        mesh['vdg'][inode, 1, icell] = tn[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 2, icell] = ux[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 3, icell] = uy[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 4, icell] = uz[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 5, icell] = bx[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 6, icell] = by[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 7, icell] = bz[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 8, icell] = ex[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 9, icell] = ey[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 10, icell] = ez[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 11, icell] = prod_op[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 12, icell] = prod_o2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 13, icell] = prod_np[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 14, icell] = prod_n2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 15, icell] = prod_nop[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 16, icell] = loss_op[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 17, icell] = loss_o2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 18, icell] = loss_np[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 19, icell] = loss_n2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 20, icell] = loss_nop[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 21, icell] = nu_op[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 22, icell] = nu_o2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 23, icell] = nu_np[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 24, icell] = nu_n2p[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 25, icell] = nu_nop[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 26, icell] = kappa_i[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 27, icell] = kappa_e[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 28, icell] = heat_i[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 29, icell] = heat_e[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 30, icell] = qei[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 31, icell] = qin[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['vdg'][inode, 32, icell] = qen[altidx[inode, icell], nodeidx[inode, icell]]

        mesh['udg'][inode, 0, icell] = numpy.log(op[altidx[inode, icell], nodeidx[inode, icell]])
        mesh['udg'][inode, 1, icell] = numpy.log(o2p[altidx[inode, icell], nodeidx[inode, icell]])
        mesh['udg'][inode, 2, icell] = numpy.log(np[altidx[inode, icell], nodeidx[inode, icell]])
        mesh['udg'][inode, 3, icell] = numpy.log(n2p[altidx[inode, icell], nodeidx[inode, icell]])
        mesh['udg'][inode, 4, icell] = numpy.log(nop[altidx[inode, icell], nodeidx[inode, icell]])
        mesh['udg'][inode, 5, icell] = vx[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 6, icell] = vy[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 7, icell] = vz[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 8, icell] = ti[altidx[inode, icell], nodeidx[inode, icell]]
        mesh['udg'][inode, 9, icell] = te[altidx[inode, icell], nodeidx[inode, icell]]

xpe,telem,xpf,tface,perm = Preprocessing.masternodes(pde['porder'],pde['nd'],1)

gpe, gwe = Preprocessing.gaussnodes(pde['pgauss'],pde['nd'],pde['elemtype'])
gpe = numpy.array(gpe,float)
shapeg = Preprocessing.mkshape(pde['porder'],xpe,gpe,1)

for d in range(0,pde['nd']+1):
    shapeg[:,:,d] = shapeg[:,:,d].transpose()

# calculate b.grad(B)
B1 = mesh['vdg'][:,5,:]
B2 = mesh['vdg'][:,6,:]
B3 = mesh['vdg'][:,7,:]
Bmag = numpy.sqrt(B1**2 + B2**2 + B3**2)
b1 = B1/Bmag
b2 = B2/Bmag
b3 = B3/Bmag
gradB = Preprocessing.gradu(shapeg[:,:,1:4], mesh['dgnodes'], Bmag)
mesh['vdg'][:,0,:] = b1*gradB[:,0,:] + b2*gradB[:,1,:] + b3*gradB[:,2,:]

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
