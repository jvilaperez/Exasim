import os
import netCDF4

cdir = os.getcwd()
ii = cdir.find("Exasim")
exec(open(cdir[0:(ii+6)] + "/install/setpath.py").read())

import Mesh

Re = 6378e3
hL = 100e3
hT = 500e3

m = 16*1.67e-27
kB = 1.38e-23
g0 = 9.81*Re*Re/(Re + hL)**2
T0 = 200
H0 = kB*T0/(m*g0)

R0 = (Re + hL)/H0
R1 = (Re + hT)/H0

_, _, dgnodes = Mesh.cubesphere(2, R0, R1, 16, 10)

data = netCDF4.Dataset(filename='coord.nc', mode='w')
data.createDimension(dimname='node_per_cell', size=dgnodes.shape[0])
data.createDimension(dimname='cell', size=dgnodes.shape[2])
x = data.createVariable(varname='x', datatype='f8', dimensions=('node_per_cell', 'cell'))
y = data.createVariable(varname='y', datatype='f8', dimensions=('node_per_cell', 'cell'))
z = data.createVariable(varname='z', datatype='f8', dimensions=('node_per_cell', 'cell'))
x[:] = dgnodes[:, 0, :]
y[:] = dgnodes[:, 1, :]
z[:] = dgnodes[:, 2, :]
data.close()
