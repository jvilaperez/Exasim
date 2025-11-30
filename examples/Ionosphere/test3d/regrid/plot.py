from netCDF4 import Dataset
from matplotlib.pyplot import figure


varname = 'N'
timeidx = 0
altidx = 0

data = Dataset(filename='data.nc')
print(data['alt'][altidx])
lon = data['lon'][:].filled()
lat = data['lat'][:].filled()
variable = data[varname][timeidx, altidx, :].filled()
data.close()

fig = figure()
axes = fig.add_subplot()
pcolor = axes.tripcolor(lon, lat, variable)
fig.colorbar(mappable=pcolor, ax=axes)
fig.savefig(fname=varname+'.png', bbox_inches='tight')
