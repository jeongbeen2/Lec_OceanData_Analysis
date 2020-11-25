import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

data = nc.Dataset('EVAP.201907.nc')

lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
evaps = data.variables['e'][:]

evap0 = evaps[0]
evap = -evap0*100*31

levels = np.arange(-4,24,4)
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0, urcrnrlon=360)
m.drawcoastlines()
draw = plt.contourf(lon,lat,evap,levels,cmap="jet", extend='both')
plt.colorbar(draw,orientation="horizontal",label="[cm]")