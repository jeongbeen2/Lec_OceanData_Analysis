import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.stats as stats



data = nc.Dataset("T2m_ERA5_1979_2018_lowR.nc", 'r' )
lon = data.variables['lon'][:] 
lat = data.variables['lat'][:] 
time = data.variables['time'][:] 
T2m = data.variables['t2m'][:,:,:]

T2m.shape
cosarray = np.zeros((len(lat),len(lon)))

for x in range(0,len(lat)):
    cosarray[x,:] = np.cos(lat[x]*np.pi/180)

wgt_T2m_year0 = np.zeros(T2m.shape)
for t in range(0, len(T2m)):
    wgt_T2m_year0[t,:,:] = T2m[t,:,:]*cosarray
    
wgt_T2m_year = np.sum(wgt_T2m_year0, 2)
avg_T2m_year = np.sum(wgt_T2m_year, 1)/np.sum(cosarray) - 273

x = np.arange(1979, 2019) # 2018 + 1
r = stats.linregress(x, avg_T2m_year)
beta = r.slope
alpha = r.intercept

y = beta*x + alpha

fig_1 = plt.plot(x,avg_T2m_year, 'go--')
#fig_2 = plt.plot(x, avg_T2m_year, 'go--',x,y,'r')


