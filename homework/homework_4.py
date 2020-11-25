import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.stats as stats

#### Import monthly evaporation data 
EVAPdir='EVAPintp/' # directory
PRECdir='PRECintp/' # directory
styear=1979 #start year
edyear=2019 #end year

common_data = nc.Dataset('EVAPintp/EVAP.201801.nc')
lat = common_data.variables['lat'][:]
lon = common_data.variables['lon'][:]
time = common_data.variables['time'][:]

evap_year = np.zeros( [ (edyear-styear+1), len(lat), len(lon)] )
prec_year = np.zeros( [ (edyear-styear+1), len(lat), len(lon)] )

nday_list1 = np.array( [31,28,31,30,31,30,31,31,30,31,30,31] )
nday_list2 = np.array( [31,29,31,30,31,30,31,31,30,31,30,31] )

for i in range(styear,edyear+1):
    if i%4 != 0:
        nday = nday_list1[:]
    else:
        nday = nday_list2[:]
    for j in range(1,13): # loop for months 
        year = str(i)
        if j < 10:
            month = "0"+str(j)
        else:
            month = str(j)
            
        evap_filename = EVAPdir+'EVAP.'+year+month+'.nc' 
        data = nc.Dataset(evap_filename) 
        evap_month = data.variables['e'][:,:] 
        evap_year[ i-styear, :, : ] += evap_month[0,:,:]*nday[j-1]
        
        prec_filename = PRECdir+'PREC.'+year+month+'.nc' 
        data = nc.Dataset(prec_filename) 
        prec_month = data.variables['tp'][:,:] 
        prec_year[ i-styear, :, : ] += prec_month[0,:,:]*nday[j-1]

#### Cosine weighting & area average  
cosarray=np.zeros((len(lat),len(lon))) 

for x in range(0, len(lat)):
	cosarray[x,:]=np.cos(lat[x]*np.pi/180)

wgt_evap_year0 = np.zeros(evap_year.shape)
wgt_prec_year0 = np.zeros(prec_year.shape)
for t in range(0,len(evap_year)):
    wgt_evap_year0[t,:,:] = -evap_year[t,:,:]*cosarray
    wgt_prec_year0[t,:,:] = prec_year[t,:,:]*cosarray

wgt_evap_year = np.sum(wgt_evap_year0,2)
wgt_prec_year = np.sum(wgt_prec_year0,2)
avg_evap_year = np.sum(wgt_evap_year,1)/np.sum(cosarray)
avg_prec_year = np.sum(wgt_prec_year,1)/np.sum(cosarray)


r = stats.linregress(avg_evap_year,avg_prec_year)
alpha = r.intercept
beta = r.slope

x = np.arange(1,1.1,0.01)
y = beta*x + alpha
plt.figure(figsize=(7,7))
plt.plot(avg_evap_year,avg_prec_year,'go',x,y,'r')
















