# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 22:00:32 2020

@author: current
"""
# load necessary modules
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import math

EVAPdir='EVAPintp/' # directory
data1 = nc.Dataset('EVAPintp/EVAP.201801.nc')
lat = data1.variables['lat'][:]
lon = data1.variables['lon'][:]
time = data1.variables['time'][:]

#### Import monthly evaporation data 
EVAPdir='EVAPintp/' # directory
styear=1979 #start year
edyear=2019 #end year

evap_year = np.zeros([(edyear-styear+1), len(lat), len(lon)])
nday_list1=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
nday_list2=np.array([31,29,31,30,31,30,31,31,30,31,30,31])
# list 
list_evaps = []
for i in range(styear,edyear+1):
    ## checking leap year
    if i%4 != 0:
        nday = nday_list1[:]
    else:
        nday = nday_list2[:]
	
    ## Save annual evaporation (m) in evap_year 
    for j in range(1,13): # loop for months 
        year = str(i)
        month = str(j).zfill(2) # two digit string
        filename = EVAPdir+'EVAP.'+year+month+'.nc' #file name
        data = nc.Dataset(filename) #read EVAP file
        evap_month = data.variables['e'][:,:] #[m/day]
        evap_year[i-styear,:,:] += evap_month[0,:,:]*nday[j-1]
        # list
        list_evaps.append(evap_month[0,:,:]*nday[j-1])
        
    print(i)

#### Cosine weighting & area average  
cosarray=np.zeros((len(lat),len(lon))) 

for x in range(0, len(lon)):
	cosarray[:,x]=np.cos(lat*np.pi/180)

wgt_evap_year0 = np.zeros(evap_year.shape)
for t in range(0,len(evap_year)):
    wgt_evap_year0[t,:,:] = evap_year[t,:,:]*cosarray

wgt_evap_year1 = np.sum(wgt_evap_year0,2)
avg_evap_year = np.sum(wgt_evap_year1,1)/np.sum(cosarray)

plt.plot(range(styear,edyear+1), -avg_evap_year)


