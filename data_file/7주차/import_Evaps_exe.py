# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 22:00:32 2020

@author: current
"""

# load necessary modules
import numpy as np
import netCDF4 as nc

print(nc.Dataset('D:/ERA5/EVAPintp/EVAP.201801.nc'))

#### Import monthly evaporation data 
EVAPdir='D:/ERA5/EVAPintp/' # directory
styear=1979 #start year
edyear=2019 #end year

evaps=np.zeros([12*(edyear-styear+1), 37, 72])

for i in range(styear,edyear+1):
    for j in range(1,13): # loop for months 
        year = 
        month = 
        filename =  
        data = nc.Dataset(filename) #read EVAP file
        evap_month = data.variables['e'][:,:] #[m/day]
        evaps[? ? ?] = evap_month
        
print(i)