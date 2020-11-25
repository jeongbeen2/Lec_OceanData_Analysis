# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 22:00:32 2020

@author: current
"""

# load necessary modules
import numpy as np
import netCDF4 as nc


#### Import monthly evaporation data 
EVAPdir='EVAPintp/' # directory
styear=1979 #start year
edyear=2019 #end year

evaps=np.zeros([12*(edyear-styear+1), 37, 72])


for i in range(styear,edyear+1):
    for j in range(1,13): # loop for months 
        year = str(i)
        if j < 10:
            month = "0"+str(j)
        else:
            month = str(j)
        filename = f"{EVAPdir}EVAP.{year}{month}.nc"
        data = nc.Dataset(filename) #read EVAP file
        evap_month = data.variables['e'][:,:] #[m/day]
        evaps[(i-styear)*12+j-1,:,:] = evap_month

print(i)

