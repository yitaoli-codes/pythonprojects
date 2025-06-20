# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 17:05:30 2023

@author: Mika
"""

#Load packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
#from netCDF4 import Dataset
import matplotlib as matplotlib
import argparse
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import scipy.optimize
from scipy import optimize
import os
import copy
#from func_brokenStick import broken_stick
#from func_brokenStick import BroSti_constright


# import data
data_dir = "/Users/Mika/Documents/Privat/Studium/Master/Lectures/Applied_Glaciology/Assignment_DATA-20231004/data_part1/" #TODO: change to your location
data_name = ("dem_fiescher_1980.asc", "dem_fiescher_2009.asc", "glacier_mask_fiescher_1980.asc", "glacier_mask_fiescher_2009.asc", "climate-reports-tables-homogenized_GRH.txt")

 # Digital Elevation Models (DEMs) of the glacier surface, providing elevation information on a regular grid; array for lat and lon
demfiescher1980 = np.loadtxt(data_dir+data_name[0], skiprows=6,) #skips first 6 rows where informations are contained
demfiescher2009 = np.loadtxt(data_dir+data_name[1], skiprows=6,)
# binary information telling whether a given grid cell is covered by ice (1) or not (0)
glacmask_fiescher1980 = np.loadtxt(data_dir+data_name[2], skiprows=6,)
glacmask_fiescher2009 = np.loadtxt(data_dir+data_name[3], skiprows=6,)


### -------------------------------------------------------------------
# =============================================================================
# # Question 1
# =============================================================================

changeinaltitude = glacmask_fiescher1980[:] #create empty array

totalsum = 0
for i in range(len(changeinaltitude[0])): # range for columns
    for k in range(len(changeinaltitude)): # range for rows
        if glacmask_fiescher1980[k,i] != 0: #array[row, column]
            changeinaltitude[k,i] = demfiescher1980[k,i]-demfiescher2009[k,i]
            totalsum = totalsum +changeinaltitude[k,i] # if there is a glacier, the difference between 1980 and 2009 is added to the total amount of icechange; 50*50 because cellsize is 50 -> I want it per cubic meter, therefore i add it for every square meter
            #TODO: is this right like this???
        else: changeinaltitude[k,i] = np.nan
    # assumption: if there was no ice anymore in 2009, it was just the glacier that reduced the elevation, so no erosion or rockfall when everything melted
#np.savetxt(data_dir+'massdifference_glacier_19802009.txt', changeinaltitude) #TODO: if you want to save the change in altitude from 1980 to 2009, unmark this line

glacmask_fiescher1980 = np.loadtxt(data_dir+data_name[2], skiprows=6,)
glacmask_fiescher2009 = np.loadtxt(data_dir+data_name[3], skiprows=6,)
gridcell1980 = gridcell2009 = 0 # to see the amount of gridcells
for i in range(len(glacmask_fiescher2009[0])): # range for columns
    for k in range(len(changeinaltitude)): # range for rows
        if glacmask_fiescher1980[k,i] == 1: #array[row, column]
            gridcell1980 = gridcell1980+1
        if glacmask_fiescher2009[k,i] == 1: #array[row, column]
            gridcell2009 = gridcell2009+1

# Calculate into water equivalent
icedensity = 850 #kg/m^3
waterdensity = 1000 #kg/m^3
icetowater = icedensity/waterdensity
numberofyears = 2009 - 1980
gridarea = 50*50 #m^2 
glaciergrids = np.mean([gridcell1980, gridcell2009])
meanarea = glaciergrids*gridarea #m^2

waterequiv_per_year = -totalsum/numberofyears/glaciergrids*icetowater # unit: m w.e./y
print("amount of m w.e. lost between 1980 and 2009",round(waterequiv_per_year*numberofyears, 2))


# =============================================================================
# #Question 2 - Data from Météo suisse Grimsel Hospiz
# =============================================================================
weatherdata = copy.copy(np.loadtxt(data_dir+data_name[4], skiprows=28 ))



# =============================================================================
# #Question 3
# =============================================================================


#mean glacier altitude method

# #initialising variables 
# glacaltitude = glacmask_fiescher1980[:]
# meanaltitude = 0
# gridnum = 0
# for i in range(len(glacaltitude[0])): # range for columns
#     for k in range(len(glacaltitude)): # range for rows
#         if glacmask_fiescher1980[k,i] != 0 and glacmask_fiescher2009[k,i] != 0: #array[row, column]
#             glacaltitude[k,i] = (demfiescher1980[k,i]+demfiescher2009[k,i])/2
#             meanaltitude+= (demfiescher1980[k,i]+demfiescher2009[k,i])/2
#             gridnum += 1
#         else: glacaltitude[k,i] = np.nan
# meanaltitude = meanaltitude/gridnum
#print("mean altitude of glacier is",meanaltitude)

# #mean temperature change
# altchange = meanaltitude-1980 #from data, in m
# tempchange = -altchange/1000*6.5

# #corrected temperatures at mean altitude
# for i in range(len(weatherdata)):
#     weatherdata[i,2]=weatherdata[i,2]+tempchange
    
# weatherdata8909 = copy.copy(weatherdata[585:933,:]) #get the data from oct 1980 to sept 2009

# A=0
# C=0
# for i in range(len(weatherdata8909[:,0])):
#     if weatherdata8909[i,2]>= 0:
#         A += weatherdata8909[i,2]
#     if weatherdata8909[i,2]<= 2:
#         C += weatherdata8909[i,3]/1000 #transform C from mm precipitation to m
# DDF=(C/numberofyears-waterequiv_per_year)*numberofyears/A/30.42*1000
# print(DDF)



###grid dependent DDF calibration
#more efficient altitude storage(?)
glacaltitudestorage = np.array([], dtype=np.float64)
for i in range(len(glacmask_fiescher2009[0])): # range for columns
    for k in range(len(glacmask_fiescher2009)): # range for rows
        if glacmask_fiescher1980[k,i] != 0: #array[row, column]
            glacaltitudestorage = np.append(glacaltitudestorage,demfiescher1980[k,i])

weatherdata = np.loadtxt(data_dir+data_name[4], skiprows=28 )
weatherdata8909 = copy.copy(weatherdata[585:933,:])
A = 0
C = 0
for i in range(len(weatherdata8909[:,0])):
    for k in range(len(glacaltitudestorage)): 
        Temp=weatherdata8909[i,2]-(glacaltitudestorage[k]-1980)/1000*6.5
        if Temp>= 0:
            A += Temp
        if Temp<= 2:
            C += weatherdata8909[i,3]/1000 #transform C from mm precipitation to m
DDF=(C/glaciergrids/numberofyears-waterequiv_per_year)*numberofyears*glaciergrids/A/30.42*1000


print("the calibrated DDF is",round(DDF,2),"mm w.e./C/d") 



# =============================================================================
# #Question 4
# =============================================================================

#get new data series
weatherdata3522 = copy.copy(weatherdata[45:1089,:])
 
B3522= np.array([],dtype=np.float64)
numberofyears=2022-1935
DDFfix=DDF*30.42/1000
A = 0
C = 0
for i in range(len(weatherdata3522)):
    for k in range(len(glacaltitudestorage)): 
        Temp=weatherdata3522[i,2]-(glacaltitudestorage[k]-1980)/1000*6.5
        if Temp>= 0:
            A += Temp*DDFfix
        if Temp<= 2:
            C += weatherdata3522[i,3]/1000 #transform C from mm precipitation to m
 
        
cumubalance = (C-A)/glaciergrids
averbalance = cumubalance/numberofyears
print("cumulative glacier mass balance is ",round(cumubalance, 2),"m w.e.")
print("average glacier mass balance is ",round(averbalance, 3),"m w.e./y")


for i in range(len(weatherdata3522)):
    A = 0
    C = 0
    if weatherdata3522[i,2]>= 0:
        A += weatherdata3522[i,2]*DDFfix
    if weatherdata3522[i,2]<= 2:
        C += weatherdata3522[i,3]/1000 #transform C from mm precipitation to m
    B3522=np.append(B3522,(C-A))

reshaped_B3522 = B3522.reshape(numberofyears,-1)
summed_B3522 = np.sum(reshaped_B3522, axis=1)
summed_B3522 /= 12

plt.figure()
plt.plot(range(len(weatherdata3522)),B3522)

plt.figure()
plt.plot((list(set(weatherdata3522[0:-12,0]))),summed_B3522)
plt.xlabel('Year')
plt.ylabel('Balance [m w.e./year]')

plt.savefig(data_dir+'glacierbalance.png', dpi=300)

# total_sum = np.sum(B9222)
# print(total_sum)

# =============================================================================
# Question 5
# =============================================================================


#1 Volume of glacier in 2009
c=0.03
gamma = 1.36
Area09 = gridcell2009*gridarea/1000000 #km^2
Area80 = gridcell1980*gridarea/1000000 #km^2
Volume09 = c*(Area09**gamma)
Volume80 = c*(Area80**gamma)

print("volume of glacier in 2009 is ",round(Volume09,2),"km^3")


#2 time till half disapears
depth09 = Volume09/(Area09)*1000*icetowater #m
depth80 = Volume80/(Area80)*1000*icetowater #m

### mass balance since 1991:
weatherdata3522= copy.copy(weatherdata[729:1089,:])
 
B9222= np.array([],dtype=np.float64)
numberofyears=2022-1992
DDFfix=DDF*30.42/1000
A = 0
C = 0
for i in range(len(weatherdata3522)):
    for k in range(len(glacaltitudestorage)): 
        Temp=weatherdata3522[i,2]-(glacaltitudestorage[k]-1980)/1000*6.5
        if Temp>= 0:
            A += Temp*DDFfix
        if Temp<= 2:
            C += weatherdata3522[i,3]/1000 #transform C from mm precipitation to m
 
        
cumubalance = (C-A)/glaciergrids
averbalance = cumubalance/numberofyears
print("cumulative glacier mass balance is ",round(cumubalance, 2),"m w.e.")
print("average glacier mass balance is ",round(averbalance, 2),"m w.e./y")


for i in range(len(weatherdata3522)):
    A = 0
    C = 0
    if weatherdata3522[i,2]>= 0:
        A += weatherdata3522[i,2]*DDFfix
    if weatherdata3522[i,2]<= 2:
        C += weatherdata3522[i,3]/1000 #transform C from mm precipitation to m
    B9222=np.append(B9222,(C-A))


# volume is in ice, so need to multiply with icedensity
glacierdensity = 900
n_year =-(Volume09*10**9/2)/glacierdensity/(averbalance)/Area09*10**-6*waterdensity
print("Remaining years " + str(round(n_year,2)))
yearend = 2009+n_year
print("Year with half volume " + str(round(yearend)))


