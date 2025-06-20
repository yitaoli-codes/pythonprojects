# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:30:02 2023

@author: Mika
"""

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
import math
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd




### -------------------------------------------------------------------
# =============================================================================
# # Question 1
# =============================================================================
x1 = 69577.37
y1 = 16321.42
z1 = 3850.29
x2 = 69488.03 
y2 = 16119.59 
z2 = 3840.78


x = x2-x1
y = y2-y1
z = z2-z1

d_horizontal = np.sqrt(x**2+y**2)

v_abs = np.sqrt(x**2+y**2) 

n=3
A = 2.4E-24
p_ice=916.7 #kg/m3
H = 536 #m
g = 9.81 #m/s2
tanalpha = abs(z/d_horizontal)

t_d = p_ice*g*H*tanalpha
print(f'Tao_d : {t_d*10**-3} kPa')


### version 1 --------------------------------------------
print(f'Absolute speed: {round(v_abs,2)} m/a')
v_dS = (2*A/(n+1))*(t_d)**n*H *365*24*3600
print(f'Surface deformation speed: {round(v_dS,2)} m/a')
#v_d = (n+1)/(n+2)*v_dS
#print(f'average deformation speed: {round(v_d,2)} m/a')
#v_sliding = v_abs - v_d
print()
#print('Result:')
#print(f'-> Sliding speed: {round(v_sliding,2)} m/a')
v_sliding = v_abs-v_dS
print('Result:')
print(f'-> Sliding speed: {round(v_sliding,2)} m/a')


### version 2-----------------------------------------------
b =1
Pw =  0#?   ###TODO: what is the value of subglacial water pressure? could not find this in script; asked chatGPT, but value seems wrong
k = 0.012
N = p_ice*g*H - Pw

v_b_sliding = k*t_d**n*N**-b # k*t_d**n*N**-b
v_b_sliding = v_b_sliding/(60*60*24*365)
print()
print('Result version 2:')
print(f'Sliding speed: {round(v_b_sliding,2)} m/a') #too high!!!



### ---------------------------
### Question 2
### ----------------------------
# Define the glacier velocity model
def glacier_velocity_model(t, v0, a, tc, m):
    return v0 + a * (tc - t) ** (-m)

data = pd.read_fwf(r"C:/Users/maike/OneDrive/Desktop/ETHZ/HeSe_23/Applied Glaciology/Assignment 2/maybeunstable_velocity_data.txt", skiprows=6, header=None, names=["time", "speed"])

data["time"] = data["time"] - 5586  # Starting time from 0 
t_data = data["time"]
v_data = data["speed"]

# Estimation of v0 during the time of stable velocity (until day 30)
v0 = data['speed'][:127].mean()

# Perform nonlinear least squares regression
p0 = [v0, 1, 5636, 1]  # initial values for v0, a, tc and m
popt, pcov = curve_fit(glacier_velocity_model, t_data, v_data, p0=p0, maxfev=10000)

# Extract the optimized parameters
opt_v0, opt_a, opt_tc, opt_m = popt

# Calculate the estimated time of failure (tc)
estimated_tc = t_data + ((v_data - opt_v0) / opt_a) ** (1 / (-opt_m))

# Plotting the estimated time of failure
plt.figure(dpi=300) #better quality of the figures
plt.scatter(t_data, v_data, label='given data points', color='green', s=3)

# Generate the regression line using the optimized parameters
regression_line = glacier_velocity_model(t_data, opt_v0, opt_a, opt_tc, opt_m)

# Plot the regression line
plt.plot(t_data, regression_line, label='regression line', color='red', linewidth=1)

# Add labels and title
plt.xlabel('time since beginning of the observed time span (in d)')
plt.ylabel('ice flow speed (in m/d)')
plt.title('glacier velocity datapoints with regression line')

# Show legend and plot
plt.legend()
plt.show()

# Output the estimated time of failure
print("Estimated Time of Failure (tc):", estimated_tc)
print("Optimized value for a:", opt_a)
print("Optimized value for m:", opt_m)
print("Optimized value for v0:", opt_v0)

#Calculate tc using t1 and v1 (random t and v combination from the dataset)
tc = 0 + ((0.653  - v0) / opt_a) ** (1 / (-opt_m)) 

#print results
print("Calculated value for tc:", tc)






### ---------------------------
### Question 3
### ----------------------------
# import data
data_dir = "/Users/Mika/Documents/Privat/Studium/Master/Lectures/Applied_Glaciology/Assignment2/" #TODO: change to your location
data_name = ("lake_dem_2004.asc", "lake_dem_2015.asc", "lake_dem_2023.asc")

 # Digital Elevation Models (DEMs) of the glavier lakes, providing elevation information on a regular grid; array for lat and lon
lake_dem_2004 = np.loadtxt(data_dir+data_name[0], skiprows=6,).reshape((53, 47)) #skips first 6 rows where informations are contained
lake_dem_2015 = np.loadtxt(data_dir+data_name[1], skiprows=6,).reshape((53, 47))
lake_dem_2023 = np.loadtxt(data_dir+data_name[2], skiprows=6,).reshape((53, 47))

glac_ext = copy.copy(lake_dem_2015-lake_dem_2023)

xll=1987.5
yll=9987.5

xspillway1=2100
yspillway1=10971
xspillway2=2150
yspillway2=10950

#distance from lower left corner. 
spillway_row1=math.floor((yspillway1-yll)/25)
spillway_col1=math.floor((xspillway1-xll)/25)

spillway_row2=math.floor((yspillway2-yll)/25)
spillway_col2=math.floor((xspillway2-xll)/25)


cutoff1=lake_dem_2004[-spillway_row1, spillway_col1]
Vol04=0
lake_dem_2004_mask = copy.copy(lake_dem_2004)
for i in range(len(lake_dem_2004[0])):  # Loop over columns
    for k in range(len(lake_dem_2004)):  # Loop over rows
        
        if lake_dem_2004[k, i] <= cutoff1:
            Vol04+=(cutoff1-lake_dem_2004[k, i])*25*25#m
            lake_dem_2004_mask[k, i] = 1
        else :
            lake_dem_2004_mask[k, i] = 0

cutoff2=lake_dem_2015[-spillway_row2, spillway_col2]
Vol15=0
lake_dem_2015_mask = copy.copy(lake_dem_2015)
for i in range(len(lake_dem_2015[0])):  # Loop over columns
    for k in range(len(lake_dem_2015)):  # Loop over rows
        
        if lake_dem_2015[k, i] <= cutoff2:
            Vol15+=(cutoff2-lake_dem_2015[k, i])*25*25#m
            lake_dem_2015_mask[k, i] = 1
        else :
            lake_dem_2015_mask[k, i] = 0
            

cutoff2=lake_dem_2023[-spillway_row2, spillway_col2]
Vol23=0
lake_dem_2023_mask = copy.copy(lake_dem_2023)
for i in range(len(lake_dem_2023[0])):  # Loop over columns
    for k in range(len(lake_dem_2023)):  # Loop over rows
        
        if lake_dem_2023[k, i] <= cutoff2 :
            Vol23+=(cutoff2-lake_dem_2023[k, i])*25*25#m
            lake_dem_2023_mask[k, i] = 1
        else :
            lake_dem_2023_mask[k, i] = 0

glac_ext_mask = copy.copy(glac_ext)
for i in range(len(glac_ext[0])):  # Loop over columns
    for k in range(len(glac_ext)):  # Loop over rows
        
        if glac_ext[k, i] >=10:
            glac_ext_mask[k, i] = 1
        else :
            glac_ext_mask[k, i] = 0
            glac_ext[k,i] = 0

### -------------
# year 2004
data=copy.copy(lake_dem_2004)
# Create a meshgrid for X and Y coordinates
x_llcorner = 1987.5
y_llcorner = 9987.5
cellsize = 25.0
x = np.arange(x_llcorner, x_llcorner + cellsize * data.shape[1], cellsize)
y = np.arange(y_llcorner, y_llcorner + cellsize * data.shape[0], cellsize)
X, Y = np.meshgrid(x, y)
'''
# Plot contour lines
data = np.flipud(data) #to correct it from starting from the bottom to the top
plt.figure(dpi=300)
plt.contour(X, Y, data, cmap="coolwarm")
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('Contour plot 2004')
plt.show()
'''
plt.figure(dpi=300)
heatmap = plt.imshow(data, cmap=cm.coolwarm, extent=(X.min(), X.max(), Y.min(), Y.max()))
plt.colorbar()#label='Data Values')
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('2004')


plt.figure(dpi=300)
plt.contour(X, Y, np.flipud(glac_ext_mask), cmap="Blues")
plt.contourf(X, Y, np.flipud(lake_dem_2004_mask), cmap="binary")
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('Glacier extent and Lake Basin 2004')
plt.show()


### -------------
# year 2015
data=copy.copy(lake_dem_2015)
# Create a meshgrid for X and Y coordinates
x_llcorner = 1987.5
y_llcorner = 9987.5
cellsize = 25.0
x = np.arange(x_llcorner, x_llcorner + cellsize * data.shape[1], cellsize)
y = np.arange(y_llcorner, y_llcorner + cellsize * data.shape[0], cellsize)
X, Y = np.meshgrid(x, y)

# Plot contour lines
'''
data = np.flipud(data) #to correct it from starting from the bottom to the top
plt.figure(dpi=300)
plt.contour(X, Y, data, cmap="coolwarm")
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('Contour plot 2015')
plt.show()
'''
plt.figure(dpi=300)
heatmap = plt.imshow(data, cmap=cm.coolwarm, extent=(X.min(), X.max(), Y.min(), Y.max()))
plt.colorbar()#label='Data Values')
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('2015')


plt.figure(dpi=300)
plt.contour(X, Y, np.flipud(glac_ext_mask), cmap="Blues")
plt.contourf(X, Y, np.flipud(lake_dem_2015_mask), cmap="binary")
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('Glacier extent and Lake Basin 2015')
plt.show()

### -------------
# year 2023
data=copy.copy(lake_dem_2023)
# Create a meshgrid for X and Y coordinates
x_llcorner = 1987.5
y_llcorner = 9987.5
cellsize = 25.0
x = np.arange(x_llcorner, x_llcorner + cellsize * data.shape[1], cellsize)
y = np.arange(y_llcorner, y_llcorner + cellsize * data.shape[0], cellsize)
X, Y = np.meshgrid(x, y)
'''
# Plot contour lines
#data = np.flipud(data) #to correct it from starting from the bottom to the top
plt.figure(dpi=300)
plt.contour(X, Y, data, cmap="coolwarm")
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('Contour plot 2023')
plt.show()
'''
plt.figure(dpi=300)
heatmap = plt.imshow(data, cmap=cm.coolwarm, extent=(X.min(), X.max(), Y.min(), Y.max()))
plt.colorbar()#label='Data Values')
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('2023')


# plt.contour(X, Y, np.flipud(lake_dem_2015_mask))
# light_blue = '#87CEEB'  # This is a hex color code for light blue
# custom_cmap = ListedColormap([light_blue])
plt.figure(dpi=300)
plt.contour(X, Y, np.flipud(glac_ext_mask), cmap="Blues")
plt.contourf(X, Y, np.flipud(lake_dem_2023_mask), cmap="binary")
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('Glacier extent and Lake Basin 2023')
plt.show()



# =============================================================================
# Question 3.2
# =============================================================================

Q=np.array([50,58,64])
V=np.array([Vol04,Vol15,Vol23])

k1=Q/(V**(2/3))
print(k1)

Qlog=np.log(Q)
Vlog=np.log(V)

def error(logk):
    temp = (Qlog-2/3*Vlog-logk)**2
    result = sum(temp)
    return result

initial_params = [1e-3]
result = minimize(error, initial_params, method='BFGS')
k_CM = np.exp(result.x)
print(k_CM)



coefficients = np.polyfit(Vlog, Qlog, 1)
slope, intercept = coefficients

Qlog_pred1 = 2/3* Vlog + np.log(k_CM)
Qlog_pred2 = slope* Vlog + intercept

# Plot the data points
plt.figure(dpi=300)
plt.scatter(Vlog, Qlog, label='DEM Data Points')

# Plot the linear model
plt.plot(Vlog, Qlog_pred1, color='red', label='Linear Model with alpha=2/3')
plt.plot(Vlog, Qlog_pred2, color='lightcoral', label='Linear Model with alpha=0.52')

# Add labels and legend
plt.xlabel('log Volume of lake')
plt.ylabel('Log maximum Flow of spill')
plt.legend()

# Show the plot
plt.show()

# =============================================================================
# Question 3.3
# =============================================================================

elev_change_rate = copy.copy(glac_ext)/8
print("maximum elevation change rate",np.max(elev_change_rate),"m/year")
plt.figure(dpi=300)
#CS=plt.contour(X, Y, np.flipud(elev_change_rate), cmap=cm.coolwarm)
#plt.xlabel('easting')
#plt.ylabel('northing')
#plt.title('elevation change rate between 2015 and 2023')

plt.figure(dpi=300)
heatmap = plt.imshow(elev_change_rate, cmap=cm.coolwarm, extent=(X.min(), X.max(), Y.min(), Y.max()))
plt.colorbar(label='Elevation Change Rate (m/a)')
plt.xlabel('easting')
plt.ylabel('northing')
plt.title('Elevation Change Rate between 2015 and 2023')

# Optionally, add contour lines for reference
#contour_lines = plt.contour(X, Y, np.flipud(elev_change_rate), colors='black', linestyles='dashed')
#plt.clabel(contour_lines, inline=True, fontsize=8)


elev_change_inlake = 0
inlake_cellcount = 0
elev_change_sanslake = 0
sanslake_cellcount = 0
for i in range(len(glac_ext[0])):  # Loop over columns
    for k in range(len(glac_ext)):  # Loop over rows
        if lake_dem_2023_mask[k,i]== 1:
            elev_change_inlake+=glac_ext[k,i]/8
            inlake_cellcount+=1
        if glac_ext_mask[k,i] == 1 and lake_dem_2023_mask[k,i] == 0:
            elev_change_sanslake+=glac_ext[k,i]/8
            sanslake_cellcount+=1

elev_change_inlake/=inlake_cellcount
print("average lake basin lowering rate is ", elev_change_inlake)
elev_change_sanslake/=sanslake_cellcount
print("average glacier lowering rate is ", elev_change_sanslake)

lake_depth_change = elev_change_inlake-elev_change_sanslake
Vol_change_rate = np.count_nonzero(lake_dem_2023_mask == 1)*lake_depth_change*25*25 #Vol change per year in m3
print("average volume increase rate is",Vol_change_rate)

Vol24 = Vol23 + Vol_change_rate
Vol27 = Vol23 + 4*Vol_change_rate
Vol32 = Vol23 + 9*Vol_change_rate

Qmax24 = k_CM*(Vol24**(2/3))
Qmax27 = k_CM*(Vol27**(2/3))
Qmax32 = k_CM*(Vol32**(2/3))


