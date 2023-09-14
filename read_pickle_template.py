#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 08:21:18 2023

@author: sschanz
"""
"""
Test unpickling of scripts and provide basic script to use the pickle files
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


#### ---- LOAD THE PICKLE FILE ---- #####
# Pickle file is a dictionary of dictionaries:
    # first level are the different variables to save, user-specified
    # second level is the variable array for each saved timestep

dict1D_variables = pd.read_pickle("./oneDvariables.pkl")

#### ----- EXTRACT VARIABLES ----- #####
# Each 1st level dict entry as it's own variable
locals().update(dict1D_variables)

# Variables are now saved as pandas Series, where the label for each entry corresponds to the timestep they were saved at

# example accessing information - valley width at 88000 years in the model:
Wv_at_88000 = Wv_save[88000]

#### ------ DO WHATEVER ANALYSIS YOU WANT! ----- ###

# example: I have Wv, dzb, eg, sed_depth, and z saved every 1000 years on a 2000m grid

time = np.arange(1, 201)*1000
x = np.arange(0, 200, 2)*1000

## ---- FIGURE 1: plot long profile over time:
fig, ax = plt.subplots(1,1)
ax.set_ylabel('Elevation (m)')
ax.set_xlabel('Distance downstream (km)')

# colors
color_bedrock = plt.cm.Greys(np.linspace(0,1,len(time)))
color_ice = plt.cm.Blues(np.linspace(0,1,len(time)))
color_sed = plt.cm.Oranges(np.linspace(0,1,len(time)))

i = 0
for ts in time:
    ax.plot(x/1000, z_save[ts]+sed_depth_save[ts]+HICE_save[ts], color = color_ice[i], lw = 1)
    ax.plot(x/1000, z_save[ts]+sed_depth_save[ts], color = color_sed[i], lw = 1)
    ax.plot(x/1000, z_save[ts], color = color_bedrock[i], lw = 2)
    i += 1

# create colorbar legend for alpha/time
cmap = mpl.cm.Greys
norm = mpl.colors.Normalize(vmin = time[0], vmax = time[-1])
fig.colorbar(plt.cm.ScalarMappable(norm = norm, cmap = cmap), ax = ax, label = 'Time')

# create legend for ice, sed, and bedrock    



## ---- FIGURE 2: change in valley width over time
fig2, ax2 = plt.subplots(1,1)
ax2.set_ylabel('Valley width (m)')
ax2.set_xlabel('Distance downstream (km)')

i = 0
color_of_lines = plt.cm.viridis(np.linspace(0,1,len(time)))
for ts in time:
    ax2.plot(x/1000, Wv_save[ts], lw = 1.5, color = color_of_lines[i]) 
    i += 1
# create colorbar legend for time/linecolor    
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin = time[0], vmax = time[-1])
fig2.colorbar(plt.cm.ScalarMappable(norm = norm, cmap = cmap), ax = ax2, label = 'Time')


## ---- FIGURE 3: Erosion rates through time
fig3, ax3 = plt.subplots(1,1)
ax3.set_ylabel('Erosion rate (mm/yr)')
ax3.set_xlabel('Distance downstream (km)')

i = 0
for ts in time:
    erosion_rate_total = (dzbforsave_save[ts] + egforsave_save[ts]) # these variables are summed up over the save interval, so divide by 1000 to get in m/yr, then multiply by 1000 to get mm/yr
    ax3.plot(x/1000, erosion_rate_total, color = color_bedrock[i], lw = 1.5, ls = '--')
    ax3.plot(x/1000, egforsave_save[ts], color = color_ice[i], lw = 1.5) 
    ax3.plot(x/1000, dzbforsave_save[ts], color = color_sed[i], lw = 1.5)
    i += 1

# create colorbar legend for time/linecolor    
cmap = mpl.cm.Greys
norm = mpl.colors.Normalize(vmin = time[0], vmax = time[-1])
fig3.colorbar(plt.cm.ScalarMappable(norm = norm, cmap = cmap), ax = ax3, label = 'Time')
