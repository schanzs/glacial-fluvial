#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 12:01:28 2024

@author: sschanz
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import os
import pandas as pd

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

"""
This script creates four subplots showing longitudinal profile development through time, and calls out the headwater and proglacial regions.

Playing around with normalization scheme where 1 = average glacial extent; below one is within the average glaciated zone and above 1 is in the proglacial zone.
"""

fig, ax = plt.subplots(3,2, sharex = True, figsize = (10,6), layout = 'constrained')
#fig2, ax2 = plt.subplots(3,1, sharex = True, figsize = (6,6), layout = 'constrained')
time = np.arange(1, 801, 1)*1000
time2 = np.array((100, 200, 300, 400, 500, 600, 700, 800))*1000
x = np.arange(0, 200, 2) * 1000
extent = np.zeros((3,8))
# linecolors = cm.terrain(np.linspace(0,1,800))
label = ('A', 'B', 'C')
d_glacial = np.zeros(3)

os.chdir('..') # get out of figures folder
os.chdir('data') # get into the data folder

i = 0

### --- COLORS ----
cmap = cm.cividis
bounds = [0, 100, 200, 300, 400, 500, 600, 700, 800]
norm = colors.BoundaryNorm(bounds, cmap.N)

for grainsize, attr in zip(('0.25', '0.1', '0.1'), ('1e-05', '1e-05', '5e-05')):
    
    # create the filename: midU, midA, glaciation on
    filename =  ('0.0005_%s_%s_True' % (grainsize, attr))
    
    dict1D_variables = pd.read_pickle("./%s.pkl" % filename)
    locals().update(dict1D_variables)
    
    j = 0
    for year2 in time2:
        ice = HICE_save[year2]
        extent[i,j] = np.where(ice==0)[0][0]*2
        j += 1
        #plot in loop so I don't need to store data
        
    mean_extent = np.mean(extent[i,:])
    glacial_km = 0.5 * mean_extent
    glacial_node = int(glacial_km / 2)
    d_glacial[i] = float(grainsize) * np.exp(-glacial_km*1000*float(attr))
    
    k = 0
    for year in time:
        valley_width = Wv_save[year]
        bed_z = z_save[year]
        lcolor = cm.cividis(norm(k))
        if year == time[-1]:
            lcolor = 'black'
        ax[i][0].plot(x/1000/mean_extent, bed_z, color = lcolor)
        ax[i][1].plot(x/1000/mean_extent, valley_width, color = lcolor)
        k += 1
    # ax[i].plot(x/1000, bed_z_initial, color = 'black', lw = 2)    
    ax[i][0].text(0.95, 0.95, '%s) D0 = %s, avg extent = %s' % (label[i], d_glacial[i], mean_extent), transform = ax[i][0].transAxes, verticalalignment = 'top', ha = 'right')
    
    i += 1

os.chdir('..') # get out of data folder
os.chdir('figure_scripts')
os.chdir('figures')


### ------- CUSTOMIZE PLOTS ----------
for axis in ax:
    for i in range(2):
        axis[i].set_xlim([0, 3])
        axis[i].tick_params(width = 1, direction = 'inout', length = 8)
    plt.setp(axis[i].spines.values(), linewidth = 1)

# for axis2 in ax2:
#     axis2.set_xlim([0,3])
#     axis2.tick_params(width = 1, direction = 'inout', length = 8)
#     plt.setp(axis2.spines.values(), linewidth = 1) 
    
fig.colorbar(cm.ScalarMappable(cmap = 'cividis', norm = colors.BoundaryNorm(bounds, cmap.N)),
             ax=ax[2], orientation='horizontal', label='Time (ky)')

ax[2,0].set_xlabel('Normalized distance, x*', weight = 'bold', fontsize = 14)
ax[2,1].set_xlabel('Normalized distance, x*', weight = 'bold', fontsize = 14)
ax[1,1].set_ylabel('Width (m)', weight = 'bold', fontsize = 14)
ax[1,0].set_ylabel('Elevation (m)', weight = 'bold', fontsize = 14)
# ax[0].set_xlim([-4, 75])


# ### --- SAVE PLOTS -----
plt.savefig('longprofiles_bedrock.svg')