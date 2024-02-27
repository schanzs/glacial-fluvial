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
"""

fig, ax = plt.subplots(3,1, sharex = True, figsize = (6,10), layout = 'constrained')
# fig2, ax2 = plt.subplots(4,1, sharex = True, figsize = (6,10), layout = 'constrained')
time = np.arange(1, 801, 1)*1000
x = np.arange(0, 200, 2) * 1000
extent = np.zeros((4,8))
linecolors = cm.Spectral(np.linspace(0,1,800))
D0 = np.array((0.5, 0.25, 0.1))
label = ('A', 'B', 'C')

os.chdir('..') # get out of figures folder
os.chdir('data') # get into the data folder

i = 0

for grainsize in ('0.5', '0.25', '0.1'):
    
    # create the filename: midU, midA, glaciation on
    filename =  ('0.0005_%s_5e-05_True' % grainsize)
    
    dict1D_variables = pd.read_pickle("./%s.pkl" % filename)
    locals().update(dict1D_variables)
    
    k = 0
    j = 0
    for year in time:
        bed_z = z_save[year]
        if year%100000 == 0:
            ice = HICE_save[year]
            extent[i,j] = np.where(ice==0)[0][0]*2
            j += 1
        #plot in loop so I don't need to store data
        
        ax[i].plot(x/1000, bed_z, color = linecolors[k])
        k += 1
    # ax[i].plot(x/1000, bed_z_initial, color = 'black', lw = 2)    
    ax[i].text(0.95, 0.95, '%s) D0 = %s, avg extent = %s' % (label[i], D0[i], np.mean(extent, axis=1)[i]), transform = ax[i].transAxes, verticalalignment = 'top', ha = 'right')
    
    i += 1

os.chdir('..') # get out of data folder
os.chdir('figure_scripts')
os.chdir('figures')


### ------- CUSTOMIZE PLOTS ----------
for axis in ax:
    axis.tick_params(width = 2, direction = 'inout', length = 8)
    plt.setp(axis.spines.values(), linewidth = 2)
        
fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(0, 800), cmap='Spectral'),
             ax=ax[2], orientation='horizontal', label='Time (ky)')

ax[2].set_xlabel('Distance downstream (km)', weight = 'bold', fontsize = 14)
ax[1].set_ylabel('Elevation (m)', weight = 'bold', fontsize = 14)
ax[0].set_xlim([-4, 75])


# ### --- SAVE PLOTS -----
plt.savefig('longprofiles_bedrock.svg')