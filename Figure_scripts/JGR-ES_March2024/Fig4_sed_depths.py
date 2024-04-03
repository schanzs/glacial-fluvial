#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:17:57 2023

@author: sschanz
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl

cmap = mpl.cm.winter
bounds = [0, 5, 10, 25, 50, 100]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

"""

1) Sediment depth versus time at two locations. Headwaters and a proglacial zone at 50 km.

                                                                                                                                      """
## CONSTANTS
                                                                                                                                         
time = np.arange(1, 801, 1)*1000
x = np.arange(1, 200, 2) * 1000

## ----------  LOOP THROUGH AND GET DATA ####

sed_us = np.zeros((45, 800))
sed_ds = np.zeros((45, 800))
average_glacial_extent = np.zeros(45)

uplift_folders = ('0.001', '0.0005', '0.0001')
D0_folders = ('0.5', '0.25', '0.1', '0.05', '0.01')
a_folders = ('0.0001', '5e-05', '1e-05')

os.chdir('..') # get out of figure_scripts folder
os.chdir('data') # get into data


i = 0
for u_folder in uplift_folders:
   
    for d_folder in D0_folders:
        
        for a_folder in a_folders:
            
            filename = ('%s_%s_%s_True' % (u_folder, d_folder, a_folder))
            
            glacial_extent = np.zeros(8)
            
            dict1D_variables = pd.read_pickle("./%s.pkl" % filename)
            locals().update(dict1D_variables)
            
            k = 0
            j = 0
            
            for year in time:
                sed_us[i,k] = sed_depth_save[year][0]
                sed_ds[i,k] = sed_depth_save[year][25]
                
                if year%100000==0:
                    ice = HICE_save[year]
                    glacial_extent[j] = np.where(ice==0)[0][0]*2
                    j += 1
                k += 1
            
            average_glacial_extent[i] = np.mean(glacial_extent)
            i += 1
            
fig, ax = plt.subplots(1,2, figsize = (8,4))

bump = np.zeros(45)
for i in range(45):
    glacial_extent = average_glacial_extent[i]
    if glacial_extent >= bounds[4]:
        bump[i] = 4.5
    elif glacial_extent >= bounds[3]:
        bump[i] = 4
    elif glacial_extent >= bounds[2]:
        bump[i] = 2
    elif glacial_extent >= bounds[1]:
        bump[i] = 1
    elif glacial_extent > bounds[0]:
        bump[i] = 0.5
    else:
        bump[i] = 0
        

for i in range(45):
    
    lcolor = mpl.cm.winter(norm(average_glacial_extent[i]))
    alphavalue = 1
    
    if average_glacial_extent[i] < 50:
        if average_glacial_extent[i] == 0:
            lcolor = 'goldenrod'
            alphavalue = 0.3

        ax[1].plot(time/1000, sed_ds[i,:]+bump[i]*50, lw = 2, color = lcolor, alpha = alphavalue)

    ax[0].plot(time/1000, sed_us[i,:]+bump[i]*1, lw = 2, color = lcolor, alpha = alphavalue)
    
    
### ----- CUSTOMIZE PLOTS ---- ##
os.chdir('..') # get out of data folder
os.chdir('figure_scripts')
os.chdir('figures')

for k in np.array((0, 0.5, 1, 2, 4, 4.5)):
    ax[0].hlines(0+k*1, 0, 800, ls = '--', color = 'grey', )
    ax[1].hlines(0 + k*50, 0, 800, ls = '--', color = 'grey')
    
for i in range(2):
    ax[i].tick_params(width = 2, direction = 'inout', length = 8)
    plt.setp(ax[i].spines.values(), linewidth = 2)

ax[0].set_ylabel('Sediment depth (m)', weight = 'bold', fontsize = 14)
fig.text(0.5, 0.04, 'Time (ky)', fontsize = 14, weight = 'bold', ha = 'center')


sm = cm.ScalarMappable(cmap = 'winter', norm = mpl.colors.BoundaryNorm(bounds, cmap.N)) 
fig.colorbar(sm, ax = ax, location = 'right', label = 'Mean glacial extent (km)')

plt.savefig('sed_depth.svg')