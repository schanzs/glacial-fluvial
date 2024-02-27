#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:51:37 2023

@author: sschanz
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import cm


plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

### ---- CONSTANTS ---- ###
time = np.arange(1, 801, 1)*1000
fig, ax = plt.subplots(1,1, figsize = (6,4))
fig2, ax2 = plt.subplots(1,1, figsize = (6,4))

linecolors = cm.viridis(np.linspace(0, 1, 50))

### ---- LOAD DATA ---- ###
os.chdir('Fig5_erosion_rate_by_glacial_extent_data')

glacial_extent = np.genfromtxt('glacial_extent_save.txt')
uplift_rate = np.genfromtxt('uplift_rate_save.txt')
erosion_0km = np.genfromtxt('erosion_0km_save.txt')
erosion_50km = np.genfromtxt('erosion_50km_save.txt')

os.chdir('..')
os.chdir('figures')

### --- NORMALIZE ---- ###

# normalization is: (erosion) / uplift * 100% showing the percent change of steady state expected erosion rates. If normalized rate = 100%, then erosion matches steady state. If 200%, then it is 2x the steady state rate. If 5%, then it is 5% of the steady state rate.

# since erosion is recorded over 1000 years and represents meters of erosion over 1000 years, the uplift rate also needs to be in meters of erosion over 1000 years (same as mm/yr). Loaded in m/yr though.

erosion_0km_normalized = np.zeros((3*5*3, 800))
erosion_50km_normalized = np.zeros((3*5*3, 800))

average_extent = np.mean(glacial_extent, axis = 1)
for i in range(3*5*3):
    erosion_0km_normalized[i,:] = (erosion_0km[i,:]) / (uplift_rate[i]*1000) * 100
    erosion_0km_normalized[i,:][erosion_0km[i,:] == 0] = np.nan
    erosion_50km_normalized[i,:] = (erosion_50km[i,:] ) / (uplift_rate[i]*1000) * 100
    erosion_50km_normalized[i,:][erosion_50km[i,:] == 0] = np.nan
    
    ### ---- ADD DATA TO PLOTS --- ###
    if average_extent[i] < 50:
        
        if average_extent[i] == 0:
            lcolor = 'goldenrod'
            alphavalue = 0.3
        else: 
            lcolor = linecolors[int(average_extent[i])]
            alphavalue = 1
        ax2.plot(time/1000, erosion_50km_normalized[i,:], color = lcolor, lw = 2, alpha = alphavalue)
    else:
        lcolor = 'yellow'
    ax.plot(time/1000, erosion_0km_normalized[i,:], color = lcolor, lw = 2, alpha = alphavalue)
        


### --- CUSTOMIZE PLOTS ---- ###
    
for axis in (ax, ax2):
    axis.set_yscale('log')
    axis.axhline(100, ls = '--', color = 'grey')
    axis.set_xlabel('Model time (ky)', weight = 'bold', fontsize = 14)
    axis.set_ylabel('Normalized erosion rate \n (% of steady state rate)', weight = 'bold', fontsize = 14)
    axis.tick_params(width = 2, direction = 'inout', length = 8)
    plt.setp(axis.spines.values(), linewidth = 2)
    
sm = cm.ScalarMappable(cmap = 'viridis', norm = plt.Normalize(vmin = 0, vmax = 50)) 

fig2.colorbar(sm, ax = ax, orientation = 'horizontal', label = 'Average glacial extent (km)')
fig.colorbar(sm, ax = ax2, orientation = 'horizontal', label = 'Average glacial extent (km)')

ax2.set_ylim([0.05, 3500])
# tight layout at end:
fig.tight_layout()
fig2.tight_layout() 

plt.savefig('norm_erosion_50km.svg')
plt.close(fig2)
plt.savefig('norm_erosion_0km.svg')
plt.close(fig)