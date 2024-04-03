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
import matplotlib as mpl


plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

### ---- CONSTANTS ---- ###
time = np.arange(1, 801, 1)*1000
fig, ax = plt.subplots(1,1, figsize = (6,6))
fig2, ax2 = plt.subplots(1,1, figsize = (6,6))

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

average_extent = np.mean(glacial_extent, axis = 1).astype(int)

uplifts = np.tile(uplift_rate*1000, (800,1)).transpose()
erosion_0km_normalized = erosion_0km - uplifts
erosion_50km_normalized = erosion_50km - uplifts

cmap = mpl.cm.winter
bounds = [0, 5, 10, 25, 50, 100]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

bump = np.zeros(45)
for i in range(45):
    glacial_extent = average_extent[i]
    if glacial_extent >= bounds[4]:
        bump[i] = 4.5
    elif glacial_extent >= bounds[3]:
        bump[i] = 3.5
    elif glacial_extent >= bounds[2]:
        bump[i] = 2.5
    elif glacial_extent >= bounds[1]:
        bump[i] = 1.5
    elif glacial_extent > bounds[0]:
        bump[i] = 0.5
    else:
        bump[i] = 0
    
for i in range(3*5*3):

    ### ---- ADD DATA TO PLOTS --- ###

    lcolor = mpl.cm.winter(norm(average_extent[i]))
    alphavalue = 1
    
    if average_extent[i] == 0:
        lcolor = 'goldenrod'
        alphavalue = 0.3    
    if average_extent[i] < 50:
        ax2.plot(time/1000, erosion_50km_normalized[i,:]+bump[i]*2, color = lcolor, lw = 2, alpha = alphavalue)

    ax.plot(time/1000, erosion_0km_normalized[i,:]+bump[i]*30, color = lcolor, lw = 2, alpha = alphavalue)
        

### --- CUSTOMIZE PLOTS ---- ###

for k in range(6):
    if k > 0:
        k -= 0.5
    ax2.hlines(0+k*2, 0, 800, ls = '--', color = 'grey', )
    ax.hlines(0 + k*30, 0, 800, ls = '--', color = 'grey')
    
for axis in (ax, ax2):   
    axis.set_xlabel('Model time (ky)', weight = 'bold', fontsize = 14)
    axis.set_ylabel('Glacial - nonglacial erosion rate (mm/yr)', weight = 'bold', fontsize = 14) #Normalized erosion rate \n (% of steady state rate)', weight = 'bold', fontsize = 14)
    axis.tick_params(width = 2, direction = 'inout', length = 8)
    plt.setp(axis.spines.values(), linewidth = 2)
  

sm = cm.ScalarMappable(cmap = 'winter', norm = mpl.colors.BoundaryNorm(bounds, cmap.N)) 

fig2.colorbar(sm, ax = ax, orientation = 'horizontal', label = 'Average glacial extent (km)')
fig.colorbar(sm, ax = ax2, orientation = 'horizontal', label = 'Average glacial extent (km)')

# tight layout at end:
fig.tight_layout()
fig2.tight_layout() 

plt.savefig('norm_erosion_50km.svg')
plt.close(fig2)
plt.savefig('norm_erosion_0km.svg')
plt.close(fig)