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
fig, ax = plt.subplots(2,1, figsize = (6,6))
fig2, ax2 = plt.subplots(2,1, figsize = (6,6))

### ---- LOAD DATA ---- ###
os.chdir('Fig5_erosion_rate_by_glacial_extent_data')

glacial_extent = np.genfromtxt('glacial_extent_save.txt')
uplift_rate = np.genfromtxt('uplift_rate_save.txt')
erosion_gl = np.genfromtxt('erosion_gl_save.txt')
erosion_pg = np.genfromtxt('erosion_pg_save.txt')

erosion_ug_gl = np.genfromtxt('erosion_ug_gl_save.txt')
erosion_ug_pg = np.genfromtxt('erosion_ug_pg_save.txt')


d_glacial = np.genfromtxt('d_glacial.txt')
d_proglacial = np.genfromtxt('d_proglacial.txt')

os.chdir('..')
os.chdir('figures')

### --- NORMALIZE ---- ###

# normalization is: simulations / fluvial-only simulations showing the percent change of steady state expected erosion rates. If normalized rate = 100%, then erosion matches steady state. If 200%, then it is 2x the steady state rate. If 5%, then it is 5% of the steady state rate.

# since erosion is recorded over 1000 years and represents meters of erosion over 1000 years, the uplift rate also needs to be in meters of erosion over 1000 years (same as mm/yr). Loaded in m/yr though.

erosion_gl_normalized = np.zeros((3*5*3, 800))
erosion_pg_normalized = np.zeros((3*5*3, 800))

average_extent = np.mean(glacial_extent, axis = 1).astype(int)

# uplifts = np.tile(uplift_rate*1000, (800,1)).transpose()
erosion_gl_normalized = erosion_gl/erosion_ug_gl
erosion_pg_normalized = erosion_pg/erosion_ug_pg

cmap = mpl.cm.viridis
bounds = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
for i in range(3*5*3):

    ### ---- ADD DATA TO PLOTS --- ###
    lcolor = mpl.cm.viridis(norm(d_glacial[i]))
    lcolor2 = mpl.cm.viridis(norm(d_proglacial[i]))
    alphavalue = 0.75
    
    if average_extent[i]==0:
        lcolor = 'goldenrod'
        lcolor2 = 'goldenrod'
    
    ### --- STACK BY GRAIN SIZE --- ###
    if d_glacial[i] > 0.095:
        stack = 0
    else:
        stack = 1
        
    if d_proglacial[i] > 0.10:
        stack2 = 0
    else:
        stack2 = 1
        
    
    ax2[stack].plot(time/1000, erosion_pg_normalized[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)
    ax[stack].plot(time/1000, erosion_gl_normalized[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)
        

### --- CUSTOMIZE PLOTS ---- ###

for axis in (ax, ax2):  
    for i in range(2):    
        axis[i].set_ylabel('Normalized erosion, E*', weight = 'bold', fontsize = 14) #Normalized erosion rate \n (% of steady state rate)', weight = 'bold', fontsize = 14)
        axis[i].tick_params(width = 1, direction = 'inout', length = 8)
        plt.setp(axis[i].spines.values(), linewidth = 1)
    axis[1].set_xlabel('Model time (ky)', weight = 'bold', fontsize = 14)
 

sm = cm.ScalarMappable(cmap = 'viridis', norm = mpl.colors.BoundaryNorm(bounds, cmap.N)) 

fig2.colorbar(sm, ax = ax[1], orientation = 'horizontal', label = 'Median grain size (m)')
fig.colorbar(sm, ax = ax2[1], orientation = 'horizontal', label = 'Median grain size (m)')

# tight layout at end:
fig.tight_layout()
fig2.tight_layout() 

plt.savefig('norm_erosion_pg.svg')
plt.close(fig2)
plt.savefig('norm_erosion_gl.svg')
plt.close(fig)