#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:32:32 2023

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
time2 = np.arange(1, 3001, 1)*1000
fig, ax = plt.subplots(2,1, figsize = (6,6))
fig2, ax2 = plt.subplots(2,1, figsize = (6,6))


### ---- LOAD DATA ---- ###
working_dir = os.getcwd()
os.chdir('Fig5_erosion_rate_by_glacial_extent_data')

glacial_extent = np.genfromtxt('glacial_extent_save.txt')
uplift_rate = np.genfromtxt('uplift_rate_save.txt')
erosion_gl = np.genfromtxt('erosion_gl_save.txt')[:,::-1]
erosion_pg = np.genfromtxt('erosion_pg_save.txt')[:,::-1]
d_glacial = np.genfromtxt('d_glacial.txt')
d_proglacial = np.genfromtxt('d_proglacial.txt')

os.chdir('..') # gets out of Fig5 data
os.chdir('..') # gets out of figure_scripts
os.chdir('..') # gets out of January2024_revision

os.chdir('May2021_runs')
os.chdir('July2021_followup')
os.chdir('3Ma_glaciations')

average_extent_3ma = np.genfromtxt('avg_glacial_extent_save.txt')
erosion_0km_3ma = np.genfromtxt('erosion_gl_save.txt')[:,::-1]
erosion_50km_3ma = np.genfromtxt('erosion_pg_save.txt')[:,::-1]

os.chdir(working_dir) # back to the figure_scripts
os.chdir('figures')


### --- AVERAGING ---- ###
average_extent = np.mean(glacial_extent, axis = 1)
avg_erates_gl = np.zeros((3*5*3, 800))
avg_erates_pg = np.zeros((3*5*3, 800))

avg_erates_gl_normalized = np.zeros((3*5*3, 800))
avg_erates_pg_normalized = np.zeros((3*5*3, 800))

avg_erates_gl_3ma = np.zeros((5, 3000))
avg_erates_pg_3ma = np.zeros((5, 3000))

avg_erates_gl_normalized_3ma = np.zeros((5, 3000))
avg_erates_pg_normalized_3ma = np.zeros((5, 3000))

### ---- PLOT LOCATIONS --- ###
cmap = mpl.cm.viridis
bounds = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        
### --- NORMALIZE ---- ###

# normalization is: (erosion) / uplift * 100% showing the percent change of steady state expected erosion rates. If normalized rate = 100%, then erosion matches steady state. If 200%, then it is 2x the steady state rate. If 5%, then it is 5% of the steady state rate.

# since erosion is recorded over 1000 years and represents meters of erosion over 1000 years, the uplift rate also needs to be in meters of erosion over 1000 years (same as mm/yr). Loaded in m/yr though.
d0_3ma = np.array((0.01, 0.05, 0.1, 0.25, 0.5))
d_glacial_3ma = d0_3ma * np.exp(-0.01/1000*average_extent_3ma*1000)
for i in range(5):
    for k in range(3000):
        year = time2[k]
        e_totals0_3ma = np.sum(erosion_0km_3ma[i,:][time2 <= year]) # meters
        e_totals50_3ma = np.sum(erosion_50km_3ma[i,:][time2 <= year]) # meters
        avg_erates_gl_3ma[i,k] = e_totals0_3ma/year # m/yr
        avg_erates_pg_3ma[i,k] = e_totals50_3ma/year # m/yr
        
    avg_erates_gl_normalized_3ma[i,:] = (avg_erates_gl_3ma[i,:]) / 0.0005
    avg_erates_pg_normalized_3ma[i,:] = (avg_erates_pg_3ma[i,:]) / 0.0005


    ### ---- ADD DATA TO PLOTS --- ###
    
    lcolor = mpl.cm.viridis(norm(d_glacial_3ma[i]))
    alphavalue = 0.75
    
    if average_extent[i] == 0:
        lcolor = 'goldenrod'
        alphavalue = 0.3    

    if d_glacial_3ma[i] > 0.095:
        stack = 0
    else:
        stack = 1
    ax2[stack].plot(time2, avg_erates_pg_normalized_3ma[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)

    ax[stack].plot(time2, avg_erates_gl_normalized_3ma[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)

for i in range(3*5*3):
    
    for k in range(800):
        year = time[k]
        e_totals0 = np.sum(erosion_gl[i,:][time <= year]) # meters
        e_totals50 = np.sum(erosion_pg[i,:][time <= year]) # meters
        avg_erates_gl[i,k] = e_totals0/year # m/yr
        avg_erates_pg[i,k] = e_totals50/year # m/yr
        
    avg_erates_gl_normalized[i,:] = (avg_erates_gl[i,:]) / (uplift_rate[i])
    avg_erates_pg_normalized[i,:] = (avg_erates_pg[i,:]) / uplift_rate[i]

    
    ### ---- ADD DATA TO PLOTS --- ###
    lcolor = mpl.cm.viridis(norm(d_glacial[i]))
    lcolor2 = mpl.cm.viridis(norm(d_proglacial[i]))
    alphavalue = 0.75
    
    if average_extent[i] == 0:
        lcolor = 'goldenrod'
        lcolor2 = 'goldenrod'
        alphavalue = 0.3    

    ## set the stacking
    if d_glacial[i] > 0.095:
        stack = 0
    else:
        stack = 1
        
    # if d_proglacial[i] > 0.10:
    #     stack2 = 0
    # else:
    #     stack2 = 1
        
    ax2[stack].plot(time, avg_erates_pg_normalized[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)

    ax[stack].plot(time, avg_erates_gl_normalized[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)


### --- CUSTOMIZE PLOTS ---- ###
  
for axis in (ax, ax2):
    for i in range(2):
        axis[i].set_xscale('log')
        axis[i].axhline(1, ls = '--', color = 'grey')
        axis[i].set_ylabel('Normalized erosion rate E*', weight = 'bold', fontsize = 14)
        axis[i].tick_params(width = 1, direction = 'inout', length = 8)
        plt.setp(axis[i].spines.values(), linewidth = 1)
    axis[1].set_xlabel('Averaging time (ky)', weight = 'bold', fontsize = 14)
sm = cm.ScalarMappable(cmap = 'viridis', norm = mpl.colors.BoundaryNorm(bounds, cmap.N)) 

fig2.colorbar(sm, ax = ax[1], orientation = 'horizontal', label = 'Median grain size (m)')
fig.colorbar(sm, ax = ax2[1], orientation = 'horizontal', label = 'Median grain size (m)')

#ax2.set_ylim([0, 400])
# tight layout at end:
fig.tight_layout()
fig2.tight_layout() 

plt.savefig('average_erosion_rates_pg.svg')
plt.close(fig2)
plt.savefig('average_erosion_rates_gl.svg')