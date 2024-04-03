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
fig, ax = plt.subplots(1,1, figsize = (6,6))
fig2, ax2 = plt.subplots(1,1, figsize = (6,6))


### ---- LOAD DATA ---- ###
working_dir = os.getcwd()
os.chdir('Fig5_erosion_rate_by_glacial_extent_data')

glacial_extent = np.genfromtxt('glacial_extent_save.txt')
uplift_rate = np.genfromtxt('uplift_rate_save.txt')
erosion_0km = np.genfromtxt('erosion_0km_save.txt')[:,::-1]
erosion_50km = np.genfromtxt('erosion_50km_save.txt')[:,::-1]

os.chdir('..') # gets out of Fig5 data
os.chdir('..') # gets out of figure_scripts
os.chdir('..') # gets out of January2024_revision

os.chdir('May2021_runs')
os.chdir('July2021_followup')
os.chdir('3Ma_glaciations')

average_extent_3ma = np.genfromtxt('avg_glacial_extent_save.txt')
erosion_0km_3ma = np.genfromtxt('erosion_0km_save.txt')[:,::-1]
erosion_50km_3ma = np.genfromtxt('erosion_50km_save.txt')[:,::-1]

os.chdir(working_dir) # back to the figure_scripts
os.chdir('figures')


### --- AVERAGING ---- ###
average_extent = np.mean(glacial_extent, axis = 1)
avg_erates_0km = np.zeros((3*5*3, 800))
avg_erates_50km = np.zeros((3*5*3, 800))

avg_erates_0km_normalized = np.zeros((3*5*3, 800))
avg_erates_50km_normalized = np.zeros((3*5*3, 800))

avg_erates_0km_3ma = np.zeros((5, 3000))
avg_erates_50km_3ma = np.zeros((5, 3000))

avg_erates_0km_normalized_3ma = np.zeros((5, 3000))
avg_erates_50km_normalized_3ma = np.zeros((5, 3000))

### ---- PLOT LOCATIONS --- ###
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
        
### --- NORMALIZE ---- ###

# normalization is: (erosion) / uplift * 100% showing the percent change of steady state expected erosion rates. If normalized rate = 100%, then erosion matches steady state. If 200%, then it is 2x the steady state rate. If 5%, then it is 5% of the steady state rate.

# since erosion is recorded over 1000 years and represents meters of erosion over 1000 years, the uplift rate also needs to be in meters of erosion over 1000 years (same as mm/yr). Loaded in m/yr though.

for i in range(5):
    for k in range(3000):
        year = time2[k]
        e_totals0_3ma = np.sum(erosion_0km_3ma[i,:][time2 <= year]) # meters
        e_totals50_3ma = np.sum(erosion_50km_3ma[i,:][time2 <= year]) # meters
        avg_erates_0km_3ma[i,k] = e_totals0_3ma/year # m/yr
        avg_erates_50km_3ma[i,k] = e_totals50_3ma/year # m/yr
        
    avg_erates_0km_normalized_3ma[i,:] = (avg_erates_0km_3ma[i,:]) - 0.0005
    avg_erates_50km_normalized_3ma[i,:] = (avg_erates_50km_3ma[i,:]) - 0.0005


    ### ---- ADD DATA TO PLOTS --- ###
    
    lcolor = mpl.cm.winter(norm(average_extent[i]))
    alphavalue = 1
    
    if average_extent[i] == 0:
        lcolor = 'goldenrod'
        alphavalue = 0.3    
    if average_extent[i] < 50:
        ax2.plot(time2, avg_erates_50km_normalized_3ma[i,:]+bump[i]*0.0015, color = lcolor, lw = 2, alpha = alphavalue)

    ax.plot(time2, avg_erates_0km_normalized_3ma[i,:]+bump[i]*0.015, color = lcolor, lw = 2, alpha = alphavalue)

for i in range(3*5*3):
    
    for k in range(800):
        year = time[k]
        e_totals0 = np.sum(erosion_0km[i,:][time <= year]) # meters
        e_totals50 = np.sum(erosion_50km[i,:][time <= year]) # meters
        avg_erates_0km[i,k] = e_totals0/year # m/yr
        avg_erates_50km[i,k] = e_totals50/year # m/yr
        
    avg_erates_0km_normalized[i,:] = (avg_erates_0km[i,:]) - (uplift_rate[i])
    avg_erates_50km_normalized[i,:] = (avg_erates_50km[i,:]) - uplift_rate[i]

    
    ### ---- ADD DATA TO PLOTS --- ###
    lcolor = mpl.cm.winter(norm(average_extent[i]))
    alphavalue = 1
    
    if average_extent[i] == 0:
        lcolor = 'goldenrod'
        alphavalue = 0.3    
    if average_extent[i] < 50:
        ax2.plot(time, avg_erates_50km_normalized[i,:]+bump[i]*0.0015, color = lcolor, lw = 2, alpha = alphavalue)

    ax.plot(time, avg_erates_0km_normalized[i,:]+bump[i]*0.015, color = lcolor, lw = 2, alpha = alphavalue)


### --- CUSTOMIZE PLOTS ---- ###
for k in range(6):
    if k > 0:
        k -= 0.5
    ax2.hlines(0+k*0.0015, 0, 900000, ls = '--', color = 'grey', )
    ax.hlines(0 + k*0.015, 0, 900000, ls = '--', color = 'grey')
    
for axis in (ax, ax2):
    axis.set_xscale('log')
    # axis.axhline(100, ls = '--', color = 'grey')
    axis.set_xlabel('Averaging time (ky)', weight = 'bold', fontsize = 14)
    axis.set_ylabel('Glacial - nonglacial erosion rate (m/yr)', weight = 'bold', fontsize = 14)
    axis.tick_params(width = 2, direction = 'inout', length = 8)
    plt.setp(axis.spines.values(), linewidth = 2)
    
sm = cm.ScalarMappable(cmap = 'winter', norm = mpl.colors.BoundaryNorm(bounds, cmap.N)) 

fig2.colorbar(sm, ax = ax, orientation = 'horizontal', label = 'Mean glacial extent (km)')
fig.colorbar(sm, ax = ax2, orientation = 'horizontal', label = 'Mean glacial extent (km)')

#ax2.set_ylim([0, 400])
# tight layout at end:
fig.tight_layout()
fig2.tight_layout() 

plt.savefig('average_erosion_rates_50km.svg')
plt.close(fig2)
plt.savefig('average_erosion_rates_0km.svg')