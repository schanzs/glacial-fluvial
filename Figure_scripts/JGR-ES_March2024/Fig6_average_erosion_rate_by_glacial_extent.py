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


plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

### ---- CONSTANTS ---- ###
time = np.arange(1, 801, 1)*1000
time2 = np.arange(1, 3001, 1)*1000
fig, ax = plt.subplots(1,1, figsize = (6,4))
fig2, ax2 = plt.subplots(1,1, figsize = (6,4))

linecolors = cm.viridis(np.linspace(0, 1, 50))

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
        
    avg_erates_0km_normalized_3ma[i,:] = (avg_erates_0km_3ma[i,:]) / (0.0005) * 100
    avg_erates_0km_normalized_3ma[i,:][avg_erates_0km_3ma[i,:] == 0] = np.nan
    avg_erates_50km_normalized_3ma[i,:] = (avg_erates_50km_3ma[i,:]) / 0.0005 * 100
    avg_erates_50km_normalized_3ma[i,:][avg_erates_50km_3ma[i,:] == 0] = np.nan

    ### ---- ADD DATA TO PLOTS --- ###
    if average_extent_3ma[i] < 50:
        if average_extent_3ma[i] == 0:
            lcolor = 'goldenrod'
            alphavalue = 0.3
        else:
            lcolor = linecolors[int(average_extent_3ma[i])]
            alphavalue = 1

        ax2.plot(time2, avg_erates_50km_normalized_3ma[i,:], color = lcolor, lw = 2, alpha = alphavalue)
    else:
        lcolor = 'yellow'
    ax.plot(time2, avg_erates_0km_normalized_3ma[i,:], color = lcolor, lw = 2, alpha = alphavalue)


for i in range(3*5*3):
    
    for k in range(800):
        year = time[k]
        e_totals0 = np.sum(erosion_0km[i,:][time <= year]) # meters
        e_totals50 = np.sum(erosion_50km[i,:][time <= year]) # meters
        avg_erates_0km[i,k] = e_totals0/year # m/yr
        avg_erates_50km[i,k] = e_totals50/year # m/yr
        
    avg_erates_0km_normalized[i,:] = (avg_erates_0km[i,:]) / (uplift_rate[i]) * 100
    avg_erates_0km_normalized[i,:][avg_erates_0km[i,:] == 0] = np.nan
    avg_erates_50km_normalized[i,:] = (avg_erates_50km[i,:]) / uplift_rate[i] * 100
    avg_erates_50km_normalized[i,:][avg_erates_50km[i,:] == 0] = np.nan
    
    ### ---- ADD DATA TO PLOTS --- ###
    if average_extent[i] < 50:
        if average_extent[i] == 0:
            lcolor = 'goldenrod'
            alphavalue = 0.3
        else:
            lcolor = linecolors[int(average_extent[i])]
            alphavalue = 1

        ax2.plot(time, avg_erates_50km_normalized[i,:], color = lcolor, lw = 2, alpha = alphavalue)
    else:
        lcolor = 'yellow'
    ax.plot(time, avg_erates_0km_normalized[i,:], color = lcolor, lw = 2, alpha = alphavalue)


### --- CUSTOMIZE PLOTS ---- ###
    
for axis in (ax, ax2):
    axis.set_yscale('log')
    axis.set_xscale('log')
    axis.axhline(100, ls = '--', color = 'grey')
    axis.set_xlabel('Averaging time (ky)', weight = 'bold', fontsize = 14)
    axis.set_ylabel('Normalized erosion rate \n (% of steady state rate)', weight = 'bold', fontsize = 14)
    axis.tick_params(width = 2, direction = 'inout', length = 8)
    plt.setp(axis.spines.values(), linewidth = 2)
    
sm = cm.ScalarMappable(cmap = 'viridis', norm = plt.Normalize(vmin = 0, vmax = 50)) 

fig2.colorbar(sm, ax = ax, orientation = 'horizontal', label = 'Average glacial extent (km)')
fig.colorbar(sm, ax = ax2, orientation = 'horizontal', label = 'Average glacial extent (km)')

ax2.set_ylim([0, 400])
# tight layout at end:
fig.tight_layout()
fig2.tight_layout() 

plt.savefig('average_erosion_rates_50km.svg')
plt.close(fig2)
plt.savefig('average_erosion_rates_0km.svg')