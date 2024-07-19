#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 09:10:35 2024

@author: sschanz
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl

cmap = mpl.cm.viridis
bounds = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

"""
This script examines the sediment depth as a function of mean glacial extent and plots based on the glacial cycle
"""

## CONSTANTS
                                                                                                                                         
time = np.arange(1, 801, 1)*1000
time2 = np.array((100, 200, 300, 400, 500, 600, 700, 800))*1000
x = np.arange(1, 200, 2) * 1000

## ----------  LOOP THROUGH AND GET DATA ####

sed_us = np.zeros((45, 800))
sed_ds = np.zeros((45, 800))
average_glacial_extent = np.zeros(45)
stdev_glacial_extent = np.zeros(45)

uplift_folders = ('0.001', '0.0005', '0.0001')
D0_folders = ('0.5', '0.25', '0.1', '0.05', '0.01')
a_folders = ('0.0001', '5e-05', '1e-05')
multiplier = np.ones(45)

os.chdir('..') # get out of figure_scripts folder
os.chdir('data') # get into data


i = 0

d_glacial = np.zeros(45)
d_proglacial = np.zeros(45)

for u_folder in uplift_folders:
   
    for d_folder in D0_folders:
        
        for a_folder in a_folders:
            
            filename = ('%s_%s_%s_True' % (u_folder, d_folder, a_folder))
            
            if u_folder == '0.001':
                if d_folder == '0.5':
                    multiplier[i] = 10
                elif d_folder == '0.25':
                    multiplier[i] = 10
            elif u_folder == '0.0005':
                if d_folder == '0.5':
                    if a_folder == '0.0001':
                        multiplier[i] = 10
                    elif a_folder == '5e-5':
                        multiplier[i] = 10
                elif d_folder == '0.25':
                    if a_folder == '0.0001':
                        multiplier[i] = 10
                    elif a_folder == '5e-5':
                        multipler[i] = 10
            
            
            glacial_extent = np.zeros(8)
            
            dict1D_variables = pd.read_pickle("./%s.pkl" % filename)
            locals().update(dict1D_variables)
            
            k = 0
            j = 0
            
            sed_size_0 = float(d_folder)
            sed_size_proglacial = float(d_folder)*np.exp(-float(a_folder)*50000)
            
            for year2 in time2:
                ice = HICE_save[year2]
                glacial_extent[j] = np.where(ice==0)[0][0]*2
                j += 1   
            
            average_extent = np.mean(glacial_extent)
            stdev_extent = np.std(glacial_extent)
            
            ### FIND THE GLACIAL NODE
            glacial_km = average_extent * 0.5
            glacial_node = int(glacial_km / 2)
            d_glacial[i] = float(d_folder)*np.exp(-glacial_km*1000*float(a_folder))
            
            ### FIND THE PROGLACIAL NODE
            proglacial_km = average_extent * 2
            if proglacial_km > 198:
                proglacial_km = 198
            proglacial_node = int(proglacial_km / 2)
            d_proglacial[i] = float(d_folder)*np.exp(-proglacial_km*1000*float(a_folder))
            
            for year in time:
                sed_us[i,k] = sed_depth_save[year][glacial_node] #multiplier[i]
                sed_ds[i,k] = sed_depth_save[year][proglacial_node] #multiplier[i]
                k += 1
            
            average_glacial_extent[i] = average_extent
            stdev_glacial_extent[i] = stdev_extent
            i += 1
            
fig, ax = plt.subplots(2,1,figsize = (6,6))
fig2, ax2 = plt.subplots(2,1, figsize = (6,6))

for i in range(45):
    
    lcolor = mpl.cm.viridis(norm(d_glacial[i]))
    #lcolor2 = mpl.cm.viridis(norm(d_proglacial[i]))
    alphavalue = 0.75
    
    if average_glacial_extent[i] == 0:
        lcolor = 'goldenrod'
    
    if d_glacial[i] > 0.095:
        stack = 0
    else:
        stack = 1
    # normalized sediment depth within each glacial cycle
    # for j in np.arange(1,9):
    #     ll = (j-1)*100
    #     ul = j*100
    #     sed_ds_norm = sed_ds[i,ll:ul]/np.max(sed_ds[i,ll:ul])
    #     sed_us_norm = sed_us[i,ll:ul]/np.max(sed_us[i,ll:ul])
    #     tplot = np.arange(1, 101, 1)
        
        # stack them by glacial cycle?
        
    ax[stack].plot(time/1000, sed_us[i,:], lw = 1.5, color = lcolor, alpha = alphavalue)
    ax2[stack].plot(time/1000, sed_ds[i,:], lw = 1.5, color = lcolor, alpha = alphavalue)  

  
    
### ----- CUSTOMIZE PLOTS ---- ##
os.chdir('..') # get out of data folder
os.chdir('figure_scripts')
os.chdir('figures')

# for k in np.array((0, 0.5, 1, 2, 4, 4.5)):
#     ax[0].hlines(0+k*1, 0, 800, ls = '--', color = 'grey', )
#     ax[1].hlines(0 + k*50, 0, 800, ls = '--', color = 'grey')

for axis in (ax, ax2):
    for i in range(2):    
        axis[i].tick_params(width = 1, direction = 'inout', length = 8)
        plt.setp(axis[i].spines.values(), linewidth = 1)
        axis[i].set_ylabel('Sediment depth (m)', weight = 'bold', fontsize = 14)
    axis[1].set_xlabel('Model time (ky)', weight = 'bold', fontsize = 14)

sm = cm.ScalarMappable(cmap = 'viridis', norm = mpl.colors.BoundaryNorm(bounds, cmap.N)) 
fig.colorbar(sm, ax = ax[1], orientation = 'horizontal', label = 'Median grain size (m)')
fig2.colorbar(sm, ax = ax2[1], orientation = 'horizontal', label = 'Median grain size (m)')

# plt.savefig('sed_depth_pg.svg')
# plt.close(fig2)
# plt.savefig('sed_depth_gl.svg')
# plt.close(fig)

# plt.savefig('sed_depth_2.svg')
# np.savetxt('d_glacial.txt', d_glacial)
# np.savetxt('d_proglacial.txt', d_proglacial)