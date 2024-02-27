#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 17:09:55 2023

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
fig2, ax2 = plt.subplots(1,1, figsize = (4,4))

linecolors = cm.viridis(np.linspace(0, 1, 100))

### ---- LOAD DATA ---- ###
os.chdir('Fig7_sed_flux_by_glacial_extent_data')

average_extent = np.genfromtxt('average_extent_save.txt')
sed_flux_unglac = np.genfromtxt('flux_out_unglac_save.txt')[:,::-1]
sed_flux_glac = np.genfromtxt('flux_out_glac_save.txt')[:,::-1]

os.chdir('..')
os.chdir('figures')

### --- AVERAGING ---- ###

avg_flux_gl = np.zeros((3*5*3, 800))
avg_flux_ug = np.zeros((3*5*3, 800))

avg_flux_gl_normalized = np.zeros((3*5*3, 800))
flux_gl_normalized = np.zeros((3*5*3, 800))

### --- NORMALIZE ---- ###

# normalization is: glaciated sed flux / unglaciated sed flux * 100% showing the percent change of steady state expected erosion rates. If normalized rate = 100%, then erosion matches steady state. If 200%, then it is 2x the steady state rate. If 5%, then it is 5% of the steady state rate.

for i in range(3*5*3):
    
    for k in range(800):
        year = time[k]
        flux_totals_gl = np.sum(sed_flux_glac[i,:][time <= year]) # meters
        flux_totals_ug = np.sum(sed_flux_unglac[i,:][time <= year]) # meters
        avg_flux_gl[i,k] = flux_totals_gl/year # m/yr
        avg_flux_ug[i,k] = flux_totals_ug/year # m/yr
    
    avg_flux_gl_normalized[i,:] = (avg_flux_gl[i,:]) / (avg_flux_ug[i,:]) * 100
    avg_flux_gl_normalized[i,:][avg_flux_ug[i,:] == 0] = np.nan

    flux_gl_normalized[i,:] = sed_flux_glac[i,:]/sed_flux_unglac[i,:] * 100
    flux_gl_normalized[i,:][sed_flux_unglac[i,:] == 0] = np.nan
    
    if np.isnan(sed_flux_glac[i,-1])==True:
        avg_flux_gl_normalized[i,:] = np.nan
        flux_gl_normalized[i,:] = np.nan
    
    if avg_flux_gl_normalized[i,-1] > 110:
        zord = 1
        
        # model iterations 35 and 38, or the 36th and 39th runs. 
        
        """
        uplift_folders = ('highU', 'midU', 'lowU')
D0_folders = ('highestD', 'higherD', 'highD', 'midD', 'lowD')
a_folders = ('highA', 'midA', 'lowA')
        So the 36th and 38th of these would be:
            1-15 is highU
            16-30 is midU
            31-45 is lowU
                31-33 is highest D where the lowA model is nan
                34-36 is higher D where 36 is lowA
                37-39 is high D where 39 is lowA
                40-42 but all are unglaciated
                
            So exceptions occur when uplift is low and attrition is low, giving large grain sizes for a fairly low profile.
        """
    else:
        zord = 20
        
        
    ### ---- ADD DATA TO PLOTS --- ###
    if average_extent[i] < 100:
        if average_extent[i] == 0:
            lcolor = 'goldenrod'
            alphavalue = 0.3
        else:
            lcolor = linecolors[int(average_extent[i])]
            alphavalue = 1

        ax.plot(time/1000, flux_gl_normalized[i,:], color = lcolor, lw = 2, alpha = alphavalue, zorder = zord)
        ax2.plot(time, avg_flux_gl_normalized[i,:], color = lcolor, lw = 2, alpha = alphavalue, zorder = zord)

### --- CUSTOMIZE PLOTS ---- ###
    
for axis in (ax, ax2):
    axis.set_yscale('log')
    axis.axhline(100, ls = '--', color = 'grey')
    # axis.set_xlabel('Averaging time (ky)', weight = 'bold', fontsize = 14)
    axis.set_ylabel('Normalized flux out \n (% of non-glacial)', weight = 'bold', fontsize = 14)
    axis.tick_params(width = 2, direction = 'inout', length = 8)
    plt.setp(axis.spines.values(), linewidth = 2)

ax2.set_xscale('log')    
sm = cm.ScalarMappable(cmap = 'viridis', norm = plt.Normalize(vmin = 0, vmax = 100)) 

ax2.set_xlabel('Averaging time', weight = 'bold', fontsize = 14)
ax.set_xlabel('Model time (ky)', weight = 'bold', fontsize = 14)
 
fig2.colorbar(sm, ax = ax, orientation = 'horizontal', label = 'Average glacial extent (km)')
fig.colorbar(sm, ax = ax2, orientation = 'horizontal', label = 'Average glacial extent (km)')

# tight layout at end:
fig.tight_layout()
fig2.tight_layout()

plt.savefig("average_rates_norm_flux_out.svg")
plt.close(fig2)
plt.savefig("norm_sedflux_time.svg")