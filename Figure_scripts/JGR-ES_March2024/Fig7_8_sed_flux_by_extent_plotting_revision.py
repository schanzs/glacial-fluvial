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
os.chdir('Fig7_sed_flux_by_glacial_extent_data')

average_extent = np.genfromtxt('average_extent_save.txt')
sed_flux_unglac = np.genfromtxt('flux_out_unglac_save.txt')[:,::-1]
sed_flux_glac = np.genfromtxt('flux_out_glac_save.txt')[:,::-1]

os.chdir('..')

os.chdir('Fig5_erosion_rate_by_glacial_extent_data')
d_glacial = np.genfromtxt('d_glacial.txt')
d_proglacial = np.genfromtxt('d_proglacial.txt')

os.chdir('..') # gets out of Fig5 data
os.chdir('..') # gets out of figure_scripts
os.chdir('..') # gets out of January2024_revision

os.chdir('May2021_runs')
os.chdir('July2021_followup')
os.chdir('3Ma_glaciations')

average_extent_3ma = np.genfromtxt('avg_glacial_extent_save.txt')
sed_flux_unglac_3ma = np.genfromtxt('sedflux_out_ug_save.txt')[:,::-1]
sed_flux_glac_3ma = np.genfromtxt('sedflux_out_save.txt')[:,::-1]

os.chdir(working_dir) # back to the figure_scripts
os.chdir('figures')

### --- AVERAGING ---- ###

avg_flux_gl = np.zeros((3*5*3, 800))
avg_flux_ug = np.zeros((3*5*3, 800))

avg_flux_gl_normalized = np.zeros((3*5*3, 800))
flux_gl_normalized = np.zeros((3*5*3, 800))

avg_flux_gl_3ma = np.zeros((5, 3000))
avg_flux_ug_3ma = np.zeros((5, 3000))

avg_flux_gl_normalized_3ma = np.zeros((5, 3000))
avg_flux_ug_normalized_3ma = np.zeros((5, 3000))

### --- COLORS ---- ###
cmap = mpl.cm.viridis
bounds = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


        
### --- NORMALIZE ---- ###

# normalization is: glaciated sed flux / unglaciated sed flux * 100% showing the percent change of steady state expected erosion rates. If normalized rate = 100%, then erosion matches steady state. If 200%, then it is 2x the steady state rate. If 5%, then it is 5% of the steady state rate.

# d0_3ma = np.array((0.01, 0.05, 0.1, 0.25, 0.5))
# d_glacial_3ma = d0_3ma * np.exp(-0.01/1000*average_extent_3ma*1000)
# for i in range(5):
#     for k in range(3000):
#         year = time2[k]
#         flux_totals_gl_3ma = np.sum(sed_flux_glac_3ma[i,:][time2 <= year]) # meters
#         flux_totals_ug_3ma = np.sum(sed_flux_unglac_3ma[i,:][time2 <= year]) # meters
#         avg_flux_gl_3ma[i,k] = flux_totals_gl_3ma/year # m/yr
#         avg_flux_ug_3ma[i,k] = flux_totals_ug_3ma/year # m/yr
        
#     avg_flux_gl_normalized_3ma[i,:] = (avg_flux_gl_3ma[i,:]) / avg_flux_ug_3ma[i,:]


#     ### ---- ADD DATA TO PLOTS --- ###
    
#     lcolor = mpl.cm.viridis(norm(d_glacial_3ma[i]))
#     alphavalue = 0.75
    
#     if average_extent_3ma[i] == 0:
#         lcolor = 'goldenrod'
#         alphavalue = 0.3    

#     if d_glacial_3ma[i] > 0.095:
#         stack = 0
#     else:
#         stack = 1
#     ax2[stack].plot(time2/1000, avg_flux_gl_normalized_3ma[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)

#     ax[stack].plot(time2/1000, avg_flux_gl_normalized_3ma[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)
    
for i in range(3*5*3):
    lcolor = mpl.cm.viridis(norm(d_glacial[i]))
    alphavalue = 0.75

    # if average_extent[i] == 0:
    #     lcolor = 'goldenrod'
        
    for k in range(800):
        year = time[k]
        flux_totals_gl = np.sum(sed_flux_glac[i,:][time <= year]) # meters
        flux_totals_ug = np.sum(sed_flux_unglac[i,:][time <= year]) # meters
        avg_flux_gl[i,k] = flux_totals_gl/year # m/yr
        avg_flux_ug[i,k] = flux_totals_ug/year # m/yr
    
    avg_flux_gl_normalized[i,:] = (avg_flux_gl[i,:]) / avg_flux_ug[i,:]
  

    flux_gl_normalized[i,:] = (sed_flux_glac[i,:]) / sed_flux_unglac[i,:]

    
    if np.isnan(sed_flux_glac[i,-1])==True:
        avg_flux_gl_normalized[i,:] = np.nan
        flux_gl_normalized[i,:] = np.nan
    
#     if avg_flux_gl_normalized[i,-1] > 110:
#         zord = 1
        
#         # model iterations 35 and 38, or the 36th and 39th runs. 
        
#         """
#         uplift_folders = ('highU', 'midU', 'lowU')
# D0_folders = ('highestD', 'higherD', 'highD', 'midD', 'lowD')
# a_folders = ('highA', 'midA', 'lowA')
#         So the 36th and 38th of these would be:
#             0-14 is highU
#             15-29 is midU
#             30-44 is lowU
#                 30-32 is highest D where the lowA model is nan
#                 33-35 is higher D where 36 is lowA
#                 36-38 is high D where 38 is lowA
#                 39-44 but all are unglaciated
                
#             So exceptions occur when uplift is low and attrition is low, giving large grain sizes for a fairly low profile.
#         """
#     else:
#         zord = 20

    ### ---- STACK BY GRAIN SIZE -----
    if d_glacial[i] > 0.095:
        stack = 0
    else:
        stack = 1
        
    ax2[stack].plot(time/1000, avg_flux_gl_normalized[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)

    ax[stack].plot(time/1000, flux_gl_normalized[i,:], color = lcolor, lw = 1.5, alpha = alphavalue)
    
### --- CUSTOMIZE PLOTS ---- ###

    
for axis in (ax, ax2):
    # axis.set_yscale('log')
    # axis.axhline(100, ls = '--', color = 'grey')
    # axis.set_xlabel('Averaging time (ky)', weight = 'bold', fontsize = 14)
    for i in range(2):
        
        axis[i].tick_params(width = 1, direction = 'inout', length = 8)
        plt.setp(axis[i].spines.values(), linewidth = 1)
        axis[i].set_ylabel('Normalized flux, Qs*', weight = 'bold', fontsize = 14)
ax2[0].set_xscale('log') 
sm = cm.ScalarMappable(cmap = 'viridis', norm = mpl.colors.BoundaryNorm(bounds, cmap.N)) 

ax2[1].set_xlabel('Averaging time', weight = 'bold', fontsize = 14)
ax[1].set_xlabel('Model time (ky)', weight = 'bold', fontsize = 14)
 
fig2.colorbar(sm, ax = ax[1], orientation = 'horizontal', label = 'Median grain size (m)')
fig.colorbar(sm, ax = ax2[1], orientation = 'horizontal', label = 'Median grain size (m)')

# tight layout at end:
fig.tight_layout()
fig2.tight_layout()

plt.savefig("average_rates_norm_flux_out_revision.svg")
plt.close(fig2)
plt.savefig("norm_sedflux_time_revision.svg")