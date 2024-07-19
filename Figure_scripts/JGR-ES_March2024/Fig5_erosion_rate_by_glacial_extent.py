#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 15:29:37 2023

@author: sschanz
"""

import numpy as np
import os
import pandas as pd


### ---- SET UP CONSTANTS --- ###
x = np.arange(0, 200, 2)*1000
time = np.arange(1, 801, 1)*1000
time2 = np.array((100, 200, 300, 400, 500, 600, 700, 800))*1000

### ----- SET UP SAVE VARIABLES ---- ###
glacial_extent_save = np.zeros((3*5*3, 8))
erosion_0km_save = np.zeros((3*5*3, 800))
erosion_50km_save = np.zeros((3*5*3, 800))
erosion_ug_0km_save = np.zeros((3*5*3, 800))
erosion_ug_50km_save = np.zeros((3*5*3, 800))
uplift_rate_save = np.zeros((3*5*3))

### ----- LOOP THROUGH ---- ###

uplift_folders = ('0.001', '0.0005', '0.0001')
D0_folders = ('0.5', '0.25', '0.1', '0.05', '0.01')
a_folders = ('0.0001', '5e-05', '1e-05')
multiplier = np.ones(45)
j = 0

os.chdir('..') # get out of the figure_scripts folder
os.chdir('data') # into the data folder


for u_folder in uplift_folders:

    for d_folder in D0_folders:

        for a_folder in a_folders:
            filename_ug = ('%s_%s_%s_False' % (u_folder, d_folder, a_folder))
            filename = ('%s_%s_%s_True' % (u_folder, d_folder, a_folder))
            
            if u_folder == '0.001':
                if d_folder == '0.5':
                    multiplier[j] = 10
                elif d_folder == '0.25':
                    multiplier[j] = 10
            elif u_folder == '0.0005':
                if d_folder == '0.5':
                    if a_folder == '0.0001':
                        multiplier[j] = 10
                    elif a_folder == '5e-5':
                        multiplier[j] = 10
                elif d_folder == '0.25':
                    if a_folder == '0.0001':
                        multiplier[j] = 10
                    elif a_folder == '5e-5':
                        multiplier[j] = 10
            
            
            dict1D_variables = pd.read_pickle("./%s.pkl" % filename)
            locals().update(dict1D_variables)                
            
            i = 0
            k = 0
            for year2 in time2:
                ice = HICE_save[year2]
                glacial_extent_save[j, k] = np.where(ice == 0)[0][0]*2
                k += 1
            
            average_extent = np.mean(glacial_extent_save[j,:])
            
            ### FIND THE GLACIAL NODE
            glacial_km = average_extent * 0.5
            glacial_node = int(glacial_km / 2)
            #d_glacial[i] = float(d_folder)*np.exp(-glacial_km*1000*float(a_folder))
            
            ### FIND THE PROGLACIAL NODE
            proglacial_km = average_extent * 2
            if proglacial_km > 198:
                proglacial_km = 198
            proglacial_node = int(proglacial_km / 2)
            #d_proglacial[i] = float(d_folder)*np.exp(-proglacial_km*1000*float(a_folder))
                 
            for year in time:
                erosion_0km_save[j, i] = egforsave_save[year][glacial_node]/multiplier[j] + dzbforsave_save[year][glacial_node]/multiplier[j]
                erosion_50km_save[j, i] = egforsave_save[year][proglacial_node]/multiplier[j] + dzbforsave_save[year][proglacial_node]/multiplier[j]
                i += 1
            
            uplift_rate_save[j] = u_folder
            
            dict1D_variables = pd.read_pickle("./%s.pkl" % filename_ug)
            locals().update(dict1D_variables)
            
            i = 0
            for year in time:
                erosion_ug_0km_save[j, i] = dzbforsave_save[year][glacial_node] 
                erosion_ug_50km_save[j, i] = dzbforsave_save[year][proglacial_node]
                i += 1
                                      
            j += 1
            

### ---- SAVE VARIABLES AS TXT FILES ----- ###

os.chdir('..') # get out of the data folder
os.chdir('figure_scripts')
os.chdir('Fig5_erosion_rate_by_glacial_extent_data')

np.savetxt('glacial_extent_save.txt', glacial_extent_save)
np.savetxt('uplift_rate_save.txt', uplift_rate_save)
np.savetxt('erosion_gl_save.txt', erosion_0km_save)
np.savetxt('erosion_pg_save.txt', erosion_50km_save)
np.savetxt('erosion_ug_gl_save.txt', erosion_ug_0km_save)
np.savetxt('erosion_ug_pg_save.txt', erosion_ug_50km_save)
       
                