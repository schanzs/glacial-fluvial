#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 17:01:04 2023

@author: sschanz
"""

import numpy as np
import os
import pandas as pd


### ---- SET UP CONSTANTS --- ###
x = np.arange(0, 200, 2)*1000
time = np.arange(1, 801, 1)*1000


### ----- SET UP SAVE VARIABLES ---- ###
average_extent_save = np.zeros(3*5*3)
flux_out_glac_save = np.zeros((3*5*3, 800))
flux_out_unglac_save = np.zeros((3*5*3, 800))

### ----- LOOP THROUGH ---- ###

uplift_folders = ('0.001', '0.0005', '0.0001')
D0_folders = ('0.5', '0.25', '0.1', '0.05', '0.01')
a_folders = ('0.0001', '5e-05', '1e-05')

os.chdir('..') # get out of figure_scripts folder
os.chdir('data') # get into data
j = 0
for u_folder in uplift_folders:

    for d_folder in D0_folders:

        for a_folder in a_folders:

            filename_ug = ('%s_%s_%s_False' % (u_folder, d_folder, a_folder))
            filename = ('%s_%s_%s_True' % (u_folder, d_folder, a_folder))
            
            flux_out_glac = np.zeros(800)
            flux_out_unglac = np.zeros(800)
            glacial_extent = np.zeros(8)
            
            dict1D_variables = pd.read_pickle("./%s.pkl" % filename)
            locals().update(dict1D_variables)
            
            i = 0
            k = 0
            for year in time:
                flux_out_glac[i] = sedoutforsave_save[year][-1] 
                
                if year%100000 == 0:
                    ice = HICE_save[year]
                    glacial_extent[k] = np.where(ice == 0)[0][0]*2
                    k += 1
                    
                i += 1
            
            avg_glacial_extent = np.mean(glacial_extent)

            dict1D_variables = pd.read_pickle("./%s.pkl" % filename_ug)
            locals().update(dict1D_variables)
                        
            i = 0
            for year in time:
                flux_out_unglac[i] = sedoutforsave_save[year][-1]# units of kg/yr
                i += 1
                       
            ## SAVE VARIABLES
            average_extent_save[j] = avg_glacial_extent
            flux_out_unglac_save[j,:] = flux_out_unglac
            flux_out_glac_save[j,:] = flux_out_glac
            
            j += 1
            


### ---- SAVE VARIABLES AS TXT FILES ----- ###
os.chdir('..') # get out of data
os.chdir('figure_scripts')
os.chdir('Fig7_sed_flux_by_glacial_extent_data')

np.savetxt('average_extent_save.txt', average_extent_save)
np.savetxt('flux_out_glac_save.txt', flux_out_glac_save)
np.savetxt('flux_out_unglac_save.txt', flux_out_unglac_save)

os.chdir('..')
               
                