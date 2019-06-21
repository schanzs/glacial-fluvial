#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:45:55 2019

@author: saschanz
"""

""" Glacier - fluvial model based on model written by B Yanites. Bedrock erodes from either stream power or saltation-abrasion (Sklar and Dietrich, 2004, 2006). Sediment is moved via Meyer-Peter Mueller equation, modified to Wong and Parker, 2006. Glacial erosion is applied following MacGregor et al., 2000 (Geology) model, with updates to the abrasion model. """

import numpy as np
from glacial_fluvial_functions import stream
import matplotlib.pyplot as plt

""" PARAMETERS TO TEST from inputs.py
backgroundU 
D0 
a
erosion_depth_threshold
glacial_discharge_sw
glacial_sed_supply_sw
isostacy
dt
dt_g
dt_i
dx
nodes
initial_slope
initial_sed_depth
"""
from inputs import backgroundU, D0, a, erosion_depth_threshold, glacial_discharge_sw, glacial_sed_supply_sw, isostacy_sw, dt, dt_g, dt_i, dx, nodes, initial_slope, initial_sed_depth


### TIME
dt = int(dt)
glacierrepeats = int(dt/dt_g)
spinuptime = 5000000
analysistime = 800000
totaltime = spinuptime + analysistime


## INSTANTIATE MODEL & GET CONSTANTS
river = stream(dx, nodes, initial_slope, initial_sed_depth, glacial_sed_supply_sw, glacial_discharge_sw)
river.get_basin_geometry()
river.get_sediment_size(D0, a)

k10 = int(10000/dt)
dz_b_save = np.ones((nodes, k10)) * backgroundU

fig,[ax1, ax2, ax3] = plt.subplots(3, 1)

### RUN MODEL SPIN UP ##################################
time = 0
dz_b_foricalc = np.zeros(nodes)
while time in range (0, spinuptime):
    dz_s, dz_b, dz_w = river.run_one_fluvial(dz_b_save, backgroundU, erosion_depth_threshold, dt)
 
    if isostacy_sw == True:
        dz_b_foricalc += dz_b
        if time%dt_i == 0:
            (dz_b_i) = river.isostacy(dz_b_foricalc, np.zeros(nodes), dt_i)
            river.z += dz_b_i
            dz_b_foricalc[:] = 0
    
    # Tracking variables:       
    dz_b_save[:,int((time/dt)%k10)] = -dz_b/dt    
    river.z += backgroundU*dt
    river.z[-1] = 0
    river.sed_depth[-1] = 0
    time += dt

        
### CALCULATE ELA PARAMETERS BASED ON FRACTION OF PROFILE
zrange = np.max(river.z) - np.min(river.z)
averageELA = 1/3 * zrange
amplitude = 1/6 * zrange


### RUN MODEL WITH GLACIERS #############################
HICE_prior = 1 * river.HICE[:]
dzbforsave = np.zeros(nodes)
egforsave = np.zeros(nodes)
while time in range(spinuptime, totaltime):

        
    #### GLACIAL EROSION ######################################################
    analysis_time = time-spinuptime
    
    Eg_total, ELA = river.run_one_glacial(analysis_time, averageELA, amplitude, dt_g, dt)
    dz_b_foricalc += Eg_total

    
    #### FLUVIAL EROSION ######################################################
    river.get_sediment_size(D0, a)
    dz_s, dz_b, dz_w = river.run_one_fluvial(dz_b_save, backgroundU, erosion_depth_threshold, dt)   
    dz_b_foricalc += dz_b


    #### ISOSTATIC ADJUSTMENTS ################################################
    if isostacy_sw == True:
        if time%dt_i == 0:
            dHICE = river.HICE - HICE_prior
            (dz_b_i) = river.isostacy(dz_b_foricalc, dHICE, dt_i)
            river.z += dz_b_i
            HICE_prior[:] = river.HICE[:]
            dz_b_foricalc[:] = 0

    # Tracking variables:
    dz_b_save[:,int((time/dt)%k10)] = -dz_b/dt
    dzbforsave += -dz_b/dt
    egforsave += Eg_total/dt
    river.z += backgroundU*dt
    river.z[-1] = 0
    river.sed_depth[-1] = 0
       
    #### PLOTTING #############################################################
    if analysis_time%100000 == 0:
        ax3.clear()
    if analysis_time%10000 == 0:
        ax1.clear()
        ax2.clear()
        ax1.plot(river.x/1000, river.z, 'k-', river.x/1000, river.z+river.sed_depth, '--r', river.x/1000, river.z+river.sed_depth+river.HICE, 'b-.')
        ax1.legend(['bedrock', 'sediment', 'ice'])
        ax1.set_title('Model year %s' % analysis_time)
        ax2.plot(river.x/1000, river.sed_depth, river.x/1000, river.HICE)
        ax2.legend(['sediment', 'ice'])
        ax2.set_xlabel('Downstream distance (km)')
        ax2.set_ylabel('Height (m)')
        ax3.plot((analysis_time)%100000, ELA, 'o', color = 'grey')
        ax3.set_ylabel('ELA elevation')
        ax3.set_xlim([0, 100000])
        ax3.set_ylim([1000, 3000])
        plt.show()
        plt.pause(0.0000001)
        print ('Calculated model year %s' % analysis_time)

        fig.savefig('time_%06d.png' % analysis_time)
        
    ##### SAVE VARIABLES ######################################################
    if analysis_time%1000 == 0:
        # save variables now
        year = analysis_time
        np.savetxt('z_%06d.txt' % year, river.z)
        np.savetxt('sed_%06d.txt' % year, river.sed_depth)
        np.savetxt('ice_%06d.txt' % year, river.HICE)
        np.savetxt('dzb_%06d.txt' % year, dzbforsave)
        np.savetxt('eg_%06d.txt' % year, egforsave)
        np.savetxt('wv_%06d.txt' % year, river.Wv)
        dzbforsave[:] = 0
        egforsave[:] = 0

    time += dt