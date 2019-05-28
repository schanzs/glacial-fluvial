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
import os

## PARAMETERS TO TEST
backgroundU = 0.5/1000 # in m/yr
D0 = 200/1000 # in m
a = 0.02/1000
erosion_depth_threshold = 0.2
glacial_discharge = True
glacial_sed_supply = True
isostacy = True


### TIME
dt = 5.
dt_g = 0.01
dt_i = 10
glacierrepeats = int(dt/dt_g)
spinuptime = int(5000/dt)
analysistime = int(800000/dt)
totaltime = spinuptime + analysistime


### MODEL SET UP
dx = 2000
nodes = 100
initial_slope = 0.02
initial_sed_depth = 0.1

## INSTANTIATE MODEL & GET CONSTANTS
river = stream(dx, nodes, initial_slope, initial_sed_depth)
river.get_basin_geometry()
river.get_sediment_size(D0, a)
river.get_fall_velocity()

k10 = int(10000/dt)
dz_b_save = np.ones((nodes, k10)) * backgroundU

fig,[ax1, ax2, ax3] = plt.subplots(3, 1)

### RUN MODEL SPIN UP ##################################
time = 0
dz_b_foricalc = np.zeros(nodes)
while time in range (0, spinuptime):
    dz_s, dz_b, dz_w = river.run_one_fluvial(dz_b_save, backgroundU, dt)
    river.sed_depth -= dz_s
    river.z += dz_b
 
    if isostacy == True:
        dz_b_foricalc += dz_b
        if time%dt_i == 0:
            (dz_b_i) = river.isostacy(dz_b_foricalc, np.zeros(nodes), dt_i)
            river.z += dz_b_i
            dz_b_foricalc[:] = 0
    
    # Tracking variables:       
    dz_b_save[:,time%k10] = -dz_b/dt    
    river.z += backgroundU*dt
    river.z[-1] = 0
    river.sed_depth[-1] = 0
    time += 1


### RUN MODEL WITH GLACIERS #############################
HICE_prior = 1 * river.HICE[:]
dzbforsave = np.zeros(nodes)
while time in range(spinuptime, totaltime):
    #### FLUVIAL EROSION ######################################################
    dz_s, dz_b, dz_w = river.run_one_fluvial(dz_b_save, backgroundU, dt)
    river.sed_depth -= dz_s
    river.z += dz_b    
    dz_b_foricalc += dz_b

    #### GLACIAL EROSION ######################################################
    time2 = 0
    Qw_melt = 0
    glacial_sed_supply = 0
    
    # calculate new ELA    
    ELA = river.get_ELA(time, dt, averageELA = 2000, amplitude = 1000, period = 100000, shape = 'sawtooth')
       
    while time2 in range(0, glacierrepeats):
        # calculate annual ELA, with a 10 C annual variation
        ELA2 = 800 * np.sin(2*np.pi*time2/glacierrepeats) + ELA
        (Eg, Qw_melt_new, dHdt, new_glacial_sed_supply) = river.run_one_glacial(ELA2, dt_g)
        
        if glacial_discharge == True:
            Qw_melt += Qw_melt_new
        if glacial_sed_supply == True:
            glacial_sed_supply += new_glacial_sed_supply
            
        time2 += 1
        dz_b_foricalc += Eg

    #### ISOSTATIC ADJUSTMENTS ################################################
    if isostacy == True:
        if time%dt_i == 0:
            dHICE = river.HICE - HICE_prior
            (dz_b_i) = river.isostacy(dz_b_foricalc, np.zeros(nodes), dt_i)
            river.z += dz_b_i
            HICE_prior[:] = river.HICE[:]
            dz_b_foricalc[:] = 0

    # Tracking variables:
    dz_b_save[:,time%k10] = -dz_b/dt
    dzbforsave += -dz_b/dt
    river.z += backgroundU*dt
    river.z[-1] = 0
    river.sed_depth[-1] = 0
       
    #### PLOTTING #############################################################
    if ((time+1)*dt)%100000 == 0:
        ax3.clear()
    if ((time+1)*dt)%1000 == 0:
        ax1.clear()
        ax2.clear()
        ax1.plot(river.x/1000, river.z, 'k-', river.x/1000, river.z+river.sed_depth, '--r', river.x/1000, river.z+river.sed_depth+river.HICE, 'b-.')
        ax1.legend(['bedrock', 'sediment', 'ice'])
        ax1.set_title('Model year %s' % ((time+1-spinuptime)*dt))
        ax2.plot(river.x/1000, river.sed_depth, river.x/1000, river.HICE)
        ax2.legend(['sediment', 'ice'])
        ax2.set_xlabel('Downstream distance (km)')
        ax2.set_ylabel('Height (m)')
        ax3.plot((time*dt)%100000, ELA, 'o', color = 'grey')
        ax3.set_ylabel('ELA elevation')
        ax3.set_xlim([0, 100000])
        ax3.set_ylim([1000, 3000])
        plt.show()
        plt.pause(0.0000001)
        print ('Calculated model year %s' % ((time+1)*dt))

        fig.savefig('time_%06d.png' % int((time+1)*dt))
        
    ##### SAVE VARIABLES ######################################################
    if ((time+1)*dt)%1000 == 0:
        # save variables now
        year = int((time+1-spinuptime)*dt)
        np.savetxt('z_%06d.txt' % year, river.z)
        np.savetxt('sed_%06d.txt' % year, river.sed_depth)
        np.savetxt('ice_%06d.txt' % year, river.HICE)
        np.savetxt('dzb_%06d.txt' % year, dzbforsave)
        dzbforsave[:] = 0

    time += 1         