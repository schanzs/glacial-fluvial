#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:45:55 2019

@author: saschanz
"""

""" 
This code is part 2 of 2 that simulates coupled glacial fluvial erosion in a 1-d finite difference grid. This script applies glacial erosion to a pre-determined landscape. To make the pre-landscape, you can upload elevation, sediment depth, and erosion rate data or use the spin_up_run.py script to simulate a landscape.

In this script, glacial erosion, bedrock erosion, sediment transport, and isostatic rebound are imposed. The glacial erosion model follows MacGregor et al, 2000 (Geology) with updates for an abrasion-erosion model of Iverson (2012). Bedrock erosion is from shear stress derivation of the stream power equation, with a sediment cover effect using a depth threshold. Sediment is moved via a Meyer-Peter and Mueller equation, modified by Wong and Parker (2006). Sediment supply is based on the 10 ky average of erosion rates. Isostacy is calculated using the methods outlined in Pollard and DeConto (2012), following Huybrechts and de Wolde (1999). 

Inputs: river geometry from spin_up_run.py (or text files of elevation, sediment depth, and erosion rate), and input parameters from inputs.py. 

Outputs: a folder of text files saving the bedrock elevation, sediment depth, ice depth, valley width, glacial erosion, and fluvial erosion every 1000 years. 

Model orginally written by B Yanites, modified and updated by S Schanz

Questions? Contact:
    Sarah Schanz
    sschanz@coloradocollege.edu
    (719)389-6513
"""

import numpy as np
import os
from glacial_fluvial_functions import stream
from inputs import *

""" 
PARAMETERS imported from inputs.py are:

    MODEL VARIABLES
        backgroundU : m/yr
        D0          : m

    MODEL TYPE ON/OFF SWITCHES
        glacial_discharge_sw    : boolean
        glacial_sed_supply_sw   : boolean
        isostacy                : boolean
        sediment_transport      : boolean
        
    TIME
        dt      : fluvial erosion timestep, years
        dt_g    : glacial erosion timestep, years - 0.01 at least
        dt_i    : isostacy calculation recurrence, years
        
    MODEL SET UP
        dx                          : node spacing, meters
        nodes                       : number of nodes
        initial_slope               : initial slope of profile, overriden in block model
        initial_sed_depth           : thickness of bedload, meters
        a                           : sediment attrition rate
        erosion_depth_threshold     : sediment cover that halts erosion, meters
"""


### CREATE NEW FOLDER TO SAVE OUTPUTS
path = os.getcwd()
new_dir = path + '/analysis'
os.mkdir(new_dir)


## INSTANTIATE MODEL
river = stream(dx, nodes, initial_slope, initial_sed_depth, glacial_sed_supply_sw, glacial_discharge_sw, sediment_transport)
river.get_basin_geometry()
river.get_sediment_size(D0, a)
k10 = int(10000/dt)   #parameter to update erosion rates for the 10ky average for sed supply

### IMPORT MODEL SPIN UP ##################################
river.z = np.genfromtxt('z_endspin.txt') # text file with 1 column and n rows, where n = nodes
river.sed_depth = np.genfromtxt('sed_endspin.txt') # text file with 1 column and n rows, where n = nodes
dz_b_save = np.genfromtxt('dzb_endspin.txt') # text file with k10 columns and n rows, where n = nodes, and k10 = 10,000 / dt - records 10,000 year history of fluvial erosion

os.chdir('analysis')
          
### SET ELA PARAMETERS
averageELA = 2500 # meters
amplitude = 1000 # meters

### RUN MODEL WITH GLACIERS #############################
dt = int(dt)
run_time = 800000 #years
time=0

# create empty variables
HICE_prior = 1 * river.HICE[:]
dz_b_foricalc = np.zeros(nodes)
dzbforsave = np.zeros(nodes)
egforsave = np.zeros(nodes)

while time in range(run_time):
    time += dt
    #### GLACIAL EROSION ######################################################    
    Eg_total, ELA = river.run_one_glacial(time, averageELA, amplitude, dt_g, dt)
    dz_b_foricalc += Eg_total

    #### FLUVIAL EROSION ######################################################
    river.get_sediment_size(D0, a) # updated based on glacial extent
    if sediment_transport == True:
        dz_s, dz_b, dz_w = river.run_one_fluvial(dz_b_save, backgroundU, erosion_depth_threshold, dt)   
    else:
        dz_b = river.run_one_fluvial_nosed(backgroundU, dt)
    dz_b_foricalc += dz_b

    #### ISOSTATIC ADJUSTMENTS ################################################
    if isostacy_sw == True:
        if time%dt_i == 0:
            dHICE = river.HICE - HICE_prior
            (dz_b_i) = river.calc_isostacy(dz_b_foricalc, dHICE, dt_i)
            river.z += dz_b_i
            HICE_prior[:] = river.HICE[:]
            dz_b_foricalc[:] = 0

    #### UPDATE EROSION FOR SEDIMENT SUPPLY
    dz_b_save[:,int((time/dt)%k10)] = -dz_b/dt
    
    #### UPDATE VARIABLES FOR OUTPUT
    dzbforsave += -dz_b/dt
    egforsave += Eg_total/dt
    
    #### IMPOSE UPLIFT FOR NEXT ROUND
    river.z += backgroundU*dt
    river.z[-1] = 0
    river.sed_depth[-1] = 0
       
    ##### SAVE VARIABLES ######################################################
    if time%1000 == 0:    # save variables every 1000 years
        year = time
        np.savetxt('z_%06d.txt' % year, river.z)
        np.savetxt('sed_%06d.txt' % year, river.sed_depth)
        np.savetxt('ice_%06d.txt' % year, river.HICE)
        np.savetxt('dzb_%06d.txt' % year, dzbforsave)
        np.savetxt('eg_%06d.txt' % year, egforsave)
        np.savetxt('wv_%06d.txt' % year, river.Wv)
        
        # clear the variable placeholders
        dzbforsave[:] = 0
        egforsave[:] = 0
