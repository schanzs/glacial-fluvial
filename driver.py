#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:45:55 2019

@author: saschanz
"""

""" 
This code is one of two that simulates a coupled glacial fluvial system in a 1-d finite difference grid. This script drives the spin up and analysis runs, to be used in Dakota for a parameter space exploration.

Bedrock erodes from stream power using shear stress and modified by the ratio of sediment supply to transport, ala Sklar and Dietrich 2004. Sediment is moved via Meyer-Peter Mueller equation, modified to Wong and Parker, 2006. 

Variables/inputs: Uplift rate, glacial erosion factor, sediment size, and sediment attrition.

Output: text files for bedrock elevation, sediment depth, ice depth, sediment flux, and bedrock erosion at the end of spinup and every 1000 years during the analysis runs. Saves to a designated folder.

Questions? Contact:
    Sarah Schanz
    sschanz@coloradocollege.edu
    (719)389-6513
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from glacial_fluvial_functions import stream
import sys

""" 
PARAMETERS imported from inputs.py are:
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
glacial_discharge_sw = True
glacial_sed_supply_sw = True
isostacy_sw = True
sediment_transport = True

dt = 1
dt_g = 0.1
dt_i = 10.

dx = 2000
nodes = 100
initial_slope = 0.02
initial_sed_depth = 0.1
erosion_depth_threshold = 0.2 

#################################
#                               #
#    STEP 1: Use                #
#     input files               #
#    to prep model run &        #
#     save files                #
#                               #
#################################

backgroundU = float(sys.argv[1])
D0 = float(sys.argv[2])
a = float(sys.argv[3])
switch_placeholder = 'True'
glaciated_sw = (switch_placeholder == sys.argv[4])
C2 = 0.0002

savefolder_name = ('%s_%s_%s_%s' % (backgroundU, D0, a, glaciated_sw))


### SPECIFY VARIABLES TO SAVE #####################
variables_to_save1D = ('z', 'Wv', 'sed_depth', 'HICE', 'dzbforsave', 'egforsave', 'sedoutforsave')

dict1D_save = {}
for variable in variables_to_save1D:
    variable_name = variable + '_save'
    dict1D_save[variable_name] = {}
    
##################################
#                                #
#     STEP 2: Run the model      #
#                                #
##################################

## INSTANTIATE MODEL
river = stream(dx, nodes, initial_slope, initial_sed_depth, glacial_sed_supply_sw, glacial_discharge_sw, sediment_transport)
river.get_basin_geometry()
river.get_sediment_size(D0, a)
river.C2 = C2

## SET UP EMPTY VARIABLES FOR SEDIMENT SUPPLY
k10 = int(10000/dt)
dz_b_save = np.zeros((nodes, k10))

## SET PLATEAU TOPOGRAPHY
river.z[:-1] = 3000

### RUN MODEL SPIN UP ##################################
time = 0
dt = int(dt)
dz_b_foricalc = np.zeros(nodes)

while time < 5000000:
    time += dt
    
    ### FLUVIAL EROSION
    if sediment_transport == True:
        dz_s, dz_b, dz_w = river.run_one_fluvial(dz_b_save, backgroundU, erosion_depth_threshold, dt, erosion_type = 'SklarDietrich')
    else:
        dz_b = river.run_one_fluvial_nosed(backgroundU, dt)
    dz_b_foricalc += dz_b
    
    ### ISOSTASY
    if isostacy_sw == True:
        dz_b_foricalc += dz_b
        if time%dt_i == 0:
            (dz_b_i) = river.calc_isostacy(dz_b_foricalc, np.zeros(nodes), dt_i)
            river.z += dz_b_i
            dz_b_foricalc[:] = 0
    
    ### UPDATE DZ_B FOR SEDIMENT SUPPLY       
    dz_b_save[:,int((time/dt)%k10)] = -dz_b/dt    #positive values are incision
    
    ### UPDATE FOR UPLIFT
    river.z += backgroundU*dt
    river.z[-1] = 0
    river.sed_depth[-1] = 0
          
    if time%1e6 == 0:
        print(time)
   
### SAVE OUTPUTS ################################
np.savetxt('z_endspin.txt', river.z)
np.savetxt('sed_endspin.txt', river.sed_depth)
np.savetxt('dzb_endspin.txt', dz_b_save)

### PLOT END SPIN UP FIGURE AND SAVE ############
plt.figure()
plt.plot(river.x, river.z, river.x, river.z+river.sed_depth)
plt.title('uplift: %s and D0: %s' % (backgroundU, D0))
plt.savefig('endspin_profile.png')


### RUN MODEL ANALYSIS ####################################
time2 = 0
### SET ELA PARAMETERS
averageELA = 2500 # meters
amplitude = 1000 # meters

#### CREATE EMPTY VARIABLES ############3
HICE_prior = 1 * river.HICE[:]
river.dzbforsave = np.zeros(nodes)
river.egforsave = np.zeros(nodes)
river.sedoutforsave = np.zeros(nodes)


while time2<800000:
    time2 += dt
    
    #### GLACIAL EROSION ######################################################    
    if glaciated_sw == True:
        Eg_total, ELA = river.run_one_glacial(time2, averageELA, amplitude, dt_g, dt)   
    else:
        Eg_total = np.zeros(nodes)
    dz_b_foricalc += Eg_total

    #### FLUVIAL EROSION ######################################################
    river.get_sediment_size(D0, a) # updated based on glacial extent
    if sediment_transport == True:
        dz_s, dz_b, dz_w = river.run_one_fluvial(dz_b_save, backgroundU, erosion_depth_threshold, dt, erosion_type = 'SklarDietrich')   
    else:
        dz_b = river.run_one_fluvial_nosed(backgroundU, dt)
    dz_b_foricalc += dz_b

    #### ISOSTATIC ADJUSTMENTS ################################################
    if isostacy_sw == True:
        if time2%dt_i == 0:
            dHICE = river.HICE - HICE_prior
            (dz_b_i) = river.calc_isostacy(dz_b_foricalc, dHICE, dt_i)
            river.z += dz_b_i
            HICE_prior[:] = river.HICE[:]
            dz_b_foricalc[:] = 0

    #### UPDATE EROSION FOR SEDIMENT SUPPLY
    dz_b_save[:,int((time2/dt)%k10)] = -dz_b/dt
    
    #### UPDATE VARIABLES FOR OUTPUT
    river.dzbforsave += -dz_b/dt
    river.egforsave += Eg_total/dt
    river.sedoutforsave += river.Qs_down[:-1]/dt
    
    #### IMPOSE UPLIFT FOR NEXT ROUND
    river.z += backgroundU*dt
    river.z[-1] = 0
    river.sed_depth[-1] = 0
    
    ##### SAVE VARIABLES ######################################################
    if time2%1000 == 0:    # save variables every 1000 years
        for item in variables_to_save1D:
            var_name = item + '_save'
            dict1D_save[var_name][time2] = getattr(river, item)*1
            
        # clear the variable placeholders
        river.dzbforsave *= 0
        river.egforsave *= 0
        river.sedoutforsave *= 0
        
##### SAVE FINAL VARIABLES TO FOLDER #######
df1D = pd.DataFrame(dict1D_save)
df1D.to_pickle("./%s.pkl" % savefolder_name)


