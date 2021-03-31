#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:45:55 2019

@author: saschanz
"""

""" 
This code is one of two that simulates a coupled glacial fluvial system in a 1-d finite difference grid. This script creates the spin up landscape, starting from a 3000m plateau. Default is 5 My of fluvial erosion. To continue the model and enact glaciation, use analysis_run.py.

Bedrock erodes from stream power using shear stress. Sediment, if turned on, is moved via Meyer-Peter Mueller equation, modified to Wong and Parker, 2006. 

Input: model parameters such as node spacing, nodes, uplift rate, grain size, etc are set in inputs.py

Output: text files for bedrock elevation, sediment depth, and erosion rate at the end of the spinup and saved plot of the long profile at the end of spinup. Saves to active folder.

Questions? Contact:
    Sarah Schanz
    sschanz@coloradocollege.edu
    (719)389-6513
"""

import numpy as np
import matplotlib.pyplot as plt
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

## INSTANTIATE MODEL
river = stream(dx, nodes, initial_slope, initial_sed_depth, glacial_sed_supply_sw, glacial_discharge_sw, sediment_transport)
river.get_basin_geometry()
river.get_sediment_size(D0, a)

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
        dz_s, dz_b, dz_w = river.run_one_fluvial(dz_b_save, backgroundU, erosion_depth_threshold, dt, erosion_type = 'Turowski')
    else:
        dz_b = river.run_one_fluvial_nosed(backgroundU, dt)
 
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

print('Spin up time = %s' % (time))
