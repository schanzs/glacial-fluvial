#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:49:50 2019

@author: saschanz
"""

""" 
STREAM FUNCTIONS
This script runs the erosion for the glacial fluvial model through the creation of the object 'stream' and its associated functions.

The object is created in either spin_up_run.py or analysis_run.py

Functions are arranged as:
    __init__
    get_(various) - alphabetical, for getting channel geometry-related variables
    calc_(various) - alphabetical, for longer calculations of channel change/erosion
    run_one_(various) - alphabetical, overarching functions used to run one timestep

Questions? contact:
    Sarah Schanz
    sschanz@coloradocollege.edu
    (719)389-6513
"""

import numpy as np
from scipy.special import kei


class stream(object):
    g = 9.81            # gravitational acceleration, meters per second squared
    rhow = 1000         # density of water, cubic meters per second
    rhoi = 917          # density of ice, cubic meters per second
    rhos = 2650         # density of sediment, cubic meters per second
    rhom = 3300         # density of asthenosphere, cubic meters per second
    beta = 0.2          # fraction of eroded hillslope material that is bedload
    R = rhos/rhow - 1   # specific gravity of sediment
    tauc = 0.0495       # critical shear stress
    n = 0.04            # Manning's roughness
    frac_yr_transport = 0.05  # fraction of the year experiencing erosion
    L = 0.2             # lambda, porosity of bedload
    kf = -5e-5        # vertical erosion coefficient
    # kf = -3e-3   # from Yanites, 2018
    kw = -7.3*kf        # lateral bedrock erosion coefficient, 7.3 times vertical from Finnegan and Dietrich, 2012
    A = 2.1*10e-18      # Arrhenius constant
    C1 = 0.0012         # glacial sliding coefficient from MacGregor et al (2000)
    C2 = 0.00125        # glacial erosion coefficient from Iverson (2012) for a Gaussian distribution of fractures in heterogenous rock
    alpha = 0.01        # meters/ (year meters) from Oerlemans, 1984

    
    def __init__(self, dx, nodes, initial_slope, initial_sed_depth, glacial_sed_supply_sw, glacial_discharge_sw, sediment_transport):
        """
        Initalization for object of class stream
        ----------
        dx : integer
            meter spacing between nodes.
        nodes : integer
            number of model nodes
        initial_slope : float
            slope of initial river, often overriden if elevation specified later.
        initial_sed_depth : float
            thickness of sediment on river.
        glacial_sed_supply_sw : boolean
            turns transfer of glacially eroded bedload on or off.
        glacial_discharge_sw : boolean
            turns transfer of glacial melt to river on or off.
        sediment_transport : boolean
            turns sediment transport on or off.

        Returns
        -------
        None.

        """
        self.nodes = nodes
        self.dx = dx
        self.initial_slope = initial_slope
        
        # use input model set up to make distance, elevation, sed depth:
        self.x = np.arange(0, self.nodes*self.dx, self.dx)
        self.z = initial_slope * (max(self.x)-self.x)
        self.z[-1] = 0
        self.sed_depth = np.ones(nodes)*initial_sed_depth

        # set up empties for later
        self.HICE = np.zeros(nodes)         # height of ice
        self.Qs_down = np.zeros(nodes+1)    # sediment flux to downstream node
        self.H = np.zeros(nodes)            # height of water
        self.Qw_change = 0                  # water discharge from glacial melt
        self.Qw = np.zeros(nodes)           # water discharge total
        self.glacial_sed_supply = 0         # sediment flux from glacial erosion
        
        # set up on/off switches
        self.glacial_sed_supply_sw = glacial_sed_supply_sw
        self.glacial_discharge_sw = glacial_discharge_sw
    
    def get_basin_geometry(self):
        """
        Calculates drainage area, discharge, and channel width based on power-law relationships to distance.
        
        Ah      : drainage area, square meters
        Qwi     : water discharge from power-law, not included glacial melt, cubic meters per second
        W       : channel width, meters
        Wv      : valley width, initially the same as channel but widens from glacial erosion, meters

        Returns
        -------
        None.

        """
        Hc = 1
        He = 1.8598
        self.Ah = Hc * self.x**He + 100000
        
        kQ = 10**-7
        eQ = 1
        self.Qwi = kQ * self.Ah**eQ
        
        kw = 5.
        eW = 0.5
        self.W = kw*self.Qwi**eW
        self.Wv = 1*self.W[:]
         
       
    def get_sediment_in(self, dz_b, dt):
        """
        Calculates the amount of sediment supplied from hillslopes. Assumes hillslopes erode at the same rate as the river incises, with that incision averaged over 10,000 years. Calculates Qs_hills, the sediment supply from hillslopes.

        Parameters
        ----------
        dz_b : 2D array of floats
            Erosion rates at each node for the last 10 ky.
        dt : integer
            timestep.

        Returns
        -------
        None.

        """
        dA = np.hstack((self.Ah[0], np.diff(self.Ah)))
        self.Qs_hills = self.rhos * self.beta * dA * dz_b * dt
            
         
    def get_sediment_size(self, D0, a):
        """
        Calculates the size of sediment in the channel using Sternberg's law and the attrition rate a set in the input file. The initial median grain size is set in the input file, and occurs at the first fluvial node - this is updated every timestep to account for glacial extent. Assumption is that glacial erosion occurs through quarrying (main erosion process in glacial erosion function) and this provides most of the bedload. Minimum grain size is set to 0.005 meters (5 mm). Creates self.D, an array of grain sizes.

        Parameters
        ----------
        D0 : float
            initial median grain size, meters
        a : float
            attrition rate.

        Returns
        -------
        None.

        """
        D = D0 * np.exp(-self.x*a)
        D[np.where(D<0.005)] = 0.005
        
        # shift D distribution for glacier (assumes glacier spits out D0 grain size):
        first_fluvial = np.min(np.where(self.HICE == 0))
        length_fluvial = self.nodes - first_fluvial
        self.D = np.zeros(self.nodes)
        self.D[first_fluvial:] = D[:length_fluvial]
    
        
    def get_slope(self):
        """
        Calculates channel slope in forward space using the topographic (bedrock + sediment) slope
        
        slope : streambed slope, m/m
        
        Returns
        -------
        None.

        """
        topo = self.sed_depth + self.z
        slope = np.diff(topo)/np.diff(self.x) # nodes-1 length, i-(i-1)
        self.slope = np.append([slope], [0])
        self.slope[-1] = 0
        self.slope[self.slope>0] = 0
        
        
    def get_taub(self):
        """
        Calculates basal shear stress in the channel, using the small angle approximation that tan(theta) = slope. Uses hydraulic radius instead of water depth. Slope is negative, so multipled by -1. Shear stress is zero where a glacier is present and where slope is 0.

        Returns
        -------
        None.

        """
        self.taub = self.rhow * self.g * (self.W * self.H / (2*self.H + self.W)) * -1 * self.slope
        self.taub[self.HICE > 0] = 0
        self.taub[self.slope == 0] = 0  
        
        
    def get_water_height(self, Qw):
        """
        Calculates the height of water in the channel using Manning's equation

        Parameters
        ----------
        Qw : 1D array of floats
            water discharge. Often self.Qw, but may be modified as self.Qwi

        Returns
        -------
        None.

        """
        self.H[self.slope < 0] = (self.n * Qw[self.slope<0]/self.W[self.slope<0])**.6 * (-1*self.slope[self.slope<0])**(-0.3)


    def calc_bank_erosion(self, dt):
        """
        Calculates lateral erosion of the bank using shear stress and a lateral erodibility coefficient. If the valley has been widened previously by glaciation, then a different lateral erodibility is used, assuming alluvial banks are preferentially eroded over bedrock. Channel width remains the same but eroded material becomes part of the sediment flux downstream; this allows glacial debris to be mobilized during interglacials.

        Parameters
        ----------
        dt : integer
            timestep, years

        Returns
        -------
        dz_w : 1D array of floats
            change in valley width, meters

        """
       
        dz_w = self.kw * self.taub * dt
        dz_w[self.Wv > self.W] = 10*dz_w[self.Wv > self.W]
        dz_w[self.HICE > 0] = 0
        
        return dz_w
    
    
    def calc_bedrock_erosion(self, dt, depth_threshold):
        """
        Calculates bedrock erosion based on shear stress derivation of stream power equation. erosion is zero where sediment cover is greater than the depth threshold, and where ice is present.

        Parameters
        ----------
        dt : integer
            timestep, years
        depth_threshold : float
            depth of sediment that protects bed from erosion, meters

        Returns
        -------
        dz_b : 1D array of floats
            change in bedrock height, negative = erosion, meters

        """

        dz_b = self.kf * self.taub * dt
     
        dz_b[self.sed_depth > depth_threshold] = 0
        dz_b[self.HICE > 0] = 0
        
        return dz_b
   
    def calc_bedrock_erosion_sklardietrich(self, dt):
        """
        Calculates bedrock erosion based on saltation-abrasion via Sklar and Dietrich (2004).

        Parameters
        ----------
        dt : integer
            timestep, years
        
        Returns
        -------
        dz_b : 1D array of floats
            change in bedrock height, negative = erosion, meters

        """
        self.F = np.zeros(self.nodes)
        # erosion factors from Sklar and Dietrich (2004)
        validF = np.where(self.Qs_cap != 0)
        self.F[validF] = 1-self.Qs_n[validF]/(np.pi*10**7*dt*self.Qs_cap[validF])
        self.F[self.F < 0] = 0
        
        dz_b = self.kf * self.taub * self.F * dt
        dz_b[self.HICE > 0] = 0
   
        return dz_b
        
    def calc_bedrock_erosion_turowski(self, dt):
        """
        Calculates bedrock erosion based on a cover effect defined by Turowski et al, 2007 in which erosion is modulated by a factor exp(-alpha*Qs/Qt) - assuming alpha = 1

        Parameters
        ----------
        dt : integer
            timestep, years
        
        Returns
        -------
        dz_b : 1D array of floats
            change in bedrock height, negative = erosion, meters

        """
        
        self.Ra = np.zeros(self.nodes)
        Qs_Qcap = np.zeros(self.nodes)
        valid_ratio = np.where(self.Qs_cap != 0)
        Qs_Qcap[valid_ratio] = self.Qs_n[valid_ratio]/(np.pi*10**7 * dt * self.Qs_cap[valid_ratio])
        
        if ("self.old_Ra" in globals()) == False:
            # first time step, must use steady state
            self.Ra = np.exp(-Qs_Qcap)
            
        else:
            # not the first timestep
            
            for k in range(self.nodes):
                if np.pi*10**7 * dt * self.Qs_cap[k] > self.Qs_n[k]:
                # no storage
                    self.Ra[k] = np.exp(-Qs_Qcap[k])
                elif self.Qs_cap[k] == 0:
                    self.Ra[k] = np.exp(-Qs_Qcap[k])
                else:
                # storage over supply
                    self.Ra[k] = -self.old_Ra[k] * (1-(Qs_Qcap[k] - self.old_Qs_Qcap[k]))


        dz_b = self.kf * self.taub * self.Ra * dt
        dz_b[self.HICE > 0] = 0
        
        self.old_Ra = 1*self.Ra
        self.old_Qs_Qcap = 1*Qs_Qcap
        
        return dz_b   
    
    def calc_bedrock_erosion_shobe(self, h_star, dt):
        """
        Calculates bedrock erosion based on a cover effect defined by exp(-h/h*) where h is sediment depth and h* is a threshold sediment depth, as explained in Shobe et al, 2017, SPACE model

        Parameters
        ----------
        dt : integer
            timestep, years
        
        Returns
        -------
        dz_b : 1D array of floats
            change in bedrock height, negative = erosion, meters

        """

        self.h_ratio = self.sed_depth/h_star
        dz_b = self.kf * self.taub * np.exp(-self.h_ratio) * dt
        dz_b[self.HICE > 0] = 0
   
        return dz_b   

    def calc_ELA(self, time, averageELA, amplitude, period = 100000, shape = 'sawtooth'):
        """
        Function calculates the ELA for each fluvial time step.
        For reference, a 1000 m amplitude is approximately a 6 degree mean annual temperature shift (Porter, 1989) noted in Sierra Nevada glaciers (Phillips et al, 1996). A shift from 3000 to 2000 m asl ELA is used by Egholm et al (2009) in their model of glacial erosion.

        Parameters
        ----------
        time : integer
            model year, year
        averageELA : float
            average elevation of ELA, meters
        amplitude : float
            range of elevations of ELA, meters
        period : integer, optional
            Periodicity of glacial cycles. The default is 100000.
        shape : text, optional
            Shape of ELA variations in glacial cycle. Options are 'sine' and 'sawtooth'. The default is 'sawtooth'.

        Returns
        -------
        ELA : float
            Elevation in meters of the ELA for the particular timestep given in 'time'

        """
    
        # ELA OPTION 1: sinusoidal, equal time coming in and out of glacials/interglacials    
        if shape == 'sine':
            ELA = -amplitude/2*np.cos(2*np.pi*time/period)+averageELA 

        # ELA OPTION 2: sawtooth, 10k rapid warming and 90k slow cooling - starts at lowest temperature/ELA
        if shape == 'sawtooth':
            cycle_i = (time)%period
            if cycle_i <= 10000: # start with rapid warming 
                ELA = amplitude/10000*cycle_i + (averageELA-amplitude/2)
            else: # slow cooling over the other 90,000 years
                ELA = -amplitude/90000*cycle_i + (averageELA+amplitude/2 + amplitude/9)
                
        return ELA
    
    def calc_ELA_vostok(self, yrBP, deltaC, analysistime, ELA_pinedale, ELA_now, dt):
        """ Calculate the ELA on an annual basis 
        """
        # trim to analysis time:
        deltaC = deltaC[yrBP<analysistime]
        yrBP = yrBP[yrBP<analysistime]
        
        # find linear relationship between deltaC and ELA
        deltaC_pinedale = np.min(deltaC)
        intercept = ELA_now
        if int(deltaC_pinedale) == 0:
            ELA = np.ones(len(deltaC))*ELA_now
        else:
            slope = (ELA_pinedale - ELA_now) / (deltaC_pinedale - 0)
            ELA = intercept + deltaC*slope
        
        ELAyr = analysistime - yrBP
        
        # use spline interpolation to get new ELA/ELAyr pairs:
        years = np.arange(0, analysistime, dt)
        newELA = np.interp(years, ELAyr[::-1], ELA[::-1])
        
        return newELA
    
    def calc_glacial_erosion(self, ELA2, dt_g):
        """
        Longitudinal glacier evolution following Macgregor et al., 2000 (Geology) formulation, with updates to the quarrying rule from Iversion (2012). Most calculations take place between nodes, hence creation of parallel set of variables representing model variables between the nodes.

        Parameters
        ----------
        ELA2 : float
            elevation of the ELA for timestep dt_g (options in run_one_glacial to adapt for interannual variability), meters
        dt_g : float
            glacier timesteps, should be ~0.01 to capture interannual variability, years

        Returns
        -------
        Eg : 1D array of floats
            glacial erosion of bedrock, meters
        dHdt : 1D array of floats
            change in ice thickness, meters

        """
        # pre-allocate empty variables
        q = np.zeros(self.nodes)       # ice flux

        # instanstiate ice thickness at first node if not there but above ela (otherwise can't grow)
        if self.z[0] > ELA2 and self.HICE[0] == 0:
            self.HICE[0] = 1.
            
        # calculate ice geometry
        topoICE = self.z+self.HICE   # elevation of ice    
        slopeICE = np.append([(topoICE[:-1] - topoICE[1:])/self.dx], [0])  # slope of ice surface
        slopeICE[slopeICE<0] = 0  # negative slopes become zero to prevent instabilities and errors
        
        # update ice thickness based on annual mass balance from eq 19 in Oerlemans, 1984 
        b = self.alpha*(topoICE-ELA2)   # meters
        b[np.where(b>2)] = 2  # limit to 2 meters positive mass balance change
        b[np.where(self.HICE==0)] = 0  # no mass balance where ice is not present

        #add mass balance to ice thickness (ice2)
        self.HICE += dt_g*b
        self.HICE[self.HICE<0] = 0  # another check to make sure negative ice thicknesses aren't created
        
        #glacier cross sectional geometry - between nodes for flux calcs
        H_ICE = np.append([self.HICE[0]], [(self.HICE[:-1] + self.HICE[1:])/2]) # height of ice between nodes for fluxes
        W_bottom = self.Wv  #width at base of ice
        trap_angle = 30 #angle of trapezoidal cross section, from MacGregor et al 2000
        W_top = W_bottom + 2*(H_ICE*np.tan(np.radians(trap_angle)))
        W_avg = (W_bottom+W_top)/2 # average width of glacier at each node
        
        # basal shear stress is modified by cross sectional shape factor F
        F = np.zeros(self.nodes)
        F[H_ICE != 0] = W_top[H_ICE != 0]/(2*H_ICE[H_ICE != 0])
        taubICE = F[:-1]*self.rhoi * self.g * H_ICE[:-1] * slopeICE[:-1]
        
        # glacier water table, after Oerlemans, 1984
        Hw = 0.2*H_ICE # water depth is 20% of ice depth by default
        Hw[H_ICE>100] = H_ICE[H_ICE>100]-75 # unless ice is greater than 100m thick, in which case it is 75 meters below glacier.

        # Ne is effective stress at the bed, imposing a water level that is 75 m less than local ice thickness (after Oerlemans, 1984)
        Ne = self.rhoi * self.g * H_ICE[1:] - self.rhow * self.g * Hw[1:]
        Ne[Ne<0] = 0
        
        # Basal sliding rule in m/yr
        Us = np.zeros(len(taubICE)) #placeholder, since Ne can be zero
        Us[Ne>1e4] = self.C1*taubICE[Ne>1e4]**2/Ne[Ne>1e4] # use only Ne values >1e4 - any lower causes instabilities in ice height to rapidly form and are probably unrealistic
        
        # Ice deformation velocity in m/yr
        Ud = 2/5*self.A*H_ICE[:-1]*taubICE**3
        
        # erosion rule of quarrying from Iverson, 2012 (in this case all eroded material is bedload), using a Gaussian distribution of fracture spacing in a heterogenous bedrock:
        erosion = np.append([self.C2 * Us**0.46], [0])
        
        # erosion checks for unusual situations
        erosion[np.isnan(erosion)] = 0
        erosion[np.isinf(erosion)] = 0
        erosion[np.where(erosion<0)] = 0
        
        """UPDATE WATER - ICE THICKNESS AND WATER DISCHARGE"""
        # ice flux between nodes (q) and updates to glacier height
        q[0] = 0
        q[1:] = W_avg[:-1] * H_ICE[:-1] * (Us + Ud)
        dqdx = (q[1:] - q[:-1])/self.dx
        dHdt = np.append([-1* dqdx/W_avg[:-1]], [0])
        dHdt[np.isnan(dHdt)] = 0
        self.HICE += dHdt*dt_g
        self.HICE[np.where(self.HICE<0)] = 0

        # give or take water from the river. The total glacier mass balance gives the water added/subtracted from Qw:
        if self.glacial_discharge_sw == True:
            mass_change = np.sum(b * W_avg * self.dx)
            water_change = self.rhoi/self.rhow * mass_change * -1  #change the sign because a negative mass_change is a positive Qw
            Qw_change = water_change/(np.pi*10**7) * dt_g
            self.Qw_change += Qw_change 
        
        """UPDATE BEDROCK - EROSION AND SEDIMENT SUPPLY"""
        # Erode the bedrock vertically and laterally
        self.z -= erosion*dt_g
        erosion_valley = erosion * np.cos(np.deg2rad(60)) # based on a trapezoid angle of 30, find horizontal component of erosion into a 30 degree valley slope
        self.Wv += erosion_valley*dt_g
        self.Wv[np.where(self.Wv>500)] = 500.     # cap the glacial width at 500 m following MacGregor et al 2000
        erosion_valley[np.where(self.Wv == 500)] = 0
        
        # Give sediment back to the river. Take mean bedrock erosion (assuming anything eroded off the bed is entrained in the ice and melted out at terminus) plus any sediment heights under now glaciated areas (this latter value should be zero unless glacier is growing). Assume all eroded material is bedload since equation is for plucking/quarrying.
        
        Eg = erosion*dt_g
        Ev = erosion_valley*dt_g
        mean_E = np.mean(Eg+Ev)
        if self.glacial_sed_supply_sw == True:
            # sum the amount of fluvial sediment pushed ('swept') downstream:
            broomvalue = sum(np.where(self.HICE>0, self.sed_depth*self.Wv*self.dx, 0)) * (1-self.L)
            
            # convert to be in kg and add in the eroded material, also converted to kg:
            if np.sum(self.HICE) == 0:
                glacial_sed_supply = 0
            else:
                glacial_sed_supply = self.rhos * (mean_E * self.dx * np.mean(self.Wv[self.HICE > 0]) + broomvalue * (1-self.L))
            
            self.glacial_sed_supply += glacial_sed_supply  # update this value - since dt_g is smaller than dt, we need to do this several times before adding back to the fluvial model.
            
        # update to get rid of sediment where there is ice, since it has now been swept downstream
        self.sed_depth[self.HICE>0] = 0

        return (Eg, dHdt)
    
    
    def calc_isostacy(self, dz, dHICE, dt_i):
        """
        Calculate the isostatic adjustment from ice loads and bedrock erosion. Based on Hubrechts and deWolde (1999) with explanations in Pollard and DeConto (2012) and help from Trevor Hillebrand.

        Parameters
        ----------
        dz : 1D array, float
            change in bedrock thickness from river and glacier erosion, meters
        dHICE : 1D array, float
            change in ice height, meters.
        dt_i : float
            time step of the isostacy model, years

        Returns
        -------
        dz_b_i : 1D array, float
            change in bedrock elevation from isostacy adjustment

        """

        # STEP 1: Find isostatic equilibrium elevations
        z_before = self.z - dz
        # changes in ice and bedrock:
        dsurf = (self.rhoi/self.rhom - 1)*dHICE + (self.rhos/self.rhom - 1)*dz
        z_eq = dsurf - dHICE + z_before
        
        
        # STEP 2: calculate how much isostatic adjustment occurs in timeframe
        
        #define constants:
        L = 132000.     # radius of relative stiffness, assumes 117 km thick crust, in meters. Equals (D/(rhom*g))**1/4
        D = 10**25      # flexural rigidity in N m
        theta = 3000.   # timescale of lithospheric response in years
        dz_b_i = np.zeros((len(self.x),))  # empty array to be filled
        
        # calculate displacements
        Pi = self.Wv * self.dx * dHICE * self.rhoi * self.g       # volumetric displacement of ice
        Pb = self.W * self.dx * dz * self.rhos * self.g      # volumetric displacement of bedrock
        for k in range(0, self.nodes):
            # for each node, calculate the displacement from all the other nodes:
            r = np.absolute(self.x-self.x[k])
            Wp = (Pi + Pb)*L**2/(2*np.pi*D) * kei(r/L)
            Wb = np.sum(Wp)
            dz_b_i[k] = -1/theta * (self.z[k]-z_eq[k]+Wb) * dt_i

        return dz_b_i  
    
    
    def calc_sediment_erosion(self, dt):
        """
        Calculates sediment transport using a modified Meyer-Peter and Mueller equation in Wong and Parker (2006)

        Parameters
        ----------
        dt : integer
            timestep, years

        Returns
        -------
        dz_s : 1d array of floats
            change in sediment thickness, positive = erosion

        """
        
        # calculate shear stress products, making exceptions for where ice is present
        tau_star = np.zeros(self.nodes)
        tau_star[self.HICE == 0] = self.taub[self.HICE == 0]/((self.rhos-self.rhow)*self.g*self.D[self.HICE == 0])
        tau_diff = tau_star - self.tauc
        tau_diff[np.where(tau_diff<0)] = 0
        
        # calculate sediment supply and capacity, set at 0 where ice is present
        self.Qs_n = self.Qs_hills + self.Qs_down[:-1]
        self.Qs_cap = self.frac_yr_transport * self.rhos * ((self.rhos-self.rhow)/self.rhow * self.g * self.D**3) ** (0.5) * 3.97 * ((tau_diff)**(3/2)) * self.W
        self.Qs_n[np.where(self.HICE>0)] = 0
        self.Qs_cap[np.where(self.HICE > 0)] = 0
        
        # calculate change in elevation based on sediment capacity and supply difference; sediment eroded comes from the channel, and sediment deposited is spread over valley
        mass_change = dt * self.Qs_cap * np.pi*10**7 - self.Qs_n
        dz_s = (1/(1-self.L)) * mass_change/(self.Wv * self.rhos * self.dx)
        erosion = np.where(mass_change>0)[0]
        dz_s[erosion] = (1/(1-self.L)) * mass_change[erosion]/(self.W[erosion] * self.rhos * self.dx)
        
        # calculate sediment passed to next node
        Qs_pass = mass_change
        erodes_all = dz_s - self.sed_depth
        eroding_all = np.where((erodes_all > 0) & (self.sed_depth > 0))[0]
        dz_s[eroding_all] = self.sed_depth[eroding_all]
        Qs_pass[eroding_all] = self.Qs_n[eroding_all]
        dz_s[np.where(dz_s > self.sed_depth)] = self.sed_depth[np.where(dz_s > self.sed_depth)]        

        # now add in what has been eroded or deposited to get the total sed downstream and update Qs_down for next timestep
        self.Qs_down[1:] = Qs_pass * 1
        self.Qs_down[self.Qs_down<0] = 0
        
        # update the sediment height:
        self.sed_depth -= dz_s
        self.sed_depth[-1] = 0
                    
        return dz_s
      

    def run_one_fluvial(self, dz_b_save, backgroundU, depth_threshold, dt, erosion_type = 'depth_threshold'):
        """
        Runs one time step of fluvial erosion with sediment transport and a depth threshold.

        Parameters
        ----------
        dz_b_save : 2D array, float
            10 ky saved history of erosion rates, meters per year.
        backgroundU : float
            uplift rate, meters per year.
        depth_threshold : float
            alluvium depth that prevents erosion, meters.
        dt : integer
            timestep, years
        erosion_type: text
            Specifies the alluvial cover effect model to be used
            'depth threshold' is binary yes or no erosion
            'SklarDietrich' uses 1-Qs/Qt from their 2004 paper
            'Turowski' uses exp(Qs/Qt) from their 2007 paper where alpha = 1
            'Shobe' uses exp(h/h*) from their 2017 paper where h* = depth_threshold, set at 1 meter

        Returns
        -------
        dz_s : 1D array, floats
            change in sediment height, meters. Positive is erosion
        dz_b : 1D array, floats
            change in bedrock height, meters. Negative is erosion.
        dz_w : 1D array, floats
            lateral movement of channel, meters. Positive is erosion

        """

        sed_supply = np.mean(dz_b_save, axis=1)
        self.get_sediment_in(sed_supply, dt)
        
        # add glacial inputs from glacial erosion
        first_fluvial = np.min(np.where(self.HICE == 0))
        self.Qs_hills[first_fluvial] += self.glacial_sed_supply/self.Wv[first_fluvial]*self.W[first_fluvial] # spread glacial pushed sed across valley and only transport that in channel. 
        self.Qw[:first_fluvial] = 0
        self.Qw[first_fluvial:] = self.Qw_change + self.Qwi[first_fluvial:]
        
        # update river geometries to prep for transport calculations
        self.get_slope()
        self.get_water_height(self.Qw)
        self.get_taub()
        
        # calculate changes to sediment depth, bedrock erosion, and bank erosion
        dz_s = self.calc_sediment_erosion(dt) # updates self.sed_depth w/in function
        if erosion_type == 'SklarDietrich': 
            dz_b = self.calc_bedrock_erosion_sklardietrich(dt)
        elif erosion_type == 'Turowski':
            dz_b = self.calc_bedrock_erosion_turowski(dt)
        elif erosion_type == 'Shobe':
            dz_b = self.calc_bedrock_erosion_shobe(depth_threshold, dt)
        else:
            dz_b = self.calc_bedrock_erosion(dt, depth_threshold)
        
        self.z += dz_b
        dz_w = self.calc_bank_erosion(dt)
        
        ## convert bank dz to values and update sediment supply
        bank_erosion = dz_w * self.H * self.dx * self.rhos
        self.Qs_down[1:] += bank_erosion

        return dz_s, dz_b, dz_w



    def run_one_fluvial_nosed(self, backgroundU, dt):
        """
        Runs one timestep of fluvial erosion without sediment transport or a depth threshold.

        Parameters
        ----------
        backgroundU : float
            rock uplift rate, meters per year.
        dt : integer
            timestep in years.

        Returns
        -------
        dz_b : 1D array of floats
            bedrock erosion, negative is erosion.

        """

        # add glacial inputs from glacial erosion
        first_fluvial = np.min(np.where(self.HICE == 0))
        self.Qw[:first_fluvial] = 0
        self.Qw[first_fluvial:] = self.Qw_change + self.Qwi[first_fluvial:]

        # update river geometries to prep for transport calculations
        self.get_slope()
        self.get_water_height(self.Qw)
        self.get_taub()
        
        # calculate change to bedrock erosion
        depth_threshold = 0.2 # no sediment anyways, so value doesn't matter
        dz_b = self.calc_bedrock_erosion(dt, depth_threshold)
        self.z += dz_b

        return dz_b  

     
        
    
    def run_one_glacial(self, analysis_time, averageELA, amplitude, dt_g, dt):
        """
        Runs one time step of glacial erosion

        Parameters
        ----------
        analysis_time : integer
            year in the model run, year
        averageELA : integer
            elevation of average ELA, meters
        amplitude : integer
            range of ELA, meters
        dt_g : float
            glacier timestep, 0.01 or less, years
        dt : integer
            model timestep

        Returns
        -------
        Eg_total : 1D array of floats
            bedrock erosion from glaciation for one model timestep, meters
        ELA : float
            elevation of ELA, meters

        """
        # set empty variables to track. time2 represents time within the larger model timestep
        time2 = 0
        Eg_total = np.zeros(self.nodes)
        self.Qw_change = 0
        self.glacial_sed_supply = 0
        
        ELA = self.calc_ELA(analysis_time, averageELA, amplitude, period = 100000, shape = 'sawtooth')
            
        while time2 < dt:
            Eg, dHdt = self.calc_glacial_erosion(ELA, dt_g)
            Eg_total += Eg
            time2 += dt_g
            
        return Eg_total, ELA
    

     

