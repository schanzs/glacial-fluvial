#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:49:50 2019

@author: saschanz
"""

""" stream functions that accompany glacial-fluvial model """
import numpy as np
from scipy.optimize import fsolve
from scipy.special import kei


class stream(object):
    g = 9.81
    rhow = 1000
    rhoi = 917
    rhos = 2650
    rhom = 3300
    beta = 0.2
    R = rhos/rhow - 1
    tauc = 0.0495
    n = 0.04
    frac_yr_transport = 0.05
    L = 0.2
    kf = -10**-5
    kw = -7.3*kf
    A = 2.1*10e-18  #Arrhenius constant
    C1 = 0.0012 # sliding coefficient
    C2 = 0.000001 # erosion coefficient
    alpha = 0.01 # meters/ (year meters) from Oerlemans, 1984
    
    dt_min = 0.1
    
    def __init__(self, dx, nodes, initial_slope, initial_sed_depth):
        self.nodes = nodes
        self.dx = dx
        self.initial_slope = initial_slope
        
        self.x = np.arange(0, nodes*dx, dx)
        self.z = initial_slope * (max(self.x)-self.x)
        self.z[-1] = 0
        
        self.sed_depth = np.ones(nodes)*initial_sed_depth

        
        # set up empties
        self.HICE = np.zeros(nodes)
        self.Qs_down = np.zeros(nodes+1)
        self.H = np.zeros(nodes)
    
    def get_basin_geometry(self):
        Hc = 1
        He = 1.8598
        self.Ah = Hc * self.x**He + 0.1
        
        kQ = 10**-7
        eQ = 1
        self.Qwi = kQ * self.Ah**eQ
        
        kw = 5.
        eW = 0.5
        self.W = kw*self.Qwi**eW
        self.Wv = 1*self.W[:]
        
    def get_slope(self):
        topo = self.sed_depth + self.z
        self.slope = np.append([np.diff(topo)/np.diff(self.x)], [0])
        self.slope[np.where(self.slope>0)] = 0
      
        
        
    def get_water_height(self, Qw):
        self.H[:-1] = (self.n * Qw[:-1]/self.W[:-1])**.6 * (-1*self.slope[:-1])**(-0.3)
        self.H[-1] = self.H[-2]
        self.H[np.where(self.slope == 0)] = 0
    
   
    def get_sediment_size(self, D0, a):
        self.D = D0 * np.exp(-self.x*a)
        self.D[np.where(self.D<0.005)] = 0.005
    
    
    
    def get_fall_velocity(self):
        visc = 0.001/1000
        Dstar = (self.rhos-self.rhow)*self.g*self.D**3/(self.rhow*visc**2)
        wstar = np.zeros((self.nodes))
        # from Dietrich, 1982:
        lgW=-3.76715+(1.92944*np.log10(Dstar))-(0.09815*((np.log10(Dstar))**2))-(.00575*((np.log10(Dstar))**3))+(.00056*((np.log10(Dstar))**4))
        wstar[np.where(Dstar>0.05)] = 10**lgW[np.where(Dstar>0.05)]
        wstar[np.where(Dstar<=0.05)] = 1.71e-4*self.D[np.where(Dstar<=0.05)]**2
        self.wf = (wstar*(self.rhos-self.rhow)*self.g*visc/self.rhow)**(1/3)
        
        
        
    def get_sediment_in(self, dz_b, dt):
        dA = np.hstack((self.Ah[0], np.diff(self.Ah)))
        self.Qs_hills = self.rhos * self.beta * dA * dz_b * dt
      
        
        
    def get_taub(self):
        self.taub = self.rhow * self.g * (self.W * self.H / (2*self.H + self.W)) * -1 * self.slope
        self.taub[np.where(self.HICE > 0)] = 0
        self.taub[np.where(self.slope == 0)] = 0
        
     
    def calc_sediment_erosion(self, dt):
        """ Calculates sediment transport based on a modified Meyer Peter and Mueller equation in Wong and Parker (2006)
        Parameters
        ---------
        dt = time step (years)
        """
        
        tau_star = self.taub/((self.rhos-self.rhow)*self.g*self.D)
        tau_diff = tau_star - self.tauc
        tau_diff[np.where(tau_diff<0)] = 0
        self.Qs_n = self.Qs_hills + self.Qs_down[:-1]
        self.Qs_cap = self.frac_yr_transport * self.rhos * ((self.rhos-self.rhow)/self.rhow * self.g * self.D**3) ** (0.5) * 3.97 * ((tau_diff)**(3/2)) * self.W
        
        self.Qs_n[np.where(self.HICE>0)] = 0
        self.Qs_cap[np.where(self.HICE > 0)] = 0
        
        mass_change = dt * self.Qs_cap * np.pi*10**7 - self.Qs_n
        dz_s = (1-(1-self.L)) * mass_change/(self.Wv * self.rhos * self.dx)
        erosion = np.where(mass_change>0)[0]
        dz_s[erosion] = (1-(1-self.L)) * mass_change[erosion]/(self.W[erosion] * self.rhos * self.dx)
        
        Qs_pass = mass_change
        erodes_all = dz_s - self.sed_depth
        eroding_all = np.where((erodes_all > 0) & (self.sed_depth > 0))[0]
        dz_s[eroding_all] = self.sed_depth[eroding_all]
        Qs_pass[eroding_all] = self.Qs_n[eroding_all]
        dz_s[np.where(dz_s > self.sed_depth)] = self.sed_depth[np.where(dz_s > self.sed_depth)]        

        # now add in what has been eroded or deposited to get the total sed downstream
        # change Qs_down
        self.Qs_down[1:] = (self.rhos*dz_s * self.Wv * self.dx) + Qs_pass
        self.Qs_down[erosion+1] = (self.rhos * dz_s[erosion] * self.W[erosion] * self.dx) + Qs_pass[erosion]
        self.Qs_down[np.where(self.Qs_down<0)] = 0

        return dz_s


      
    def calc_bedrock_erosion(self, dt, erosion_type=1, depth_threshold=0.2):
        """Calculates bedrock erosion based on stream power. Placeholder for saltation abrasion is put in. 
        Parameters
        ---------
        dt = model timestep (year)
        erosion_type = type of erosion model, 1 = stream power, 2 = saltation_abrasion
        depth_threshold = depth of sediment to protect bed from erosion (meters)
        """
        
        def stream_power(dt):
            return self.kf * self.taub * dt
        
        def saltation_abrasion(dt):
            return dz_b
        
        if erosion_type == 1:
            dz_b = stream_power(dt)
        if erosion_type == 2:
            dz_b = saltation_abrasion(dt)
        
        dz_b[np.where(self.sed_depth > depth_threshold)] = 0
        dz_b[np.where(self.HICE > 0)] = 0
        
        return dz_b
   

    def calc_bank_erosion(self, dt):
        """ Calculates the stream bank erosion using stream power and a lateral erodibility coefficient. Keeps track of whether valley is widened previously by glaciation to use a different erodibility.
        Parameters
        ---------
        dt = time step (year)
        """
        
        dz_w = self.kw * self.taub * dt
        dz_w[np.where(self.Wv > self.W)] = 10*dz_w[np.where(self.Wv > self.W)]
        
        dz_w[np.where(self.HICE > 0)] = 0
        
        return dz_w



    def isostacy(self, dz, dHICE, dt_i):
        """ Calculate the isostatic adjustment from ice loads and bedrock erosion. Based on Hubrechts and deWolde (1999) with explanations in Pollard and DeConto (2012) and help from Trevor Hillebrand.
        Parameters
        ----------
        dz = erosion of bedrock by river (m)
        dHICE = change in ice thickness (m)
        dt_i = time step of the isostacy model (year)
        """

        # find equilibrium elevations
        z_before = self.z - dz
        
        # changes in ice and bedrock:
        dsurf = (self.rhoi/self.rhom - 1)*dHICE + (self.rhos/self.rhom - 1)*dz
        z_eq = dsurf - dHICE + z_before
        
        
        # define constants:
        L = 132000.     # radius of relative stiffness, assumes 117 km thick crust, in meters. Equals (D/(rhom*g))**1/4
        D = 10**25      # flexural rigidity in N m
        theta = 3000.   # timescale of lithospheric response in years
        dx = np.mean(np.diff(self.x))
        dz_b_i = np.zeros((len(self.x),))
        Pi = self.Wv * dx * dHICE * self.rhoi * self.g       # volumetric displacement of ice
        Pb = self.W * dx * dz * self.rhos * self.g      # volumetric displacement of bedrock
        for k in range(0, self.nodes):
            # for each node, calculate the displacement from all the other nodes:
            r = np.absolute(self.x-self.x[k])
            Wp = (Pi + Pb)*L**2/(2*np.pi*D) * kei(r/L)
            Wb = np.sum(Wp)
            dz_b_i[k] = -1/theta * (self.z[k]-z_eq[k]+Wb) * dt_i

        return (dz_b_i)     
    

    def get_ELA(self, time, dt, averageELA = 2000, amplitude = 500, period = 100000, shape = 'sawtooth'):
        """ Function calculates the ELA for each time step 
        Parameters:
        ---------
        time = model year (years)
        dt = time step of model (years)
        averageELA = average elevation of ELA (meters)
        amplitude = variation in ELA (meters)
        period = how often glacial cycles occur (years)
        shape = shape of ela curve - sine or sawtooth
        """
    
    # ELA OPTION 1: sinusoidal, equal time in glacial and interglacial    
    # calculate ELA elevation based on time variation sinusoid. 1000 m amplitude is approximately a 6 degree mean annual temperature shift (Porter, 1989), noted in Sierra Nevada glaciers (Phillips et al, 1996).
        if shape == 'sine':
            ELA = -amplitude/2*np.cos(2*np.pi*time*dt/period)+averageELA 

    # ELA OPTION 2: sawtooth, 10k rapid warming and 90k slow cooling
    # total change in temperature is 6 degrees based on Porter, 1989, corresponding to 1000 m amplitude in ELA.
        if shape == 'sawtooth':
            cycle_i = (time*dt)%period
            if cycle_i <= 10000: # start with rapid warming 
                ELA = amplitude/10000*cycle_i + (averageELA-amplitude/2) #create linear equation with slope of 1000m/10,000 years and a value of 2000 at 5,000 years
            else:
                ELA = -amplitude/90000*cycle_i + (averageELA+amplitude/2 + amplitude/9)
                
        return ELA
        

    def run_one_glacial(self, ELA2, dt_g):
        """ Longitudinal glacier evolution following Macgregor et al., 2000 (Geology) formulation, with updates to the quarrying rule from Iversion (2012).
        
        Parameters:
        ----------
        ELA2 = ELA incorporating both the annual and glacial cycle temperature swings
        dt_g = time step for glaciers, in years
        
        """
        xcells = self.nodes

        # pre-allocate
        q = np.zeros(xcells)
        H_ICE=np.zeros(xcells)
        
        topoICE = self.z+self.HICE

        if self.z[0] > ELA2 and self.HICE[0] == 0:
            self.HICE[0] = 1.
        #annual mass balance from eq 19 in Oerlemans, 1984 - this is the snowfall and melt contributions?
        b = self.alpha*(topoICE-ELA2)
        b[np.where(b>2)] = 2
        b[np.where(self.HICE==0)] = 0

        #add mass balance to ice thickness (ice2)
        self.HICE += dt_g*b
        self.HICE[np.where(self.HICE<0)] = 0
        
        #glacier geometry
        # slopeICE[0] = (topoICE[0] - topoICE[1])/self.dx
        H_ICE[0] = self.HICE[0]
        H_ICE[1:] = (self.HICE[:-1] + self.HICE[1:])/2
        W_bottom = self.Wv  #width at base of ice
        trap_angle = 30 #angle of trapezoidal cross section
        W_top = W_bottom + 2*(H_ICE/np.tan(np.radians(trap_angle)))
        W_avg = (W_bottom+W_top)/2
        
        slopeICE = np.zeros_like(topoICE)
        slopeICE[:-1] = (topoICE[:-1] - topoICE[1:])/self.dx
        
        # basal shear stress is modified by cross sectional shape factor F
        F = W_top/(2*H_ICE)
        F[np.where(H_ICE==0)] = 0
        taubICE = F[:-1]*self.rhoi * self.g * H_ICE[:-1] * slopeICE[:-1]
        
        # glacier water table, after Oerlemans, 1984
        Hw = H_ICE - 75
        
        # Ne is effective stress at the bed, imposing a water level that is 75 m less than local ice thickness (after Oerlemans, 1984)
        Ne = self.rhoi * self.g * H_ICE[1:] - self.rhow * self.g * 0.2 * Hw[1:]
        
        # Basal sliding rule in m/yr
        Us = self.C1*taubICE**2/Ne
        
        # Ice deformation velocity in m/yr
        Ud = 2/5*self.A*H_ICE[:-1]*taubICE**3
        
        # erosion rule, proportional to sliding velocity. No erosion at last node. Lots of checks for where there is no ice and Us is either Inf or NaN. 
        #erosion = np.append([self.C2 * Us], [0])
        # erosion rule of quarrying from Iverson, 2012 (in this case all eroded material is bedload), using a Gaussian distribution of fracture spacing in a homogeneous bedrock:
        erosion = np.append([self.C2 * Us**1.86], [0])
        erosion[np.isnan(erosion)] = 0
        erosion[np.isinf(erosion)] = 0
        erosion[np.where(erosion<0)] = 0
        
        
        # ice flux at each node
        q[0] = 0
        q[1:] = W_avg[:-1] * H_ICE[:-1] * (Us + Ud)
        dqdx = (q[1:] - q[:-1])/self.dx
        dHdt = np.append([b[:-1] - dqdx/W_avg[:-1]], [0])
        dHdt[np.isnan(dHdt)] = 0
        self.HICE = self.HICE + dHdt*dt_g
        self.HICE[np.where(self.HICE<0)] = 0

        # give water to the river. If glacier is losing mass, convert the mass lost in each time step to meltwater released.
        mass_loss = np.sum(b[np.where(b<0)] * W_avg[np.where(b<0)] * self.dx)
        if mass_loss < 0:
            # convert to water from ice
            water_released = self.rhoi/self.rhow * mass_loss  * -1
            # in cubic meters per year: convert to cm/s and only for timestep
            Qw_melt = water_released/(np.pi*10**7) * dt_g
        else:
            Qw_melt = 0

        
        #Erode the bedrock vertically and laterally
        self.z -= erosion*dt_g
        erosion_valley = erosion * np.cos(np.deg2rad(60)) # based on a trapezoid angle of 30, find   horizontal component of erosion into a 30 degree valley slope
        self.Wv += erosion_valley*dt_g
        self.Wv[np.where(self.Wv>500)] = 500.     # cap the glacial width at 500 m
        erosion_valley[np.where(self.Wv == 500)] = 0
        
        # Give sediment back to the river. Take mean bedrock erosion (assuming anything eroded off the bed is entrained in the ice and melted out at terminus) plus any sediment heights under now glaciated areas (this latter value should be zero unless growing). Assume all eroded material is bedload since equation is for plucking/quarrying.
        Eg = erosion*dt_g
        Ev = erosion_valley*dt_g
        mean_E = np.mean(Eg+Ev)
        broomvalue = sum(np.where(self.HICE>0, self.sed_depth*self.Wv*self.dx, 0)) * (1-self.L)
        # needs to be in kg:
        if np.sum(self.HICE) == 0:
            glacial_sed_supply = 0
        else:
            glacial_sed_supply = self.rhos * (mean_E * self.dx * np.mean(self.Wv[np.where(self.HICE > 0)]) + broomvalue * (1-self.L))
        
        # get rid of sediment where there is ice
        self.sed_depth[np.where(self.HICE>0)] = 0

        
        # totalice[i] = np.sum(ice2)*self.dx*W_avg
        # glac_sed = Eg*np.sum(np.diff(A))

        
        return (Eg, Qw_melt, dHdt, glacial_sed_supply)

    def run_one_fluvial(self, dz_b_save, backgroundU, dt):
        """ runs one time step of fluvial erosion 
        Parameters
        ----------
        dz_b_save = last 10ky of bedrock erosion, in m/yr
        backgroundU = uplift rate, in m/yr
        dt = model timestep, years
        """
        
        sed_supply = np.mean(dz_b_save, axis=1)
        sed_supply[np.where(sed_supply<backgroundU)] = backgroundU
        self.get_sediment_in(sed_supply, dt)
        
        self.get_slope()
        self.get_water_height(self.Qwi)
        self.get_taub()
        
        dz_s = self.calc_sediment_erosion(dt)
        dz_b = self.calc_bedrock_erosion(dt, erosion_type=1, depth_threshold=0.2)
        dz_w = self.calc_bank_erosion(dt)
        
        ## convert dz to values
        bank_erosion = dz_w * self.H * self.dx * self.rhos
        self.Qs_down[1:] += bank_erosion

        
        return dz_s, dz_b, dz_w
        
        