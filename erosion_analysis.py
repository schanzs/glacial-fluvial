#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 10:12:29 2019

@author: saschanz
"""

"""
plot dz_b for the glacial-fluvial runs
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


# MODEL TYPE?
ELA_type = 'sawtooth'


# TIME AND SPACE
x = np.arange(0, 150000, 1000)
time = np.arange(1, 800, 1)*1000
dt = 1000

# CALCULATE ELA
amplitude = 1000; period = 100000; averageELA = 2000
time2 = np.arange(0, 100000, 1000)
if ELA_type == 'sine':
    ELA = -amplitude/2*np.cos(2*np.pi*time2/period)+averageELA 

if ELA_type == 'sawtooth':
    cycle_i = (time2)%period
    ELA = -amplitude/90000*cycle_i + (averageELA+amplitude/2 + amplitude/9)
    ELA[np.where(cycle_i <=10000)] = amplitude/10000*cycle_i[np.where(cycle_i <= 10000)] + (averageELA-amplitude/2)

# CALCULATE REFERENCE CONCAVITY
x = np.arange(0, 200000, 2000)
Hc = 1
He = 1.8598
DA = Hc * x**He + 0.1

z_in = np.genfromtxt('z_000000.txt')
sed_in = np.genfromtxt('sed_000000.txt')
slope = -1*np.diff(z_in+sed_in)/np.diff(x)
ref_concavity, ksn = np.polyfit(np.log10(DA[:-1]), np.log10(slope), 1)
ksn_init = slope/(DA[:-1]**ref_concavity)

# CALCULATE EROSION  
z = np.empty((len(x), len(time)))
sed = np.empty_like(z)
ice = np.empty_like(z)
dzb = np.empty_like(z)
ksn = np.empty((len(x)-1, len(time)))
ksn_diff = np.empty_like(ksn)
eg = np.empty_like(z)
wv = np.empty_like(z)

for year in time:
    i = int(year/dt-1)
    z[:,i] = np.genfromtxt('z_%06d.txt' % (year))
    sed[:,i] = np.genfromtxt('sed_%06d.txt' % (year))
    ice[:,i] = np.genfromtxt('ice_%06d.txt' % (year))
    dzb[:,i] = np.genfromtxt('dzb_%06d.txt' % (year))/dt
    eg[:,i] = np.genfromtxt('eg_%06d.txt' % (year))/dt
    wv[:,i] = np.genfromtxt('wv_%06d.txt' % year)
    slope = -1*np.diff(z[:,i]+sed[:,i])/np.diff(x)
    ksn[:,i] = slope/(DA[:-1]**ref_concavity)
    ksn_diff[:,i] = ksn_init - slope/(DA[:-1]**ref_concavity)
    
#dzb = np.diff(z, axis=1)/dt # m/ky
dzs = np.diff(sed, axis=1)/dt # m/ky
dzi = np.diff(ice, axis = 1)/dt # m/ky
dzb[(dzb <= 0)] = np.nan
eg[(eg <= 0)] = np.nan
dzi[(dzi == 0)] = np.nan
dzs[(dzs == 0)] = np.nan
sed[sed < 0.2] = np.nan
ice[ice == 0] = np.nan
ksn_diff[np.abs(ksn_diff)<0.5] = np.nan

# PLOT
tin, xin = np.meshgrid(time[:], x[:60]/1000)
fig, ax = plt.subplots(4,1, sharex=True)
minor_ticksx = np.arange(0, 800, 50)*1000
minor_ticksy = np.arange(0, 61, 10)
plotting = (dzb,
            eg,
            sed,
            ice,
            wv,
            ksn_diff)

colormap = ('pink',
            'pink',
        'pink',
        'pink',
        'pink',
        'BrBG')

vminv = (np.nanmin(plotting[0]),
        np.nanmin(plotting[1]),
        np.nanmin(plotting[2]),
        np.nanmin(plotting[3]),
        np.nanmin(plotting[4]),
        -1*np.nanmax(np.abs(plotting[5])))

vmaxv = (np.nanmax(plotting[0]),
        np.nanmax(plotting[1]),
        np.nanmax(plotting[2]),
        np.nanmax(plotting[3]),
        np.nanmax(plotting[4]),
        np.nanmax(np.abs(plotting[5])))

for axis, i in zip(ax, np.arange(0, 6)):
    name = axis.pcolormesh(tin, xin, plotting[i][:60,:], cmap = colormap[i], vmax = vmaxv[i], vmin = vminv[i])
    axis.set_xticks(minor_ticksx, minor=True)
    axis.set_yticks(minor_ticksy, minor=True)
    axis.grid(which='minor', alpha = 0.5, linestyle=':', color='darkgray')
    axis.grid(which='major', alpha = 0.7, linestyle = '-.', color='darkgray') 
    fig.colorbar(name, ax=axis)
    
ax[5].set_xlabel('Time (yr)')
ax[0].set_title('Bedrock erosion (m/ky)')
ax[1].set_title('Glacial erosion (m/ky)')
ax[2].set_title('Sediment thickness (m)')
ax[3].set_title('Ice thickness (m)')
ax[4].set_title('Valley width (m)')
ax[5].set_title('initial Ksn - Ksn with ref concavity of: %s' % np.round(ref_concavity, 2))
