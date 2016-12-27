#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 22:45:55 2016

@author: throop
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
import astropy.table as table
from astropy.coordinates import SkyCoord
import math
import hbt
import spiceypy as sp
import astropy.units as u
import astropy.constants as c
    
file_tm = '/Users/throop/git/NH_rings/kernels_nh_pluto_mu69.tm'

# Read HD into a table. Its positions are 1950. Coords are radians.

t = read_hd_all()  

# Initialize 

sp.furnsh(file_tm)

utc_ca = '2019 1 jan 03:09:00'
et_ca  = sp.utc2et(utc_ca)
fov_lorri = 0.3*hbt.d2r  # Radians

hour = 3600
day  = 24 * hour

plot_tick_every = 360  # Plot a time-tick every __ seconds

num_dt = 2000 # Number of timesteps

mag_limit = 8

et_start = et_ca + 0.3*hour
et_end   = et_ca + 96*hour

et = hbt.frange(et_start, et_end, num_dt)

dt = et[1] - et[0] # Timestep

# Precess stars from 1900 to 2000

print("Precessing to J2000")

mx = sp.pxform('B1950', 'J2000', et[0])

ra_1950 = t['RA']
dec_1950 = t['Dec']

num_stars = np.size(ra_1950)

ra_2000 = []
dec_2000 = []
pt_1950 = []

# Loop over every star and precess it individually. This routine is a real bottleneck.

for i in range(num_stars):
    
  pt_1900 = sp.radrec(1, ra_1950[i], dec_1950[i])
  pt_1950 = sp.mxv(mx, pt_1900)
  pt_2000 = sp.mxv(mx, pt_1950) # SPICE will not precess 1900 -> 2000. So we apply 1950 -> 2000, do it 2x.
  d, ra_2000_i, dec_2000_i = sp.recrad(pt_2000)
  dec_2000.append(dec_2000_i)
  ra_2000.append(ra_2000_i)  
  
t['RA_2000'] = ra_2000
t['Dec_2000'] = dec_2000
  
ra_kbo = []
dec_kbo = []
dist_kbo = []

color_kbo = 'pink'
color_lorri = 'lightgreen'
color_roche = 'lightblue'

# Look up NH position

print("Looking up NH position...")

for et_i in et:

    (st,lt) = sp.spkezr('MU69', et_i, 'J2000', 'LT+S', 'New Horizons')
    (dist_i, ra_i, dec_i) = sp.recrad(st[0:3])
    
    ra_kbo.append(ra_i)    # MU69 RA/Dec, in radians
    dec_kbo.append(dec_i)
    dist_kbo.append(dist_i)
#    dist_kbo.append(dist_i)

ra_kbo = np.array(ra_kbo)
dec_kbo = np.array(dec_kbo)
ra = t['RA_2000']*hbt.r2d           # Look up stellar coordinates - be sure to use J2000 
dec = t['Dec_2000']*hbt.r2d
dist_kbo = np.array(dist_kbo)*u.km  # NH-MU69 distance

# Compute the size of the Hill sphere

rho_kbo =  2.5*u.gram/u.cm**3
r_kbo   = 16.5 * u.km
m_kbo   = 4/3 * math.pi * r_kbo**3 * rho_kbo
a_kbo = 43*u.au

a_hill = (a_kbo * (m_kbo / c.M_sun/3)**(1/3)).to('km')

angle_hill = (a_hill / dist_kbo).to('').value # Convert to radians, and drop the units
angle_hill = np.clip(angle_hill, -math.pi, math.pi) # Clip it (same as IDL > or < operator)

angle_kbo  = (r_kbo / dist_kbo).to('').value # Radians

a_roche = 2.5 * r_kbo
angle_roche = (a_roche / dist_kbo).to('').value # Radians

xlim = hbt.mm(ra_kbo*hbt.r2d)   # Plotting limits, in deg
ylim = hbt.mm(dec_kbo*hbt.r2d)

pad_ra_deg  = 2 # Padding at edge of plot, RA = x dir. Degrees.
pad_dec_deg = 2 # Padding at edge of plot, Dec = y dir. Degrees.

# Pad the edges fof this a bit
xlim = np.array(xlim) + pad_ra_deg * np.array([-1,1])
ylim = np.array(ylim) + pad_dec_deg * np.array([-1,1])

is_good =  np.logical_and(
                            t['Ptg'] < mag_limit,
               
                            np.logical_and( 
                                np.logical_and(xlim[0] < ra, 
                                               ra < xlim[1]),     
                                np.logical_and(ylim[0] < dec,
                                               dec < ylim[1])  )    )

hbt.figsize((15,10))
hbt.set_fontsize(15)

#==============================================================================
# Make the plot
#==============================================================================

# Draw the main trajectory

plt.plot(ra_kbo*hbt.r2d, dec_kbo*hbt.r2d)

# Draw the LORRI FOVs -- 0.3 x 0.3 deg square

dec_kbo_plus_lorri  = dec_kbo + (fov_lorri/2) * np.sqrt(2)  # All radians
dec_kbo_minus_lorri = dec_kbo - (fov_lorri/2) * np.sqrt(2)

#plt.plot(ra_kbo*hbt.r2d, dec_kbo_plus_lorri*hbt.r2d,  linestyle = '--', color=color_lorri)
#plt.plot(ra_kbo*hbt.r2d, dec_kbo_minus_lorri*hbt.r2d, linestyle = '--', color=color_lorri)

plt.fill_between(ra_kbo*hbt.r2d, (dec_kbo_plus_lorri)*hbt.r2d, (dec_kbo_minus_lorri)*hbt.r2d, 
                 color=color_lorri, alpha=0.5, label = 'LORRI FOV')


# Draw the Hill sphere

plt.plot(ra_kbo*hbt.r2d, (dec_kbo - angle_hill)*hbt.r2d, linestyle = '--')
plt.plot(ra_kbo*hbt.r2d, (dec_kbo + angle_hill)*hbt.r2d, linestyle = '--')

# Draw the Roche radius

#plt.plot(ra_kbo*hbt.r2d, (dec_kbo + angle_roche)*hbt.r2d, linestyle = '--', color = color_roche)
#plt.plot(ra_kbo*hbt.r2d, (dec_kbo - angle_roche)*hbt.r2d, linestyle = '--', color = color_roche)

plt.fill_between(ra_kbo*hbt.r2d, (dec_kbo+angle_roche)*hbt.r2d, (dec_kbo-angle_roche)*hbt.r2d, 
                 color=color_roche, alpha=1, label = 'MU69 Roche radius')


# Draw the MU69 radius

#plt.plot(ra_kbo*hbt.r2d, (dec_kbo + angle_kbo)*hbt.r2d, linestyle = '--', color=color_kbo)
#plt.plot(ra_kbo*hbt.r2d, (dec_kbo - angle_kbo)*hbt.r2d, linestyle = '--', color=color_kbo)

plt.fill_between(ra_kbo*hbt.r2d, (dec_kbo+angle_kbo)*hbt.r2d, (dec_kbo-angle_kbo)*hbt.r2d, 
                 color=color_kbo, alpha=1, label = 'MU69 radius')

# Draw '+' symbols every hour

for i,et_i in enumerate(et):
    
    if (np.mod(et_i - et[0], plot_tick_every)) < dt:
        if (i == 0): # Pass the info for this to plt.label()... but only for the first call!
            kwargs = {'label' : 'MU69 position'}
        else:
            kwargs = {}
        plt.plot(ra_kbo[i]*hbt.r2d, dec_kbo[i]*hbt.r2d, marker='+', 
                 markersize=20, color='black', **kwargs )

# Plot the stars

plt.plot(t['RA_2000'][is_good]*hbt.r2d, t['Dec_2000'][is_good]*hbt.r2d, linestyle='none', marker='.', 
         label = 'HD star')

# Label the stars

for i in range(np.size(is_good)):
    if (is_good[i]):
        plt.text(t['RA_2000'][i]*hbt.r2d, t['Dec_2000'][i]*hbt.r2d, 
            '  {}  {:.1f} {}'.format(t['Type'][i], t['Ptg'][i], t['HD'][i]),
#            '  {}  {:.1f}'.format(t['Type'][i], t['Ptg'][i]),

            fontsize = 8)
        
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
plt.title('MU69 Outbound, HD Stars V < {}. Ticks every {} min. K+{}h .. K+{}h'.format(
          mag_limit, plot_tick_every/60, (et_start - et_ca)/hour, (et_end-et_ca)/hour))
plt.legend(framealpha=0.5)
plt.show()

#plt.set_xlim(xlim)

#==============================================================================
# Do some test cases
#==============================================================================
