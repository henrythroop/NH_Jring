#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 10:07:07 2017

@author: throop
"""

# This program makes a plot of the mosaic for MU69 rings observations.
# This is not trivial because the target distance is changing rapidly!
# This plot goes into the Wiki, for planning the MU69 K+30 min 
# outbound rings obs, 
# https://www.spaceops.swri.org/nh/wiki/index.php/KBO/MT1.3-G4-1.2
#
# HBT 2-Apr-2017

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy
from   astropy.table import Table
import astropy.table as table
from   astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import math
import hbt
import spiceypy as sp
import astropy.units as u
import astropy.constants as c
import pickle
from   astroquery.vo_conesearch import conesearch
import astroquery
from   astroquery.vizier import Vizier 
from matplotlib.patches import Ellipse
import os.path # For expanduser

#==============================================================================
# Define kernels and files
#==============================================================================

file_tm = '/Users/throop/git/NH_rings/kernels_nh_pluto_mu69.tm'       # C/A time 11:00:00 (??) - older 
file_tm = '/Users/throop/git/NH_rings/kernels_nh_pluto_mu69_tcm22.tm' # C/A time 07:00:00 - newer version

file_hd_pickle = '/Users/throop/git/NH_rings/cat_hd.pkl'

#==============================================================================
# Initialize settings, define constants, start SPICE, etc.
#==============================================================================

sp.furnsh(file_tm) # Start up SPICE

hour      = 3600
day       = 24 * hour
minute    = 60

name_target = 'MU69'
name_observer = 'New Horizons'

DO_INVERT_SAUCER = True # This is John Spencer's terminology. He thinks my plot looks like a teacup on a saucer.
                        # During MT call, we decided to flip it around, so now the saucer will be inverted.

ang_per_pix = (5.7 / 5024) * hbt.d2r  # Angular resolution, radians, MVIC framing

width_fov_rad = 5.7 * hbt.d2r
height_fov_rad = 0.15 * hbt.d2r

dt_slew   = 30 # Time between slews
exptime   = 10  # Exposure time, MVIC.
exp_per_footprint = 3   # Number of total exposures per pointing. ('3' for 2+1.)
 
frac_overlap = 0.92   # Fractional overlap between footprints

fs        = 15   # General font size

utc_ca    = '2019 1 Jan 07:00:00'

n_footprints = 27

et_ca     = sp.utc2et(utc_ca)

et_start  = et_ca    + 30*minute
et_end    = et_start + 30*minute

radius_image = 2500*u.km
radius_kbo   = 20*u.km 
radius_ring  = 1000*u.km

radius_ring_km = radius_ring.to('km').value
radius_ring_km = np.array([150, 500, 1250])

#==============================================================================
# Set up the times for each exposure
#==============================================================================

n_exp = n_footprints * exp_per_footprint # Total number of exposures

index_exp = hbt.frange(0,n_exp-1)        
index_footprint = np.trunc(index_exp/3).astype('int')

# Define 'et' as an array, one per exposure, of the time that exposure is taken

et = et_start + (index_exp * exptime) + (index_footprint * dt_slew)

dist_kbo = np.zeros(n_exp)

for i,et_i in enumerate(et):
    (st, junk) = sp.spkezr(name_target, et_i, 'J2000', 'LT+S', name_observer)
    dist_kbo[i] = sp.vnorm(st[0:2])

dist_kbo = dist_kbo * u.km

width_fov_km = (dist_kbo * width_fov_rad).to('km').value    # Width, in km
height_fov_km = (dist_kbo * height_fov_rad).to('km').value  # Height, in km
      
#==============================================================================
# Make the plot
#==============================================================================
#%%
fig, ax = plt.subplots()

hbt.set_fontsize(20)
hbt.figsize((10,10))
lim_image_km = (-radius_image.to('km').value, radius_image.to('km').value)

ax.plot([0], [0])
ax.set_xlim(lim_image_km)
ax.set_ylim(lim_image_km)
ax.set_aspect('equal')

# Draw MU69

dx_kbo_deg = 1.   # Uncertainty in MU69 position, halfwidth, at K+30 min. From GV.
dx_kbo_km = dist_kbo[0].to('km').value * dx_kbo_deg*hbt.d2r

plt.errorbar(0, 0, xerr = dx_kbo_km, marker = 'o', linestyle = 'none', label = 'MU69', markersize=5, color='red')

# Draw the rings

# Draw the MU69 positional uncertainty, which is 1° in each direction, approx.

xy = (0,0)

for radius_ring_km_i in radius_ring_km:
    ell = matplotlib.patches.Ellipse(xy=xy, width=radius_ring_km_i*2, height=radius_ring_km_i*2, angle = 0, alpha=0.5, 
                                     facecolor = 'none', edgecolor = 'grey', linewidth=1)
    ax.add_patch(ell)

# Loop over and draw the FOVs!

y_last = -height_fov_km[0] * 6.5  # Place the starting Y-pos of the scan
index_footprint_last = 0

for i in range(n_exp):
    x = -width_fov_km[i]/2 
    y = y_last

    if (i == 60):              # After this many images, do a one-time slew to the bottom half of MU69
        if (DO_INVERT_SAUCER): 
            y = y - 1600           # Size of slew in 'inverted' case 
        else:  
            y = y - 2500           # Size of slew in 'non-inverted' case
      
    if (index_footprint[i] > index_footprint_last):
        if (DO_INVERT_SAUCER) and (i >= 60):
            y -= height_fov_km[i] * frac_overlap
        else:
            y += height_fov_km[i] * frac_overlap

    rect = matplotlib.patches.Rectangle(xy=(x,y), width=width_fov_km[i], height=height_fov_km[i], angle=0, 
                                        edgecolor='green', facecolor='none', 
                                        linewidth=0.2)    
    y_last = y
    index_footprint_last = index_footprint[i]

    ax.add_patch(rect)

    if (np.mod(i,3) == 1): 
        plt.text(x, y+25, index_footprint[i], fontsize=10)
        
# Make a legend

rect = matplotlib.patches.Rectangle(xy=(10000,10000), width=width_fov_km[i], height=height_fov_km[i], angle=0, 
                                        edgecolor='green', facecolor='none', 
                                        linewidth=2)          

ell = matplotlib.patches.Ellipse(xy=(10000,10000), width=radius_ring_km_i, height=radius_ring_km_i, angle = 0, alpha=0.5, 
                                     facecolor = 'none', edgecolor = 'grey', linewidth=2)
    
rect.set_label('MVIC FOV, numbered by footprint')
ax.add_patch(rect)

ell.set_label('MU69 rings at {}, {}, {} km'.format(radius_ring_km[0], radius_ring_km[1], radius_ring_km[2]))
ax.add_patch(ell)

ax.legend(loc = 'upper right')

t_start_rel_str = 'K + {:.0f} min'.format( (et_start - et_ca)/minute )

ax.set_xlabel('Projected Distance [km]')
ax.set_ylabel('Projected Y distance [km]')
ax.set_title('MU69 Outbound Rings, ' + t_start_rel_str )

# Write the image to disk

dir_out = os.path.expanduser('~') + '/git/NH_rings/out/'
file_out = 'mu69_rings_mosaic_mvic_plot' + ('_inverted' if DO_INVERT_SAUCER else '') + '.png'

plt.savefig(dir_out + file_out)
print("Wrote: " + dir_out + file_out)
    
plt.show()

#%%

