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
import pickle
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs

    
file_tm = '/Users/throop/git/NH_rings/kernels_nh_pluto_mu69.tm'
file_hd_pickle = '/Users/throop/git/NH_rings/cat_hd.pkl'

#==============================================================================
# Initialize settings
#==============================================================================

sp.furnsh(file_tm) # Start up SPICE

hour      = 3600
day       = 24 * hour

utc_ca    = '2019 1 jan 03:09:00' # MU69 C/A time. I got this from GV.
et_ca     = sp.utc2et(utc_ca)
dt_tof    = 180                   # Time-of-flight uncertainty (halfwidth, seconds)

hbt.figsize((15,10))
mag_limit = 8 # Plot stars brighter than this.
plot_tick_every = 360  # Plot a time-tick every __ seconds

DO_PLOT_TOF = False               # In general, plotting TOF is not necessary -- it
                                  # is a translation in position along track,
                                  # and changes the occultation times but not positions.

DO_PLOT_LORRI = True              # Plot limits of LORRI FOV?
                            
# Define the start and end time for the plot. This is the main control to use.

case = 4.5

if (case == 1): # Outbound, 0.3 .. 96 hour
    et_start  = et_ca + 0.3*hour
    et_end    = et_ca + 96*hour
    pad_ra_deg  = 1 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 1
    DO_LABEL_HD = True # Plot star IDs on chart
    DO_PLOT_HD = True
    hbt.figsize((15,10))
    mag_limit = 8 # Plot stars brighter than this.

if (case == 2): # Inbound, -96h .. -0.1h
    et_start  = et_ca -96*hour
    et_end    = et_ca -0.1*hour
    pad_ra_deg  = 2 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 2
    DO_PLOT_HD  = True
    DO_LABEL_HD = True
    hbt.figsize((25,10))
    mag_limit = 8 # Plot stars brighter than this.

if (case == 3): # Outbound, 2 .. 200 hour
    et_start  = et_ca + 2*hour
    et_end    = et_ca + 200*hour
    pad_ra_deg  = 0.5 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 0.5
    DO_PLOT_HD = True
    DO_LABEL_HD = True # Plot star IDs on chart
    hbt.figsize((15,10))
    mag_limit = 10 # Plot stars brighter than this.

if (case == 4): # Outbound, 2d .. 10 *** This one is key!! It has one good occultation.
    et_start  = et_ca + 2*day
    et_end    = et_ca + 20*day
    pad_ra_deg  = 0.01 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 0.01
    DO_PLOT_HD  = False
    DO_LABEL_HD = DO_PLOT_HD # Plot star IDs on chart
    DO_PLOT_USNO   = True
    DO_LABEL_USNO  = DO_PLOT_USNO
    plot_tick_every = 7200*5
    hbt.figsize((15,10))
    mag_limit = 12 # Plot stars brighter than this. HD is probably complete to 10 or so.
    DO_PLOT_LORRI = False

if (case == 4.5): # Outbound, zoom on the good one!
    et_start  = et_ca + 3.0*day
    et_end    = et_ca + 20*day
    pad_ra_deg  = 0.005 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 0.005
    DO_PLOT_HD  = False
    DO_LABEL_HD = DO_PLOT_HD # Plot star IDs on chart
    DO_PLOT_USNO   = True
    DO_LABEL_USNO  = DO_PLOT_USNO
    plot_tick_every = 7200*5
    hbt.figsize((15,10))
    mag_limit = 12 # Plot stars brighter than this. HD is probably complete to 10 or so.
    DO_PLOT_LORRI = False
    
if (case == 5): # Inbound, 2d .. 10
    et_start  = et_ca - 20*day
    et_end    = et_ca - 2*day
    pad_ra_deg  = 0.01 # Add additional padding at edge of plots, RA = x dir. Degrees.
    pad_dec_deg = 0.01
    DO_PLOT_HD  = False
    DO_LABEL_HD = DO_PLOT_HD # Plot star IDs on chart
    DO_PLOT_USNO   = True
    DO_LABEL_USNO  = DO_PLOT_USNO
    plot_tick_every = 7200*5
    hbt.figsize((15,10))
    mag_limit = 12 # Plot stars brighter than this. HD is probably complete to 10 or so.
    DO_PLOT_LORRI = False

    
fov_lorri = 0.3*hbt.d2r  # Radians. LORRI width.

num_dt = 4000 # Number of timesteps

# Define plot colors

color_kbo_radius = 'pink'
color_kbo_center = 'red'
color_lorri      = 'lightgreen'
color_roche      = 'lightblue'
color_stars      = 'green'
color_tof        = 'blue'

hbt.set_fontsize(15)
    
try:
    lun = open(file_hd_pickle, 'rb')
    print("Loading file: " + file_hd_pickle)
    hd = pickle.load(lun)
    lun.close()
    
except IOError: # If Pickle file not found, then go ahead and read the catalog from scratch
    
    # Read HD into a table. Its positions are 1900. Coords are radians.

    print ("Reading HD catalog...")
    
    hd = read_hd_all()  
    num_stars = np.size(ra_1950)

    # Precess stars from 1900 to 2000
    
    print("Precessing to J2000")
    
    mx = sp.pxform('B1950', 'J2000', et[0]) # SPICE knows about 1950, but not 1900!
    
    ra_1950  = t['RA']
    dec_1950 = t['Dec']
        
    ra_2000  = []
    dec_2000 = []
    pt_1950  = []
    
    # Loop over every star and precess it individually. This loop is a real bottleneck.
    
    for i in range(num_stars):
        
        pt_1900 = sp.radrec(1, ra_1950[i], dec_1950[i])
        pt_1950 = sp.mxv(mx, pt_1900)
        pt_2000 = sp.mxv(mx, pt_1950) # SPICE will not precess 1900 -> 2000. So we apply 1950 -> 2000, do it 2x.
        d, ra_2000_i, dec_2000_i = sp.recrad(pt_2000)
        dec_2000.append(dec_2000_i)
        ra_2000.append(ra_2000_i)  
      
    hd['RA_2000']  = ra_2000
    hd['Dec_2000'] = dec_2000
    
    # Now save this pickle file, so we never have to run this routine again
    
    lun = open(file_hd_pickle, 'wb')
    pickle.dump(hd, lun)
    lun.close()
    print("Wrote: " + file_hd_pickle)

# Use the virtual observatory to get some UCAC2 stellar positions
# crval is a two-element array of [RA, Dec], in degrees
# Use conesearch.list_catalogs() to get list

encounter_phase = 'Inbound' if (et_start < et_ca) else 'Outbound'

# Set up the ET array

et = hbt.frange(et_start, et_end, num_dt)

index_asymptote = 0 if (encounter_phase == 'Outbound') else 1

# We make one call to SPICE to look up the asymptote position
# Most SPICE calls are later, but we need to do this one now, ahead of order.

(st,lt) = sp.spkezr('MU69', et[index_asymptote], 'J2000', 'LT', 'New Horizons')
(junk, ra_asymptote, dec_asymptote) = sp.recrad(st[0:3])   
crval = np.array([ra_asymptote, dec_asymptote]) * hbt.r2d
                  
radius_search = 0.1 # degrees # We only use these catalog for very fine searches, so narrow is OK.

DO_PLOT_GSC = False # Goes to about v=12

if DO_PLOT_GSC:                 
    name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
    stars = conesearch.conesearch(crval, radius_search, cache=False, catalog_db = name_cat)
    ra_stars  = np.array(stars.array['RAJ2000'])*hbt.d2r # Convert to radians
    dec_stars = np.array(stars.array['DEJ2000'])*hbt.d2r # Convert to radians
    mag_stars = np.array(stars.array['Pmag'])

if DO_PLOT_USNO:                 
    name_cat = u'The USNO-A2.0 Catalogue (Monet+ 1998) 1'
    stars = conesearch.conesearch(crval, radius_search, cache=False, catalog_db = name_cat)
    ra_stars  = np.array(stars.array['RAJ2000'])*hbt.d2r # Convert to radians
    id_stars = np.array(stars.array['USNO-A2.0']) # ID
    dec_stars = np.array(stars.array['DEJ2000'])*hbt.d2r # Convert to radians
    mag_b_stars = np.array(stars.array['Bmag'])
    mag_r_stars = np.array(stars.array['Rmag'])

    usno = Table([id_stars, ra_stars, dec_stars, mag_b_stars, mag_r_stars], 
                 names = ['ID', 'RA_2000', 'Dec_2000', 'Bmag', 'Rmag'])

#==============================================================================
# Look up NH position (ie, KBO position)
#==============================================================================

print("Looking up NH position...")


dt = et[1] - et[0] # Timestep

ra_kbo   = []
dec_kbo  = []
dist_kbo = []

ra_kbo_tof_minus = []
ra_kbo_tof_plus  = []
dec_kbo_tof_minus = []
dec_kbo_tof_plus  = []

for et_i in et:

    # Nominal position
    (st,lt) = sp.spkezr('MU69', et_i, 'J2000', 'LT', 'New Horizons')
    (dist_i, ra_i, dec_i) = sp.recrad(st[0:3])   
    ra_kbo.append(ra_i)    # MU69 RA/Dec, in radians
    dec_kbo.append(dec_i)
    dist_kbo.append(dist_i)

    # TOF uncertainty -- negative 
    (st,lt) = sp.spkezr('MU69', et_i-dt_tof, 'J2000', 'LT+S', 'New Horizons')
    (dist_i, ra_i, dec_i) = sp.recrad(st[0:3]) 
    ra_kbo_tof_minus.append(ra_i)    # MU69 RA/Dec, in radians
    dec_kbo_tof_minus.append(dec_i) 

    # TOF uncertainty -- positive 
    (st,lt) = sp.spkezr('MU69', et_i+dt_tof, 'J2000', 'LT+S', 'New Horizons')
    (dist_i, ra_i, dec_i) = sp.recrad(st[0:3]) 
    ra_kbo_tof_plus.append(ra_i)
    dec_kbo_tof_plus.append(dec_i) 
  
# Convert all of these from lists to NP arrays
    
dec_kbo           = np.array(dec_kbo)
ra_kbo            = np.array(ra_kbo)
dec_kbo_tof_plus  = np.array(dec_kbo_tof_plus)
dec_kbo_tof_minus = np.array(dec_kbo_tof_minus)
ra_kbo_tof_plus   = np.array(ra_kbo_tof_plus)
ra_kbo_tof_minus  = np.array(ra_kbo_tof_minus)

ra  = hd['RA_2000']*hbt.r2d           # Look up stellar coordinates - be sure to use J2000 
dec = hd['Dec_2000']*hbt.r2d
dist_kbo = np.array(dist_kbo)*u.km  # NH-MU69 distance

# Compute the size of the Hill sphere

rho_kbo = 2.5*u.gram/u.cm**3
r_kbo   = 16.5 * u.km
m_kbo   = 4/3 * math.pi * r_kbo**3 * rho_kbo
a_kbo   = 43*u.au

a_hill = (a_kbo * (m_kbo / c.M_sun/3)**(1/3)).to('km')

# Compute angular size of Hill sphere
# In general the Hill sphere is just too large to plot... LORRI covers it only a couple days after C/A

angle_hill = (a_hill / dist_kbo).to('').value # Convert to radians, and drop the units
angle_hill = np.clip(angle_hill, -math.pi, math.pi) # Clip it (same as IDL > or < operator)

# Compute angular size of KBO itself

angle_kbo  = (r_kbo / dist_kbo).to('').value # Radians

# Compute angular size of Roche limit

a_roche     = 2.5 * r_kbo
angle_roche = (a_roche / dist_kbo).to('').value # Radians

# Set the x and y limits for the plot

xlim = hbt.mm(ra_kbo*hbt.r2d)   # Plotting limits, in deg
ylim = hbt.mm(dec_kbo*hbt.r2d)

xlim = np.array(xlim) + pad_ra_deg * np.array([-1,1]) # Pad the plot by the specified amount
ylim = np.array(ylim) + pad_dec_deg * np.array([-1,1])



# Filter the stars. Set a flag for each one we want to plot, based on mag and position.

is_good =  np.logical_and(
                            hd['Ptg'] < mag_limit,
               
                            np.logical_and( 
                                np.logical_and(xlim[0] < ra, 
                                               ra < xlim[1]),     
                                np.logical_and(ylim[0] < dec,
                                               dec < ylim[1])  )    )

#==============================================================================
# Make the plot
#==============================================================================

# Draw the main trajectory

DO_PLOT_TRAJECTORY = True

if DO_PLOT_TRAJECTORY:
    plt.plot(ra_kbo*hbt.r2d, dec_kbo*hbt.r2d, color=color_kbo_center)

# Plot the time-of-flight uncertainty.
# Concl: TOF error does not change the position of any of these occultations. It just changes the time 
# at which we see them. This wasn't obvious to me, but these two curves lay over each other identically.
    
if DO_PLOT_TOF:
    plt.plot(ra_kbo_tof_minus*hbt.r2d, dec_kbo_tof_minus*hbt.r2d, color=color_tof, marker = '+', ms=5)
    plt.plot(ra_kbo_tof_plus*hbt.r2d,  dec_kbo_tof_plus*hbt.r2d,  color=color_tof, marker = '+', ms=10)
    
# Draw the LORRI FOVs -- 0.3 x 0.3 deg square

if DO_PLOT_LORRI:
    dec_kbo_plus_lorri  = dec_kbo + (fov_lorri/2) * np.sqrt(2)  # All radians
    dec_kbo_minus_lorri = dec_kbo - (fov_lorri/2) * np.sqrt(2)
    
    
    plt.fill_between(ra_kbo*hbt.r2d, (dec_kbo_plus_lorri)*hbt.r2d, (dec_kbo_minus_lorri)*hbt.r2d, 
                     color=color_lorri, alpha=0.5, label = 'LORRI FOV')

# Draw the Hill sphere

plt.plot(ra_kbo*hbt.r2d, (dec_kbo - angle_hill)*hbt.r2d, linestyle = '--')
plt.plot(ra_kbo*hbt.r2d, (dec_kbo + angle_hill)*hbt.r2d, linestyle = '--')

# Draw the Roche radius

plt.fill_between(ra_kbo*hbt.r2d, (dec_kbo+angle_roche)*hbt.r2d, (dec_kbo-angle_roche)*hbt.r2d, 
                 color=color_roche, alpha=0.5, label = 'MU69 Roche radius')

# Draw the MU69 radius

plt.fill_between(ra_kbo*hbt.r2d, (dec_kbo+angle_kbo)*hbt.r2d, (dec_kbo-angle_kbo)*hbt.r2d, 
                 color=color_kbo_radius, alpha=1, label = 'MU69 radius')

# Draw '+' symbols every hour

for i,et_i in enumerate(et):
    
    if (np.mod(et_i - et[0], plot_tick_every)) < dt:
        if (i == 0): # Pass the info for this to plt.label()... but only for the first call!
            kwargs = {'label' : 'MU69 position'}
        else:
            kwargs = {}
        plt.plot(ra_kbo[i]*hbt.r2d, dec_kbo[i]*hbt.r2d, marker='+', 
                 markersize=20, color='black', **kwargs )

# Plot the HD stars
if DO_PLOT_HD:
    
    plt.plot(hd['RA_2000'][is_good]*hbt.r2d, hd['Dec_2000'][is_good]*hbt.r2d, linestyle='none', marker='.', 
         label = 'HD, V < {}'.format(mag_limit), color=color_stars)

# Label the HD stars

    if (DO_LABEL_HD): 

        for i in range(np.size(is_good)):
            if (is_good[i]):
                    string_hd = hd['ID'][i]
            else:
                string_hd = ''
            
            plt.text(hd['RA_2000'][i]*hbt.r2d, hd['Dec_2000'][i]*hbt.r2d, 
                '  {}  {:.1f}  {}'.format(hd['Type'][i], hd['Ptg'][i], string_hd),
                fontsize = 8, clip_on = True)

# Plot the USNO stars
#plt.plot(usno['RA_2000'][is_good]*hbt.r2d, hd['Dec_2000'][is_good]*hbt.r2d, linestyle='none', marker='.', 
#         label = 'HD star', color=color_stars)

if DO_PLOT_USNO:
    plt.plot(usno['RA_2000']*hbt.r2d, usno['Dec_2000']*hbt.r2d, linestyle='none', marker='.', 
         label = 'USNO', color=color_stars)

if DO_LABEL_USNO:
    for i in range(np.size(usno)):
        plt.text(usno['RA_2000'][i]*hbt.r2d, usno['Dec_2000'][i]*hbt.r2d, 
            '  B={:.1f}, R={:.1f}  {}'.format(usno['Bmag'][i], usno['Rmag'][i], usno['ID'][i].decode("utf-8")),
            fontsize = 8, clip_on = True)
           
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
    
plt.title('MU69 {}, ticks every {:.0f} min, K{:+}h .. K{:+}h'.format(
          encounter_phase, plot_tick_every/60, (et_start - et_ca)/hour, (et_end-et_ca)/hour))
plt.legend(framealpha=0.8, loc = 'lower right')
plt.show()

#==============================================================================
# Do some test cases
#==============================================================================
