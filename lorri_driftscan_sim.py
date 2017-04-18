#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:44:27 2017

@author: throop
"""

# lorri_driftscan_sim -- simulate a LORRI driftscan of a given image
# This is for HBT's NH MU69 MT KBO/MT1.3-G4-4.2
# https://www.spaceops.swri.org/nh/wiki/index.php/KBO/MT1.3-G4-4.2
# HBT 24-Jan-2017

# General python imports

import pdb
import glob
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
import os.path
import os
import subprocess

import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
    

# HBT imports

import hbt

dir  = '/Users/throop/data/NH_misc/'
file = dir + 'lor_0235025618_0x633_sci_1.fit' # NGC3532 - open cluster - not a good choice - 0.4 sec 4X4
#file = dir + 'lor_0235132638_0x633_sci_1.fit'
#file = dir + 'lor_0268266239_0x630_sci_1.fit' # 0.1 sec opnav 1X1
#
ang_deadband_deg = 0.03
hbt.figsize((8,8))   # Image size
stretch_percent = 99
ang_rate_nominal = 30e-6        # Angular drift rate, radians/sec
length_sec = 29.9  # output exposure time
rate_drift = 30e-6 # rad/sec. For 30 urad/sec, that is 5 1x1 pixels/sec, or 150 pixels in a 30-sec exposure.
                   # That's pretty good -- i.e., most images will have zero thruster hits in them.
num_steps = 500 # Add this many frames to create output frame


hdulist = fits.open(file)
im = hdulist['PRIMARY'].data

stretch = astropy.visualization.PercentileInterval(stretch_percent)

plt.set_cmap('Greys_r')
header = hdulist['PRIMARY'].header
exptime = header['EXPTIME']
plt.imshow(stretch(im))

IS_4X4 = (header['SFORMAT'] == '4X4')  # Is the image we loaded in 4X4 format?

# Get the angular pixel size
         
ang_fov = 0.3*u.degree.to('radian') # radians in LORRI FOV

if IS_4X4:
    ang_pix = ang_fov / 256     # Radians per pixe, 4X4
else:
    ang_pix = ang_fov / 1024    # Radians per pixel, 1X1    

ang_deadband = ang_deadband_deg * hbt.d2r  # Deadband size


drift_per_sec_pix = rate_drift / ang_pix

num_sum = num_steps   # corect?? XXX I had to add this line; not sure if it is right or not.

dt = length_sec / num_sum  # Timestep

t = hbt.frange(0, length_sec, num_sum)

# Need to make a function which 

ang_x = np.zeros(num_steps)     # Angle, in radians, at each timestep
ang_y = np.zeros(num_steps)
ang_rate_x = np.zeros(num_steps) # Drift rate, in radians/second, at each timestep
ang_rate_y = np.zeros(num_steps)

# Set the position and drift rate at initial timestep

ang_x[0] = 0
ang_y[0] = 0
ang_rate_x[:] = ang_rate_nominal
ang_rate_y[:] = ang_rate_nominal/2

for i in range(num_steps):
    ang_x[i] = ang_x[i-1] + ang_rate_x[i]*dt
    ang_y[i] = ang_y[i-1] + ang_rate_y[i]*dt
    
    if (abs(ang_x[i]) > ang_deadband):
        ang_rate_x[i:-1] *= -1
        
    if (abs(ang_y[i]) > ang_deadband):
        ang_rate_y[i:-1] *= -1          

im_out = im.copy()
        
for i in range(num_steps):
    dx_pix =  int(ang_x[i] / ang_pix)
    dy_pix =  int(ang_y[i] / ang_pix)
    im_i = np.roll(im,   dx_pix, 0)  # Shift it
    im_i = np.roll(im_i, dy_pix, 1)  # Sum it
    im_out = im_i + im_out
#    print ("Rolled {}, {}".format(dx_pix, dy_pix))

# Convert to 4x4 if it is not.

#if not(IS_4X4):
#    im_out = im_out[::4,::4]  # Downsample it to 4x4

                    
ax = plt.imshow((im_out))
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)

plt.title('LORRI driftscan {} sec, {} urad/sec, deadband {} deg'.format(length_sec, rate_drift*1e6, ang_deadband_deg))

file_out = dir + 'out/lorri_driftscan_t{}_rate{}_db{}.png'.format(length_sec, rate_drift*1e6, ang_deadband_deg)
plt.savefig(file_out)
print('Wrote: ' + file_out)
plt.show()

#==============================================================================
# Now, do an unrelated calculation to determine the SNR that we'd get for a driftscan observation
# if it happened to go across a ring.
#==============================================================================

#%%
file_tm_dayside = '/Users/throop/git/NH_rings/kernels_nh_pluto_mu69_tcm22.tm' # Sort of hacked this together

utc_ca = '2019 1 Jan 07:00:00'

sp.furnsh(file_tm_dayside)

#Get the typical 

et_ca = sp.utc2et(utc_ca)

hour = 3600
minute = 60

dt = -12 * hour   # What is the offset of our observation, from KBO C/A. K-1h, K+2d, etc.

ddt = 1*minute   # Time offset for my calculation of velocity

width_pix_rad_4x4 = 4 * (0.3*hbt.d2r / 1024)  # LORRI 4x4 pixel size, in radians
width_pix_rad_1x1 =     (0.3*hbt.d2r / 1024)  # LORRI 4x4 pixel size, in radians

et_0 = et_ca + dt

(state_0, lt_0) = sp.spkezr('MU69', et_0,            'J2000', 'LT+S', 'New Horizons')
(state_1, lt_1) = sp.spkezr('MU69', et_0 + ddt, 'J2000', 'LT+S', 'New Horizons')

(junk, ra_0, dec_0) = sp.recrad(state_0[0:3])
(junk, ra_1, dec_1) = sp.recrad(state_1[0:3])

omega_kbo = sp.vsep(state_0[0:3], state_1[0:3]) / ddt  # Radians per second of sky motion that the KBO has, from NH

dist_kbo = sp.vnorm(state_0[0:3])  # Distance to KBO, in km

# Calculate the shadow velocity of the KBO

v_shadow_kbo = omega_kbo * dist_kbo  # km/sec of the shadow

# Calculate the time resolution of the LORRI driftscan.
# That is, how long does it take for LORRI to drift one 4x4 LORRI pixel?
# [A: Basically one sec: 0.681 sec, technically.]

dt_lorri_driftscan = width_pix_rad_4x4 / rate_drift

# Calculate Roche radius 

r_kbo       = 16.5 * u.km
r_roche     = 2.5 * r_kbo
r_roche_km  = r_roche.to('km').value

res_occ_kbo = v_shadow_kbo / dt_lorri_driftscan  # Resolution, in km

print("-----")
print("T =  K{:-.2f} h, dist =  {:.2f} km".format(dt/hour, dist_kbo))
print("KBO Shadow velocity = {:.2f} km/s".format(v_shadow_kbo))
print("LORRI drift rate (deadband) = {:.2e} rad/sec = {:.2f} 4x4 pix/sec.".format(
        rate_drift, rate_drift / width_pix_rad_4x4))
print("LORRI time resolution in a driftscan (ie, time to move one pixel) = {:.2f} sec.".format(dt_lorri_driftscan))
print("Time for KBO to move one LORRI 4x4 pixel = {:.2f} sec".format( 
        width_pix_rad_4x4 * dist_kbo / v_shadow_kbo))

# Compare speeds of KBO and of LORRI drift

print("LORRI is drifting {:.2f} x as fast as KBO".format(rate_drift / omega_kbo))

# Get our final answer, for the spatial resolution we can observe with driftscan

# XXX I don't think this one is correct below

print("Final spatial resolution of a driftscan observation at K{:-.2f} h = {:.2f} km".format(
        dt/hour, v_shadow_kbo / dt_lorri_driftscan))

print("At this same time, the 4x4 imaging resolution of LORRI is: {:.2f} km/pix".format(
        dist_kbo * width_pix_rad_4x4))
print("r_roche = {:.2f} km = {:.2f} LORRI 4x4 pix.".format(r_roche_km, 
      r_roche_km /dist_kbo / width_pix_rad_4x4))
      