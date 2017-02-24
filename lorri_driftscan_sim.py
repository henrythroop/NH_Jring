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

