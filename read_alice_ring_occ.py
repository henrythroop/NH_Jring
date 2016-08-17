# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:58:26 2016

READ_ALICE_RING_OCC
@author: throop
"""

# General python imports

import pdb
import glob
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.
from   subprocess import call
import warnings
import pdb
import os.path
import os
import subprocess

import astropy
from   astropy.io import fits
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import cspice
import skimage
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils import daofind
import wcsaxes
import time
from scipy.interpolate import griddata


import imreg_dft as ird
import re # Regexp
import pickle # For load/save

# Imports for Tk

import Tkinter
import ttk
import tkMessageBox
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

# HBT imports

import hbt

##########
# Select which of the sequences we want to read in
##########
# (NB: The spellling / capitalization is inconsistent in SAPNAME vs. VISITNAM. I have standardized it here.)

#sequence 	= 'O_RING_OC3'
sequence 	= 'O_RING_OC2'

fs           	= 15		# Font size

dir_images = '/Users/throop/Data/NH_Alice_Ring/' + sequence + '/data/pluto/level2/ali/all'
file_tm = '/Users/throop/gv/dev/gv_kernels_new_horizons.txt'

cspice.furnsh(file_tm)

#########
# Read the Alice data
#########

file_list = glob.glob(dir_images + '/*fit')

#file_list = file_list[0:10]
met_all = []  # A long array with a list of all of the timestamps
count_rate_all = []

for file in file_list:
    data = hbt.read_alice(file)
#    print "."
    startmet = data['header_spect']['STARTMET']
    spcutcal = data['header_spect']['SPCUTCAL']	# Mid-obs time

    visitnam = data['header_spect']['VISITNAM']
    sapname   = data['header_spect']['SAPNAME']

    dt       = data['header_count_rate']['SAMPLINT']
    count_rate_i = data['count_rate']              # Extract count rate from file

    met_i = startmet + dt * hbt.frange(0, np.shape(count_rate_i)[0]-1)   # Create the MET by interpolation

    print "Read " + os.path.basename(file) + ", MET " + repr(startmet) + ' = ' + hbt.met2utc(startmet) + \
      ', N = ' + repr(len(count_rate_i)) + ' samples' + \
      ', duration = ' + hbt.trunc(dt * len(count_rate_i),3) + ' sec'
#       print "  " + visitnam + " " + spcutcal

    met_all.append(met_i)
    count_rate_all.append(count_rate_i)

count_rate  = np.array([item for sublist in count_rate_all for item in sublist])  # Flatten the count rate (from 2D, to 1D)
count_rate  = np.array(count_rate, dtype=float)					  # Convert to float. Otherwise get wraparound.
met         = np.array([item for sublist in met_all for item in sublist])         # Flatten the MET array  (from 2D, to 1D)

# Compute UTC and ET for the initial timestep

utc_start= hbt.met2utc(np.min(met))
et_start = cspice.utc2et(utc_start)

# Compute MET and ET for all timesteps. Both ET and MET are in seconds, but their offset is different.

et       = met - met[0] + et_start    # ET is now fully populated and correct
t        = et - et[0]                 # Seconds since start of observation

num_dt   = np.size(et)

# Computed a smoothed count rate. _s indicates smoothed. Array is cropped at edges too.

binning      = 3000		# Smoothing. 25000 is too much (shows opposite trend!). 5000 and 1000 look roughly similar.
                            # To be meaningful, the binning timescale must be less than the deadband timescale (~20-30 sec RT).
                            # At binning=3000, it's 12 sec... so that is the longest we'd really want to go.
                            
count_rate_s = hbt.smooth_boxcar(count_rate, binning)[binning:-binning]

# Compute truncated versions of the time arrays, just in case they are useful

t_s          = t[binning:-binning]
et_s         = et[binning:-binning]
met_s        = met[binning:-binning]

##########
# Compute the angle from the star to Alice boresight, for every timestep.
##########

# NB: For FSS, I got the FSS-Sun angle directly from Gabe -- I didn't get it from SPICE.

# Get vector to star. 67 Ori = HR 2159 = HD 41753, a V=4.2 B3V. RA=91.89, Dec=14.77.

ra_star  = 91.89 * hbt.d2r
dec_star = 14.77 * hbt.d2r

ra       = np.zeros(num_dt)
dec      = np.zeros(num_dt)

vec_star_j2k = cspice.radrec(1., ra_star, dec_star)

name_fov = 'NH_ALICE_AIRGLOW'

vec_bsight_alice = (-1, 0, 0)  # -X defines Alice SOC FOV

vsep = np.zeros(np.size(et))

for i,et_i in enumerate(et):
  mx = cspice.pxform(name_fov, 'J2000', et[i])
  vec_alice_j2k = cspice.mxvg(mx, vec_bsight_alice)		# ICY did not have mxvg(), but pdstools does.
  vsep[i] = cspice.vsep(vec_star_j2k, vec_alice_j2k)   # Angular separation, in radians

  (junk, ra[i], dec[i]) = cspice.recrad(vec_alice_j2k)

#plt.plot(t, count_rate, marker = '.', linestyle='none')

##########
# Calculate statistics
##########

# Sigma_s = S / SNR = (mean) / (sqrt(n * mean) / (n * mean))
# This is the expected stdev for binned data (i.e., 63% of data should be within this amount of the mean, etc.)

sigma_s = np.mean(count_rate) * np.sqrt(np.mean(count_rate * binning)) / (np.mean(count_rate * binning))

##########
# Make plots
##########

plt.rcParams['figure.figsize'] = 15,5

# Plot of count rate vs. time

plt.plot(t_s, count_rate_s, marker = '.', linestyle='none', ms=0.1)
plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + ' = ' + repr(dt * binning) + ' sec' + \
                     ', 1$\sigma$ = ' + hbt.trunc(sigma_s,4), fontsize=fs)

plt.errorbar(100, np.mean(count_rate_s) - 4 * sigma_s, xerr=binning*dt/2, yerr=None, label='Binning Width') 
			# X 'error bar' -- show the bin width
plt.errorbar(300, np.mean(count_rate_s) - 4 * sigma_s, xerr=None, yerr=sigma_s/2, label='1$\sigma$ shot noise') 
			# Y 'error bar' -- show the binned shot noise error
plt.legend()
plt.ylabel('Count Rate', fontsize=fs)
plt.xlabel('Time since ' + utc_start + ' [sec]', fontsize=fs)
plt.show()

# Plot of pointing vs. time

#   plt.plot(et - et[0], (vsep % 0.02) * hbt.r2d - 0.65 - ((et-et[0] - 1900) / 1000))
plt.plot(et - et[0], (vsep % 0.02) * hbt.r2d)
#   plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + ' = ' + \
#  repr(dt * binning) + ' sec', fontsize=fs)
plt.xlabel('Seconds since ' + utc_start, fontsize=fs)
plt.ylabel('Degrees from star', fontsize=fs)
#plt.xlim(0,500)
plt.show()

##########
# Do another plot of pointing v. time
#########

#plt.plot(t_s, count_rate_s, marker = '.', linestyle='none', ms=0.1)
#plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + ' = ' + \
# repr(dt * binning) + ' sec', fontsize=fs)
#plt.ylabel('Count Rate', fontsize=fs)
#plt.xlim((0,210))
#plt.show()
#
#plt.plot(et - et[0], (vsep % 0.02) * hbt.r2d)
#plt.xlim((0,210))
#plt.ylim((0,0.04))
##   plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + \
#  ' = ' + repr(dt * binning) + ' sec', fontsize=fs)
#plt.xlabel('Seconds', fontsize=fs)
#plt.ylabel('Degrees from star', fontsize=fs)
#plt.show()

##########
# Make a plot of motion thru the deadband
##########

plt.rcParams['figure.figsize'] = 5,5
plt.plot(ra*hbt.r2d, dec*hbt.r2d, linestyle='none', marker='.', ms=1)
plt.title(sequence + ', start = ' + utc_start)
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
ax = plt.gca()
ax.ticklabel_format(useOffset=False)  # Turn off the 'offset' that matplotlib can use
plt.show()

### Testing for pointing

utc_test = '2015 jul 16 21:25:00'  # First timestep of OC3 movie
et_test  = cspice.utc2et(utc_test)

mx = cspice.pxform('NH_ALICE_AIRGLOW', 'J2000', et_test)
vec_alice_j2k = cspice.mxvg(mx, vec_bsight_alice)		# ICY did not have mxvg(), but pdstools does.
vsep_test = cspice.vsep(vec_star_j2k, vec_alice_j2k)   # Angular separation, in radians

(junk, ra_test, dec_test) = cspice.recrad(vec_alice_j2k)

## Now do some correlation between position and pointing, to see if there is anything I should unwrap there.
#
# Concl: there is quite a bit of trend across the detector. So, a lot of the variation I'm seeing
# in DN value is probably not due to statistics, but due to sensitivity across the detector.

# Now I need to make an XY scatter plot of all these data points. I really have no idea how to do this.

ra_s = ra[binning:-binning]    # Same as RA, but cropped, so the edges are missing. Has same # elements as count_rate_s (smoothed)
dec_s = dec[binning:-binning]

plt.rcParams['figure.figsize'] = 10,10
plt.plot(ra_s, count_rate_s, linestyle='none', marker='.', ms=0.1)
plt.xlabel('RA [deg]', fontsize=fs)
plt.title(sequence, fontsize=fs*1.5)
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.ylabel('DN (binned ' + repr(binning) + ')', fontsize=fs)
plt.show()

plt.plot(dec_s, count_rate_s, linestyle='none', marker='.', ms=0.1)

num_dx = 100
num_dy = num_dx

ra_arr  = np.linspace(np.min(ra),  np.max(ra),  num_dx)
dec_arr = np.linspace(np.min(dec), np.max(dec), num_dy)

count_rate_s_arr = griddata((ra_s, dec_s), count_rate_s, (ra_arr[None,:], dec_arr[:,None]), method='cubic')

plt.imshow(count_rate_s_arr, interpolation='none', vmin=16.2,vmax=16.45)
plt.imshow(count_rate_s_arr, interpolation='none', vmin=16.4,vmax=16.8)

plt.title(sequence + ', Mean DN vs position, binning=' + repr(binning), fontsize=fs)
plt.xlabel('RA bin', fontsize=fs)
plt.ylabel('Dec bin', fontsize=fs)
plt.colorbar()
plt.show()

# Same, but unbinned

#count_rate_arr = griddata((ra, dec), count_rate, (ra_arr[None,:], dec_arr[:,None]), method='cubic')
#
#plt.imshow(count_rate_arr, interpolation='none', vmin=15, vmax=17)
#plt.title(sequence + ', Mean DN vs position, binning=none', fontsize=fs)
#plt.xlabel('RA bin', fontsize=fs)
#plt.ylabel('Dec bin', fontsize=fs)
#plt.colorbar()
#plt.show()


# 


  