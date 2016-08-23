# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Test for reading Alice files.



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

from mpl_toolkits.axes_grid1 import host_subplot # For adding a second axis to a plot
import mpl_toolkits.axisartist as AA             # For adding a second axis to a plot

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

sequence 	= 'O_RING_OC3'
#sequence 	= 'O_RING_OC2'

fs           	= 15		# Font size

dir_images = '/Users/throop/Data/NH_Alice_Ring/' + sequence + '/data/pluto/level2/ali/all'
file_tm = '/Users/throop/gv/dev/gv_kernels_new_horizons.txt'

cspice.furnsh(file_tm)

#==============================================================================
# Read the Alice data
#==============================================================================

file_list = glob.glob(dir_images + '/*fit')

met_all = []               # A long array with a list of all of the timestamps, one per 4 ms (i.e., at 250 hz)
count_rate_fits_all = []   # The count rate as read from the COUNT_RATE extension directly
count_rate_all = []        # Count rate computed from the PIXEL_LIST_TABLE. Should match that in COUNT_RATE extension
count_rate_target_all = [] # Count rate for the target only, extracted by spatially filtering the PIXEL_LIST_TABLE

d_target_summed = np.zeros((5,540))

for i,file in enumerate(file_list):
    
    hdulist = fits.open(file)
    d = hdulist['PRIMARY'].data # Units of this are float, but I'm not sure what they are. I would prefer raw counts.
    d_target = d[13:18, 370:911]
    d_target_summed += d_target
    p = hdulist['PIXEL_LIST_TABLE'].data
    count_rate_fits_i = hdulist['COUNT_RATE'].data
    num_samples = hdulist['COUNT_RATE'].header['NAXIS1'] # Number of samples in this file
    dt          = hdulist['COUNT_RATE'].header['SAMPLINT']  # Count rate sampling interval [sec]
    
    bins = hbt.frange(0, num_samples) # Get a list of all of the timestep bins, inclusive, for this file.
                                        # Use '+1' so we create the histogram upper size bin.
    
    # Now downselect the pixel list for just the photons in the proper X and Y position on the detector
    
    is_good = (p['Y_INDEX'] < 18) & (p['Y_INDEX'] >= 13) & (p['X_INDEX'] >= 370) & (p['X_INDEX'] < 910)

    # Now we have a list of all of the good pixels. For each of these, now we want to grab its timestep.

    timesteps_good = p['TIMESTEP'][is_good]
    timesteps_all  = p['TIMESTEP']

# Now count how many photons are in each timestep bin. I have defined those timestep bins up above.

    (count_rate_target_i, junk) = np.histogram(timesteps_good, bins)
    (count_rate_i, junk)        = np.histogram(timesteps_all,  bins) 
    met_i = hdulist['PRIMARY'].header['STARTMET'] + dt * np.array(range(num_samples))
          
#    print "File " + os.path.basename(file) + ' : ' + repr(np.sum(count_rate_target_i)) + ' / ' + repr(np.sum(count_rate_i))

# Now append these into the output lists

    count_rate_fits_all.append(count_rate_fits_i)
    count_rate_all.append(count_rate_i)
    count_rate_target_all.append(count_rate_target_i)
    met_all.append(met_i)
    
    hdulist.close()

count_rate_fits  = np.array([item for sublist in count_rate_fits_all for item in sublist])  # Flatten the count rate (from 2D, to 1D)
count_rate_fits  = np.array(count_rate_fits, dtype=float)					          # Convert to float. Otherwise get wraparound.

count_rate  = np.array([item for sublist in count_rate_all for item in sublist])
count_rate  = np.array(count_rate, dtype=float)					  

count_rate_target  = np.array([item for sublist in count_rate_target_all for item in sublist]) 
count_rate_target  = np.array(count_rate_target, dtype=float)					  

met         = np.array([item for sublist in met_all for item in sublist])         # Flatten the MET array  (from 2D, to 1D)
met         = np.array(met, dtype=float)

#==============================================================================
# Done reading Alice data
#==============================================================================

plt.plot(count_rate_target[0:100]) # drawstyle='steps'  <-- to make histogram style
plt.plot(count_rate[0:100])
plt.plot(count_rate_fits[0:100]+1)

plt.plot(met, count_rate_target, linestyle='none', ms=0.5, marker='.')

plt.plot(count_rate[0:100], drawstyle='steps')
plt.plot(count_rate[0:100]+1, drawstyle='steps')

count_rate_target_s = hbt.smooth_boxcar(count_rate_target, 10000)
count_rate_s        = hbt.smooth_boxcar(count_rate, 10000)

plt.rcParams['figure.figsize'] = 15,5

plt.plot(count_rate_target_s)
plt.ylim((13.5,13.8))
plt.plot(count_rate_s - count_rate_target_s + 11.1)
#plt.plot(count_rate_s)
plt.show()

# Compare to Steffl's OC3 plot

plt.plot(hbt.smooth_boxcar(count_rate_target, 2500)*250)
plt.ylim((3300,3500))



met         = np.array([item for sublist in met_all for item in sublist])         # Flatten the MET array  (from 2D, to 1D)
duration    = np.array(duration_all)

quit
    
plt.imshow(d_target_summed, aspect = 100, interpolation = 'none')
plt.show()

file = file_list[0]
hdulist = fits.open(file)


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
    duration_i = dt * len(count_rate_i)    

    
    print "Read " + os.path.basename(file) + ", MET " + repr(startmet) + ' = ' + hbt.met2utc(startmet) + \
      ', N = ' + repr(len(count_rate_i)) + ' samples' + \
      ', duration = ' + hbt.trunc(duration_i,3) + ' sec'
#       print "  " + visitnam + " " + spcutcal

    met_all.append(met_i)
    count_rate_all.append(count_rate_i)
    duration_all.append(duration_i)

quit