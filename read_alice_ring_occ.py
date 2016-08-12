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

####

sequence 	= 'O_RING_Occ3'

#   sequence 	= 'O_Ring_Occ2'

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
    dt = data['header_count_rate']['SAMPLINT']
    count_rate_i = data['count_rate']              # Extract count rate from file

    met_i = startmet + dt * hbt.frange(0, np.shape(count_rate_i)[0]-1)   # Create the MET by interpolation

    print "Read file " + os.path.basename(file) + ", MET " + repr(startmet) + ' = ' + hbt.met2utc(startmet) + \
      ', N = ' + repr(len(count_rate_i)) + ' samples' + \
      ', duration = ' + repr(dt * len(count_rate_i)) + ' sec'
    
    met_all.append(met_i)
    count_rate_all.append(count_rate_i)

count_rate  = np.array([item for sublist in count_rate_all for item in sublist])  # Flatten the count rate (from 2D, to 1D)
met         = np.array([item for sublist in met_all for item in sublist])         # Flatten the MET array  (from 2D, to 1D)

# Compute ET for each timestep

et_start = cspice.utc2et(hbt.met2utc(startmet))
et       = met - met[0] + et_start    # ET is now fully populated and correct

t        = et - et[0]                 # Seconds since start of observation

# Computed a smoothed count rate. _s indicates smoothed. Array is cropped at edges too.

binning      = 25000		# Smoothing 
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

vec_star_j2k = cspice.radrec(1., ra_star, dec_star)

name_fov = 'NH_ALICE_AIRGLOW'

vec_bsight_alice = (-1, 0, 0)  # -X defines Alice SOC FOV

vsep = np.zeros(np.size(et))

for i,et_i in enumerate(et):
  mx = cspice.pxform(name_fov, 'J2000', et[i])
  vec_alice_j2k = cspice.mxvg(mx, vec_bsight_alice)		# ICY did not have mxvg(), but pdstools does.
  vsep[i] = cspice.vsep(vec_star_j2k, vec_alice_j2k)   # Angular separation, in radians

#plt.plot(t, count_rate, marker = '.', linestyle='none')

##########
# Make plots
##########

plt.rcParams['figure.figsize'] = 15,5

# Plot of count rate vs. time

plt.plot(t_s, count_rate_s, marker = '.', linestyle='none', ms=0.1)
plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + ' = ' + repr(dt * binning) + ' sec', fontsize=fs)
plt.ylabel('Count Rate', fontsize=fs)
plt.xlabel('t since start [sec]', fontsize=fs)
plt.show()

# Plot of pointing vs. time

#   plt.plot(et - et[0], (vsep % 0.02) * hbt.r2d - 0.65 - ((et-et[0] - 1900) / 1000))
plt.plot(et - et[0], (vsep % 0.02) * hbt.r2d)
#   plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + ' = ' + repr(dt * binning) + ' sec', fontsize=fs)
plt.xlabel('Seconds', fontsize=fs)
plt.ylabel('Degrees from star', fontsize=fs)
plt.show()

##########
# Do another plot of pointing v. time
#########

plt.plot(t_s, count_rate_s, marker = '.', linestyle='none', ms=0.1)
plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + ' = ' + repr(dt * binning) + ' sec', fontsize=fs)
plt.ylabel('Count Rate', fontsize=fs)
plt.xlim((0,210))
plt.show()

plt.plot(et - et[0], (vsep % 0.02) * hbt.r2d)
plt.xlim((0,210))
plt.ylim((0,0.04))
#   plt.title(sequence + ', dt = ' + repr(dt) + ' sec, binned x ' + repr(binning) + ' = ' + repr(dt * binning) + ' sec', fontsize=fs)
plt.xlabel('Seconds', fontsize=fs)
plt.ylabel('Degrees from star', fontsize=fs)
plt.show()

