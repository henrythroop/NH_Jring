#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:21:34 2016

@author: throop
"""

# Now that we have generated the radial profiles in nh_jring_gui.py, read and plot the results.

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
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
import wcsaxes
import time
from scipy.interpolate import griddata

import re # Regexp
import pickle # For load/save


# HBT imports

import hbt

hbt.figsize((12,8))

dir_input = '/Users/throop/data/NH_Jring/out/'
files_input = glob.glob(dir_input + '*_export.pkl')

num_radius         = 300
num_azimuth        = 300
    
numfiles           = np.size(files_input)
ang_phase          = np.zeros((numfiles))
ang_elev           = np.zeros((numfiles))
ew                 = np.zeros((numfiles))
profile_radius_dn  = np.zeros((numfiles, num_radius))
profile_azimuth_dn = np.zeros((numfiles, num_azimuth))
et                 = np.zeros((numfiles))
index_image        = np.zeros((numfiles))
index_group        = np.zeros((numfiles))

image              = np.zeros((numfiles, num_radius, num_azimuth))

dist_inner = 128000
dist_outer = 132000

# Read in all the files

for i,file in enumerate(files_input):
    lun = open(file, 'rb')
    vals = pickle.load(lun)
    lun.close

    (image[i,:,:], et[i], radius, azimuth, profile_radius_dn[i,:], profile_azimuth_dn[i,:], \
       ang_elev[i], ang_phase[i], 
       index_image[i], index_group[i]) = vals

    
# Cross-correlate these signals and shift them radially in order to align them using radial profiles

shift = np.zeros((numfiles)).astype(int)
profile_radius_dn_roll = profile_radius_dn.copy()

correl = np.zeros((numfiles, 41)) 
for i in range(numfiles):
    for j,dx in enumerate(hbt.frange(-20, 20).astype(int)):
        correl[i,j] = np.correlate(profile_radius_dn[0,:], np.roll(profile_radius_dn[i,:],dx))

    shift[i] = (hbt.wheremax(correl[i,:]))
    profile_radius_dn_roll[i,:] = np.roll(profile_radius_dn[i,:],shift[i])

# Shift each image vertically

bin_inner_vnorm = hbt.x2bin(131000, radius)
bin_outer_vnorm = hbt.x2bin(133000, radius)

for i in range(numfiles):
    profile_radius_dn_roll[i,:] -= np.mean(profile_radius_dn_roll[i,bin_inner_vnorm:bin_outer_vnorm])
    
# Make a plot of all the radial profiles, now aligned both vertically and horizontally

dy = 3

hbt.figsize((12,8))

for i in range(numfiles):    
    plt.plot(radius, profile_radius_dn_roll[i,:] + i * dy)

    plt.vlines(dist_inner, -10, 30)
plt.vlines(dist_outer, -10, 30)
plt.show()

# Calculate EW

bin_inner = hbt.x2bin(129000, radius)
bin_outer = hbt.x2bin(131000, radius)
    
for i in range(numfiles):
    ew[i] = np.sum(profile_radius_dn_roll[i,bin_inner:bin_outer])
    
# Make a plot of EW vs. phase

hbt.figsize((5,5))
plt.plot(ang_phase*hbt.r2d, ew, marker = 'o', linestyle='none')
plt.xlabel('Phase Angle [deg]') 
plt.ylabel('EW [DN * km]')
plt.show()



# Make a plot of all of the images, ganged

stretch = astropy.visualization.PercentileInterval(1)  # PI(90) scales array to 5th .. 95th %ile. 

hbt.figsize((12,8))

hbt.figsize((5,2*numfiles))      # x, y
for i in range(numfiles):
    plt.subplot(numfiles,1,i+1)  # y, x, n
    plt.imshow(image[i,:,:], aspect=0.3, vmin=-10,vmax=20)
plt.show()    

