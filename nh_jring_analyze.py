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

hbt.figsize((20,20))

dir_input = '/Users/throop/data/NH_Jring/out/'
files_input = glob.glob(dir_input + '*_export.pkl')
    
numfiles = np.size(files_input)
ang_phase = np.zeros((numfiles))
ang_elev  = np.zeros((numfiles))
ew        = np.zeros((numfiles))

profile_radius_dn = np.zeros((numfiles, 300))
profile_azimuth_dn = np.zeros((numfiles, 300))

dist_inner = 128000
dist_outer = 130000

for i,file in enumerate(files_input):
    lun = open(file, 'rb')
    vals = pickle.load(lun)
    lun.close

    (image, radius, azimuth, profile_radius_dn[i,:], profile_azimuth_dn[i,:], ang_elev[i], ang_phase[i]) = vals

    plt.plot(radius, profile_radius_dn[i,:] + i*3)
    
    bin_inner = hbt.x2bin(dist_inner, radius)
    bin_outer = hbt.x2bin(dist_outer, radius)
    
    ew[i] = np.sum(profile_radius_dn[i,bin_inner:bin_outer])

    print "File {}: bin {} .. {}: EW = {:.3f}".format(i, bin_inner, bin_outer, ew[i])

plt.vlines(dist_inner, -10, 30)
plt.vlines(dist_outer, -10, 30)

plt.show()

hbt.figsize((5,5))
plt.plot(ang_phase*hbt.r2d, ew)
plt.xlabel('Angle [deg]') 
plt.ylabel('EW [DN * km]')
plt.show()



# Now cross-correlate these signals

shift = np.zeros((numfiles)).astype(int)
profile_radius_dn_roll = profile_radius_dn.copy()

correl = np.zeros((numfiles, 41)) 
for i in range(numfiles):
    for j,dx in enumerate(hbt.frange(-20, 20).astype(int)):
        correl[i,j] = np.correlate(profile_radius_dn[0,:], np.roll(profile_radius_dn[i,:],dx))

    shift[i] = (hbt.wheremax(correl[i,:]))
    profile_radius_dn_roll[i,:] = np.roll(profile_radius_dn[i,:],shift[i])

hbt.figsize((20,20))

for i in range(numfiles):    
    plt.plot(radius, profile_radius_dn_roll[i,:])
plt.show()

hbt.figsize((5,5)) 
plt.plot(radius, np.sum(profile_radius_dn_roll,axis=0))
    