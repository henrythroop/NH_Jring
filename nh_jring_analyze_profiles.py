#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:51:59 2017

Program to analyze rings data

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
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling

import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

import re # Regexp
import pickle # For load/save

from   matplotlib.figure import Figure

# HBT imports

import hbt

# Local NH rings imports

from  nh_jring_mask_from_objectlist import nh_jring_mask_from_objectlist

from nh_jring_mask_from_objectlist             import nh_jring_mask_from_objectlist
from nh_jring_unwrap_ring_image                import nh_jring_unwrap_ring_image
from nh_jring_extract_profile_from_unwrapped   import nh_jring_extract_profile_from_unwrapped   

file_pickle = 'nh_jring_read_params_571.pkl' # Filename to read to get filenames, etc.
dir_out    = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
    
lun = open(dir_out + file_pickle, 'rb')
t = pickle.load(lun)
lun.close()

# Process the group names. Some of this is duplicated logic -- depends on how we want to use it.

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

# Get a list of all of the exported analysis files

files_analysis = glob.glob(dir_out + '*_analysis.pkl')

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

profile_radius_dn_arr = []
profile_azimuth_dn_arr = []
exptime_arr = []
ang_phase_arr = []
azimuth_arr = []

for file in files_analysis:
    lun = open(file, 'rb')
    vals = pickle.load(lun)
    lun.close()

    (image_unwrapped,     # Unwrapped image itself
                mask_unwrapped,      # Boolean mask
                radius,              # Axis values for the unwrapped image
                azimuth,             # Axis values for the unwrapped image 
                profile_radius_dn,   # Radial profile (several, in a dictionary)
                profile_azimuth_dn,  # Az profile (several, in a dictionary)
                range_of_azimuth,
                range_of_radius,
                exptime,             # Exposure time
                et,                  # ET
                ang_elev,            # Elevation angle above ring
                ang_phase,           # Phase angle (mean to rings -- not to planet center)
                bg_method,
                bg_argument,
                index_image,         # Index of image
                index_group) = vals  # Index of image group

    profile_radius_dn_arr.append(profile_radius_dn)
    profile_azimuth_dn_arr.append(profile_azimuth_dn)
    exptime_arr.append(exptime)
    ang_phase_arr.append(ang_phase)
    azimuth_arr.append(azimuth) 

#==============================================================================
# Now put these into arrays (not lists). Ideally we should put these into an astropy table (so we can sort, etc.)
#==============================================================================

ang_phase_arr = np.array(ang_phase_arr)
azimuth_arr   = np.array(azimuth_arr)
    
#==============================================================================
# Make a consolidated plot of radial profile
#==============================================================================

hbt.figsize((10,5))
dy = 1.5 # Vertical offset between lines

profile_radius_sum = profile_radius_dn_arr[0]['core']*0

for i,profile_radius in enumerate(profile_radius_dn_arr):
    plt.plot(radius, (i * dy*2) + profile_radius['core'], color='green', alpha=0.15)
    plt.text(132000, (i*dy*2 - 3), ang_phase_arr[i] * hbt.r2d)  # Label the phase angle here
    profile_radius_sum += profile_radius['core']
    
plt.xlabel('Radial Distance [km]')
plt.ylabel('DN')
plt.show()

#==============================================================================
# Make a consolidated plot of azimuthal profile
#==============================================================================

hbt.figsize((10,2))

for i,profile_azimuth in enumerate(profile_azimuth_dn_arr):
    plt.plot(azimuth_arr[i,:]*hbt.r2d, profile_azimuth['net'], color='green', alpha=0.15)
plt.xlabel('Azimuth [deg]')
plt.ylabel('DN')
plt.show()
       