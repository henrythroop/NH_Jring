#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 13:26:44 2019

@author: throop
"""

# Make a radial profile from the backplanes directly, not from the unwrapped images
# The idea here is that this will be higher resolution and higher sensitivity, since
# all of the hi-res pixels haven't been dumped into one when unwrapping.
#
# This does not negate the work of the GUI program. I still use that to output the processed FITS files,
# and figure out the light scattering, etc.
#
# HBT 25-Jan-2018

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

from   scipy import signal, fftpack

import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy import units as u           # Units library
from   astropy import constants as c
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

from   astropy.convolution import Box1DKernel, Gaussian1DKernel, convolve
from   pymiecoated import Mie

from   scipy.stats import linregress

import re # Regexp
import pickle # For load/save

from   matplotlib.figure import Figure
from   scipy import interpolate

# HBT imports

import hbt

# Local NH rings imports

from  nh_jring_mask_from_objectlist            import nh_jring_mask_from_objectlist

from  nh_jring_mask_from_objectlist             import nh_jring_mask_from_objectlist
from  nh_jring_unwrap_ring_image                import nh_jring_unwrap_ring_image

from  nh_jring_extract_profile_from_unwrapped   import nh_jring_extract_profile_from_unwrapped

from  scatter_mie_ensemble                      import scatter_mie_ensemble
from  area_between_line_curve                   import area_between_line_curve

# =============================================================================
# Start main program.
# This is just a test, to eventually be made into a function for use elsewhere.
# =============================================================================

# We want to load the images of 5/1-6. Then load their backplanes. Then make a profile from each one.

file_pickle = 'nh_jring_read_params_571.pkl'     # Filename to read to get filenames, etc.
dir_out     = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
    
# self.dir_out = dir_out

# self.is_flattened = False

lun = open(dir_out + file_pickle, 'rb')
t = pickle.load(lun)                        # Self.t is the *entire* table for all J-ring obs, not subset.
lun.close()

# Process the group names. Some of this is duplicated logic -- depends on how we want to use it.

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

index_group = 5
index_images = hbt.frange(1,6)

# Extract the one group we are looking at 
# NB: This creates an entirely new table -- *not* a view into the original table. 
# If we modify values here, we need to explicitly write the changes back to the original.

groupmask = t['Desc'] == groups[index_group]
t_group = t[groupmask]  # 
    
num_images_group = np.size(t_group)

num_bins_radius  = 300
num_bins_azimuth = 300

bins_radius = hbt.frange(126_000, 131_000, num_bins_radius)
bins_azimuth = hbt.frange(0, 2*math.pi, num_bins_azimuth)

A_METIS    = 127979.8*u.km        # Orbital distance, in km. From SCW07
A_ADRASTEA = 128980.5*u.km        # Orbital distance, in km. From SCW07

frac_azimuth = 1  # What fraction of all the azimuth angles do we use? 0.5 → "Use central 50% of az angles", etc.
profiles = []
 
# Loop over all the images
 
for i in index_images:

    # Create the proper filenames
    
    file_orig = t_group[i]['Filename']
    file_processed = os.path.join(dir_out, os.path.basename(file_orig)).replace('_opnav.fit', '_processed.fit')
    file_backplanes = file_processed.replace('_processed.fit', '_opnav_planes.pkl')
    
    # Load the processed image
    
    print(f'Loading image {index_group}/{i}: {file_processed}')

    hdu = fits.open(file_processed)
    im_processed = hdu['PRIMARY'].data
    hdu.close
    
    # Load the backplanes
    
    lun = open(file_backplanes, 'rb')
    planes = pickle.load(lun)
    lun.close()
    plane_radius = planes['Radius_eq']
    plane_azimuth = planes['Longitude_eq']
  
    # Now extract the profile
    
    profile = 0. * bins_radius

    median_az = np.median(plane_azimuth)
    range_az = np.max(plane_azimuth) - np.min(plane_azimuth)
    is_good_azimuth = ((plane_azimuth > (median_az - (frac_azimuth/2)*range_az)) & \
                       (plane_azimuth < (median_az + (frac_azimuth/2)*range_az)))
    
    for j in range(num_bins_radius-1):
        is_good_radius = (plane_radius > bins_radius[j-1]) & (plane_radius < bins_radius[j])
        is_good = is_good_radius & is_good_azimuth
        val = hbt.nanmedian(im_processed[is_good])
        profile[j] = val
    
    # Save the profile
    
    profiles.append(profile)

# Make a plot of all the profiles
    
for i,index in enumerate(index_images):
    plt.plot(bins_radius/1000, profiles[i]+i)
    

# Add a plot of the merged profile

profile_merged = np.nanmedian(np.array(profiles),axis=0)

plt.plot(bins_radius/1000, profile_merged + i + 2, linewidth=3)
plt.title(f'Profiles from backplanes, {a}, frac={frac_azimuth}')
plt.xlabel('Radius [1000 km]')
plt.ylabel('DN')

plt.axvline(x=A_METIS.to('km').value/1000)
plt.axvline(x=A_ADRASTEA.to('km').value/1000)

plt.show()

#     
    # For each one: load the processed image, and load the backplane
    
    # Define a few constants within the method. We're going to keep using these in many places,
    # so best to set them up here. As per PEP-8, constants should be ALL_CAPS.
    
