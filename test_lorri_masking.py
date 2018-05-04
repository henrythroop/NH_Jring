#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Nov 17 10:55:04 2017

TEST_LORRI_MASKING.PY

This program reads in one LORRI Jupiter-ring image. It then applies several different masks, and
extract a radial profile using each one.

The goals here are:
    o To compare the sensitivity between the amount of stray light, and the radial profile.
    o To investigate the importance of applying the sfit() *after* doing masking, vs. before.
    
Results:
    o In general, it is better to over-mask than under-mask. Here, I am throwing away 90% of the ring
      pixels, and that results in a much better profile. I have remove anything that is questionable.    

    o 

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
import spiceypy as py
from   astropy.wcs import WCS
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
#import wcsaxes
import time
from   scipy.interpolate import griddata
from   astropy.stats import sigma_clip

import re # Regexp
import pickle # For load/save

from skimage.io import imread, imsave

from nh_jring_extract_profile_from_unwrapped import nh_jring_extract_profile_from_unwrapped

from nh_jring_unwrap_ring_image import nh_jring_unwrap_ring_image

# HBT imports

import hbt

plt.set_cmap('Greys_r')

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'

dir_jring = '/Users/throop/data/NH_Jring/'
dir_out = dir_jring + 'out/'
dir_backplanes = dir_out

filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in

index_group = 7  # Jupiter ring phase curve
index_image = 32 # Image number within the group
        
# Load pickle file

#if os.path.exists(self.filename_save):

print("Loading file: {}".format(filename_save))

#    self.load(verbose=False) # This will load self.t
#    t = self.t
    
lun = open(dir_out + filename_save, 'rb')
t = pickle.load(lun)
lun.close()
        
# Extract the one group we are looking at 
# NB: This creates an entirely new table -- *not* a view into the original table. 
# If we modify values here, we need to explicitly write the changes back to the original.

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

groupmask = t['Desc'] == groups[index_group]
t_group = t[groupmask]

num_images_group = np.size(t_group)     

file = t_group[index_image]['Filename']

# Load the backplanes

dir_backplanes = '/Users/throop/data/NH_Jring/out/'
    
file_backplane = dir_backplanes + t_group['Shortname'][index_image].replace('.fit', '_planes.pkl')
				
lun = open(file_backplane, 'rb')
planes = pickle.load(lun)
lun.close()
       
# Read the image

hdulist = fits.open(file)
data = hdulist['PRIMARY'].data

# Clip the image, removing extrema

#data_clip = sigma_clip(data, sigma=3, iters=10)

stretch = astropy.visualization.PercentileInterval(90)  # PI(90) scales array to 5th .. 95th %ile. 

styles = ['', 'aggressive', 'nominal', 'ring_broad_only', 'ring_narrow_only']

# Set up plotting parameters

n_cols = 3
n_rows = len(styles)

imagenum = 1

hbt.figsize((20,30))

for style in styles:
    
    if (style):
        file_mask = dir_jring + '/masks/mask_7_32-35_{}.png'.format(style)
        mask = (imread(file_mask)[:,:,0])   # Flatten to 2D
        mask = (mask > 0)

    else:
        mask = np.ones(np.shape(data))
        mask = (mask > 0)

    # Remove the sfit from the masked image.
    # We have two different ways of doing this -- a flag chooses.
    # Concl: sfit() should be subtracted based on the masked data.
    
    DO_SFIT_TO_MASKED = False
    
    if (DO_SFIT_TO_MASKED):      
        out = hbt.remove_sfit(data, mask = mask)      # Remove sfit using just masked data
    else:
        out = hbt.remove_sfit(data)                   # Remove sfit using *all* data
        
    # Calculate the min and max values. But do this of the data, where it is masked True.
    # And, exclude outliers (like CRs) when we do this.
    # I am doing this carefully, and it looks good in the end. 
    # There is no need to do an a second stretch on top of this.
    
    out_masked = out[mask]
    out_masked_clipped = sigma_clip(out_masked, sigma_lower=1000, sigma_upper=3, iters=10)
    
    mm = hbt.mm(out_masked_clipped)

    stretch_mm = astropy.visualization.ManualInterval(vmin=mm[0], vmax=mm[1])
                
    # Set the defaults for unwrapping the images
    dx = dy = 0
    
    r_ring_inner = 114000   # Follow same limits as in Throop 2004 J-ring paper fig. 7
    r_ring_outer = 135000

    num_bins_azimuth = 300    # 500 is OK. 1000 is too many -- we get bins ~0 pixels
    num_bins_radius  = 500
    
    image_processed = out
    mask_stray = mask
    mask_objects = np.ones(np.shape(mask_stray), dtype=bool)
    range_of_azimuth = 1
    
    # Do the unwrapping
    
    (image_unwrapped,                                     # NB: This has a lot of NaNs in it.  
     mask_stray_unwrapped, mask_objects_unwrapped,        # The masks, unwrapped. True = good pixel
     radius_unwrapped, azimuth_unwrapped                  # The bins which define the unwrapped coordinates
    ) = \
     nh_jring_unwrap_ring_image(image_processed, 
                                 num_bins_radius, 
                                 (r_ring_inner, r_ring_outer),
                                 num_bins_azimuth, 
                                 planes, 
                                 dx=-dx, dy=-dy,
                                 mask_objects=mask_objects,
                                 mask_stray=mask_stray)

    # Take the radial profile
    
    profile_radius = nh_jring_extract_profile_from_unwrapped(
                                          image_unwrapped, 
                                          radius_unwrapped,   # Array defining bins of radius
                                          azimuth_unwrapped,  # Array defining bins of azimuth
                                          range_of_azimuth,        # ie, range of az used for rad profile
                                          'radius',
                                          mask_unwrapped = mask_stray_unwrapped)   

    # Now plot everything

    # Plot the image, using these properly calculated robust min and maxes, in the masked area only.
    
    ax = plt.subplot(n_rows, n_cols, imagenum);     imagenum += 1
    plt.imshow(stretch_mm(out*mask))
    plt.title(style)

    ax = plt.subplot(n_rows, n_cols, imagenum);     imagenum += 1
    plt.imshow(stretch_mm(image_unwrapped * mask_stray_unwrapped))   
    ax.set_aspect(3/5.)

    plt.title('Mask = {}'.format(file_mask.split('/')[-1]))
    
#    ax.set_aspect('equal')

    plt.subplot(n_rows, n_cols, imagenum);     imagenum += 1
    plt.plot(radius_unwrapped, profile_radius)
    plt.title(file.split('/')[-1])
    plt.xlim((119000, 135000))
    plt.ylim((-3,10))
    plt.xlabel('Radius [km]')
    plt.ylabel('DN')
    
plt.show()
    