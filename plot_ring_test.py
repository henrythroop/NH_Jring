#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 14:56:53 2017

This is a demo routine that incorporates all of my ring image reduction functions.
  o Load an image
  o Read the backplanes
  o Read WCS
  o Plot ring
  o Extract unwrapped ring image
  
All of these are in modular functions.
  
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

from   astropy.utils import data
from   astropy.wcs import WCS

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

import cProfile # For profiling

from   matplotlib.figure import Figure

# HBT imports

import hbt

from nh_jring_mask_from_objectlist             import nh_jring_mask_from_objectlist
from nh_jring_unwrap_ring_image                import nh_jring_unwrap_ring_image
from nh_jring_extract_profiles_from_unwrapped  import nh_jring_extract_profiles_from_unwrapped

dir_image = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'
file_tm    = "/Users/throop/git/NH_rings/kernels_nh_jupiter.tm"  # SPICE metakernel

#file_image = 'lor_0034961819_0x630_sci_1_opnav.fit' # noise, can't see much

file_image = 'lor_0034962025_0x630_sci_1_opnav.fit' # Metis on ansa.

#file_image + 'lor_0034962461_0x630_sci_1_opnav.fit' # Ring, but not on ansa -- not a good test

#file_image = 'lor_0034616523_0x630_sci_1_opnav.fit' # Adrastea in middle right. dt =40 works. Good test.
#file_image = 'lor_0034618323_0x630_sci_1_opnav.fit' # Adrastea in on ansa
#file_image = 'lor_0034620123_0x630_sci_1_opnav.fit' # Adrastea in middle left. Fits dt=40 better than dt=0.
#                                                   # dt=0 has cross 5 pix down right. 

#file_image = 'lor_0035103963_0x633_sci_1_opnav.fit' # 4x4, gossamer, Amal in middle. Bad navigation - not useful.
#file_image = 'lor_0035117161_0x633_sci_1_opnav.fit' # 4x4, gossamer. Bad navigation - not useful.

#file_image = 'lor_0034715044_0x630_sci_1_opnav.fit' # 1x1, Metis on LHS. Fits dt = 60. Cross is 10 pix down from o.
#file_image = 'lor_0034627923_0x630_sci_1_opnav.fit' # 1x1, bright sat Amalthea passing by. Can't tell dt from it.

#file_image = 'lor_0035078888_0x633_sci_1_opnav.fit' # Not sure -- I think ring image off edge?

file_image = 'lor_0035079398_0x633_sci_1_opnav.fit'

# Initialize SPICE

sp.furnsh(file_tm)

# Set up the paths and filenames

dir_out    = '/Users/throop/Data/NH_Jring/out/'

dir_backplanes = dir_out

file_objects = file_image.replace('.fit', '_objects.txt')

path_objects = dir_out + file_objects
    
dir_images = '/Users/throop/Data/NH_Jring/data/jupiter/level2/lor/all/'

file_backplane = file_image.replace('.fit', '_planes.pkl')

file_short = file_image.replace('.fit', '').replace('_sci', '').replace('_opnav', '')[0:-8]

#==============================================================================
# Read image and all associated files
#==============================================================================

# Read the image and header

power_sfit = 5

im     = hbt.read_lorri(dir_images + file_image)
im     = hbt.lorri_destripe(im)
im     = hbt.remove_sfit(im, power_sfit)

header = hbt.get_image_header(dir_images + file_image)

# Load the WCS coords (so we can overlay the ring)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    w = WCS(dir_images + file_image)

# Process headaer

dx_pix = header['NAXIS1']
dy_pix = header['NAXIS2']
et     = header['SPCSCET']

# Load the mask

mask = nh_jring_mask_from_objectlist(path_objects,do_plot=True)    

# Load the backplane

lun = open(dir_backplanes + file_backplane, 'rb')
planes = pickle.load(lun)
lun.close()

im_masked = im.copy()
im_masked[mask == True] = np.median(im)

im_masked_t = im.copy()
im_masked_t[np.transpose(mask) == True] = np.median(im)

im_halfmasked = im.copy() + 10 * mask

#==============================================================================
# Find location of ring points
#==============================================================================

a_ring_outer_km = 129300
a_ring_inner_km = 122000

x_ring1, y_ring1 = hbt.get_pos_ring(et, name_body='Jupiter', radius=a_ring_inner_km, units='pixels', abcorr='LT', wcs=w)
x_ring2, y_ring2 = hbt.get_pos_ring(et, name_body='Jupiter', radius=a_ring_outer_km, units='pixels', abcorr='LT', wcs=w)

#==============================================================================
# Do the ring extraction
#==============================================================================

(numrad, rj_array) = sp.bodvrd('JUPITER', 'RADII', 3) # 71492 km
rj = rj_array[0]
r_ring_inner = 126000 # * rj   # Follow same limits as in Throop 2004 J-ring paper fig. 7
r_ring_outer = 132000 # rj

num_bins_azimuth = 300    # 500 is OK. 1000 is too many -- we get bins ~0 pixels
num_bins_radius  = 300

limits_radius = (r_ring_inner, r_ring_outer) # Distances in km

# Define an offset in X and Y for the ring. This is the residual navigation error, which we want to apply here.

dx = 0
dy = 0

# Do the unwrapping. If there are no ring points, then we catch the error.

try:
    (im_unwrapped, mask_unwrapped, bins_radius, bins_azimuth) = nh_jring_unwrap_ring_image(im, 
                                                                   num_bins_radius, limits_radius,
                                                                   num_bins_azimuth, 
                                                                   planes, dx=dx, dy=dy, mask=mask)
    IS_UNWRAPPED = True
    
except ValueError:
    print("No valid ring points found.")
    IS_UNWRAPPED = False
      
#==============================================================================
# Make a plot!
#==============================================================================

# Set up all plotting parameters

color_phot  = 'red'            # Color for stars found photometrically
color_cat   = 'lightgreen'     # Color for stars in catalog  
color_sat   = 'yellow'
marker_sat  = '+'
marker_ring = 'None'
color_ring  = 'blue'
alpha_ring  = 0.6

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile.     
    
DO_PLOT_RING_INNER = True
DO_PLOT_RING_OUTER = True

hbt.figsize((10,10))
plt.set_cmap('Greys_r')

fig, ax = plt.subplots(1, 3, figsize=(20, 10))  

# Plot the original images themselves

ax[0].imshow(stretch(im))
ax[1].imshow(stretch(im_halfmasked))
ax[2].imshow(stretch(np.roll(np.roll(im, dy, 0), dx, 1)))

ax[0].contour(planes['Radius_eq'])

# Set the plot parameters, limits, etc.

titles = ['Raw', 'Merged', 'Raw, Rolled dx={}, dy={}'.format(dx,dy)]

for i in [0,1,2]:  # Loop over all the images

    ax[i].set_xlim([0,dx_pix-1])  
    ax[i].set_ylim([dx_pix-1,0])
    ax[i].set_title(titles[i])

# Plot the rings

    if (DO_PLOT_RING_OUTER):
        ax[i].plot(x_ring2, y_ring2, marker=marker_ring, color = color_ring, ls = '--',
                      alpha = alpha_ring, label = 'Ring outer')
    
    if (DO_PLOT_RING_INNER):
        ax[i].plot(x_ring1, y_ring1, marker=marker_ring, color='green', ls = '-', \
            ms=8, alpha = alpha_ring, label='Ring inner')
    
    ax[i].legend(loc = 'upper left') 

plt.show()

# Plot the unwrapped images

if (IS_UNWRAPPED):
    
    extent = [bins_azimuth[0], bins_azimuth[-1], bins_radius[0], bins_radius[-1]]
    f      = (np.max(bins_radius) - np.min(bins_radius)) / (np.max(bins_azimuth) - np.min(bins_azimuth))
    aspect = 0.5 * 1/f  # Set the aspect ratio
    
    fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    
    ax[0].imshow(stretch(im_unwrapped),                       extent=extent, aspect=aspect, origin='lower')
    ax[1].imshow(stretch(im_unwrapped + 10 * mask_unwrapped), extent=extent, aspect=aspect, origin='lower')
    
    for i in [0,1]:
        ax[i].set_xlim([bins_azimuth[0], bins_azimuth[-1]])
        ax[i].set_ylim([bins_radius[0],  bins_radius[-1]])
        ax[i].set_title(titles[i])
    
        ax[i].set_xlabel('Azimuth [radians]')
        ax[i].set_ylabel('Radius [$R_J$]')
    
    plt.show()

#==============================================================================
# Extract radial and azimuthal profiles from the unwrapped images
#==============================================================================

    (profile_radius, profile_azimuth) \
                      = nh_jring_extract_profiles_from_unwrapped(im_unwrapped, bins_radius, bins_azimuth, 
                                              0.2, # radius_out -- ie, fraction of radius used for az profile
                                              0.5, # azimuth_out -- ie, fraction of az used for radial profile
                                              mask_unwrapped=mask_unwrapped)

    print('No valid ring points found')
    
# Plot the radial and azimuthal profiles

    fig, ax = plt.subplots(2, 1, figsize=(10, 10))
    ax[0].plot(bins_azimuth, profile_azimuth)
    ax[0].set_title('Azimuthal Profile, {}'.format(file_short))
    ax[0].set_xlabel('Azimuth [radians]')
    ax[0].set_ylabel('DN')
    
    ax[1].plot(bins_radius, profile_radius)
    ax[1].set_title('Radial Profile, {}'.format(file_short))
    ax[1].set_xlabel('Radius [km]')
    ax[1].set_ylabel('DN')
    plt.show()