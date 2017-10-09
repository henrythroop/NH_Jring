#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 22:00:03 2017

@author: throop
"""

# This program just does a bunch of misc testing on the various images taken in Sep-2017.
# Histograms, statistics, stacking, subtraction, navigation, etc. 
# My goal was to try to find anything interesting. 
# MU69 wasn't visible. These were images of the right field, so we could prep our exposures
# for the real hazard search.

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
from   scipy.stats import mode
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)
import imreg_dft as ird                    # Image translation

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt

from image_stack import image_stack    # I am not sure how to properly do this. Module name and filename are the same.
                                       # If I just do 'import image_stack', then I have to call 
                                       # as image_stack.image_stack().

from image_stack import get_radial_profile_circular

def nh_find_simulated_rings_lorri():

# =============================================================================
# Now go thru the synthetic ring images. 
#    Load and stack the synthetic implanted images.
#    Load and stack the original 'raw' frames
#    Difference them, and see if we can find a ring in there.
# =============================================================================

    dir_porter = '/Users/throop/Dropbox/Data/NH_KEM_Hazard/Porter_Sep17/'
    dir_synthetic = '/Users/throop/Dropbox/Data/NH_KEM_Hazard/synthetic/'
    
    do_subpixel = False  # Flag: Do we use sub-pixel shifting when doing the flattening? 
                         # It is slower and in theory better, but in reality makes a trivial difference.

    # Start up SPICE
    
    file_kernel = 'kernels_kem.tm'
    sp.furnsh(file_kernel)
    
    # Load the images into a table
    
    images_raw = image_stack(dir_porter)
    images_syn = image_stack(dir_synthetic, do_force=False)
    
    stretch = astropy.visualization.PercentileInterval(95)
    plt.set_cmap('Greys_r')

    # =============================================================================
    # If desired, do a one-time routine for the synthetic images:
    #  extract the I/F and ring size from the filenames, and append that to the table.
    # This routine should be run after creating new synthetic images (e.g., adding an I/F value) 
    # =============================================================================
    
    DO_APPEND = False
    if (DO_APPEND):

        t_syn = images_syn.t
        num_images_syn = (np.shape(t_syn))[0]

        iof_ring  = np.zeros(num_images_syn, dtype=float)
        size_ring = np.zeros(num_images_syn, dtype='U30')
    
        for i in range(num_images_syn):
            f = t_syn['filename_short'][i]
            m = re.search('ring_(.*)_iof(.*)_K', f)  # Call regexp to parse it.
            iof_ring[i] = eval(m.group(2))
            size_ring[i] = m.group(1)
            
        t_syn['size_ring'] = size_ring
        t_syn['iof_ring']  = iof_ring
        images_syn.t = t_syn
        images_syn.save()           # Save the whole pickle archive (including images and table) back to disk
    
    data_raw = images_raw.data
    data_syn = images_syn.data
    
    t_raw = images_raw.t
    t_syn = images_syn.t
    
    num_images_raw = (np.shape(t_raw))[0]
    num_images_syn = (np.shape(t_syn))[0]

    # Look up the time offset, from the image title. (Would be better to have it stored in table, but this will do.)

    match = re.search('_K(.*)d', t_syn['filename_short'][0])
    
    dt_ca = ((match.group(1)*u.day).to('s'))  # Weird: we don't need .value here. I can't explain it.

    utc_ca = '2019 1 Jan 05:33'
    et_ca  = sp.utc2et(utc_ca)
    et_obs = et_ca + dt_ca
    
    # Set the pixel scale
    
    vec,lt = sp.spkezr('2014 MU69', et_obs, 'J2000', 'LT', 'New Horizons')
    vec_sc_targ = vec[0:3]
    dist_target_km = (sp.vnorm(vec_sc_targ)*u.km).value    
    scale_pix_lorri_1x1_rad = 0.3*hbt.d2r / 1024
    scale_pix_lorri_4x4_rad = scale_pix_lorri_1x1_rad * 4
    scale_pix_km_dict = {'1X1' : scale_pix_lorri_1x1_rad * dist_target_km,
                         '4X4' : scale_pix_lorri_4x4_rad * dist_target_km}  # We are 
    
    # Create a bunch of possible image sets, based on various parameters
    
    # Indices for 'raw' images
    
    indices_sep17_raw = t_raw['et'] > sp.utc2et('15 sep 2017')  # The positon of MU69 has changed a few pixels.
                                                                # We can't blindly co-add between sep and pre-sep
    indices_jan17_raw = t_raw['et'] < sp.utc2et('1 sep 2017')
                                                        
    indices_rot0_raw  = t_raw['angle'] < 180   # One rotation angle
    indices_rot90_raw = t_raw['angle'] > 180   # The other rotation angle
    indices_10sec_raw = np.logical_and( t_raw['exptime'] < 10, t_raw['exptime'] > 5  )
    indices_20sec_raw = np.logical_and( t_raw['exptime'] < 20, t_raw['exptime'] > 10 )
    
    indices_30sec_raw = np.logical_and( t_raw['exptime'] < 30, t_raw['exptime'] > 20 )
    
    indices_1x1_raw = t_raw['naxis1'] == 1024
    indices_4x4_raw = t_raw['naxis1'] == 256
    
    indices_30sec_4x4_raw = np.logical_and(indices_4x4_raw, indices_30sec_raw) # 94
    
    # Indices for synthetic images
    
    indices_ring_small_syn = t_syn['size_ring'] == 'small'
    indices_ring_large_syn = t_syn['size_ring'] == 'large'
    
    indices_iof_1em7_syn = t_syn['iof_ring'] == 1e-7
    indices_iof_3em7_syn = t_syn['iof_ring'] == 3e-7
    indices_iof_1em6_syn = t_syn['iof_ring'] == 1e-6
    indices_iof_1em5_syn = t_syn['iof_ring'] == 1e-5
    indices_iof_1em4_syn = t_syn['iof_ring'] == 1e-4
    
    indices_small_1em7_syn = np.logical_and(indices_iof_1em7_syn, indices_ring_small_syn)
    indices_small_3em7_syn = np.logical_and(indices_iof_3em7_syn, indices_ring_small_syn)
    indices_small_1em6_syn = np.logical_and(indices_iof_1em6_syn, indices_ring_small_syn)
    indices_small_1em5_syn = np.logical_and(indices_iof_1em5_syn, indices_ring_small_syn)
    indices_small_1em4_syn = np.logical_and(indices_iof_1em4_syn, indices_ring_small_syn)
    indices_large_1em7_syn = np.logical_and(indices_iof_1em7_syn, indices_ring_large_syn)
    indices_large_3em7_syn = np.logical_and(indices_iof_3em7_syn, indices_ring_large_syn)
    indices_large_1em6_syn = np.logical_and(indices_iof_1em6_syn, indices_ring_large_syn)
    indices_large_1em5_syn = np.logical_and(indices_iof_1em5_syn, indices_ring_large_syn)
    indices_large_1em4_syn = np.logical_and(indices_iof_1em4_syn, indices_ring_large_syn)

    # Choose which indiex. ** THIS IS WHERE WE SET THE RING TO USE!!
    
    indices_raw = indices_30sec_4x4_raw.copy()   # 94 of 344
    indices_syn = indices_small_1em6_syn.copy()  # 94 of 752

    # Now take the first half of the synthetic indices, and the second half of the raw ones
    # This is to assure that we are using different images for the two stacks! Otherwise, the results are trivial.
    
    frames_max = int(np.sum(indices_raw) / 2)         # Total number of frames (94)
    
    w = np.where(indices_raw)[0]
    indices_raw[w[frames_max]:] = False          # De-activate all frames *below* frames_max
    
    w = np.where(indices_syn)[0]
    indices_syn[:w[frames_max]] = False          # De-activate all frames above frames_max
    
    # Set the indices
    
    images_raw.set_indices(indices_raw)
    images_syn.set_indices(indices_syn)
    
    # Do the flattening
        
    arr_raw = images_raw.flatten(do_subpixel=do_subpixel)
    arr_syn = images_syn.flatten(do_subpixel=do_subpixel)
    
#    arr_raw_sub = images_raw.flatten(do_subpixel=True)    
#    arr_syn_sub = images_syn.flatten(do_subpixel=True)
    
    # Extract various fields from the data table. We can look up from any of the images -- they should be all the same.
    
    t_syn       = images_syn.t  # Get the data table
    
    iof_ring    = t_syn[indices_syn]['iof_ring'][0]
    size_ring   = t_syn[indices_syn]['size_ring'][0]
    exptime     = t_syn[indices_syn]['exptime'][0]
    
    # The two flattened images need some offsetting. Do that.
    
    shift = ird.translation(arr_raw, arr_syn)['tvec']
#    shift = np.round(shift).astype('int')
    
#    arr_syn_shift = np.roll(np.roll(arr_syn, int(round(shift[0])), axis=0), int(round(shift[1])), axis=1)
    arr_syn_shift = scipy.ndimage.shift(arr_syn, shift, order=5)  # This allows sub-pixel shifts, apparently. *NO*!

#    a = arr_syn.copy()
#    a_05_05 = scipy.ndimage.shift(arr_syn, (0.5, 0.5), order=5)  # Ugh. 0.5, 0.5 and 1, 1 are *exactly* the same.
#    a_1_05 = scipy.ndimage.shift(arr_syn, (1, 0.5), order=5)
#    a_1_1 = scipy.ndimage.shift(arr_syn, (1, 1), order=5)
#    a_1_15 = scipy.ndimage.shift(arr_syn, (1, 1.5), order=5)
#    a_1_0 = scipy.ndimage.shift(arr_syn, (1, 0), order=5)
#    a_05_0 = scipy.ndimage.shift(arr_syn, (0.5, 0), order=5)
    
    arr_diff  = arr_syn_shift  - arr_raw
    
    pos = (images_raw.y_pix_mean*4, images_raw.x_pix_mean*4)
    
    # Set the binning width of the radial profiles

    binning_pix = 5
    
    # Extract the radial profiles
    
    (dist_pix_1d, profile_1d_median) = get_radial_profile_circular(arr_diff, pos, method='median', width=binning_pix)
    (dist_pix_1d, profile_1d_mean)   = get_radial_profile_circular(arr_diff, pos, method='mean', width=binning_pix)

    str_title = ('Synthetic ring - raw, I/F = {:.0e}, {}, {} x {:.1f}s'.format(
            iof_ring, size_ring, frames_max, exptime))
    
    plt.imshow(stretch(arr_diff))
    plt.title(str_title)
    plt.plot(pos[1], pos[0], marker='.', color='red')
    plt.show()
    
    # Set the scale for the effective mode of these observations. Many are taken as 4x4, but we've rebinned to 1x1
    
    if (np.shape(arr_raw)[0] == 1024):
        scale_mode = '1X1'
    else:
        scale_mode = '4X4'
    scale_pix_km = scale_pix_km_dict[scale_mode]
    
    # Make a plot of the radial profile. Don't plot the innermost bin. It is useless, since it has so few pixels in it.
    
    hbt.figsize((12,8))
    
    plt.plot(dist_pix_1d[1:] * scale_pix_km, profile_1d_median[1:], label = 'Annulus median', alpha = 0.7)
#    plt.plot(dist_pix_1d[1:] * scale_pix_km, profile_1d_mean[1:],   label = 'Mean',   alpha = 0.2)
    plt.xlabel('Distance [km]')
    plt.ylabel('DN per pixel')
    plt.title(str_title + ', binning = {}'.format(binning_pix))
    plt.xlim((0,30000))
    
    # Set the y axis range. This is really stupid. Can't matplotlib figure this out itself?
    
    ax = plt.gca()
    lims = ax.get_xlim()
    i = np.where( (dist_pix_1d * scale_pix_km > lims[0]) &  (dist_pix_1d*scale_pix_km < lims[1]) )[0]
    ax.set_ylim( profile_1d_median[i].min(), profile_1d_median[i].max() ) 
    
    plt.legend()
    plt.show()
    plt.savefig()
    
nh_find_simulated_rings_lorri()
    