#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 22:00:03 2017

@author: throop
"""

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

class image_stack:

# =============================================================================
# # This class stacks images. It takes a list of files, and it does various processing on them.
# =============================================================================
    
# =============================================================================
# Init method: load the index and all files
# Mandatory argument: the directory holding all the files
# =============================================================================
   
    def __init__(self, dir) :   
#        file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel
        file_tm = "kernels_kem.tm"  # SPICE metakernel
        
        self.dir = (dir + '/').replace('//', '/') # Make sure it terminates in exactly one /
                
        files = glob.glob(self.dir + '*/*.fits')     # Look in subdirs
        files.extend(glob.glob(self.dir + '*.fits')) # And look at root-level files
        
        num_files = len(files)

        self.file_save = self.dir + 'kem_hazard_picksar_n{}.pkl'.format(num_files)
        
# If a save file was found, then load it, and immediately return
        
        if (os.path.isfile(self.file_save)):
            self.load()
            return
        
#        files = files[0:10]
        
        stretch_percent = 95    
        stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
        
        plt.set_cmap('Greys_r')
        
        # Start up SPICE
        
        sp.furnsh(file_tm)
        
        mode     = []
        exptime  = []
        filename_short = []
        exptime  = []
        visitnam = []
        sapname  = []
        sapdesc  = []
        reqid    = []
        et       = []
        utc      = []
        target   = []
        
        # Set up a dictionary for the data itself
        
        self.data     = {}
        
        # Set up the table. "The standard method of representing Python 3 strings in numpy is via the unicode 'U' dtype."
        # Disadvantage of this is memory space, but not an issue for me.
        
        self.t = Table(  [[],              [],          [],         [],        [],       [],       [],      [],   [], 
                                                    [], [], [], [], [], [] ],
            names=('filename_short', 'exptime', 'visitname', 'sapname', 'sapdesc', 'target', 'reqid', 'et',  'utc', 
                                                    'x_pix', 'y_pix', 'ra', 'dec', 'angle', 'naxis1'),
            dtype = ('U50',           'float64', 'U50',      'U50',     'U50',     'U50',    'U50',   'float64', 'U50', 
                                                    'float64', 'float64', 'float64', 'float64', 'float64', 'float64'  ))
        
        
        for i,file in enumerate(files):
        
            # Read the data
        
            hdulist        = fits.open(file)
            arr            = hdulist[0].data
            err            = hdulist[1].data
            quality        = hdulist[2].data
            
            filename_short = file.split('/')[-1].replace('.fits', '').replace('lor_', '')\
                     .replace('_0x633_pwcs','')\
                     .replace('_0x630_pwcs','')
            exptime = hdulist[0].header['EXPTIME']
            visitnam = hdulist[0].header['VISITNAM']
            sapname = hdulist[0].header['SAPNAME']
            sapdesc = hdulist[0].header['SAPDESC']
            target = hdulist[0].header['TARGET']
            reqid = hdulist[0].header['REQID']
            et = hdulist[0].header['SPCSCET']
            utc = sp.et2utc(et, 'C', 1)
        
            print("Read {}/{} {}".format(i, len(files), filename_short))
            hdulist.close()
        
            # Read the WCS coords, and grab the location of MU69
            
            w = WCS(file)
            
            vec,lt = sp.spkezr('2014 MU69', et, 'J2000', 'LT', 'New Horizons')
            vec_sc_targ = vec[0:3]
            (junk,ra,dec) = sp.recrad(vec_sc_targ) # Get the RA / Dec of the object
            x_pix, y_pix    = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0) # Convert to pixels
            
        # Put the FITS parameters into an astropy table row
        
            self.t.add_row(
                      [filename_short, exptime, visitnam, sapname, sapdesc, target, reqid, et, utc, x_pix, y_pix,
                       hdulist[0].header['SPCBRRA'],   # RA (pre-WCS I assume, but not sure)
                       hdulist[0].header['SPCBRDEC'],  # Dec
                       hdulist[0].header['SPCEMEN'],  # EME North angle -- ie, boresight rotation angle
                       hdulist[0].header['NAXIS1']])   # NAXIS1 (and 2) -- either 256 (4x4) or 1024 (1X1)
        
        
        # Save the data into the dictionary. We would save it to the table, but that doesn't seem possible.
        # The key is the short filename
            
            self.data[filename_short] = arr
        
        # Sort by ET.
            
            self.t.sort('et')
        
    # Plot one image. Plot MU69 position on top of it to be sure that WCS is correct.
    
    # Finally, remove a few columns that we don't need, or that are wrong.
        
        self.t.remove_column('sapdesc')
        self.t.remove_column('sapname')
        self.t.remove_column('target')
        self.t.remove_column('visitname')

    # Initialize the 'indices' vector, which indicates which planes we use for flattening
    
        self.indices = np.ones(len(self.t), dtype=bool)
        
    # Return. Looks like an init method should not return anything.

# =============================================================================
# Load a pickle file from disk. 
# Running this is just faster than reading and processing the files explicitly -- it duplicates no functionality
# =============================================================================

    def load(self):
        lun = open(self.file_save, 'rb')
        (self.t, self.data, self.indices) = pickle.load(lun)
        print("Read: " + self.file_save)
        lun.close() 

# =============================================================================
# Save current state into a pickle file.
# =============================================================================

    def save(self):

        lun = open(self.file_save, 'wb')
        pickle.dump((self.t, self.data, self.indices), lun)
        lun.close()
        print("Wrote: " + self.file_save)     
        
# =============================================================================
# Print the table
# =============================================================================

    def print(self):
        self.t.pprint(max_lines = -1, max_width = -1)
        
# =============================================================================
# Set the indices to extract
# =============================================================================

    def set_indices(self, indices):
        self.indices = indices
        
# =============================================================================
# Flatten a stack of images as per the currently-set indices
# =============================================================================

    def flatten(self, method='mean'):
        
        self.num_planes = np.sum(self.indices)

        # Create an output array for all of the images to go into
        
        self.arr = np.zeros((self.num_planes, 1024, 1024))
        w = np.where(self.indices)[0]  # List of all the indices

# Get the mean offset in X and Y for this set of images
# Reference this relative to the mean for *all* the images. This assures that any stack of the 
# same imageset will be centered with any other stack of it.
        
#        self.x_pix_mean = np.mean(t['x_pix'][self.indices])  # x_pix itself is MU69 position, as per SPICE and WCS
#        self.y_pix_mean = np.mean(t['y_pix'][self.indices])  # This does seem to allow boolean indexing just fine

        self.x_pix_mean = np.mean(self.t['x_pix'])  # x_pix itself is MU69 position, as per SPICE and WCS
        self.y_pix_mean = np.mean(self.t['y_pix'])  # This does seem to allow boolean indexing just fine


        for j,w_i in enumerate(w):  # Loop over all the indices. 'j' counts from 0 to N. w_i is each individual index.
            im = self.data[self.t['filename_short'][w_i]]
            im_expand = scipy.ndimage.zoom(im,4)
            dx = self.t['x_pix'][w_i] - self.x_pix_mean
            dy = self.t['y_pix'][w_i] - self.y_pix_mean
        
            # Apply the proper roll in X and Y. What I am calling 'x' and 'y' 
            
#            if (do_subpixel):
#                im_expand = scipy.ndimage.shift(im_expand, (-dx*4, -dy*4)) # 
#            else:
            
            im_expand = np.roll(im_expand, int(round(-dy*4)), axis=0) # positive roll along axis=0 shifts downward
            im_expand = np.roll(im_expand, int(round(-dx*4)), axis=1) # positive roll along axis=1 shifts rightward
            
            self.arr[j,:,:] = im_expand                 # Stick it into place
        
        # Merge all the individual frames, using mean or median
        
        print("Flattening array with dimension {}".format(np.shape(self.arr)))
        
        if (method == 'mean'):
            self.arr_flat   = np.nanmean(self.arr,0)    # Fast
            
        if (method == 'median'):
            self.arr_flat = np.nanmedian(self.arr,0)  # Slow -- about 15x longer
            
        return self.arr_flat

#=============================================================================
# Extract a radial profile. This is the simplest possible case: circular
# ** This is NOT part of the method! It is a standalone function.        
#=============================================================================

def get_radial_profile_circular(arr, pos = (0,0), width=1, method='mean'):
  
    dx = hbt.sizex(arr) 
    dy = hbt.sizey(arr)
    
    xx, yy = np.mgrid[:dx, :dy]  # What is this syntax all about? That is weird.
                                 # A: mgrid is a generator. np.meshgrid is the normal function version.
    
    dist_pix_2d = np.sqrt( ((xx - pos[0]) ** 2) + ((yy - pos[1]) ** 2) )

    dist_pix_1d = hbt.frange(0, int(np.amax(dist_pix_2d)/width))*width
    
    profile_1d    = 0. * dist_pix_1d.copy()
    
#    profile_1d_median  = 0. * dist_pix_1d.copy()
    
    for i in range(len(dist_pix_1d)-2):

    # Identify the pixels which are at the right distance
    
        is_good = np.logical_and(dist_pix_2d >= dist_pix_1d[i],
                                 dist_pix_2d <= dist_pix_1d[i+1]) 
    
        if (method == 'mean'):
            profile_1d[i]   = np.mean(arr[is_good])
    #        profile_1d = profile_1d_mean
    
        if (method == 'median'):
            profile_1d[i] = np.median(arr[is_good])
    #        profile_1d = profile_1d_median
        
    return (dist_pix_1d, profile_1d)
