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

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt

class image_stack:

# =============================================================================
# Init method: load the index and all files
# =============================================================================

    def __init__(self) :   
        file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel
        
        dir = '/Users/throop/Dropbox/Data/NH_KEM_Hazard/'
        files = glob.glob(dir + '*/*.fits')
#        files = files[0:10]
        
        stretch_percent = 95    
        stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
        
        plt.set_cmap('Greys_r')
        
        # Start up SPICE
        
        #file_tm = 'kernels_nh_jupiter.tm'  # SPICE metakernel
        sp.furnsh(file_tm) # Commented out for testing while SPICE was recopying kernel pool.
        
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
# Stack some images as per the currently-set indices
# =============================================================================

    def stack(self, method='mean'):
        
        self.num_planes = np.sum(self.indices)

        self.arr = np.zeros((self.num_planes, 1024, 1024))
        w = np.where(self.indices)[0]  # List of all the indices

# Get the mean offset in X and Y for this set of images

        self.x_pix_mean = np.mean(t['x_pix'][self.indices])  # x_pix itself is MU69 position, as per SPICE and WCS
        self.y_pix_mean = np.mean(t['y_pix'][self.indices])  # This does seem to allow boolean indexing just fine

        for j,w_i in enumerate(w):  # Loop over all the indices. 'j' counts from 0 to N. w_i is each individual index.
            im = data[t['filename_short'][w_i]]
            im_expand = scipy.ndimage.zoom(im,4)
            dx = self.t['x_pix'][w_i] - self.x_pix_mean
            dy = self.t['y_pix'][w_i] - self.y_pix_mean
        
            # Apply the proper roll in X and Y. What I am calling 'x' and 'y' 
            
            im_expand = np.roll(im_expand, int(round(-dy*4)), axis=0) # positive roll along axis=0 shifts downward
            im_expand = np.roll(im_expand, int(round(-dx*4)), axis=1) # positive roll along axis=1 shifts rightward
            self.arr[j,:,:] = im_expand                 # Stick it into place
        
        # Merge all the individual frames, using mean or median
        
        if (method == 'mean'):
            self.arr_flat   = np.nanmean(arr,0)    # Fast
            
        if (method == 'median'):
            self.arr_flat = np.nanmedian(arr,0)  # Slow -- about 15x longer
            
        return self.arr_flat
        
#for i in range(120,150):
#
#        
#    plt.imshow(stretch(arr))
#    plt.title("{}, exptime {}".format(t['filename_short'][i], t['exptime'][i]))
#    plt.plot(x_pix, y_pix, marker = '+', color='red')
#    
#    plt.ylabel('Y pixels')
#    plt.xlabel('X pixels')
#    
#    plt.show()
#    
#    hdulist.close()

#=============================================================================
# OK! Now time to do some data analysis.
#=============================================================================


# =============================================================================
# Now do some tests!
# =============================================================================

stack = image_stack()

# Extract the data table

t = stack.t
data = stack.data  # Dictionary with *all* the images

plt.plot(t['x_pix'], t['y_pix'], marker = '+', linestyle='none')
plt.xlabel('MU69 X pos [pix]')
plt.ylabel('MU69 Y pos [pix]')
plt.show()

plt.plot(t['ra'], t['dec'], marker = '+', linestyle = 'none')
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
plt.show()


indices_sep17 = t['et'] > sp.utc2et('15 sep 2017')  # The positon of MU69 has changed a few pixels.
                                                    # We can't blindly co-add between sep and pre-sep
                                                    
indices_rot0  = t['angle'] < 180   # One rotation angle
indices_rot90 = t['angle'] > 180   # The other rotation angle
indices_10sec = np.logical_and( t['exptime'] < 10, t['exptime'] > 5  )
indices_20sec = np.logical_and( t['exptime'] < 20, t['exptime'] > 10 )
indices_30sec = np.logical_and( t['exptime'] < 30, t['exptime'] > 20 )

indices_1x1 = t['naxis1'] == 1024
indices_4x4 = t['naxis1'] == 256

indices_30sec_rot0  = np.logical_and(indices_30sec, indices_rot0)  # 94
indices_30sec_rot90 = np.logical_and(indices_30sec, indices_rot90) # 0

indices_20sec_rot0  = np.logical_and(indices_20sec, indices_rot0)  # 93
indices_20sec_rot90 = np.logical_and(indices_20sec, indices_rot90) # 0

indices_10sec_rot0  = np.logical_and(indices_10sec, indices_rot0)  # 96
indices_10sec_rot90 = np.logical_and(indices_10sec, indices_rot90) # 48

indices_30sec_rot0_sep  = np.logical_and(indices_30sec_rot0, indices_sep17)  # 94
indices_30sec_rot90_sep = np.logical_and(indices_30sec_rot90, indices_sep17) # 0

indices_20sec_rot0_sep  = np.logical_and(indices_20sec_rot0, indices_sep17)  # 93
indices_20sec_rot90_sep = np.logical_and(indices_20sec_rot90, indices_sep17) # 0

indices_10sec_rot0_sep  = np.logical_and(indices_10sec_rot0, indices_sep17)  # 48
indices_10sec_rot90_sep = np.logical_and(indices_10sec_rot90, indices_sep17) # 0




# =============================================================================
# Now set some indices, and stack them.
# =============================================================================

        
sequence = '10sec_rot0_sep'       # Which of these image sets defined above do we use?
indices  = eval('indices_' + sequence)   # Python trick: evaluate 

# Set the indices

stack.set_indices(indices)

w = np.where(indices)[0]

# Perform the stack, and get the result back

arr_flat_median = stack.stack(method = 'median')
arr_flat_mean   = stack.stack(method = 'mean')
arr             = stack.arr   # array, N x 1024 x 1024, with just the most recent stack

plt.imshow(stretch(arr_flat))
plt.title(sequence)
plt.show()

# Make a plot, showing difference between 

plt.plot(arr_flat_median[512,:], marker='.', ls='none', label = 'Summed Median')
plt.plot(arr_flat_mean[512,:],   marker='.', ls='none', label = 'Summed Mean')

plt.yscale('log')
plt.ylim((2,100))
#plt.plot(data[0,512,:], label = 'Individual', ls='none', marker='.', ms=0.5)
#plt.plot(data[10,512,:], label = 'Individual', ls='none', marker='.', ms=0.5)
#plt.plot(data[20,512,:], label = 'Individual', ls='none', marker='.', ms=0.5)

plt.legend()
plt.show()


# =============================================================================
# Now do a different stack ****
# =============================================================================

sequence = '30sec_rot0_sep'       # Which of these image sets defined above do we use?
indices  = eval('indices_' + sequence)   # Python trick: evaluate 

stack.set_indices(indices)

# Perform the stack, and get the result back

arr_flat = stack.stack()
plt.imshow(stretch(arr_flat))
plt.title(sequence)
plt.show()


#indices = indices_30sec_rot0_sep

# Create the output numpy array


    
# Set up a two-component stretch. This is because the second one (Sinh, or Sqrt) requires the data to 
# already be in range 0 .. 1.

stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
stretch2 = astropy.visualization.SqrtStretch()

plt.imshow(stretch2(stretch(arr_flat)))

plt.plot(x_pix_mean*4, y_pix_mean*4, marker = 'o', linestyle='none', color='red')
plt.title('{}, N = {}, {:0.2f} deg, {:0.2f} sec, {}'.format(
          sequence,
          num_planes, 
          t['angle'][w[0]]-180,
          t['exptime'][w[0]], 
          t['utc'][w[0]]) )
plt.show()

# Plot the original (ie, one frame, non-summed)

plt.imshow(stretch2(stretch(arr[0,:,:])))

plt.title('N = {}, {:0.2f} deg, {:0.2f} sec, {}'.format(
          sequence,  
          1, 
          t['angle'][w[0]]-180,
          t['exptime'][w[0]], 
          t['utc'][w[0]]) )
plt.show() 

# Make linear plots

row = 512

hbt.figsize((20,8))
plt.plot((arr_flat[row, :]), 
                 marker = '.', label = 'Summed, N = {}'.format(num_planes),linestyle='-', ms=1)

plt.plot((arr[0,row,:]), marker = '+', label = 'Individual frame')
plt.title(sequence)
plt.yscale('log')
plt.ylim((10,200))
plt.legend()
plt.show()

# Make a histogram plot

plt.plot()

bins = hbt.frange(0, 200) # Get a list of all of the timestep bins, inclusive, for this file.

(hist_flat, junk) = np.histogram(arr_flat, bins)

(hist_single, junk) = np.histogram(arr[0,:,:], bins)

hbt.set_fontsize(16)
plt.plot(bins[0:-1], hist_single, linestyle = '-', marker = '.', ms=1, label = 'Single')
plt.plot(bins[0:-1], hist_flat, linestyle = '-', marker = '.', ms=1, 
         label = 'Flattened, N={}'.format(num_planes))
plt.yscale('log')
plt.legend()
plt.title(sequence)
plt.ylabel('# of bins')
plt.xlabel('DN value per pixel')           
plt.show()

# Make a plot of the PSF size vs. exptime

# Take some     
# First get their indices

# Then expand them by 4x4

# Then shift them as per WCS

# Then co-add, and/or median them

# =============================================================================
# Calculate how many DN MU69 should be at encounter (K-20d, K-10d, etc.)
# Or alternatively, convert all of my DN values, to I/F values
# =============================================================================
