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
# Init method: load the index and all files
# =============================================================================

    def load(self):
        lun = open(self.file_save, 'rb')
        (self.t, self.data, self.indices) = pickle.load(lun)
        print("Read: " + self.file_save)
        lun.close() 

    def save(self):

        lun = open(self.file_save, 'wb')
        pickle.dump((self.t, self.data, self.indices), lun)
        lun.close()
        print("Wrote: " + self.file_save)        
        
    def __init__(self) :   
        file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel
        
        self.dir = '/Users/throop/Dropbox/Data/NH_KEM_Hazard/'
        
        files = glob.glob(self.dir + '*/*.fits')
        
        num_files = len(files)

        self.file_save = self.dir + 'kem_hazard_picsar_n{}.pkl'.format(num_files)
        
# If a save file was found, then load it, and immediately return
        
        if (os.path.isfile(self.file_save)):
            self.load()
            return
        
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

    # Initialize the 'indices' vector, which indicates which planes we use for flattening
    
        self.indices = np.ones(len(self.t), dtype=bool)
        
    # Return. Looks like an init method should not return anything.
        
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

        self.x_pix_mean = np.mean(t['x_pix'])  # x_pix itself is MU69 position, as per SPICE and WCS
        self.y_pix_mean = np.mean(t['y_pix'])  # This does seem to allow boolean indexing just fine


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
        
        print("Flattening array with dimension {}".format(np.shape(self.arr)))
        
        if (method == 'mean'):
            self.arr_flat   = np.nanmean(self.arr,0)    # Fast
            
        if (method == 'median'):
            self.arr_flat = np.nanmedian(self.arr,0)  # Slow -- about 15x longer
            
        return self.arr_flat

#=============================================================================
# OK! Now time to do some data analysis.
#=============================================================================

images = image_stack()

# Extract the data table

t = images.t
data = images.data  # Dictionary with *all* the images

DO_PLOT_POSITIONS = False

if (DO_PLOT_POSITIONS):
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
indices_jan17 = t['et'] < sp.utc2et('1 sep 2017')
                                                    
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

indices_10sec_rot0_jan = np.logical_and(indices_10sec_rot0, indices_jan17) # 48

# =============================================================================
# Now set some indices, and stack them.
# =============================================================================

        
sequence = '10sec_rot0_sep'       # Which of these image sets defined above do we use?
indices  = eval('indices_' + sequence)   # Python trick: evaluate 

# Set the indices

images.set_indices(indices)

w = np.where(indices)[0]

# Perform the stack, and get the result back

arr_flat_median = images.flatten(method = 'median')
arr_flat_mean   = images.flatten(method = 'mean')
arr             = images.arr   # array, N x 1024 x 1024, with just the most recent stack

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
# Now load three stacks, and compare them. For 10 sec, 20 sec, 30 sec
# =============================================================================

hist_single = []

hist_flattened = [] # Histogram from stacked frames
hist_single    = [] # Histogram from single frame 
n              = [] # Number of frames summed
arr_flattened  = []

sequences = ['30sec_rot0_sep', '20sec_rot0_sep', '10sec_rot0_sep']
dn_min = -50
dn_max = 200
bins = hbt.frange(dn_min, dn_max) # Get a list of all of the timestep bins, inclusive, for this file.
frames_max = 40  # Max frames in a stack. Some stacks have more or less -- we need to force to be the same, to compare.

colors = ['red', 'green', 'blue']

for sequence in sequences:
    
    indices  = eval('indices_' + sequence).copy()   # Python trick: evaluate a statement
    w = np.where(indices)[0]
    indices[w[frames_max]:] = False          # De-activate all frames above frames_max
    
# XXX need to limit this to N frames max, for consistency against all the stacks
    
    images.set_indices(indices)
    
# Perform the flattening, and get the result back

    arr_flat = images.flatten()
    
    arr = images.arr    # Get the raw stack, un-flattened. Only valid after calling flatten().
    
    (h, junk)  = np.histogram(arr_flat, bins)
    (h2, junk) = np.histogram(arr[0,:,:], bins)

    # Save the output from this stack calculation
    
    hist_flattened.append(h)   # Histogram of stacked images
    hist_single.append(h2)     # Histogram of 0th image in the stack
    n.append(np.sum(indices))  # Number of images in each stack
    arr_flattened.append(arr_flat) # The flattened array itself (so we can do stats on it)
    
# Make a nice plot comparing noise from the three sequences

hbt.figsize((15,10))
for i, sequence in enumerate(sequences):    
    plt.plot(bins[0:-1], hist_flattened[i], label="{}, N={}".format(sequence, n[i]), marker = '', color=colors[i],
             alpha = 0.5)
    plt.plot(bins[0:-1], hist_single[i],    label="{}, N={}".format(sequence, 1),    marker = '.', color=colors[i], 
             linestyle = 'none')
    
plt.yscale('log')
plt.legend()
plt.xlabel('DN per pixel')
plt.ylabel('# of pixels')
plt.title('Noise Histogram vs. Exposure Stacking, KEM Approach Images September 2017')          
plt.show()

# =============================================================================
# Compare noise: 10 x 30 sec  vs  30 x 10 sec. Also 30 x 30 sec
# =============================================================================

hist_single = []

hist_flattened = [] # Histogram from stacked frames
hist_single    = [] # Histogram from single frame 
n              = [] # Number of frames summed
arr_flattened  = []
exptime        = []

sequences = ['30sec_rot0_sep', '10sec_rot0_sep', '30sec_rot0_sep']
dn_min = -50
dn_max = 200
bins = hbt.frange(dn_min, dn_max) # Get a list of all of the timestep bins, inclusive, for this file.
frames_max = [10, 30, 30]  # Max frames in a stack.
scaling    = [3, 1, 1]

colors = ['red', 'green', 'blue']

for i,sequence in enumerate(sequences):
    
    indices  = eval('indices_' + sequence).copy()   # Python trick: evaluate a statement
    w = np.where(indices)[0]
    indices[w[frames_max[i]]:] = False          # De-activate all frames above frames_max
    
# XXX need to limit this to N frames max, for consistency against all the stacks
    
    images.set_indices(indices)
    
# Perform the flattening, and get the result back

    arr_flat = images.flatten()
    t = images.t
    
    arr = images.arr    # Get the raw stack, un-flattened. Only valid after calling flatten().
    
    (h, junk)  = np.histogram(arr_flat, bins)
    (h2, junk) = np.histogram(arr[0,:,:], bins)

    # Save the output from this stack calculation
    
    hist_flattened.append(h)   # Histogram of stacked images
    hist_single.append(h2)     # Histogram of 0th image in the stack
    n.append(np.sum(indices))  # Number of images in each stack
    arr_flattened.append(arr_flat) # The flattened array itself (so we can do stats on it)
    exptime.append(t['exptime'][indices][0])
    
# Make a plot comparing noise from the sequences
# Not really sure about the conversion to DN/sec here.    

hbt.figsize((15,10))
for i, sequence in enumerate(sequences):    
    plt.plot(bins[0:-1]/exptime[i], hist_flattened[i]*exptime[i], label="{}, N={}".format(sequence, n[i]), marker = '', color=colors[i],
             alpha = 0.5)
    
plt.yscale('log')
plt.legend()
plt.xlabel('DN per pixel per sec')
plt.ylabel('# of pixels')
plt.title('Noise Histogram vs. Exposure Stacking, KEM Approach Images September 2017')          
plt.show()

# =============================================================================
# Now see how well we can subtract one stack from another. For this, must use median.
# =============================================================================

hist_single = []

hist_flattened = [] # Histogram from stacked frames
hist_single    = [] # Histogram from single frame 
n              = [] # Number of frames summed
arr_flattened  = []

sequences = ['30sec_rot0_sep', '20sec_rot0_sep', '10sec_rot0_jan']
bins = hbt.frange(0, 200) # Get a list of all of the timestep bins, inclusive, for this file.

for sequence in sequences:
    
    indices  = eval('indices_' + sequence)   # Python trick: evaluate 

    images.set_indices(indices)
    
# Perform the flattening, and get the result back

    arr_flat = images.flatten(method='median')
    
    arr = images.arr    # Get the raw stack, un-flattened. Only valid after calling flatten().
    
    (h, junk)  = np.histogram(arr_flat, bins)
    (h2, junk) = np.histogram(arr[0,:,:], bins)

    hist_flattened.append(h)
    hist_single.append(h2)
    n.append(np.sum(indices))
    arr_flattened.append(arr_flat)

# Normalize them to each other -- September
    
ratio = np.nanmean(arr_flattened[0]) / np.nanmean(arr_flattened[1])
arr_flattened[1] = arr_flattened[1] * ratio

# Normalize them to each other -- September:Jan
    
ratio = np.nanmean(arr_flattened[0]) / np.nanmean(arr_flattened[2])
arr_flattened[2] = arr_flattened[2] * ratio

# Now make a plot showing how well we can subtract one stack from another. To how many DN. Sep-Sep.

diff = (arr_flattened[1] - arr_flattened[0])
(h, junk) = np.histogram(diff, bins)
h_cum = np.cumsum(h)
h_cum_frac = h_cum / np.amax(h_cum)

hbt.figsize((12,10))

plt.plot(bins[0:-1], h_cum_frac*100)
plt.yscale('linear')
plt.title("{} - {}".format(sequences[0], sequences[1]))
plt.xlabel('Per-pixel change in $\Delta$DN after subtraction of two fields')
plt.ylabel('% of pixels with < $\Delta$DN')
plt.ylim((90,100))
plt.show()

plt.imshow(stretch(diff))
plt.title("{} - {}".format(sequences[0], sequences[1]))
plt.show()

plt.imshow(stretch(arr_flattened[0]))
plt.title('{}, N = {}'.format(sequences[0], n[0]))
plt.show()
    
plt.imshow(stretch(arr_flattened[1]))
plt.title('{}, N = {}'.format(sequences[1], n[1]))
plt.show()

# Now do again for Jan-Sep

tvec = ird.translation(arr_flattened[2], arr_flattened[0])['tvec']

diff = (arr_flattened[2] - np.roll(np.roll(arr_flattened[0], int(tvec[0]), axis=0), int(tvec[1]), axis=1))

(h, junk) = np.histogram(diff[300:,:], bins)  # Exclude the portion that is not overlapped, of course!
h_cum = np.cumsum(h)
h_cum_frac = h_cum / np.amax(h_cum)


plt.plot(bins[0:-1], h_cum_frac*100)
plt.yscale('linear')
plt.title("{} - {}".format(sequences[0], sequences[2]))
plt.xlabel('Per-pixel change in $\Delta$DN after subtraction of two fields')
plt.ylabel('% of pixels with < $\Delta$DN')
plt.ylim((90,100))
plt.show()

plt.imshow(stretch(diff))
plt.title("{} - {}".format(sequences[0], sequences[2]))
plt.show()

plt.imshow(stretch(arr_flattened[0]))
plt.title('{}, N = {}'.format(sequences[0], n[0]))
plt.show()
    
plt.imshow(stretch(arr_flattened[2]))
plt.title('{}, N = {}'.format(sequences[2], n[2]))
plt.show()


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
# Calculate how many DN MU69 should be at encounter (K-20d, etc.)
# Or alternatively, convert all of my DN values, to I/F values
# =============================================================================

        # Convert DN values in array, to I/F values
        
RSOLAR_LORRI_1X1 = 221999.98  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
RSOLAR_LORRI_4X4 = 3800640.0  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)

C = arr[0] # Get the DN values of the ring. Typical value is 1 DN.

# Define the solar flux, from Hal's paper.

FSOLAR_LORRI  = 176.	     	    # We want to be sure to use LORRI value, not MVIC value!
F_solar = FSOLAR_LORRI # Flux from Hal's paper

RSOLAR = RSOLAR_LORRI_4X4

# Calculate the MU69-Sun distance, in AU (or look it up). 

r_sun_mu69 = 43.28 # Solar distance in AU, from GV.

r_nh_mu69 = 0.16671 # NH distance at K-20d, from GV. 0.16 AU = 24.9 M km.

TEXP = t['exptime'][0]   # Look up exposure time, in seconds
TEXP = 30

I = C / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. All v similar, except for spectrum assumed.

# Apply Hal's conversion formula from p. 7, to compute I/F and print it.

IoF = math.pi * I * r_sun_mu69**2 / F_solar # Equation from Hal's paper

# Now convert to 'normal I/F'

# Define mu = cos(e), where e = emission angle, and e=0 is face-on.

e = 0 # Assume ring is face-on. The sunflower orbits are. 

mu = np.cos(e)

#        mu_2D = np.transpose(np.tile(mu, (self.num_bins_radius(),1)))

# Calculate the normal I/F

IoF_normal = 4 * mu * IoF

# Save this result into 'self,' replacing the DN values with the normal I/F values

#        self.profile_radius_arr = IoF_normal

# Change the units

#        self.profile_radius_units = 'Normal I/F'

# Do the same for azimuthal profile

TEXP_2D = np.transpose(np.tile(self.exptime_arr, (self.num_bins_azimuth(),1)))
mu_2D = np.transpose(np.tile(mu, (self.num_bins_azimuth(),1)))
self.profile_azimuth_arr = self.profile_azimuth_arr / TEXP_2D / RSOLAR * math.pi * r**2 / F_solar * 4 * mu_2D

self.profile_azimuth_units = 'Normal I/F'

# Calculate how many DN is MU69 at K-10d. It is unresolved, so it is a slightly different calculation than above.

# Caclulate how many DN is I/F = 1e-6 ring. This is large spatial extent, so it should be just the inverse of above.


IoF_normal = 1e-6
TEXP = 30
IoF = IoF_normal / (4 * mu)     
I = IoF / math.pi / r_sun_mu69**2 * F_solar
C = I * TEXP * RSOLAR  # C is the value in DN
print("I/F_normal = {}, TEXP = {} s -> DN = {}".format(IoF_normal, TEXP, C))
        
#        

#;PPLUTO =  8.7249995e+15 * 1.15 = 1.0033749e+16
        
# Keywords from Hal Weaver email
#        
#;For 1x1:
#;V = -2.5 * alog10(S/texp) +18.76 + CC – AC
#;Where “S” is the raw DN integrated over aperture being used,
#;“texp” is the exposure time in seconds, “CC” is the color correction,
#;and “AC” is the aperture correction (i.e., the value need to correct
#;from the photometry aperture used to an “infinite” aperture).
# 
#;For early-type stars, e.g., O,B,A --> CC = -0.06
#;For F and G stars --> CC = 0
#;For K stars --> CC = +0.4
#;For M stars --> CC = +0.6
#;Pluto --> CC = -0.037
#;Charon --> CC = -0.014
#;Jupiter --> CC = -0.138
#;Pholus --> CC =  0.213
# 
#;Diffuse sensitivity keywords
#;Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
#;RSOLAR =  221999.98
#;RPLUTO =  214583.33
#;RPHOLUS =  270250.00
#;RCHARON =  219166.66
#;RJUPITER =  195583.33
# 
#;Point Source sensitivity keywords
#;Units are (DN/s)/(erg/cm^2/s/A)
#;PSOLAR =  9.0249987e+15
#;PPLUTO =  8.7249995e+15
#;PPHOLUS =  1.0991666e+16
#;PCHARON =  8.9083329e+15
#;xPJUPITER =  7.9524992e+15
# 
#;For 4x4:
#;V = -2.5 * alog10(S/texp) +18.91 + CC – AC
#;with the same definitions as above.
# 
#;Diffuse sensitivity keywords
#;Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
#;RSOLAR =  221999.98 * 1.07 * 16. = 3800640.0
#;RPLUTO =  214583.33 * 1.07 * 16. = 3673666.8
#;RPHOLUS =  270250.00  * 1.07 * 16. = 4626680.0
#;RCHARON =  219166.66  * 1.07 * 16. = 3752133.2
#;RJUPITER =  195583.33  * 1.07 * 16. = 3348386.8
# 
#;Point Source sensitivity keywords
#;Units are (DN/s)/(erg/cm^2/s/A)
#;PSOLAR =  9.0249987e+15 * 1.15 = 1.0378749e+16
#;PPLUTO =  8.7249995e+15 * 1.15 = 1.0033749e+16
#;PPHOLUS =  1.0991666e+16 * 1.15 = 1.2640416e+16
#;PCHARON =  8.9083329e+15 * 1.15 = 1.0244583e+16
#;PJUPITER =  7.9524992e+15 * 1.15 = 9.1453737e+15
        