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


from   matplotlib.figure import Figure

# HBT imports

import hbt

file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

dir = '/Users/throop/Dropbox/Data/NH_KEM_Hazard/'
files = glob.glob(dir + '*/*.fits')

stretch_percent = 90    
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

data     = {}

# Set up the table. "The standard method of representing Python 3 strings in numpy is via the unicode 'U' dtype."
# Disadvantage of this is memory space, but not an issue for me.

t = Table(  [[],              [],          [],         [],        [],       [],       [],      [],   [], 
                                            [], [], [], [], []],
    names=('filename_short', 'exptime', 'visitname', 'sapname', 'sapdesc', 'target', 'reqid', 'et',  'utc', 
                                            'x_pix', 'y_pix', 'ra', 'dec', 'angle'),
    dtype = ('U50',           'float64', 'U50',      'U50',     'U50',     'U50',    'U50',   'float64', 'U50', 
                                            'float64', 'float64', 'float64', 'float64', 'float64'  ))


for i,file in enumerate(files):

    # Read the data

    hdulist        = fits.open(file)
    arr            = hdulist[0].data
    err            = hdulist[1].data
    quality        = hdulist[2].data
    
    filename_short = file.split('/')[-1].replace('.fits', '').replace('lor_', '').replace('_0x633_pwcs','')
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

    t.add_row([filename_short, exptime, visitnam, sapname, sapdesc, target, reqid, et, utc, x_pix, y_pix,
               hdulist[0].header['SPCBRRA'],   # RA (pre-WCS I assume, but not sure)
               hdulist[0].header['SPCBRDEC'],  # Dec
               hdulist[0].header['SPCEMEN']])  # EME North angle -- ie, boresight rotation angle


# Save the data into the dictionary. We would save it to the table, but that doesn't seem possible.
# The key is the short filename
    
    data[filename_short] = arr
    
# Sort by ET.
    
t.sort('et')
    
# Plot one image. Plot MU69 position on top of it to be sure that WCS is correct.

# Remove a few columns that we don't need, or that are wrong.
    
t.remove_column('sapdesc')
t.remove_column('sapname')
t.remove_column('target')
t.remove_column('visitname')

# Print the table

t.pprint(max_lines = -1, max_width = -1)
 
for i in range(120,150):

        
    plt.imshow(stretch(arr))
    plt.title("{}, exptime {}".format(t['filename_short'][i], t['exptime'][i]))
    plt.plot(x_pix, y_pix, marker = '+', color='red')
    
    plt.ylabel('Y pixels')
    plt.xlabel('X pixels')
    
    plt.show()
    
    hdulist.close()


plt.plot(t['x_pix'], t['y_pix'], marker = '+', linestyle='none')
plt.xlabel('MU69 X pos [pix]')
plt.ylabel('MU69 Y pos [pix]')
plt.show()

plt.plot(t['ra'], t['dec'], marker = '+', linestyle = 'none')
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
plt.show()

# Q: can I please put a numpy array into a table?

t2 = Table(  [ [], [], [] ],
    names=('filename_short', 'data', 'reqid'),
    dtype = ('U30', 'float64', 'U30') )
t2.add_row([filename_short, arr, reqid])

t3 = Table()

#=============================================================================
# OK! Now time to do some data analysis.
#=============================================================================

indices_sep17 = t['et'] > sp.utc2et('15 sep 2017')
indices_rot0  = t['angle'] < 180   # One rotation angle
indices_rot90 = t['angle'] > 180   # The other rotation angle
indices_10sec = t['exptime'] < 10
indices_20sec = np.logical_and( t['exptime'] < 20, t['exptime'] > 10 )
indices_30sec = np.logical_and( t['exptime'] < 30, t['exptime'] > 20 )

indices_30sec_rot0  = np.logical_and(indices_30sec, indices_rot0)  # 94
indices_30sec_rot90 = np.logical_and(indices_30sec, indices_rot90) # 0

indices_20sec_rot0  = np.logical_and(indices_20sec, indices_rot0)  # 93
indices_20sec_rot90 = np.logical_and(indices_20sec, indices_rot90) # 0

indices_10sec_rot0  = np.logical_and(indices_10sec, indices_rot0)  # 96
indices_10sec_rot90 = np.logical_and(indices_10sec, indices_rot90) # 48


# Let's start by trying to co-add all of the 30-sec exposures

indices = indices_10sec_rot0.copy()

# Create the output numpy array
num_planes = np.sum(indices)
arr = np.zeros((num_planes, 1024, 1024))
w = np.where(indices)[0]  # List of all the indices
for j,w_i in enumerate(w):  # Loop over all the indices. 'j' counts from 0 to N. w_i is each individual index.
    im = data[t['filename_short'][w_i]]
    im_expand = scipy.ndimage.zoom(im,4)   
    arr[j,:,:] = im_expand                 # Stick it into place

arr_flat = np.nanmean(arr,0)
    
# First get their indices

# Then expand them by 4x4

# Then shift them as per WCS

# Then co-add, and/or median them