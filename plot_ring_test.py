#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 14:56:53 2017

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

from nh_jring_mask_from_objectlist import nh_jring_mask_from_objectlist

dir_image = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'
file_tm    = "/Users/throop/git/NH_rings/kernels_nh_jupiter.tm"  # SPICE metakernel

#file_image = 'lor_0034961819_0x630_sci_1_opnav.fit' # noise, can't see much

#file_image = 'lor_0034962025_0x630_sci_1_opnav.fit' # Metis on ansa.

#file_image + 'lor_0034962461_0x630_sci_1_opnav.fit' # Ring, but not on ansa -- not a good test

file_image = 'lor_0034616523_0x630_sci_1_opnav.fit' # Adrastea in middle right. dt =40 works. Good test.
#file_image = 'lor_0034618323_0x630_sci_1_opnav.fit' # Adrastea in on ansa
#file_image = 'lor_0034620123_0x630_sci_1_opnav.fit' # Adrastea in middle left. Fits dt=40 better than dt=0.
#                                                   # dt=0 has cross 5 pix down right. 

#file_image = 'lor_0035103963_0x633_sci_1_opnav.fit' # 4x4, gossamer, Amal in middle. Bad navigation - not useful.
#file_image = 'lor_0035117161_0x633_sci_1_opnav.fit' # 4x4, gossamer. Bad navigation - not useful.

#file_image = 'lor_0034715044_0x630_sci_1_opnav.fit' # 1x1, Metis on LHS. Fits dt = 60. Cross is 10 pix down from o.
#file_image = 'lor_0034627923_0x630_sci_1_opnav.fit' # 1x1, bright sat Amalthea passing by. Can't tell dt from it.

# Initialize SPICE

sp.furnsh(file_tm)

dir_out    = '/Users/throop/Data/NH_Jring/out/'

file_objects = file_image.replace('.fit', '_objects.txt')

path_objects = dir_out + file_objects
    
#file_image = objectfile_base.replace('_objects', '').replace('.txt', '.fit')
dir_images = '/Users/throop/Data/NH_Jring/data/jupiter/level2/lor/all/'

# Read the image and header

im     = hbt.read_lorri(dir_images + file_image)
im     = hbt.lorri_destripe(im)

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

# And make some plot to show if it matches up, or not

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile.     
    
im_masked = im.copy()
im_masked[mask == True] = np.median(im)

im_masked_t = im.copy()
im_masked_t[np.transpose(mask) == True] = np.median(im)

#==============================================================================
# Find location of ring points
#==============================================================================

a_ring_outer_km = 129300
a_ring_inner_km = 122000

x_ring1, y_ring1 = hbt.get_pos_ring(et, name_body='Jupiter', radius=a_ring_inner_km, units='pixels', abcorr='LT', wcs=w)
x_ring2, y_ring2 = hbt.get_pos_ring(et, name_body='Jupiter', radius=a_ring_outer_km, units='pixels', abcorr='LT', wcs=w)
      
#==============================================================================
# Make a plot!
#==============================================================================

color_phot = 'red'            # Color for stars found photometrically
color_cat  = 'lightgreen'     # Color for stars in catalog  
color_sat  = 'yellow'
marker_sat = '+'
marker_ring = 'None'
color_ring = 'blue'
alpha_ring = 0.3

DO_PLOT_RING_INNER = True
DO_PLOT_RING_OUTER = True

hbt.figsize((10,10))
plt.set_cmap('Greys_r')

fig, ax = plt.subplots(1, 2, figsize=(20, 10))

dx = 0
dy = 0        

for i in [0,1]:
    ax[i].set_xlim([0,dx_pix-1])  
    ax[i].set_ylim([dx_pix-1,0])

ax[0].imshow(stretch(im))
ax[1].imshow(stretch(im_masked))

# Plot the rings

if (DO_PLOT_RING_OUTER):
    ax[0].plot(x_ring2, y_ring2, marker=marker_ring, color = color_ring, ls = '--',
                  alpha = alpha_ring, label = 'Ring outer')

if (DO_PLOT_RING_INNER):
    ax[0].plot(x_ring1, y_ring1, marker=marker_ring, color='green', ls = '-', \
        ms=8, alpha = alpha_ring, label='Ring inner')

ax[0].legend(loc = 'upper left') 

plt.show()

plt.imshow(stretch(im + np.transpose(10*mask)))
  
    