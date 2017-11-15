#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 10:55:04 2016

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
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
#import wcsaxes
import time
from   scipy.interpolate import griddata
from astropy.stats import sigma_clip

import re # Regexp
import pickle # For load/save

from skimage.io import imread, imsave

# HBT imports

import hbt

plt.set_cmap('Greys_r')

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'

dir_jring = '/Users/throop/data/NH_Jring/'
dir_out = dir_jring + 'out/'

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
        
# Read the image

hdulist = fits.open(file)
data = hdulist['PRIMARY'].data

# Clip the image, removing extrema

#data_clip = sigma_clip(data, sigma=3, iters=10)

stretch = astropy.visualization.PercentileInterval(90)  # PI(90) scales array to 5th .. 95th %ile. 

styles = ['', 'aggressive', 'nominal', 'ring_broad_only', 'ring_narrow_only']

for style in styles:
    if (style):
        file_mask = dir_jring + '/masks/mask_7_32-35_{}.png'.format(style)
        mask = (imread(file_mask)[:,:,0])   # Flatten to 2D
        mask = (mask > 0)

    else:
        mask = np.ones(np.shape(data))
        mask = (mask > 0)

    # Remove the sfit from the masked image
    
    out = hbt.remove_sfit(data, mask = mask)
    
    # Calculate the min and max values. But do this of the data, where it is masked True.
    # And, exclude outliers (like CRs) when we do this.
    
    out_masked = out[mask]
    out_masked_clipped = sigma_clip(out_masked, sigma_lower=1000, sigma_upper=3, iters=10)
    
    mm = hbt.mm(out_masked_clipped)

    stretch_mm = astropy.visualization.ManualInterval(vmin=mm[0], vmax=mm[1])

    # Plot the image, using these properly calculated robust min and maxes, in the masked area only
    
    plt.imshow(stretch_mm(out*mask))
        
    plt.title(style)
    plt.show()

# Now I need to take a radial profile of each of these...
    


# OK. We want to scale it to the min and max *within the masked region*.
