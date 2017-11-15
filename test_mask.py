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

import re # Regexp
import pickle # For load/save

from   matplotlib.figure import Figure

# HBT imports

import hbt

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'

dir_out = '/Users/throop/data/NH_Jring/out/'

file = 'lor_0034612923_0x630_sci_1.fit'

filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in

index_group_default = 8 # Jupiter ring phase curve
index_image_default = 20 # Image number within the group
        
# Load pickle file

#if os.path.exists(self.filename_save):

print("Loading file: {}".format(filename_save))

#    self.load(verbose=False) # This will load self.t
#    t = self.t
    
lun = open(dir_out + filename_save, 'rb')
t = pickle.load(lun)
lun.close()
        
index_group = index_group_default
index_image  = index_image_default

# Extract the one group we are looking at 
# NB: This creates an entirely new table -- *not* a view into the original table. 
# If we modify values here, we need to explicitly write the changes back to the original.

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

groupmask = t['Desc'] == groups[index_group]
t_group = t[groupmask]

num_images_group = np.size(t_group)
        
# Read the image

hdulist = fits.open(dir_images + '/' + file)
data = hdulist['PRIMARY'].data
stretch = astropy.visualization.PercentileInterval(10)  # PI(90) scales array to 5th .. 95th %ile. 

# Create the backplane

planes = hbt.create_backplane(dir_images + '/' + file)

am = planes['Ang_Metis']

dx_total =  ( t_group['dx_opnav'][index_image] )
dy_total =  ( t_group['dy_opnav'][index_image])

(dx_total, dy_total) =(52,70)
  
am_roll = np.roll(np.roll(am, -dy_total, axis=0), -dx_total, axis=1)

# Plot it

data_stretch = stretch(hbt.remove_sfit(data,degree=5))

plt.imshow(data_stretch)
plt.show()

plt.imshow( data_stretch * (am_roll > 0.0001))
plt.show()
