#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 13:39:08 2018

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
from   astroquery.vo_conesearch import conesearch # Virtual Observatory, ie star catalogs
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

import cProfile # For profiling

# Imports for Tk

#import Tkinter # change Tkinter -> tkinter for py 2 - 3?
import tkinter
from tkinter import ttk
from tkinter import messagebox
tkinter.messagebox
#import tkMessageBox #for python2
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

from   astropy.stats import sigma_clip

# HBT imports

import hbt

plt.set_cmap('plasma')

# Now list all of the Gossamer observations. I have determined these groupings manually, based on timings.
# In general it looks like usually these are stacks of four frames at each pointing.

index_group = 6
index_image = hbt.frange(59,62)
index_image = hbt.frange(63,66)
index_image = hbt.frange(67,70)
index_image = hbt.frange(71,74)
index_image = hbt.frange(75,78)

index_image = hbt.frange(95,99)

index_image = hbt.frange(112,115)
index_image = hbt.frange(116,119)
index_image = hbt.frange(120,123)
index_image = hbt.frange(124,127)
index_image = hbt.frange(128,131)

#index_image = hbt.frange(100,102)

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters 

dir_out    = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
dir_backplanes = '/Users/throop/data/NH_Jring/out/'

lun = open(dir_out + filename_save, 'rb')
t_all = pickle.load(lun)
lun.close()

groups = astropy.table.unique(t_all, keys=(['Desc']))['Desc']

groupmask = t_all['Desc'] == groups[index_group]

t_group = t_all[groupmask]  # 

arr_sum = np.zeros((256,256))

for i in index_image:        
    t_i = t_group[i]  # Grab this, read-only, since we use it a lot.
    
    arr = hbt.read_lorri(t_i['Filename'])
    arr_sum += arr

arr_sum_s = hbt.remove_sfit(arr_sum,degree=3)
plt.imshow(stretch(arr_sum_s))
plt.show()

# =============================================================================
# Make a plot showing the center position of a lot of these frames.
# =============================================================================

file_tm = 'kernels_nh_jupiter.tm'  # SPICE metakernel

sp.furnsh(file_tm)

index_group = 6
index_images = np.array(
              list(hbt.frange(10,15)) +
              list(hbt.frange(46,54)) +
              list(hbt.frange(59,135)) )

ra_arr      = []
dec_arr     = []
dist_rj_arr = []

for index_image in index_images:
    t_i = t_group[index_image]  # Grab this, read-only, since we use it a lot.
    
    im = hbt.read_lorri(t_i['Filename'])
    
    w = WCS(t_i['Filename'])

    # Get the RA / Dec for the four corner points
    
    radec = w.wcs_pix2world([[0,0],[0,256],[256,256],[256,0],[0,0]],0)
    radec_center = w.wcs_pix2world([[128,128]],0)
    
    ra = radec[:,0]
    dec = radec[:,1]

    plt.plot(ra, dec) # Plot all the four points for this box
    
#    ra_arr.append(radec)
#    dec_arr.append(dec_deg)

    # Now, of course I don't really want RA / Dec, but projected distance (RJ), as well as height above the midplane.
    
    # So for each of these points, take a ray. 
    # Make a plane, which is normal to observer, and goes thru Jupiter. 
    # Intersect the ray, with this plane
    # Get the position in that plane (which will probably be in units of KM from Jupiter center)
    # Save that position, and plot them.
    
    # Create a plane centered at Jupiter
    
    abcorr = 'LT+S'
    
    (vec,lt) = sp.spkezr('Jupiter', t_i['ET'], 'J2000', abcorr, 'New Horizons')
    vec_nh_jup = vec[0:3]

    (vec,lt) = sp.spkezr('Jupiter', t_i['ET'], 'J2000', abcorr, 'Sun')
    vec_sun_jup = vec[0:3]  # This is position of Jupiter center, in J2000, from Sun. (Maybe should do from SS baryctr?)

    (vec,lt) = sp.spkezr('New Horizons', t_i['ET'], 'J2000', abcorr, 'Sun')
    vec_sun_nh = vec[0:3]  # This is position of NH, in J2000.

    plane = sp.nvp2pl(vec_nh_jup, vec_sun_jup)
    
    vec_nh_pixel = sp.radrec(1, radec_center[0,0], radec_center[0,1])  # This is the vector from NH, to the specific LORRI pixel
    
    # Get the intersection of the pixel ray, with the Jupiter-centered plane
    
    (nxpts, xpt) = sp.inrypl(vec_sun_nh, vec_nh_pixel, plane)
    
    dist_rj = sp.vnorm(xpt - vec_sun_jup) / 70000

    dist_rj_arr.append(dist_rj)
    
plt.xlabel('RA  [deg]')    
plt.ylabel('Dec [deg]')
plt.title(f'All Gossamer images, {index_group}/{np.amin(index_image)} .. {np.amax(index_image)}')
 
plt.show()    

# And hang on. We've already made backplanes for everything. Can I use those here? 
# I don't know. These images are obviously very edge-on.

        file_backplane = dir_backplanes + t_group['Shortname'][index_image].replace('.fit', '_planes.pkl')

        # Save the shortname associated with the current backplane. 
        # That lets us verify if the backplane for current image is indeed loaded.

        file_backplane_shortname = t_group['Shortname'][index_image]
				
        if (os.path.isfile(file_backplane)):
            lun = open(file_backplane, 'rb')
            planes = pickle.load(lun)
            lun.close()

# I could do RA/Dec relative to RA/Dec of Jupiter.
            