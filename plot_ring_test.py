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
import wcsaxes
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

# HBT imports

import hbt

dir = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'

file = dir + 'lor_0034961819_0x630_sci_1_opnav.fit' # noise, can't see much
file = dir + 'lor_0034962025_0x630_sci_1_opnav.fit' # Sat on edge
#file = dir + 'lor_0034962461_0x630_sci_1_opnav.fit' # Ring, but not on ansa -- not a good test

file = dir + 'lor_0034616523_0x630_sci_1_opnav.fit' # Adrastea in middle right
file = dir + 'lor_0034618323_0x630_sci_1_opnav.fit' # Adrastea in on ansa
file = dir + 'lor_0034620123_0x630_sci_1_opnav.fit' # Adrastea in middle left

file_tm    = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

# Initialize SPICE

sp.furnsh(file_tm)

dt = 0  # Offset in flight time along line-of-sight.

w = WCS(file)

im = hbt.read_lorri(file)
header = hbt.get_image_header(file)
et = header['SPCSCET']
utc = sp.et2utc(et, 'C', 1)

stretch_percent = 90

stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile. 

NUM_STARS_PHOT = 100  # How many stars to use from DAOPhot. For noisy images, DAO will find a lot of
                          # fake stars, so we need to crank this up higher than the # of cat stars.
NUM_STARS_CAT  = 50  # How many stars to use from star catalog

#==============================================================================
# Find the catalog stars to plot
#==============================================================================

name_cat = u'Guide Star Catalog v2 1'

radius_search_deg = 0.15
    
with data.conf.set_temp('remote_timeout', 30): # This is the very strange syntax to set a timeout delay.
                                               # The default is 3 seconds, and that times out often.
    stars_cat = conesearch.conesearch(w.wcs.crval, radius_search_deg, cache=True, catalog_db = name_cat)

ra_stars_cat  = np.array(stars_cat.array['ra'])*hbt.d2r # Convert to radians
dec_stars_cat = np.array(stars_cat.array['dec'])*hbt.d2r # Convert to radians

        
mag       = np.array(stars_cat.array['Mag'])

# Extract just the brightest catalog stars

order = np.argsort(mag)
order = np.array(order)[0:NUM_STARS_CAT]

ra_stars_cat = ra_stars_cat[order]
dec_stars_cat = dec_stars_cat[order]

radec_stars_cat        = np.transpose(np.array((ra_stars_cat, dec_stars_cat)))

# Convert RA/Dec values back into pixels

x_stars_cat,    y_stars_cat      = w.wcs_world2pix(radec_stars_cat[:,0]*hbt.r2d,   radec_stars_cat[:,1]*hbt.r2d, 0)

points_cat = np.transpose(np.array([x_stars_cat, y_stars_cat])) # Make a list of the catalog stars

#==============================================================================
# Find the DAOphot stars
#==============================================================================

points_stars_phot = hbt.find_stars(im, num=50) # Returns N x 2 aray. 0 = Row = y; 1 = Column = x.

y_stars_phot =(points_stars_phot[:,0]) # xy is correct -- see above
x_stars_phot =(points_stars_phot[:,1]) # 

#==============================================================================
# Find location of Metis, Adrastea, Thebe
#==============================================================================

name_bodies = np.array(['Metis', 'Adrastea', 'Thebe', 'Amalthea', 'Io'])        

# Look up the times. XXX This does not properly incorporate the offset in time. It should be another argument,
# and it is not just additive to et.

x_bodies_pix,  y_bodies_pix   = hbt.get_pos_bodies(et + dt, name_bodies, units='pixels', wcs=w)
x_bodies_deg,  y_bodies_deg   = hbt.get_pos_bodies(et + dt, name_bodies, units='degrees', wcs=w)
            
#==============================================================================
# Make a plot of stars: Catalog vs. Photometric
#==============================================================================

#%%
color_phot = 'red'            # Color for stars found photometrically
color_cat  = 'lightgreen'     # Color for stars in catalog  
color_sat  = 'red'
marker_sat = '+'

hbt.figsize((10,10))

plt.imshow(stretch(im))

dx = 0
dy = 0

plt.plot(x_stars_cat + dx, y_stars_cat + dy, 
         marker='o', ls='None', 
         color=color_cat, alpha = 0.25, ms=12, mew=1, label = 'Cat Stars, using _opnav')
         
#plt.plot(x_pos_star_cat, y_pos_star_cat, 
#         eval('np.' + t['y_pos_star_cat']), 
#         marker='o', ls='None', 
#         color=color_cat, alpha = 1, ms=4, mew=1, label = 'Cat Stars, raw')

plt.plot(x_stars_phot, y_stars_phot, 
         marker='o', ls='None', 
         color='none', markersize=10, mew=1, mec=color_phot, alpha = 1,
         label = 'DAOfind Stars')               


plt.xlim([0,1023])  # This is an array and not a tuple. Beats me, like so many things with mpl.
plt.ylim([1023,0])

# Plot the satellites
for i in range(np.size(name_bodies)):
    plt.plot(x_bodies_pix[i], y_bodies_pix[i], marker = marker_sat, color=color_sat, markersize=20, linestyle='none')
    plt.text(x_bodies_pix[i], y_bodies_pix[i], '   ' + name_bodies[i], clip_on=True)

ra_adras_gv = 252.16805*hbt.d2r
dec_adras_gv = -18.76189*hbt.d2r

radec_adras_gv = np.transpose(np.array((ra_adras_gv, dec_adras_gv)))
x_adras_gv, y_adras_gv = w.wcs_world2pix(radec_adras_gv[0]*hbt.r2d, radec_adras_gv[1]*hbt.r2d, 0)

plt.plot(x_adras_gv, y_adras_gv, marker = '*', color='lightblue', markersize=5)

plt.legend(loc = 'upper left') 
plt.show()

#==============================================================================
# Debugging: Print out RA/Dec of satellite positions, to compare to GV
#==============================================================================

#%%

print()
print("file = {}".format(file))
print("UTC  = {}".format(utc))
print("ET   = {}".format(et))

for i in range(np.size(name_bodies)):
    print("{:10}: RA={:.4f}, Dec={:.4f}".format(name_bodies[i], x_bodies_deg[i], y_bodies_deg[i]))
    
    