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
file_tm    = "/Users/throop/git/NH_rings/kernels_nh_jupiter.tm"  # SPICE metakernel

file = dir + 'lor_0034961819_0x630_sci_1_opnav.fit' # noise, can't see much
file = dir + 'lor_0034962025_0x630_sci_1_opnav.fit' # Metis on ansa. dt=40 way off. 
#file = dir + 'lor_0034962461_0x630_sci_1_opnav.fit' # Ring, but not on ansa -- not a good test

#file = dir + 'lor_0034616523_0x630_sci_1_opnav.fit' # Adrastea in middle right. dt =40 works.
#file = dir + 'lor_0034618323_0x630_sci_1_opnav.fit' # Adrastea in on ansa
#file = dir + 'lor_0034620123_0x630_sci_1_opnav.fit' # Adrastea in middle left. Fits dt=40 better than dt=0.
#                                                   # dt=0 has cross 5 pix down right. 

#file = dir + 'lor_0035103963_0x633_sci_1_opnav.fit' # 4x4, gossamer, Amal in middle. Bad navigation - not useful.
#file = dir + 'lor_0035117161_0x633_sci_1_opnav.fit' # 4x4, gossamer. Bad navigation - not useful.

#file = dir + 'lor_0034715044_0x630_sci_1_opnav.fit' # 1x1, Metis on LHS. Fits dt = 60. Cross is 10 pix down from o.
#file = dir + 'lor_0034627923_0x630_sci_1_opnav.fit' # 1x1, bright sat Amalthea passing by. Can't tell dt from it.

# Initialize SPICE

sp.furnsh(file_tm)

dt = 0.001  # Offset in flight time along line-of-sight. Nominally zero, but might be otherwise.

with warnings.catch_warnings():   # Block the warnings
    warnings.simplefilter("ignore")    
    w  = WCS(file)           # Look up the WCS coordinates for this frame

# Load the image and process the header

im     = hbt.read_lorri(file)
header = hbt.get_image_header(file)
et     = header['SPCSCET']
dx_pix = header['NAXIS1']
dy_pix = header['NAXIS2']
utc    = sp.et2utc(et, 'C', 1)

stretch_percent = 90

stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile. 

NUM_STARS_PHOT = 100  # How many stars to use from DAOPhot. For noisy images, DAO will find a lot of
                          # fake stars, so we need to crank this up higher than the # of cat stars.
NUM_STARS_CAT  = 50  # How many stars to use from star catalog

#==============================================================================
# Find the catalog stars to plot
#==============================================================================

# Define the star catalog.
# On tomato, for some reason name lookup does not work. But URL lookup always works.

name_cat = u'Guide Star Catalog v2 1'
url_cat = 'http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=GSC23&'

radius_search_deg = 0.15
    
with data.conf.set_temp('remote_timeout', 30): # This is the very strange syntax to set a timeout delay.
                                               # The default is 3 seconds, and that times out often.
    with warnings.catch_warnings():   # Block the warnings
        warnings.simplefilter("ignore")    
        stars_cat = conesearch.conesearch(w.wcs.crval, radius_search_deg, cache=True, catalog_db = url_cat)

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
#%%

name_bodies = np.array(['Metis', 'Adrastea', 'Thebe', 'Amalthea', 'Io'])        

# Look up the times. XXX This does not properly incorporate the offset in time. It should be another argument,
# and it is not just additive to et.

# For abcorr, the stellar aberration is a constant shift which has been applied during opnav (ie, the offset between
# catalog positions, and observed positions). We do not want to apply it again here to just the sats. So, use 'LT'
# alone.

x_bodies_pix,  y_bodies_pix   = hbt.get_pos_bodies(et, name_bodies, units='pixels',  abcorr='LT', wcs=w, dt=dt)
x_bodies_deg,  y_bodies_deg   = hbt.get_pos_bodies(et, name_bodies, units='degrees', abcorr='LT', wcs=w, dt=dt)

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
DO_PLOT_RING_INNER = True
DO_PLOT_RING_OUTER = True
marker_ring = 'None'

color_ring = 'blue'
alpha_ring = 0.3

hbt.figsize((10,10))
plt.set_cmap('Greys_r')

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


plt.xlim([0,dx_pix-1])  # This is an array and not a tuple. Beats me, like so many things with mpl.
plt.ylim([dx_pix-1,0])

# Plot the rings

if (DO_PLOT_RING_OUTER):
    plt.plot(x_ring2, y_ring2, marker=marker_ring, color = color_ring, ls = '--',
                  alpha = alpha_ring, label = 'Ring outer')

if (DO_PLOT_RING_INNER):
    plt.plot(x_ring1, y_ring1, marker=marker_ring, color='green', ls = '-', \
        ms=8, alpha = alpha_ring, label='Ring inner')

# Plot the satellites

for i in range(np.size(name_bodies)):
    plt.plot(x_bodies_pix[i], y_bodies_pix[i], marker = marker_sat, color=color_sat, markersize=20, linestyle='none')
    plt.text(x_bodies_pix[i], y_bodies_pix[i], '   ' + name_bodies[i] + ', code', clip_on=True, color=color_sat)

# Load values for Adrastea, from GV-dev local default

ra_adras_gv = 252.16805*hbt.d2r
dec_adras_gv = -18.76189*hbt.d2r

# Load values for Adrastea, from GV 15sci_rhr on server

ra_adras_gv = 252.15338*hbt.d2r
dec_adras_gv = -18.76665*hbt.d2r

radec_adras_gv = np.transpose(np.array((ra_adras_gv, dec_adras_gv)))
x_adras_gv, y_adras_gv = w.wcs_world2pix(radec_adras_gv[0]*hbt.r2d, radec_adras_gv[1]*hbt.r2d, 0)

if (0<x_adras_gv<1023) and (0 < y_adras_gv < 1023):
    plt.plot(x_adras_gv, y_adras_gv, marker = '*', color='lightblue', markersize=10, clip_on=True)
    plt.text(x_adras_gv, y_adras_gv, '                                     Adrastea, GV', color='lightblue', clip_on=True)

plt.legend(loc = 'upper left') 
plt.show()

#==============================================================================
# Debugging: Print out RA/Dec of satellite positions, to compare to GV
#==============================================================================

#%%

print()
print("file   = {}".format(file))
print("UTC    = {}".format(utc))
print("ET     = {}".format(et))
print("TOF dt = {} sec".format(dt))

for i in range(np.size(name_bodies)):
    print("{:10}: RA={:.4f}, Dec={:.4f}".format(name_bodies[i], x_bodies_deg[i], y_bodies_deg[i]))
    
    