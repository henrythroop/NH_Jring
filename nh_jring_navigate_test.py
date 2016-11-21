#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:54:05 2016

@author: throop
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 14:55:00 2016

@author: throop

"""

# -*- coding: utf-8 -*-
"""
# Program display a GUI of NJ J-ring data to navigate and extract it.

# Possible features to be added:
#  o Output images as PNG
#  o Select multiple images and display them as an arrayf
#  o Slider to scale the image brightness / contrast
#  o Subtract median? Subract neighbor image?

# Widget program created 18-May-2016 based on previous text version (same name), and salt_interact.py .

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
from scipy.interpolate import griddata

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
import warnings

# HBT imports

import hbt

# First we define any general-purpose functions, which are not part of the class/module.
# We can move these to a different file at some point.

plt.set_cmap('Greys_r')

dir_data   = '/Users/throop/Dropbox/Data/NH_Jring/data/jupiter/level2/lor/all/'
file       = dir_data + 'lor_0034765323_0x630_sci_1.fit'
file_tm    = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel
dir_out    = '/Users/throop/Dropbox/Data/NH_Jring/out/'
file_short = file.split('/')[-1]  

d2r     = hbt.d2r
r2d     = hbt.r2d

# Load the image

hdulist  = fits.open(file)
image    = hdulist['PRIMARY'].data
header   = hdulist['PRIMARY'].header

# Read the WCS coordinates

with warnings.catch_warnings():
    warnings.simplefilter("ignore")    
    w = WCS(file)                # Look up the WCS coordinates for this frame
                                 # Otherwise it gives "FITSFixedWarning: 'unitfix': 'Changed units: 'DEG' -> 'deg'"
# Read the WCS parameters
           
center  = w.wcs.crval  # degrees. # crval is a two-element array of [RA, Dec], in degrees

# Initialize SPICE

sp.furnsh(file_tm)
et = header['SPCSCET']
utc = sp.et2utc(et, 'C', 1)

# Stretch the image

stretch_percent = 90
stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales array to 5th .. 95th %ile. 

# Display it

plt.imshow(stretch(image))

# Load matching stars

DO_GSC1     = False    # Stopped working 2-Oct-2016
DO_GSC2     = True
DO_USNOA2   = False

#==============================================================================
# def navigate(file, et):        
#==============================================================================

DO_GSC1     = False    # Stopped working 2-Oct-2016
DO_GSC2     = True
DO_USNOA2   = False

if (DO_GSC1):
    name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
    radius_search = 0.15
    stars = conesearch.conesearch(w.wcs.crval, radius_search, cache=False, catalog_db = name_cat)
    ra_stars  = np.array(stars.array['RAJ2000'])*d2r # Convert to radians
    dec_stars = np.array(stars.array['DEJ2000'])*d2r # Convert to radians
#            table_stars = Table(stars.array.data)

if (DO_GSC2):
    name_cat = u'Guide Star Catalog v2 1'
    file_pickle = dir_out + file_short.replace('.fit', '') + '.stars_gsc.pkl'

    # If there is already a saved pickle file, then load from disk
    
    if os.path.isfile(file_pickle):
        
        print("Loading file: " + file_pickle)
        lun = open(file_pickle, 'rb')
        (ra_stars, dec_stars) = pickle.load(lun)
        lun.close()
            
    else:     
#            name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating

#            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
        from astropy.utils import data
        
        with data.conf.set_temp('remote_timeout', 30): # This is the very strange syntax to set a timeout delay.
                                                       # The default is 3 seconds, and that times out often.
            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
    
        ra_stars  = np.array(stars.array['ra'])*d2r # Convert to radians
        dec_stars = np.array(stars.array['dec'])*d2r # Convert to radians
    
        mag       = np.array(stars.array['Mag'])
        
        print("Stars downloaded: {}; mag = {} .. {}".format(np.size(mag), np.nanmin(mag), np.nanmax(mag)))
        print("RA = {} .. {}".format(np.nanmin(ra_stars)*r2d, np.nanmax(ra_stars)*r2d))
        
        # Now sort by magnitude, and keep the 100 brightest
        # This is because this GSC catalog is huge -- typically 2000 stars in LORRI FOV.
        # We need to reduce its size to fit in our fixed astropy table string length.
    
        num_stars_max = 100            
        order = np.argsort(mag)
        order = np.array(order)[0:num_stars_max]
    
        ra_stars = ra_stars[order]
        dec_stars = dec_stars[order]

        lun = open(file_pickle, 'wb')
        pickle.dump((ra_stars, dec_stars), lun)
        lun.close()
        print("Wrote: " + file_pickle)
        
#            table_stars = Table(stars.array.data)

if (DO_USNOA2):  
    name_cat = u'The USNO-A2.0 Catalogue (Monet+ 1998) 1' # Works but gives stars down to v=17; I want to v=13 
    stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
    table_stars = Table(stars.array.data)
    mask = table_stars['Bmag'] < 13
    table_stars_m = table_stars[mask]            

    ra_stars  = table_stars_m['RAJ2000']*d2r # Convert to radians
    dec_stars = table_stars_m['DEJ2000']*d2r # Convert to radians

# Get an array of points along the ring

ra_ring1, dec_ring1 = hbt.get_pos_ring(et, name_body='Jupiter', radius=122000, units='radec', wcs=w)
ra_ring2, dec_ring2 = hbt.get_pos_ring(et, name_body='Jupiter', radius=129000, units='radec', wcs=w)

# Return as radians

x_ring1, y_ring1    = w.wcs_world2pix(ra_ring1*r2d, dec_ring1*r2d, 0) # Convert to pixels
x_ring2, y_ring2    = w.wcs_world2pix(ra_ring2*r2d, dec_ring2*r2d, 0) # Convert to pixels

# Get position of Metis, in pixels

(vec6,lt) = sp.spkezr('Metis', et, 'J2000', 'LT+S', 'New Horizons')
(junk, ra_metis, dec_metis) = sp.recrad(vec6[0:3])

# Look up velocity of NH, for stellar aberration

abcorr = 'LT+S'
frame = 'J2000'
st,ltime = sp.spkezr('New Horizons', et, frame, abcorr, 'Sun') # Get velocity of NH 
vel_sun_nh_j2k = st[3:6]

# Correct stellar RA/Dec for stellar aberration

radec_stars        = np.transpose(np.array((ra_stars,dec_stars)))
radec_stars_abcorr = hbt.correct_stellab(radec_stars, vel_sun_nh_j2k) # Store as radians

# Convert ring RA/Dec for stellar aberration

radec_ring1        = np.transpose(np.array((ra_ring1,dec_ring1)))
radec_ring1_abcorr = hbt.correct_stellab(radec_ring1, vel_sun_nh_j2k) # radians
radec_ring2        = np.transpose(np.array((ra_ring2,dec_ring2)))
radec_ring2_abcorr = hbt.correct_stellab(radec_ring2, vel_sun_nh_j2k) # radians

# Convert RA/Dec values back into pixels

x_stars_cat,    y_stars_cat      = w.wcs_world2pix(radec_stars[:,0]*r2d,   radec_stars[:,1]*r2d, 0)

x_stars_abcorr, y_stars_abcorr   = w.wcs_world2pix(radec_stars_abcorr[:,0]*r2d, radec_stars_abcorr[:,1]*r2d, 0)
x_ring1_abcorr, y_ring1_abcorr   = w.wcs_world2pix(radec_ring1_abcorr[:,0]*r2d, radec_ring1_abcorr[:,1]*r2d, 0)
x_ring2_abcorr, y_ring2_abcorr   = w.wcs_world2pix(radec_ring2_abcorr[:,0]*r2d, radec_ring2_abcorr[:,1]*r2d, 0)

points_stars        = np.transpose((x_stars_cat, y_stars_cat))
points_stars_abcorr = np.transpose((x_stars_abcorr, y_stars_abcorr))

# Read the image file from disk

image_polyfit = hbt.read_lorri(file, frac_clip = 1.,  
                             bg_method = 'Polynomial', bg_argument = 4)
image_raw     = hbt.read_lorri(file, frac_clip = 0.9, 
                             bg_method = 'None')

# Use DAOphot to search the image for stars. It works really well.

points_phot = hbt.find_stars(image_polyfit, num=50)

y_stars_dao =(points_phot[:,0]) # XXX Strangely, I had to swap x and y from how it is in nh_jring_gui to get to work...
x_stars_dao =(points_phot[:,1])

# Now look up the shift between the photometry and the star catalog. 
# Do this by making a pair of fake images, and then looking up image registration on them.
# I call this 'opnav'. It is returned in order (y,x) because that is what imreg_dft uses, even though it is a bit weird.
#
# For this, I can use either abcorr stars or normal stars -- whatever I am going to compute the offset from.        

#points_phot = np.array(x_stars_cat, y_stars_cat)

(dy_opnav, dx_opnav) = hbt.calc_offset_points(points_phot, points_stars, np.shape(image_raw), do_plot=True)

points_dao = np.transpose(np.array([x_stars_dao, y_stars_dao]))
points_cat = np.transpose(np.array([x_stars_cat, y_stars_cat]))

dy = dy_opnav
dx = dx_opnav

# Now make a plot!

hbt.figsize((10,10))

plt.imshow(stretch(image))

x_pos_ring1 = x_ring1 # Convert from string (which can go in table) to array
y_pos_ring1 = y_ring1
x_pos_ring2 = x_ring2
y_pos_ring2 = y_ring2          

# Get the user offset position
#
#dx = t['dx_opnav'] + self.slider_offset_dx.get()
#dy = t['dy_opnav'] + self.slider_offset_dy.get()

# Plot the stars -- catalog, and DAO

plt.plot(x_stars_cat, y_stars_cat, 
         marker='o', ls='None', 
         color='lightgreen', alpha = 0.5, ms=12, mew=1, label = 'Cat Stars, Raw')
         
plt.plot(x_stars_dao, y_stars_dao, 
         marker='o', ls='None', 
         color='pink', alpha = 0.5, label = 'DAOfind Stars')               

# Get position of satellites

name_bodies = np.array(['Metis', 'Adrastea', 'Thebe', 'Amalthea', 'Io'])        

x_bodies,  y_bodies   = hbt.get_pos_bodies(et, name_bodies, units='pixels', wcs=w)
ra_bodies, dec_bodies = hbt.get_pos_bodies(et, name_bodies, units='radec', wcs=w)

# Plot satellites
  
plt.plot(x_bodies+dx, y_bodies+dy, marker = '+', color='red', markersize=20, linestyle='none')

# Plot the ring
    
DO_PLOT_RING_INNER = False
DO_PLOT_RING_OUTER = True

if (DO_PLOT_RING_OUTER):
    plt.plot(x_pos_ring2, y_pos_ring2, marker='o', color = 'blue', ls = '--',
                  label = 'Ring, OpNav only')

    plt.plot(x_pos_ring2 + dx, y_pos_ring2 + dy, marker='o', color = 'lightblue', ls = '--',
                  label = 'Ring, OpNav+User')
        
if (DO_PLOT_RING_INNER):
    plt.plot(x_pos_ring1, y_pos_ring1, marker='o', color='green', ls = '-', \
        ms=8, label='Ring, LT')

    plt.plot(x_pos_ring1 + dx, y_pos_ring1 + dy, \
        marker='o', color='purple', ls = '-', ms=8, label='Ring, LT, Shifted')

plt.legend()  # Draw legend. Might be irrel since remove() might keep it; not sure.

plt.imshow(stretch(image))
plt.show

plt.imshow(stretch(image))
plt.show()



    
##### 

# Just a scratch spot to test star search, etc.


from photutils import daofind
from astropy.stats import sigma_clipped_stats
from photutils import find_peaks

num = 200 # Number of brightest objects to keep

image_s = hbt.remove_sfit(image,4)

mean, median, std = sigma_clipped_stats(image, sigma=3.0, iters=5)
    
sources = daofind(image, fwhm=2.0, threshold=2.*std)

threshold = median + (10.0 * std) # Ten sigma
tbl = find_peaks(image, threshold, box_size=5)

sources.sort('flux')  # Sort in-place
tbl.sort('peak_value')
    
if (num > 0):  
    index_start = -num
else:
    index_start = 0
        
x_phot = np.array(sources['xcentroid'][index_start:].data)
y_phot = np.array(sources['ycentroid'][index_start:].data)


x_phot_2 = np.array(tbl['x_peak'][index_start:].data)
y_phot_2 = np.array(tbl['y_peak'][index_start:].data)

# Make a plot

plt.imshow(stretch(image_s))
plt.plot(x_phot,   y_phot,   marker='o', alpha=0.5, mec='lightgreen', mew=3, ls = 'none', markersize=12, mfc='none',
         label = 'DAOfind')
plt.plot(x_phot_2, y_phot_2, marker='o', alpha=0.5, mec='lightblue',  mew=3, ls = 'none', markersize=12, mfc='none',
         label = 'find_peaks')

plt.plot(x_stars_cat, y_stars_cat, 
         marker='o', ls='None', mfc='none', 
         mec='red', alpha = 0.5, ms=12, mew=3, label = 'Cat Stars, Raw')

plt.legend(framealpha=0.9)
plt.xlim((0,1000))
plt.ylim((1000,0))

plt.show()

x_phot_merged = np.append(x_phot, x_phot_2)
y_phot_merged = np.append(y_phot, y_phot_2)

points_phot_merged = np.transpose(np.array([x_phot_merged, y_phot_merged]))

(dy_opnav, dx_opnav) = hbt.calc_offset_points(points_phot_merged, points_stars, np.shape(image_raw), do_plot=True)


