#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:54:05 2016

@author: throop
"""

# -*- coding: utf-8 -*-
"""

@author: throop

"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 14:55:00 2016

# Short stand-alone program to debug image navigation. This program includes code taken from nh_jring_gui.pro .
# 
# Ideally, I will debug it here, and then write a standalone function, which nh_jring_gui will call.


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
from   astropy.utils import data

from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astroquery.vo_conesearch import conesearch
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
import time
from   scipy.interpolate import griddata

import re # Regexp
import pickle # For load/save

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

radius_search_deg = 0.15

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
# Get stars from star catalogs     
#==============================================================================

if (DO_GSC1):
    name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
    stars = conesearch.conesearch(w.wcs.crval, radius_search, cache=False, catalog_db = name_cat)
    ra_stars  = np.array(stars.array['RAJ2000'])*d2r # Convert to radians
    dec_stars = np.array(stars.array['DEJ2000'])*d2r # Convert to radians
#            table_stars = Table(stars.array.data)

if (DO_GSC2):
    name_cat = u'Guide Star Catalog v2 1'
    file_pickle = dir_out + file_short.replace('.fit', '') + '.stars_gsc.pkl'

#    # If there is already a saved pickle file, then load from disk
#    
#    if os.path.isfile(file_pickle):
#        
#        print("Loading file: " + file_pickle)
#        lun = open(file_pickle, 'rb')
#        (ra_stars, dec_stars) = pickle.load(lun)
#        lun.close()
#            
#    else:     
#            name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors;investigating

#            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
        
    with data.conf.set_temp('remote_timeout', 30): # This is the very strange syntax to set a timeout delay.
                                                   # The default is 3 seconds, and that times out often.
        stars = conesearch.conesearch(w.wcs.crval, radius_search_deg, cache=True, catalog_db = name_cat)

    ra_stars  = np.array(stars.array['ra'])*d2r # Convert to radians
    dec_stars = np.array(stars.array['dec'])*d2r # Convert to radians

    mag       = np.array(stars.array['Mag'])
    
    print("Stars downloaded: N = {}; mag = {:.2f} .. {:.2f}".format(np.size(mag), np.nanmin(mag), np.nanmax(mag)))
    print("RA = {:.2f} .. {:.2f}".format(np.nanmin(ra_stars)*r2d, np.nanmax(ra_stars)*r2d))
    
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

ra_stars_cat = ra_stars
dec_stars_cat = dec_stars

#==============================================================================
# Calculate the ring location
#==============================================================================

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

radec_stars_cat        = np.transpose(np.array((ra_stars_cat, dec_stars_cat)))
radec_stars_cat_abcorr = hbt.correct_stellab(radec_stars_cat, vel_sun_nh_j2k) # Store as radians

# Convert ring RA/Dec for stellar aberration

radec_ring1        = np.transpose(np.array((ra_ring1,dec_ring1)))
radec_ring1_abcorr = hbt.correct_stellab(radec_ring1, vel_sun_nh_j2k) # radians
radec_ring2        = np.transpose(np.array((ra_ring2,dec_ring2)))
radec_ring2_abcorr = hbt.correct_stellab(radec_ring2, vel_sun_nh_j2k) # radians

# Convert RA/Dec values back into pixels

x_stars_cat,    y_stars_cat      = w.wcs_world2pix(radec_stars_cat[:,0]*r2d,   radec_stars_cat[:,1]*r2d, 0)

x_stars_cat_abcorr, y_stars_cat_abcorr   = w.wcs_world2pix(radec_stars_cat_abcorr[:,0]*r2d, 
                                                           radec_stars_cat_abcorr[:,1]*r2d, 0)
x_ring1_abcorr, y_ring1_abcorr   = w.wcs_world2pix(radec_ring1_abcorr[:,0]*r2d, radec_ring1_abcorr[:,1]*r2d, 0)
x_ring2_abcorr, y_ring2_abcorr   = w.wcs_world2pix(radec_ring2_abcorr[:,0]*r2d, radec_ring2_abcorr[:,1]*r2d, 0)

points_stars_cat        = np.transpose((y_stars_cat,        x_stars_cat))  # Yes, order is supposed to be (y,x)
points_stars_cat_abcorr = np.transpose((y_stars_cat_abcorr, x_stars_cat_abcorr))

# Read the image file from disk

image_polyfit = hbt.read_lorri(file, frac_clip = 1.,  
                             bg_method = 'Polynomial', bg_argument = 4)
image_raw     = hbt.read_lorri(file, frac_clip = 0.9, 
                             bg_method = 'None')

#==============================================================================
# Use DAOphot to search the image for stars.
#==============================================================================

points_stars_phot = hbt.find_stars(image_polyfit, num=50) # Returns N x 2 aray. 0 = Row = y; 1 = Column = x.

y_stars_phot =(points_stars_phot[:,0]) # xy is correct -- see above
x_stars_phot =(points_stars_phot[:,1]) # 

#==============================================================================
# Make a plot showing the DAO stars on the image
#==============================================================================

color_usno = 'red'

DO_PLOT_DAO = True

if (DO_PLOT_DAO):
    plt.imshow(stretch(image_polyfit))
    plt.plot(x_stars_phot, y_stars_phot, linestyle='none', 
             marker='o', markersize=9, mec=color_usno, mew=1, color='none', 
             label = 'DAO photometric stars') # plot() uses x, y 
    plt.ylim((1024,0))
    plt.xlim((0,1024))
    plt.legend(loc = 'upper left')
    plt.show()
    
#==============================================================================
# Look up the shift between the photometry and the star catalog. 
# Do this by making a pair of fake images, and then looking up image registration on them.
#==============================================================================

# I call this pointing process 'opnav'. 
# It is returned in order (y,x) because that is what imreg_dft uses, even though it is a bit weird.
#
# For this, I can use either abcorr stars or normal stars -- whatever I am going to compute the offset from.        

#points_phot = np.array(x_stars_cat, y_stars_cat)

(dy_opnav, dx_opnav) = hbt.calc_offset_points(points_stars_phot, points_stars_cat, np.shape(image_raw),
    labels = ['DAO', 'GSC Catalog'], do_plot_before=True, do_plot_after=True)

points_dao = np.transpose(np.array([x_stars_dao, y_stars_dao])) # Make a list of the measured stars
points_cat = np.transpose(np.array([x_stars_cat, y_stars_cat])) # Make a list of the catalog stars

dy = dy_opnav
dx = dx_opnav

#==============================================================================
# Make a plot, showing DAO positions, catalog positions, computed ring positions, and more
#==============================================================================

hbt.figsize((10,10))

plt.imshow(stretch(image))

x_pos_ring1 = x_ring1 # Convert from string (which can go in table) to array
y_pos_ring1 = y_ring1
x_pos_ring2 = x_ring2
y_pos_ring2 = y_ring2          

# Plot the stars -- catalog, and DAO

plt.plot(x_stars_cat + dx_opnav, y_stars_cat + dy_opnav, 
         marker='o', ls='None', 
         color='lightgreen', alpha = 0.5, ms=12, mew=1, label = 'Cat Stars, adjusted')

plt.plot(x_stars_cat, y_stars_cat, 
         marker='o', ls='None', 
         color='lightgreen', alpha = 1, ms=4, mew=1, label = 'Cat Stars, raw')
         
plt.plot(x_stars_phot, y_stars_phot, 
         marker='o', ls='None', 
         color='none', markersize=10, mew=1, mec='red', alpha = 1, label = 'DAOfind Stars')               

# Get position of satellites

name_bodies = np.array(['Metis', 'Adrastea', 'Thebe', 'Amalthea', 'Io'])        

x_bodies,  y_bodies   = hbt.get_pos_bodies(et, name_bodies, units='pixels', wcs=w)
ra_bodies, dec_bodies = hbt.get_pos_bodies(et, name_bodies, units='radec', wcs=w)

# Plot satellites

DO_PLOT_SATELLITES = False

if (DO_PLOT_SATELLITES):
    plt.plot(x_bodies+dx, y_bodies+dy, marker = '+', color='red', markersize=20, linestyle='none')

# Plot the ring
    
DO_PLOT_RING_INNER = False
DO_PLOT_RING_OUTER = False

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

plt.title('After navigation, with dx = {:.1f}, dy = {:.1f}'.format(dx_opnav, dy_opnav))
plt.legend()  # Draw legend. Might be irrel since remove() might keep it; not sure.

plt.imshow(stretch(image))
plt.show

    
##### 

# Just a scratch spot to test star search, etc.

from photutils import daofind
from astropy.stats import sigma_clipped_stats
from photutils import find_peaks

DO_TEST = False

if (DO_TEST):
    
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
    
    #(dy_opnav, dx_opnav) = hbt.calc_offset_points(points_phot_merged, points_stars, np.shape(image_raw), do_plot=True)


