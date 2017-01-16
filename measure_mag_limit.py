#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 10:36:23 2017

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
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
from astropy.utils import data
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
import wcsaxes
import time
from scipy.interpolate import griddata

import re # Regexp
import pickle # For load/save

import hbt

# Determine the stellar mag limit for an NH outbound rings image.

file = '/Users/throop/Data/NH_MVIC_Ring/mvic_d305_sum_mos_v1_wcs_header.fits'

hdulist = fits.open(file)

im = hdulist['PRIMARY'].data
w = WCS(file)                  # Look up the WCS coordinates for this frame
center  = w.wcs.crval  # degrees

#==============================================================================
# Get stars from a star catalog
#==============================================================================

name_cat = u'Guide Star Catalog v2 1'
#            name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating

#            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)

with data.conf.set_temp('remote_timeout', 200): # This is the very strange syntax to set a timeout delay.
                                           # The default is 3 seconds, and that times out often.
    stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)

# Extract proper fields
   
ra_cat  = np.array(stars.array['ra'])*hbt.d2r # Convert to radians
dec_cat = np.array(stars.array['dec'])*hbt.d2r # Convert to radians
mag_cat       = np.array(stars.array['Mag'])

# Sort the stars. I don't know why stars.sort('mag') doesn't work, but it fails.

order = np.argsort(mag_cat)

ra_cat = ra_cat[order]
dec_cat = dec_cat[order]
mag_cat = mag_cat[order]

print("Stars downloaded: {}; mag = {} .. {}".format(np.size(mag_cat), np.nanmin(mag_cat), np.nanmax(mag_cat)))
print("RA = {} .. {}".format(np.nanmin(ra_cat)*hbt.r2d, np.nanmax(ra_cat)*hbt.r2d))

#==============================================================================
# Do photometry on the NH image
#==============================================================================

# Use DAOphot to search the image for stars. It works really well.
# However, it returns position only -- no magnitudes.

points_phot = hbt.find_stars(im, do_flux=True)

order = (np.argsort(points_phot[:,2]))[::-1] # Sort from brightest to faintest
points_phot = points_phot[order]

#==============================================================================
# Now make a plot 
#==============================================================================

# First plot the data image

        
# Set the color map


#plt.rc('image', cmap='Greys_r')


# Render the main LORRI frame
# *** In order to get things to plot in main plot window, 
# use self.ax1.<command>, not plt.<command>
# ax1 is an instance of Axes, and I can call it with most other methods (legend(), imshow(), plot(), etc)

stretch_percent = 99

stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile. 

hbt.figsize((30,70))
hbt.figsize((15,15))

plt.set_cmap('Greys_r')

# Plot the image itself.
               
plt.imshow(stretch(im))
            
# Plot the catalog stars

num_plot_cat = 300

x_cat, y_cat    = w.wcs_world2pix(ra_cat[0:num_plot_cat]*hbt.r2d, dec_cat[0:num_plot_cat]*hbt.r2d, 0)
plt.plot(x_cat, y_cat, marker = 'o', ms=10, mew=0, ls='None', mec = 'none', mfc = 'lightgreen', alpha=0.5)

# Plot the DAOphot stars

num_plot_phot = -1

#for point in points_phot[0:num_plot_phot]:
#    plt.plot(point[1], point[0], marker = 'o', ms=10, mew=1, ls='None', mec='pink', mfc='None', alpha=0.5)

plt.plot(points_phot[0:num_plot_phot,1], points_phot[0:num_plot_phot,0], 
         marker = 'o', ms=10, mew=1, ls='None', mec='red', mfc='None', alpha=0.5)

plt.ylim((2000,3000))
plt.show()