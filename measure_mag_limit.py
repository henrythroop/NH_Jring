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


from   astropy.stats import sigma_clipped_stats
from   photutils import daofind

import re # Regexp
import pickle # For load/save

import hbt

# Determine the stellar mag limit for an NH outbound rings image.

sequence = 'MVIC_D305'
sequence = 'MVIC_D211'
#sequence = 'MVIC_D202'
#sequence = 'LORRI_D202'
#sequence = 'LORRI_D305'

if (sequence == 'MVIC_D305'):
    file = '/Users/throop/Data/NH_MVIC_Ring/mvic_d305_sum_mos_v1_wcs_header.fits'
    ylim = (1,2500)
    fwhm = 4
    radius = 0.5 # Stellar search radius, in degrees
    dir_out = '/Users/throop/Data/NH_MVIC_Ring/out/'
    ylim_data = (2100,3100)
    exptime = 750 # Exposure time. Only used in figure caption. Total integrated time per pix, in seconds.
  
if (sequence == 'MVIC_D202'):
    file = '/Users/throop/Data/NH_MVIC_Ring/mvic_d202_mos_v1_wcs_header.fits'
    ylim = (1,2500)
    fwhm = 4
    radius = 0.5
    ylim_data = (2100,3100)
    exptime = 250

if (sequence == 'MVIC_D211'): # NB: there are a lot more artifacts here than the other sequences!
    file = '/Users/throop/Data/NH_MVIC_Ring/mvic_d211_mos_v1_wcs_header.fits'
    ylim = (1,2500)
    fwhm = 4
    radius = 0.5
    ylim_data = (600,1600)
    exptime = 250

if (sequence == 'LORRI_D202'): # NB: there are a lot more artifacts here than the other sequences!
    file = '/Users/throop/Data/NH_LORRI_Ring/ring_lor202_mos_v1_wcs_header.fits'
    ylim = (1,2500)
    fwhm = 4
    radius = 0.3
    ylim_data = (300,1300)
    exptime = 6

if (sequence == 'LORRI_D305'): # NB: there are a lot more artifacts here than the other sequences!
    file = '/Users/throop/Data/NH_LORRI_Ring/ring_lor305_mos_v3_wcs_header.fits'
    ylim = (1,2500)
    fwhm = 4
    radius = 0.3
    ylim_data = (300,1300)
    exptime = 18

    
fwhm = 4   # Works well for MVIC framing

hdulist = fits.open(file)

im = hdulist['PRIMARY'].data
w = WCS(file)                  # Look up the WCS coordinates for this frame
center  = w.wcs.crval  # degrees

et = hdulist['PRIMARY'].header['SPCSCET']

file_tm = '/Users/throop/gv/dev/gv_kernels_new_horizons.txt'

sp.furnsh(file_tm)

(st_sun,lt) = sp.spkezr('Sun', et, 'J2000', 'LT', 'New Horizons')
(st_pluto,lt) = sp.spkezr('Pluto', et, 'J2000', 'LT', 'New Horizons')

ang_sep = sp.vsep(st_sun[0:3], st_pluto[0:3])
ang_phase_deg = (math.pi - ang_sep)*hbt.r2d

# Generate the pickle filename

file_stars_pkl = dir_out + sequence + '_stars_deg' + repr(radius) + '.pkl'

#==============================================================================
# Get stars from a star catalog
#==============================================================================

try:
    lun = open(file_stars_pkl, 'rb')
    cat = pickle.load(lun)
    lun.close()
    print("Loaded: " + file_stars_pkl)
    
except FileNotFoundError:
    
    name_cat = u'Guide Star Catalog v2 1'
    #            name_cat = u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
    
    #            stars = conesearch.conesearch(w.wcs.crval, 0.3, cache=False, catalog_db = name_cat)
    
    with data.conf.set_temp('remote_timeout', 200): # This is the very strange syntax to set a timeout delay.
                                               # The default is 3 seconds, and that times out often.
        stars = conesearch.conesearch(w.wcs.crval, radius, cache=False, catalog_db = name_cat)
    
    # Extract proper fields
       
    ra_cat  = np.array(stars.array['ra'])*hbt.d2r # Convert to radians
    dec_cat = np.array(stars.array['dec'])*hbt.d2r # Convert to radians
    mag_cat       = np.array(stars.array['Mag'])
    
    # Sort the stars. I don't know why stars.sort('mag') doesn't work, but it fails.
    
    order = np.argsort(mag_cat)
    
    ra_cat = ra_cat[order]
    dec_cat = dec_cat[order]
    mag_cat = mag_cat[order]
    
    # Convert from RA Dec -> XY
    
    x_cat, y_cat    = w.wcs_world2pix(ra_cat*hbt.r2d, dec_cat*hbt.r2d, 0)
    
    # Stuff the results into an Astropy table
    
    cat = Table([ra_cat, dec_cat, x_cat, y_cat, mag_cat], names = ['RA', 'Dec', 'X', 'Y', 'Mag'])
    
    # Add two more fields to this table
    
    cat['dist_match'] = np.zeros(np.size(x_cat))-1
    cat['matched'] = np.zeros(np.size(x_cat), dtype=bool)
    
    print("Stars downloaded: {}; mag = {} .. {}".format(np.size(mag_cat), np.nanmin(mag_cat), np.nanmax(mag_cat)))
    print("RA = {} .. {}".format(np.nanmin(ra_cat)*hbt.r2d, np.nanmax(ra_cat)*hbt.r2d))

    lun = open(file_stars_pkl, 'wb')
    pickle.dump(cat, lun)
    lun.close()
    print("Wrote: " + file_stars_pkl)
        
#==============================================================================
# Do photometry on the NH image
#==============================================================================

# Use DAOphot to search the image for stars. It works really well.
# However, it returns position only -- no magnitudes.

#if (sequence == 'MVIC_D211'):
#    im = im[:,0:350]
    
mean, median, std = sigma_clipped_stats(im, sigma=3.0, iters=5)
sources = daofind(im - median, fwhm=4.0, threshold=2.*std)

x_phot = sources['xcentroid']
y_phot = sources['ycentroid']
flux   = sources['flux']

points_phot = np.transpose((y_phot, x_phot, flux)) # Create an array N x 2

# Sort them, bright to faint

order = (np.argsort(points_phot[:,2]))[::-1] # Sort from brightest to faintest
points_phot = points_phot[order]

#==============================================================================
# For each star in the catalog, go and see if there is a star we found with a center within, say, 2 pixels.
#==============================================================================

dist_match = 1  # Stars within this number of pixels are counted as a match
                # For LORRI, dist_match = 1 is fine
if ('211' in sequence):
    dist_match = 2
    
for i,star in enumerate(cat):
    x_i, y_i = star['X'], star['Y'] # Get position of this star
    dist_mutual = np.sqrt((x_i - x_phot)**2 + (y_i - y_phot)**2)
    amin = np.amin(dist_mutual)
    if (amin < dist_match):
#        print("Matched with distance {}".format(amin))
        cat['matched'][i] = True
        cat['dist_match'][i] = amin
        
#==============================================================================
# Now make a plot 
#==============================================================================

stretch_percent = 99

stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales to 5th..95th %ile. 

#hbt.figsize((30,70))

hbt.figsize((15,15))

plt.set_cmap('Greys_r')

# Plot the image itself.
               
plt.imshow(stretch(im))
            
# Plot the catalog stars

num_plot_cat = 2000

#x_cat, y_cat    = w.wcs_world2pix(cat['RA'][0:num_plot_cat]*hbt.r2d, cat['Dec'][0:num_plot_cat]*hbt.r2d, 0)

plt.plot(cat['X'][0:num_plot_cat], cat['Y'][0:num_plot_cat], 
         marker = 'o', ms=8, mew=0, ls='None', mec = 'none', mfc = 'lightgreen', alpha=0.3,
         label = 'GSC stars')

# Plot the DAOphot stars

num_plot_phot = -1  # -1 to plot them all

#for point in points_phot[0:num_plot_phot]:
#    plt.plot(point[1], point[0], marker = 'o', ms=10, mew=1, ls='None', mec='pink', mfc='None', alpha=0.5)

plt.plot(points_phot[0:num_plot_phot,1], points_phot[0:num_plot_phot,0], 
         marker = 'o', ms=10, mew=1, ls='None', mec='cyan', mfc='None', alpha=0.3,
         label = 'Photometric stars')

# Plot the catalog stars, which are flagged as a match

is_match = cat['matched'] == True

plt.plot(cat['X'][is_match[0:num_plot_cat]], cat['Y'][is_match[0:num_plot_cat]], 
         marker = 'o', ms=8, mew=1, ls='None', mec = 'red', mfc = 'none', alpha=1,
         label = 'Matched stars')


plt.ylim(ylim_data)
plt.legend(loc = 'lower right')
plt.title("{}, {:d} sec, {:.1f} deg".format(sequence, exptime, ang_phase_deg))
plt.xlabel('Column')
plt.ylabel('Row')

file_out = dir_out + sequence + '_maglimit_image.png'
plt.savefig(file_out)
print("Wrote: " + file_out)

plt.show()

#==============================================================================
# Make a histogram plot, showing fraction of mag10 matched, mag11 matched, etc.
#==============================================================================

# Plot a histogram of all stars

hbt.set_fontsize(15)
hbt.figsize((12,8))

num_bins = 40
xlim = (8,18)

(num,bins) = np.histogram(cat['Mag'], bins=num_bins, range = (5,20))
plt.plot(bins[0:-1], num, color = 'blue', label = 'GSC', marker = 'o')

(num,bins) = np.histogram(cat['Mag'][is_match], bins=num_bins, range = (5,20))
plt.plot(bins[0:-1], num, color = 'red', label='Detected', marker = 'o')
plt.yscale('log')
plt.xlabel('Magnitude')
plt.ylabel('Number')
plt.xlim(xlim)
plt.title('Stellar mag limit, {}, {:d} sec, {:.1f} deg'.format(sequence, exptime, ang_phase_deg))
plt.legend(loc='upper left')
plt.xticks(np.arange(xlim[0], xlim[1]+1, 1.0))

file_out = dir_out + sequence + '_maglimit.png'
plt.savefig(file_out)
print("Wrote: " + file_out)

plt.show()



    