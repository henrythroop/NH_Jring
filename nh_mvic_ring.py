# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 22:55:34 2016

@author: throop
"""

dir = '/Users/throop/Data/NH_MVIC_Ring/'
file = 'mvic_d305_sum_mos_v1.fits'

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
import astropy.visualization
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import cspice
import skimage
from   itertools   import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy     import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils   import daofind
import wcsaxes
import time
from scipy.interpolate import griddata

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

# Imports for Tk

import Tkinter
import ttk
import tkMessageBox
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

###

import hbt

dir  = '/Users/throop/Data/NH_MVIC_Ring' 

file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

cspice.furnsh(file_tm) # Start up SPICE

#==============================================================================
# Navigate and fix the A_RINGDEP_01 image as necessary
#==============================================================================

sequence = 'A_RINGDEP_01'
file_raw = dir + '/ringdep_mos_v1.fits'

hdulist = fits.open(file_raw)
im = hdulist['PRIMARY'].data

name_cat = 'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
crval = [91, 15]  # Array with RA, Dec of center position, in degrees

stars     = conesearch.conesearch(crval, 5, cache=False, catalog_db = name_cat) # center [deg], radius [deg]
table_stars = Table(stars.array.data)
mask = table_stars['Pmag'] < 9
table_stars_m = table_stars[mask]  

ra_stars  = table_stars_m['RAJ2000']*hbt.d2r # Convert to radians
dec_stars = table_stars_m['DEJ2000']*hbt.d2r # Convert to radians

radec_stars        = np.transpose(np.array((ra_stars,dec_stars)))
x_stars,        y_stars          = w.wcs_world2pix(radec_stars[:,0]*hbt.r2d,   radec_stars[:,1]*hbt.r2d, 0)        
points_stars        = np.transpose((y_stars, x_stars))

shape = (5000,700)
diam_kernel=20

image_1 = hbt.image_from_list_points(points_stars, shape, diam_kernel)
image_2 = hbt.image_from_list_points(points_phot,  shape, diam_kernel)
        
points_phot = find_stars(im) # Weirdly, find_stars does not return magnitudes -- only positions
            
(dy_opnav, dx_opnav) = calc_offset_points(points_phot, points_stars, np.shape(im), plot=False)

#==============================================================================
# Look at D305 image. Fix it if necessary
#==============================================================================

# '-new-image' is added by astrometry.net
# '_fixed'     is added by me to indicate that I have fixed the header (that it, added SCET field)
#files = ['mvic_d305_sum_mos_v1.fits', 'ringdep_mos_v1.fits']

#file = dir + '/mvic_d305_sum_mos_v1.fits'

sequence        = 'D305'
DO_FIX_D305     = False
DO_ANALYZE_D305 = True

file_wcs    = dir + '/mvic_d305_sum_mos_v1-new-image.fits' # Load the navigated image, with WCS    
file_fixed  = file_wcs.replace('.fits', '_fixed.fits')       # File with fixed FITS ET info
file_planes = file_fixed.replace('.fits', '_planes.pkl')# File for backplanes
    
if (DO_FIX_D305):

# Add a missing header field and write out
    
    utc = '2015::305 00:00:00'  # Set the time for D305 image.
    et = cspice.utc2et(utc)
    
    hdulist['PRIMARY'].header['SPCSCET'] =     (et, "[s past J2000] Spacecraft mid-obs time, TDB")  
    
    hdulist.writeto(file_fixed, clobber=True)
    print "Wrote new FITS file with fixed header: " + file_fixed

# Create the backplanes and write out

    planes = hbt.create_backplane(file_fixed, frame = 'IAU_PLUTO', name_target='Pluto', name_observer='New Horizons')
    
    lun = open(file_planes, 'wb')
    pickle.dump(planes, lun)
    lun.close()
    
    print "Wrote backplane file: " + file_planes    

if (DO_ANALYZE_D305):
    
# Load the original image and display it
    
file = file_wcs
hbt.figsize((6,28))

hdulist = fits.open(file)
im = hdulist['PRIMARY'].data
    
plt.set_cmap('Greys_r')
stretch = astropy.visualization.PercentileInterval(99)

lun = open(file_planes, 'rb')
planes = pickle.load(lun)
lun.close()

radius = planes['Radius_eq']

#==============================================================================
# Make a plot of the data
#==============================================================================

plt.subplot(1,3,1)
plt.imshow(stretch(im))
plt.title('D305')
plt.gca().get_xaxis().set_visible(False)

plt.subplot(1,3,2)
plt.imshow(planes['Radius_eq'], cmap='plasma')
plt.title('Radius')
plt.gca().get_xaxis().set_visible(False)

plt.subplot(1,3,3)
plt.imshow(planes['Longitude_eq'], cmap='plasma')
plt.title('Longit')
plt.gca().get_xaxis().set_visible(False)

plt.show()

#==============================================================================
# Measure the radial profile, from Pluto
#==============================================================================

nbins_radius = 100
bins_radius = hbt.frange(np.amin(radius), np.amax(radius), nbins_radius) # Set up radial bins, in km
#h = np.histogram()
flux_arr = np.zeros(nbins_radius)
for i in range(nbins_radius-1):
    is_good = np.array(radius > bins_radius[i]) & np.array(radius < bins_radius[i+1])
    flux_arr[i] = np.mean(im[is_good])

hbt.imsize((5,4))
plt.plot(bins_radius / 1000, flux_arr)
plt.xlabel('Distance from Pluto [1000 km]')
plt.show()
    
#==============================================================================
# Label Pluto, Charon, Nix, Hydra, etc.
#==============================================================================

hbt.figsize((12,112))
plt.imshow(stretch(im))

offset_x, offset_y = 100,100

name_body = ['Pluto', 'Charon', 'Nix', 'Hydra', 'Styx', 'Kerberos']
for name_body_i in name_body:
    vec,lt = cspice.spkezr(name_body_i, et, 'J2000', 'LT', 'New Horizons')
    vec_sc_targ = vec[0:3]
    (junk,ra,dec) = cspice.recrad(vec_sc_targ) # Get the RA / Dec of the object
    x_pix, y_pix    = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0) # Convert to pixels
    plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
             markeredgecolor='red', markeredgewidth=2)
    plt.text(x_pix + offset_x, y_pix + offset_y, name_body_i[0], color='red', fontsize=15)

plt.title(sequence)
plt.xlim((0,np.shape(im)[1]))
plt.ylim((0,np.shape(im)[0]))
plt.show()


#==============================================================================
# Look at the second file (larger)
#==============================================================================

#file = dir + files[1]
#
#plt.rcParams['figure.figsize'] = 16,156
#hdulist = fits.open(file)
#im = hdulist['PRIMARY'].data
#
#plt.set_cmap('Greys_r')
#stretch = astropy.visualization.PercentileInterval(99.8)
#plt.imshow(stretch(im))
#plt.show()

#==============================================================================
# Now generate a backplane
#==============================================================================

planes = hbt.create_backplane(file, frame = 'IAU_PLUTO', name_target='Pluto', name_observer='New Horizons')

file_out = dir_out + file_short    
file_out = file_out.replace('.fits', '_planes.pkl')
    
#header['SPCSCET'] =     (490257526.5058528, "[s past J2000] Spacecraft mid-obs time, TDB")   

#file_out = file.replace('.fits', '_fixed.fits')
#hdulist.writeto(file_out)

    # Write one variable to a file        

file_out = file_out.replace('.fit', '_planes.pkl')

lun = open(file_out, 'wb')
pickle.dump(plane, lun)
lun.close()

print "Wrote: " + file_out
print
    