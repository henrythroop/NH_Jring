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
# Look at D305 image. Fix it if necessary
#==============================================================================

# '-new-image' is added by astrometry.net
# '_fixed'     is added by me to indicate that I have fixed the header (that it, added SCET field)
#files = ['mvic_d305_sum_mos_v1.fits', 'ringdep_mos_v1.fits']

#file = dir + '/mvic_d305_sum_mos_v1.fits'

sequence        = 'D305'
#sequence        = 'A_RINGDEP_01'  # Departure imaging, closest MVIC image of the whole system

DO_FIX_FITS     = False
DO_ANALYZE      = True

#==============================================================================
# Initialize parameters for each of the possible MVIC mosaics we can analyze
#==============================================================================

if (sequence == 'D305'):
  file_wcs    = dir + '/mvic_d305_sum_mos_v1-new-image.fits' # Load the navigated image, with WCS  
  utc = '2015::305 00:00:00'  # Set the time 
  stretch = astropy.visualization.PercentileInterval(99)
  
if (sequence == 'A_RINGDEP_01'):
  file_wcs = dir + '/ringdep_mos_v1-new-image.fits'
  utc = '2015 Jul 15 18:50:00'  # Set the time.
  stretch = astropy.visualization.PercentileInterval(99.6)
 
file_fixed  = file_wcs.replace('.fits', '_fixed.fits')       # File with fixed FITS ET info
file_planes = file_fixed.replace('.fits', '_planes.pkl')# File for backplanes
#==============================================================================
# Repair / edit the FITS files if needed
#==============================================================================
    
if (DO_FIX_FITS):

# Add a missing header field and write out
    
    et = cspice.utc2et(utc)
    hdulist = fits.open(file_raw)
    im = hdulist['PRIMARY'].data    
    hdulist['PRIMARY'].header['SPCSCET'] =     (et, "[s past J2000] Spacecraft mid-obs time, TDB")  
    
    hdulist.writeto(file_fixed, clobber=True)
    print "Wrote new FITS file with fixed header: " + file_fixed

# Create the backplanes and write out

    planes = hbt.create_backplane(file_fixed, frame = 'IAU_PLUTO', name_target='Pluto', name_observer='New Horizons')
    
    lun = open(file_planes, 'wb')
    pickle.dump(planes, lun)
    lun.close()
    
    print "Wrote backplane file: " + file_planes    
    
#==============================================================================
# Load the image + backplane
#==============================================================================
    
file = file_fixed
hbt.figsize((6,28))

hdulist = fits.open(file)
im = hdulist['PRIMARY'].data
    
lun = open(file_planes, 'rb')
planes = pickle.load(lun)
lun.close()

radius = planes['Radius_eq']

et = hdulist['PRIMARY'].header['SPCSCET'] # Get the ET from the file
utc = cspice.et2utc(et, 'C', 0)

w = WCS(file)

#==============================================================================
# Make a plot of the image + backplanes
#==============================================================================

plt.set_cmap('Greys_r')

hbt.figsize((20,40))
plt.subplot(1,3,1)
plt.imshow(stretch(im))
plt.title(sequence)
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
# Make a plot of radial profile
#==============================================================================
    
#hbt.figsize((5,4))
#plt.plot(bins_radius / 1000, flux_arr)
#plt.xlabel('Distance from Pluto [1000 km]')
#if (sequence == 'A_RINGDEP_01'):
#    plt.ylim((-0.1,0.1))
#plt.show()
    
#==============================================================================
# Make a labled plot with Pluto, Charon, Nix, Hydra, etc.
#==============================================================================

hbt.figsize((12,112))
plt.imshow(stretch(im))

offset_x, offset_y = 100,100  # Offset PCNHSK labels by this many pixels

pos_body_pix = {}  # Create a new dictionary to save positions of each body

name_body = ['Pluto', 'Charon', 'Nix', 'Hydra', 'Styx', 'Kerberos']
for name_body_i in name_body:
    vec,lt = cspice.spkezr(name_body_i, et, 'J2000', 'LT', 'New Horizons')
    vec_sc_targ = vec[0:3]
    (junk,ra,dec) = cspice.recrad(vec_sc_targ) # Get the RA / Dec of the object
    x_pix, y_pix    = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0) # Convert to pixels
    plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
             markeredgecolor='red', markeredgewidth=2)
    plt.text(x_pix + offset_x, y_pix + offset_y, name_body_i[0], color='red', fontsize=15)

    pos_body_pix[name_body_i] = np.array([y_pix, x_pix])  # Save position of each body as a dictionary entry tuple
        
plt.title(sequence)
plt.xlim((0,np.shape(im)[1]))
plt.ylim((0,np.shape(im)[0]))
plt.show()

#plt.imshow(dist_pluto_pix)
#plt.xlim((0,np.shape(im)[1]))
#plt.ylim((0,np.shape(im)[0]))
#plt.title('Distance from Pluto')
#plt.show()
#
#plt.imshow(dist_charon_pix)
#plt.xlim((0,np.shape(im)[1]))
#plt.ylim((0,np.shape(im)[0]))
#plt.title('Distance from Charon')
#plt.show()

#==============================================================================
# Create new backplane set, with pixel distances from each body
#==============================================================================

dist_body_pix = {}

for name_body_i in name_body:
    
    (dist_body_pix_y,dist_body_pix_x) = np.meshgrid(\
                     np.array(range(np.shape(im)[1])) - int(pos_body_pix[name_body_i][1]),
                     np.array(range(np.shape(im)[0])) - int(pos_body_pix[name_body_i][0]))


    dist_body_pix[name_body_i] = np.sqrt(dist_body_pix_x**2 + dist_body_pix_y**2)

#==============================================================================
# Filter stars, CRs, and other outliers out of the image.
#==============================================================================

# Idea: look for pixels that are more that 2sigma away from mean of the image

std = np.std(im)
im_clean = im.copy()
is_outlier = im > (np.median(im) + 3*std)
im_clean[is_outlier] = np.median(im)

# Do it again

std = np.std(im_clean)
im_clean2 = im_clean.copy()
is_outlier = im_clean > (np.median(im_clean) + 3*std)
im_clean2[is_outlier] = np.median(im_clean)

#==============================================================================
# Measure the radial profile, from Pluto
#==============================================================================

nbins_radius = 100

bins_radius = hbt.frange(np.amin(radius), np.amax(radius), nbins_radius) # Set up radial bins, in km

dradius = bins_radius[1] - bins_radius[0]

bin_number_2d = np.copy(im)*0

flux_mean_arr   = np.zeros(nbins_radius)
flux_median_arr = np.zeros(nbins_radius)
flux_std_arr    = np.zeros(nbins_radius)
npix_arr        = np.zeros(nbins_radius)

flux_mean_clean_arr   = np.zeros(nbins_radius)
flux_median_clean_arr = np.zeros(nbins_radius)

for i in range(nbins_radius-1):
    is_good = np.array(radius > bins_radius[i]) \
            & np.array(radius < bins_radius[i+1]) \
            & (dist_body_pix['Pluto'] > 100) \
            & (dist_body_pix['Charon'] > 100) \
            & (im != 0)  # Mask out the zero pixels (e.g., between MVIC mosaic frames)
                         # Also need to mask out area near Pluto and near Charon
            
    bin_number_2d[is_good] = i          # For each pixel in the array, assign a bin number

    flux_mean_arr[i]         = np.mean(im[is_good])
    flux_mean_clean_arr[i]   = np.mean(im_clean2[is_good])
    flux_median_arr[i]       = np.median(im[is_good])
    flux_median_clean_arr[i] = np.median(im_clean2[is_good])
    flux_std_arr[i]          = np.std(im[is_good]) # stdev
    
    npix_arr[i] = np.sum(is_good) # Count of number of pixels in this radial bin


#==============================================================================
# Calc orbital distances (from Pluto) for each body
#==============================================================================

d_pluto_body = {}
mask_orbit = {}

if (sequence == 'D305'):
  d_radius = 4000  # Halfwidth to plot, of orbit, in km
if (sequence == 'A_RINGDEP01'):
    radius = 100
    
# Look up the distance using SPICE

for name_body_i in name_body:
  (vec_body,junk) = cspice.spkezr(name_body_i, et, 'IAU_PLUTO', 'LT', 'Pluto')
  d_pluto_body[name_body_i] = cspice.vnorm(vec_body[0:3])

# Generate a pixel mask showing the orbit of each body

  mask_orbit[name_body_i] = \
    np.array(radius > (d_pluto_body[name_body_i] - d_radius)) & \
    np.array(radius < (d_pluto_body[name_body_i] + d_radius))
    
#==============================================================================
# Make plots of radial profile
#==============================================================================

# First plot: Mean

hbt.figsize((10,10))

if (nbins_radius == 100):
    offset = 0.02
    ylim = (-0.1, 0.1)

if (nbins_radius == 1000):
    offset = 0.2
    ylim = ((-0.3, 0.3))
    
plt.subplot(2,1,1)
#plt.plot(bins_radius / 1000, flux_mean_arr, label='Mean')
plt.plot(bins_radius / 1000, flux_mean_clean_arr, label='Mean, clean')
plt.plot(bins_radius / 1000, flux_median_clean_arr + offset, label='Median, clean + offset')

plt.ylim(ylim)

plt.title(sequence + ', nbins = ' + repr(nbins_radius) + \
                     ', $\delta r$ = ' + hbt.trunc(dradius,0) + ' km')

plt.xlabel('Orbital Distance [1000 km]')
plt.ylabel('DN')
plt.legend(loc = 'bottom center')

# Plot lines for satellite orbits
for name_body_i in name_body:
  plt.vlines(d_pluto_body[name_body_i]/1000, 0,1, linestyle='--')
  plt.text(  d_pluto_body[name_body_i]/1000, 0.8, ' ' + name_body_i[0])
    
# Second plot: Median

plt.subplot(2,1,2)
plt.plot(bins_radius / 1000, flux_median_arr, label='Median + offset')
plt.plot(bins_radius / 1000, flux_median_clean_arr, label='Median Clean + offset')
plt.ylim((-0.1, 0.1))

plt.title(sequence + ', nbins = ' + repr(nbins_radius) + \
                     ', $\delta r$ = ' + hbt.trunc(dradius,0) + ' km')

# Plot lines for satellite orbits
if (sequence == 'A_RINGDEP01'):
    for name_body_i in name_body:
      plt.vlines(d_pluto_body[name_body_i]/1000, 0,1, linestyle='--')
      plt.text(  d_pluto_body[name_body_i]/1000, 0.8, ' ' + name_body_i[0])

if (sequence == 'D305'):
    rh = d_pluto_body['Hydra']  # hydra radii
    
    for rh_i in [10, 20, 30, 40, 50]:
      plt.vlines(rh_i * r_h, 0,1, linestyle='--')
      plt.text((rh_i + 2) * r_h, 0.08, ' ' + repr(rh_i) + ' $R_H$')
  
plt.legend()
plt.show()

#Make a masking array with orbit of Charon shown



#==============================================================================
# Make an image showing the satellite orbits, and labled with satellite names
#==============================================================================
if (sequence == 'A_RINGDEP_01'):    
hbt.figsize((15,15))

# Make a composite image showing data, superimposed with orbits

# Plot first subplot : Annotated
    
    plt.subplot(1,2,1)
    
    im_composite = im.copy()
    
    mask_composite = ((mask_orbit['Charon'])   & (dist_body_pix['Charon'] > 100)) + \
                     ((mask_orbit['Styx'])     & (dist_body_pix['Styx'] > 100)) + \
                     ((mask_orbit['Hydra']) & (dist_body_pix['Hydra'] > 100)) + \
                     ((mask_orbit['Kerberos']) & (dist_body_pix['Kerberos'] > 100)) + \
                     ((mask_orbit['Nix'])      & (dist_body_pix['Nix'] > 100))
                     
    im_composite[mask_composite] = 10
                     
    plt.imshow(stretch(im_composite))
    
    for name_body_i in name_body:
        y_pix, x_pix    = pos_body_pix[name_body_i]
        plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
                 markeredgecolor='red', markeredgewidth=2, markersize=30)
        plt.text(x_pix + offset_x, y_pix + offset_y, name_body_i[0], color='white', weight='bold', fontsize=15)
        
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))
    plt.title(sequence)
    
    # Plot second subplot: Un-annotated
    
    plt.subplot(1,2,2)
    plt.imshow(stretch(im))
    plt.xlim((0,np.shape(im)[1]))
    plt.ylim((0,np.shape(im)[0]))
    for name_body_i in name_body:
        y_pix, x_pix    = pos_body_pix[name_body_i]
        plt.plot(x_pix, y_pix, marker = 'o', linestyle='none', markerfacecolor='none',
                 markeredgecolor='red', markeredgewidth=2, markersize=30)
    plt.title(sequence)
    plt.show()

#####
    
if (sequence == 'D305'):
    hbt.figsize((50,50))
    im_composite = im.copy()
    
    mask_composite = (mask_orbit['Charon']) + (mask_orbit['Hydra'])

#                     (mask_orbit['Styx'])      
#                     (mask_orbit['Kerberos'])  
#                     (mask_orbit['Nix'])     
                     
    im_composite[mask_composite] = 10
    plt.imshow(stretch(im_composite))

    for name_body_i in ['Pluto', 'Hydra']:
        y_pix, x_pix    = pos_body_pix[name_body_i]

        plt.text(x_pix + offset_x, y_pix + offset_y, name_body_i[0], weight='bold', color='red', fontsize=12)
        
#==============================================================================
# XXX OLD CODE: Navigate the image using star catalog
#==============================================================================

#sequence = 'A_RINGDEP_01'
#file_raw = dir + '/ringdep_mos_v1.fits'
#
#hdulist = fits.open(file_raw)
#im = hdulist['PRIMARY'].data
#
#name_cat = 'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1' # works, but 1' errors; investigating
#crval = [91, 15]  # Array with RA, Dec of center position, in degrees
#
#stars     = conesearch.conesearch(crval, 5, cache=False, catalog_db = name_cat) # center [deg], radius [deg]
#table_stars = Table(stars.array.data)
#mask = table_stars['Pmag'] < 9
#table_stars_m = table_stars[mask]  
#
#ra_stars  = table_stars_m['RAJ2000']*hbt.d2r # Convert to radians
#dec_stars = table_stars_m['DEJ2000']*hbt.d2r # Convert to radians
#
#radec_stars        = np.transpose(np.array((ra_stars,dec_stars)))
#x_stars, y_stars   = w.wcs_world2pix(radec_stars[:,0]*hbt.r2d,   radec_stars[:,1]*hbt.r2d, 0)        
#points_stars        = np.transpose((y_stars, x_stars))  # 
#
#shape = (5000,2250)
#diam_kernel=5
#
#image_1 = hbt.image_from_list_points(points_stars, shape, diam_kernel) # star catalog image: working
#image_2 = hbt.image_from_list_points(points_phot,  shape, diam_kernel) # photometry: cut off
#        
#points_phot = hbt.find_stars(im) # Weirdly, find_stars does not return magnitudes -- only positions
#            
#(dy_opnav, dx_opnav) = hbt.calc_offset_points(points_phot, points_stars, np.shape(im), plot=False)
