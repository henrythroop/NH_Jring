#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 00:24:12 2019

@author: throop
"""

# Read the teacup images

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 11:43:11 2018

This program is to read the mosaic'ed FITS file that Tod Lauer has created of the teacup.

HBT March 2019

@author: throop
"""

import glob
import math
import os.path
import os
from astropy.table import Table
import astropy
from   astropy.io import fits
import astropy.table
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
import spiceypy as sp
from   astropy import units as u           # Units library
import pickle # For load/save

from   astropy.wcs import WCS

import scipy
import copy

# HBT imports

import hbt

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular
from   get_radial_profile_circular import get_radial_profile_circular_quadrant
from   get_radial_profile_backplane import get_radial_profile_backplane
from   get_radial_profile_backplane import get_radial_profile_backplane_quadrant
from   plot_img_wcs import plot_img_wcs
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes
from   scipy.optimize import curve_fit
from   wcs_translate_pix import wcs_translate_pix, wcs_zoom

# Load the proper kernels

file_tm = '/Users/throop/git/NH_rings/kernels_kem_prime.tm'
hbt.unload_kernels_all()
sp.furnsh(file_tm)

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

# First read the *raw* files. We read these just to get the header info (exptime, et, etc)

file_mosaic_teacup = '/Users/throop/Data/MU69/teacup/dpdeep_mos_tea_v2.fits'
dir = '/Users/throop/Data/MU69/teacup/'

files_raw = glob.glob(dir + 'mpf*.fit*')

ra_mu69_arr    = []
dec_mu69_arr   = []
ra_bsight_arr  = []
dec_bsight_arr = []
dist_arr       = []
et_arr         = []
vsep_arr       = []
file_short_arr = []
img_arr        = []
exptime_arr    = []

# Just load one file

files = [files_raw[0]]

# Loop over all the raw frames

for file in files:
    
    hdu = fits.open(file)
    header = hdu['PRIMARY'].header
    et = header['SPCSCET']
    
    img = hdu['PRIMARY'].data[0]
    hdu.close()
    
    file_short = os.path.basename(file)
    
    (st,lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    
    vec_sc_targ = st[0:3]
    
    dist = sp.vnorm(st)  # Convert the vector distance to MU69, in km

    (_, ra_mu69, dec_mu69) = sp.recrad(vec_sc_targ)
    
# Now compute the projected distance. I guess this is just the angular separation between the midpoint of 
# this frame, and MU69.
    

    # Set up WCS coord for this frame.

    # w = WCS.wcs(file)  # Load WCS -- but that does not exist for these individual MVIC frames!

    # OK, doesn't exist -- don't do that!
    # Neither Tod's mosaic nor the original FITS file has a WCS.
    
    # I could eye-ball it. I could create a SPICE kernel for an MU69 satellite (or ring) at 3000 km, in the proper
    # plane which is defined by the known MU69 pole position. I can plot that in GV. And then I can manually
    # tweak the ellipse equation to rotate and tilt an ellipse in python to match that?
    
    ra_bsight  = header['SPCBRRA']
    dec_bsight = header['SPCBRDEC']
    exptime    = header['EXPTIME']
    
    # Turn this into a vector
    
    vec_inst_bsight = sp.radrec(1, ra_bsight*hbt.d2r, dec_bsight*hbt.d2r)
    
    # Get the separation angle, in radians, from MU69
    
    vsep = sp.vsep(vec_inst_bsight, vec_sc_targ)

    # Save all the variables
        
    dec_bsight_arr.append(dec_bsight)
    ra_bsight_arr.append(ra_bsight)
    et_arr.append(et)
    vsep_arr.append(vsep)
    ra_mu69_arr.append(ra_mu69)
    dec_mu69_arr.append(dec_mu69)
    file_short_arr.append(file_short)
    dist_arr.append(dist)
    img_arr.append(img)
    exptime_arr.append(exptime)

# Print the UT of the first image
# I used this to get the roll angle orientation. Looks like we are rotated relative to the Sun. 
# Ugh this makes getting the RA/Dec of each pixel complicated.
    
ut = sp.et2utc(et_arr[0], 'C', 0)
print(f'First image taken at UT {ut}')
    
# Now put them all into a table
    
dec_bsight = np.array(dec_bsight_arr)
ra_bsight  = np.array(ra_bsight_arr)
dec_mu69   = np.array(dec_mu69_arr)*hbt.r2d
ra_mu69    = np.array(ra_mu69_arr)*hbt.r2d
vsep       = np.array(vsep_arr)
et         = np.array(et_arr)
file_short = np.array(file_short_arr)
dist       = np.array(dist_arr)
exptime    = np.array(exptime_arr)

header_raw = header
pixfov     = header_raw['PIXFOV'] # urad per pixel

# And for the image array, just stack them all next to each other in the Y direction

num_files = len(files)
img_mosaic = np.zeros((num_files * np.shape(img)[0], np.shape(img)[1]))

for i in range(len(files)):
    a_flat = hbt.remove_sfit( (img_arr[i]))
    a_flat_crop = a_flat[2:127,:]
    img_mosaic[i * 125 : i*125 + 125, :] = a_flat_crop
    
    print(f'Loaded array {i}')

# Show the image

hbt.figsize((12,10))    
plt.imshow(stretch(img_mosaic), origin='lower')

# =============================================================================
# Now load and process Tod's teacup mosaic.
# =============================================================================

# Load the FITS file

hdu = fits.open(file_mosaic_teacup)
header = hdu['PRIMARY'].header
img_lauer = hdu['PRIMARY'].data

pos_x_mu69 = 3400 # Position of MU69, from Tod
pos_y_mu69 = 450

pixscale_km = 1.7646756177340408  # We calculate this later, but for now, we define it as a constant.

# Convert all the 0.0's to NaN. Lauer makes them 0.

indices_bg = img_lauer == 0.
img_lauer[indices_bg] = np.nan

# Plot the teacup

plt.set_cmap('Greys_r')
plt.imshow(stretch(img_lauer),origin='lower')
plt.plot([pos_x_mu69],[pos_y_mu69], marker = 'o', color='red')
plt.title(f'{os.path.basename(file_mosaic_teacup)}, {pixscale_km:.3} km/pix')
plt.ylabel('y pix')
plt.xlabel('x pix')
plt.show()

# Plot the teacup, with an annulus overlaid at 500 pixels radius

plt.set_cmap('Greys_r')
plt.imshow(stretch(img_lauer),origin='lower')
plt.plot([pos_x_mu69],[pos_y_mu69], marker = 'o', color='none', mec='white', ms=130)
plt.plot([pos_x_mu69],[pos_y_mu69], marker = 'o', color='red')
plt.title(f'{os.path.basename(file_mosaic_teacup)}, {pixscale_km:.3} km/pix')
plt.ylabel('y pix')
plt.xlabel('x pix')
plt.show()

# Now, for each pixel, calc an RA / Dec for it.
# Then use that to calculate the expected impact parameter radius at MU69, given the plane.
# - Set up a plane, centered on MU69 and pointing at the RA / Dec of the pole (from Hal)
# - Calc RA / Dec for each pixel
# - Draw a vector to each pixel.
# - Calc intersection btwn that vector, andthe plane defined above.
# - Calc the distance from that point, to the center of MU69.
# - Keep track of that radius for each pixel.
# - Then take the radial profile from that.

# Take some radial profiles of it. 

binwidth_profile = 20
(radius_pix, dn_profile) = get_radial_profile_circular(img_lauer, pos=(pos_y_mu69, pos_x_mu69), 
                                                        width=binwidth_profile, method='median')

plt.plot(radius_pix, dn_profile)
plt.xlim((0,1000))
plt.show()

# Take a quadrant profile

# Set up a few different binwidths

binwidth_profile = [2,5,10,20]

binwidth_profile = [60]
radius_pix = {}
dn_profile_quad = {}

# And take radial profiles at each of them

for width in binwidth_profile:
    print(f'Generating profile for binwidth {width}...')
    (radius_pix[width], dn_profile_quad[width]) = get_radial_profile_circular_quadrant(img_lauer, 
                                                        pos=(pos_y_mu69, pos_x_mu69), 
                                                        width=width, method='median')

for width in binwidth_profile:

    hbt.figsize((18,6))
    hbt.fontsize(12)
    for i in range(4):
        plt.plot(radius_pix[width] * pixscale_km, 
                 dn_profile_quad[width][i,:],label=f'Quadrant {i+1}', lw=3, alpha=0.6)
    
    title = f'{os.path.basename(file_mosaic_teacup)}, binning = {width} pix = {(width*pixscale_km):.1f} km'
    plt.xlim((0,2000))
    plt.ylim((-0.5,0.5))
    plt.title(title)
    plt.ylabel('DN')
    plt.xlabel('Distance [km]')
    plt.legend()
    plt.show()

# For comparison, do the same on the raw frames.

# binwidth_profile = 2
# (radius_pix, dn_profile) = get_radial_profile_circular(img_mosaic, width=binwidth_profile)
# plt.plot(radius_pix, dn_profile)
# plt.xlim((0,1000))
# plt.ylim((-2, 1))
# plt.show()

# Calculate the pixel scale

pixscale = np.mean(dist_arr)

# Convert from DN to IoF:

# =============================================================================
# Convert the image from DN to I/F
# =============================================================================
    
    # Apply Hal's conversion formula from p. 7, to compute I/F and print it.

RSOLAR_MVIC  = 98313.172   #       /(DN/s)/(erg/cm^2/s/Ang/sr), Solar spectrum  # MVIC value in MVIC frame header

# RSOLAR_LORRI_1X1 = 221999.98  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
# RSOLAR_LORRI_4X4 = 3800640.0  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)

# Define the solar flux, from Hal's paper.
FSOLAR_LORRI  = 176.	     	    # We want to be sure to use LORRI value, not MVIC value!
                                 # X XXX 176 is LORRI value. What is MVIC value??
F_solar       = FSOLAR_LORRI # Flux from Hal's paper
RSOLAR        = RSOLAR_MVIC

km2au = 1 / (u.au/u.km).to('1')
    
# Calculate the MU69-Sun distance, in AU (or look it up).         

# exptime        = header['EXPTIME'] # ET for the final image stack

TEXP        = exptime[0]
(st,lt)   = sp.spkezr('MU69', et_arr[0], 'J2000', 'LT', 'New Horizons')
r_nh_mu69 = sp.vnorm(st[0:3]) * km2au # NH distance, in AU
(st,lt)   = sp.spkezr('MU69', et_arr[0], 'J2000', 'LT', 'Sun')
r_sun_mu69= sp.vnorm(st[0:3]) * km2au # NH distance, in AU
pixscale_km =  (r_nh_mu69/km2au) * pixfov / 1e6 # km per pix (assuming LORRI 4x4)

# TEXP        = header['exptime']  # Exposure time of the field frames. All are 29.967 sec.

# Now convert from DN to I/F

# I                = dn_profile / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
# profile_iof      = math.pi * I * r_sun_mu69**2 / F_solar # Equation from Hal's paper

profile_iof_quad = {}
 
for width in binwidth_profile:
    I                       = dn_profile_quad[width] / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    profile_iof_quad[width] = math.pi * I * r_sun_mu69**2 / F_solar # Equation from Hal's paper

hbt.fontsize(15)
for width in binwidth_profile:

    quadrants_to_plot = [1,2,3,4]
    # quadrants_to_plot = [1]
    
    for j in quadrants_to_plot:  # Plot all four quadrants
        plt.plot(radius_pix[width] * pixscale_km, profile_iof_quad[width][j-1,:], alpha=0.5, 
                 lw = 3, label=f'Quadrant {j}')
    
    title = f'{os.path.basename(file_mosaic_teacup)}, binning = {width} pix = {(width * pixscale_km):.1f} km'
    plt.ylim((-5e-6, 2e-5))
    # plt.xlim((0,800))
    plt.xlabel('Radius [km]')
    plt.ylabel('I/F')
    plt.title(title)
    plt.legend()
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
    plt.show()
    
    # img_superstack_mean_iof   = math.pi * I_mean   * r_sun_mu69**2 / F_solar # Equation from Hal's paper

t = Table(   [et, ra_mu69, dec_mu69, ra_bsight, dec_bsight, vsep, dist, file_short],
              names = ['ET', 'RA_mu69', 'Dec_mu69', 'RA_bsight', 'Dec_bsight', 'Angle', 'Range', 'File' ])

t['ET'].format = '.0f'
t['RA_mu69'].format = '.3f'        
t['Dec_mu69'].format = '.3f'        
t['RA_bsight'].format = '.3f'        
t['Dec_bsight'].format = '.3f'        
t['Angle'].format = '.4f'    
t['Range'].format = '.0f'
t.sort('ET')
t['#'] = range(len(et))
t.pprint(max_width=200) 

