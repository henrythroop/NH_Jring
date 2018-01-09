#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 22:00:03 2017

@author: throop
"""

#=============================================================================
#  This program just does a bunch of misc testing on the various images taken in Sep-2017.
#  Histograms, statistics, stacking, subtraction, navigation, etc. 
#  My goal was to try to find anything interesting. 
#  MU69 wasn't visible. These were images of the right field, so we could prep our exposures
#  for the real hazard search.
#=============================================================================

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
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)
import imreg_dft as ird                    # Image translation

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt

from image_stack import image_stack

#=============================================================================
# OK! Now time to do some data analysis.
#=============================================================================

stretch_percent = 95    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

dir = '/Users/throop/Dropbox/Data/NH_KEM_Hazard/Porter_Sep17/'

images = image_stack(dir)

# Extract the data table

t = images.t
data = images.data  # Dictionary with *all* the images

DO_PLOT_POSITIONS = False

if (DO_PLOT_POSITIONS):
    plt.plot(t['x_pix'], t['y_pix'], marker = '+', linestyle='none')
    plt.xlabel('MU69 X pos [pix]')
    plt.ylabel('MU69 Y pos [pix]')
    plt.show()
    
    plt.plot(t['ra'], t['dec'], marker = '+', linestyle = 'none')
    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    plt.show()


indices_sep17 = t['et'] > sp.utc2et('15 sep 2017')  # The positon of MU69 has changed a few pixels.
                                                    # We can't blindly co-add between sep and pre-sep
indices_jan17 = t['et'] < sp.utc2et('1 sep 2017')
                                                    
indices_rot0  = t['angle'] < 180   # One rotation angle
indices_rot90 = t['angle'] > 180   # The other rotation angle
indices_10sec = np.logical_and( t['exptime'] < 10, t['exptime'] > 5  )
indices_20sec = np.logical_and( t['exptime'] < 20, t['exptime'] > 10 )
indices_30sec = np.logical_and( t['exptime'] < 30, t['exptime'] > 20 )

indices_1x1 = t['naxis1'] == 1024
indices_4x4 = t['naxis1'] == 256

indices_30sec_rot0  = np.logical_and(indices_30sec, indices_rot0)  # 94
indices_30sec_rot90 = np.logical_and(indices_30sec, indices_rot90) # 0

indices_20sec_rot0  = np.logical_and(indices_20sec, indices_rot0)  # 93
indices_20sec_rot90 = np.logical_and(indices_20sec, indices_rot90) # 0

indices_10sec_rot0  = np.logical_and(indices_10sec, indices_rot0)  # 96
indices_10sec_rot90 = np.logical_and(indices_10sec, indices_rot90) # 48

indices_30sec_rot0_sep  = np.logical_and(indices_30sec_rot0, indices_sep17)  # 94
indices_30sec_rot90_sep = np.logical_and(indices_30sec_rot90, indices_sep17) # 0

indices_20sec_rot0_sep  = np.logical_and(indices_20sec_rot0, indices_sep17)  # 93
indices_20sec_rot90_sep = np.logical_and(indices_20sec_rot90, indices_sep17) # 0

indices_10sec_rot0_sep  = np.logical_and(indices_10sec_rot0, indices_sep17)  # 48
indices_10sec_rot90_sep = np.logical_and(indices_10sec_rot90, indices_sep17) # 0

indices_10sec_rot0_jan = np.logical_and(indices_10sec_rot0, indices_jan17) # 48

# =============================================================================
# Now set some indices, and stack them.
# =============================================================================

        
sequence = '10sec_rot0_sep'       # Which of these image sets defined above do we use?
indices  = eval('indices_' + sequence)   # Python trick: evaluate 

# Set the indices

images.set_indices(indices)

w = np.where(indices)[0]

# Perform the stack, and get the result back

arr_flat_median = images.flatten(method = 'median')
arr_flat_mean   = images.flatten(method = 'mean')
arr             = images.arr   # array, N x 1024 x 1024, with just the most recent stack

plt.imshow(stretch(arr_flat))
plt.title(sequence)
plt.show()

# Make a plot, showing difference between 

plt.plot(arr_flat_median[512,:], marker='.', ls='none', label = 'Summed Median')
plt.plot(arr_flat_mean[512,:],   marker='.', ls='none', label = 'Summed Mean')

plt.yscale('log')
plt.ylim((2,100))
#plt.plot(data[0,512,:], label = 'Individual', ls='none', marker='.', ms=0.5)
#plt.plot(data[10,512,:], label = 'Individual', ls='none', marker='.', ms=0.5)
#plt.plot(data[20,512,:], label = 'Individual', ls='none', marker='.', ms=0.5)

plt.legend()
plt.show()


# =============================================================================
# Now load three stacks, and compare them. For 10 sec, 20 sec, 30 sec
# =============================================================================

hist_single = []

hist_flattened = [] # Histogram from stacked frames
hist_single    = [] # Histogram from single frame 
n              = [] # Number of frames summed
arr_flattened  = []

sequences = ['30sec_rot0_sep', '20sec_rot0_sep', '10sec_rot0_sep']
dn_min = -50
dn_max = 200
bins = hbt.frange(dn_min, dn_max) # Get a list of all of the timestep bins, inclusive, for this file.
frames_max = 40  # Max frames in a stack. Some stacks have more or less -- we need to force to be the same, to compare.

colors = ['red', 'green', 'blue']

for sequence in sequences:
    
    indices  = eval('indices_' + sequence).copy()   # Python trick: evaluate a statement
    w = np.where(indices)[0]
    indices[w[frames_max]:] = False          # De-activate all frames above frames_max
    
# XXX need to limit this to N frames max, for consistency against all the stacks
    
    images.set_indices(indices)
    
# Perform the flattening, and get the result back

    arr_flat = images.flatten()
    
    arr = images.arr    # Get the raw stack, un-flattened. Only valid after calling flatten().
    
    (h, junk)  = np.histogram(arr_flat, bins)
    (h2, junk) = np.histogram(arr[0,:,:], bins)

    # Save the output from this stack calculation
    
    hist_flattened.append(h)   # Histogram of stacked images
    hist_single.append(h2)     # Histogram of 0th image in the stack
    n.append(np.sum(indices))  # Number of images in each stack
    arr_flattened.append(arr_flat) # The flattened array itself (so we can do stats on it)
    
# Make a nice plot comparing noise from the three sequences

hbt.figsize((15,10))
for i, sequence in enumerate(sequences):    
    plt.plot(bins[0:-1], hist_flattened[i], label="{}, N={}".format(sequence, n[i]), marker = '', color=colors[i],
             alpha = 0.5)
    plt.plot(bins[0:-1], hist_single[i],    label="{}, N={}".format(sequence, 1),    marker = '.', color=colors[i], 
             linestyle = 'none')
    
plt.yscale('log')
plt.legend()
plt.xlabel('DN per pixel')
plt.ylabel('# of pixels')
plt.title('Noise Histogram vs. Exposure Stacking, KEM Approach Images September 2017')          
plt.show()

# =============================================================================
# Compare noise: 10 x 30 sec  vs  30 x 10 sec. Also 30 x 30 sec
# =============================================================================

hist_single = []

hist_flattened = [] # Histogram from stacked frames
hist_single    = [] # Histogram from single frame 
n              = [] # Number of frames summed
arr_flattened  = []
exptime        = []

sequences = ['30sec_rot0_sep', '10sec_rot0_sep', '30sec_rot0_sep']
dn_min = -50
dn_max = 200
bins = hbt.frange(dn_min, dn_max) # Get a list of all of the timestep bins, inclusive, for this file.
frames_max = [10, 30, 30]  # Max frames in a stack.
scaling    = [3, 1, 1]

colors = ['red', 'green', 'blue']

for i,sequence in enumerate(sequences):
    
    indices  = eval('indices_' + sequence).copy()   # Python trick: evaluate a statement
    w = np.where(indices)[0]
    indices[w[frames_max[i]]:] = False          # De-activate all frames above frames_max
    
# XXX need to limit this to N frames max, for consistency against all the stacks
    
    images.set_indices(indices)
    
# Perform the flattening, and get the result back

    arr_flat = images.flatten()
    t = images.t
    
    arr = images.arr    # Get the raw stack, un-flattened. Only valid after calling flatten().
    
    (h, junk)  = np.histogram(arr_flat, bins)
    (h2, junk) = np.histogram(arr[0,:,:], bins)

    # Save the output from this stack calculation
    
    hist_flattened.append(h)   # Histogram of stacked images
    hist_single.append(h2)     # Histogram of 0th image in the stack
    n.append(np.sum(indices))  # Number of images in each stack
    arr_flattened.append(arr_flat) # The flattened array itself (so we can do stats on it)
    exptime.append(t['exptime'][indices][0])
    
# Make a plot comparing noise from the sequences
# Not really sure about the conversion to DN/sec here.    

hbt.figsize((15,10))
for i, sequence in enumerate(sequences):    
    plt.plot(bins[0:-1]/exptime[i], hist_flattened[i]*exptime[i], label="{}, N={}".format(sequence, n[i]), marker = '', 
             color=colors[i],
             alpha = 0.5)
    
plt.yscale('log')
plt.legend()
plt.xlabel('DN per pixel per sec')
plt.ylabel('# of pixels')
plt.title('Noise Histogram vs. Exposure Stacking, KEM Approach Images September 2017')          
plt.show()

# =============================================================================
# Now see how well we can subtract one stack from another. For this, must use median.
# =============================================================================

hist_single = []

hist_flattened = [] # Histogram from stacked frames
hist_single    = [] # Histogram from single frame 
n              = [] # Number of frames summed
arr_flattened  = []

sequences = ['30sec_rot0_sep', '20sec_rot0_sep', '10sec_rot0_jan']
bins = hbt.frange(0, 200) # Get a list of all of the timestep bins, inclusive, for this file.

for sequence in sequences:
    
    indices  = eval('indices_' + sequence)   # Python trick: evaluate 

    images.set_indices(indices)
    
# Perform the flattening, and get the result back

    arr_flat = images.flatten(method='median')
    
    arr = images.arr    # Get the raw stack, un-flattened. Only valid after calling flatten().
    
    (h, junk)  = np.histogram(arr_flat, bins)
    (h2, junk) = np.histogram(arr[0,:,:], bins)

    hist_flattened.append(h)
    hist_single.append(h2)
    n.append(np.sum(indices))
    arr_flattened.append(arr_flat)

# Normalize them to each other -- September
    
ratio = np.nanmean(arr_flattened[0]) / np.nanmean(arr_flattened[1])
arr_flattened[1] = arr_flattened[1] * ratio

# Normalize them to each other -- September:Jan
    
ratio = np.nanmean(arr_flattened[0]) / np.nanmean(arr_flattened[2])
arr_flattened[2] = arr_flattened[2] * ratio

# =============================================================================
# Now make a plot showing how well we can subtract one stack from another. To how many DN. Sep-Sep.
# =============================================================================

diff = (arr_flattened[1] - arr_flattened[0])
(h, junk) = np.histogram(diff, bins)
h_cum = np.cumsum(h)
h_cum_frac = h_cum / np.amax(h_cum)

hbt.figsize((12,10))

plt.plot(bins[0:-1], h_cum_frac*100)
plt.yscale('linear')
plt.title("{} - {}".format(sequences[0], sequences[1]))
plt.xlabel('Per-pixel change in $\Delta$DN after subtraction of two fields')
plt.ylabel('% of pixels with < $\Delta$DN')
plt.ylim((90,100))
plt.show()

plt.imshow(stretch(diff))
plt.title("{} - {}".format(sequences[0], sequences[1]))
plt.show()

plt.imshow(stretch(arr_flattened[0]))
plt.title('{}, N = {}'.format(sequences[0], n[0]))
plt.show()
    
plt.imshow(stretch(arr_flattened[1]))
plt.title('{}, N = {}'.format(sequences[1], n[1]))
plt.show()

# =============================================================================
# Make a plot showing the ring, as seen from K-22d
# =============================================================================


utc_ca = '2019 JAN 01 05:33:00'
utc_km22 = '2018 dec 12 00:00:00'

day = 86400
frame = 'J2000'
abcorr = 'LT+S'

ang_pix = 0.3 / 1024 * hbt.d2r # LORRI 1x1 pixels, in radians
et_ca = sp.utc2et(utc_ca)
et_obs = et_ca - 22*day

(vec,lt) = sp.spkezr('MU69', et_obs, frame, abcorr, 'New Horizons')
dist_km = sp.vnorm(vec[0:3]) # distance in km

scale_pix_km = ang_pix * dist_km

radius_km = np.array([3500, 10000])  # Ring radius in km

radius_pix = radius_km / dist_km / ang_pix
area_pix = math.pi * radius_pix**2  # Total area covered by each ring


# Plot the image (differential, for fun)

hbt.figsize((19,19))
plt.imshow(stretch(diff))
# Plot MU69

#plt.plot(images.x_pix_mean*4, images.y_pix_mean*4, marker='.', color='red')

# Now draw a circle at each of these radii

tau = 0.05

circle1=plt.Circle((images.x_pix_mean*4, images.y_pix_mean*4),radius_pix[1]*1,color='red', alpha=tau)
circle2=plt.Circle((images.x_pix_mean*4, images.y_pix_mean*4),radius_pix[0],color='green', alpha=tau)

plt.gcf().gca().add_artist(circle1)
plt.gcf().gca().add_artist(circle2)
#
plt.title("Simulated Image at K-22d, I/F = {}".format(tau))
    
plt.show()


# Plot as seen at k-22d

# Plot a ring: 3500 km

# Plot a ring: 10,000 km

# =============================================================================
# Now do again for Jan-Sep
# =============================================================================

tvec = ird.translation(arr_flattened[2], arr_flattened[0])['tvec']

diff = (arr_flattened[2] - np.roll(np.roll(arr_flattened[0], int(tvec[0]), axis=0), int(tvec[1]), axis=1))

(h, junk) = np.histogram(diff[300:,:], bins)  # Exclude the portion that is not overlapped, of course!
h_cum = np.cumsum(h)
h_cum_frac = h_cum / np.amax(h_cum)


plt.plot(bins[0:-1], h_cum_frac*100)
plt.yscale('linear')
plt.title("{} - {}".format(sequences[0], sequences[2]))
plt.xlabel('Per-pixel change in $\Delta$DN after subtraction of two fields')
plt.ylabel('% of pixels with < $\Delta$DN')
plt.ylim((90,100))
plt.show()

plt.imshow(stretch(diff))
plt.title("{} - {}".format(sequences[0], sequences[2]))
plt.show()

plt.imshow(stretch(arr_flattened[0]))
plt.title('{}, N = {}'.format(sequences[0], n[0]))
plt.show()
    
plt.imshow(stretch(arr_flattened[2]))
plt.title('{}, N = {}'.format(sequences[2], n[2]))
plt.show()


# Create the output numpy array

    
# Set up a two-component stretch. This is because the second one (Sinh, or Sqrt) requires the data to 
# already be in range 0 .. 1.

stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
stretch2 = astropy.visualization.SqrtStretch()

plt.imshow(stretch2(stretch(arr_flat)))

plt.plot(x_pix_mean*4, y_pix_mean*4, marker = 'o', linestyle='none', color='red')
plt.title('{}, N = {}, {:0.2f} deg, {:0.2f} sec, {}'.format(
          sequence,
          num_planes, 
          t['angle'][w[0]]-180,
          t['exptime'][w[0]], 
          t['utc'][w[0]]) )
plt.show()

# Plot the original (ie, one frame, non-summed)

plt.imshow(stretch2(stretch(arr[0,:,:])))

plt.title('N = {}, {:0.2f} deg, {:0.2f} sec, {}'.format(
          sequence,  
          1, 
          t['angle'][w[0]]-180,
          t['exptime'][w[0]], 
          t['utc'][w[0]]) )
plt.show() 

# Make linear plots

row = 512

hbt.figsize((20,8))
plt.plot((arr_flat[row, :]), 
                 marker = '.', label = 'Summed, N = {}'.format(num_planes),linestyle='-', ms=1)

plt.plot((arr[0,row,:]), marker = '+', label = 'Individual frame')
plt.title(sequence)
plt.yscale('log')
plt.ylim((10,200))
plt.legend()
plt.show()

# Make a histogram plot

plt.plot()

bins = hbt.frange(0, 200) # Get a list of all of the timestep bins, inclusive, for this file.

(hist_flat, junk) = np.histogram(arr_flat, bins)

(hist_single, junk) = np.histogram(arr[0,:,:], bins)

hbt.set_fontsize(16)
plt.plot(bins[0:-1], hist_single, linestyle = '-', marker = '.', ms=1, label = 'Single')
plt.plot(bins[0:-1], hist_flat, linestyle = '-', marker = '.', ms=1, 
         label = 'Flattened, N={}'.format(num_planes))
plt.yscale('log')
plt.legend()
plt.title(sequence)
plt.ylabel('# of bins')
plt.xlabel('DN value per pixel')           
plt.show()

# Make a plot of the PSF size vs. exptime

# Take some     
# First get their indices

# Then expand them by 4x4

# Then shift them as per WCS

# Then co-add, and/or median them

# =============================================================================
# Calculate how many DN MU69 should be at encounter (K-20d, etc.)
# Or alternatively, convert all of my DN values, to I/F values
# =============================================================================

        # Convert DN values in array, to I/F values
        
RSOLAR_LORRI_1X1 = 221999.98  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
RSOLAR_LORRI_4X4 = 3800640.0  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)

C = arr[0] # Get the DN values of the ring. Typical value is 1 DN.

# Define the solar flux, from Hal's paper.

FSOLAR_LORRI  = 176.	     	    # We want to be sure to use LORRI value, not MVIC value!
F_solar = FSOLAR_LORRI # Flux from Hal's paper

RSOLAR = RSOLAR_LORRI_4X4

# Calculate the MU69-Sun distance, in AU (or look it up). 

r_sun_mu69 = 43.28 # Solar distance in AU, from GV.

r_nh_mu69 = 0.16671 # NH distance at K-20d, from GV. 0.16 AU = 24.9 M km.

TEXP = t['exptime'][0]   # Look up exposure time, in seconds
TEXP = 30

I = C / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. All v similar, except for spectrum assumed.

# Apply Hal's conversion formula from p. 7, to compute I/F and print it.

IoF = math.pi * I * r_sun_mu69**2 / F_solar # Equation from Hal's paper

# Now convert to 'normal I/F'

# Define mu = cos(e), where e = emission angle, and e=0 is face-on.

e = 0 # Assume ring is face-on. The sunflower orbits are. 

mu = np.cos(e)

#        mu_2D = np.transpose(np.tile(mu, (self.num_bins_radius(),1)))

# Calculate the normal I/F

IoF_normal = 4 * mu * IoF

# Save this result into 'self,' replacing the DN values with the normal I/F values

#        self.profile_radius_arr = IoF_normal

# Change the units

#        self.profile_radius_units = 'Normal I/F'

# Do the same for azimuthal profile

TEXP_2D = np.transpose(np.tile(self.exptime_arr, (self.num_bins_azimuth(),1)))
mu_2D = np.transpose(np.tile(mu, (self.num_bins_azimuth(),1)))
self.profile_azimuth_arr = self.profile_azimuth_arr / TEXP_2D / RSOLAR * math.pi * r**2 / F_solar * 4 * mu_2D

self.profile_azimuth_units = 'Normal I/F'

# Calculate how many DN is MU69 at K-10d. It is unresolved, so it is a slightly different calculation than above.

# Caclulate how many DN is I/F = 1e-6 ring. This is large spatial extent, so it should be just the inverse of above.


IoF_normal = 1e-6
TEXP = 30
IoF = IoF_normal / (4 * mu)     
I = IoF / math.pi / r_sun_mu69**2 * F_solar
C = I * TEXP * RSOLAR  # C is the value in DN
print("I/F_normal = {}, TEXP = {} s -> DN = {}".format(IoF_normal, TEXP, C))


# =============================================================================
# If desired, do a one-time routine for the synthetic images:
#  extract the I/F and ring size from the filenames, and append that to the table.
# =============================================================================

DO_APPEND = False
if (DO_APPEND):

    iof_ring = np.zeros(num_images, dtype=float)
    size_ring = np.zeros(num_images, dtype='U30')

    for i in range(num_images):
        f = t['filename_short'][i]
        m = re.search('ring_(.*)_iof(.*)_K', f)  # Call regexp to parse it.
        iof_ring[i] = eval(m.group(2))
        size_ring[i] = m.group(1)
        
    t['size_ring'] = size_ring
    t['iof_ring']  = iof_ring
    images.t = t
    images.save()           # Save it back to disk
        
# =============================================================================
# Now go thru the synthetic ring images. 
#    Load and stack the synthetic implanted images.
#    Load and stack the original 'raw' frames
#    Difference them, and see if we can find a ring in there.
# =============================================================================

dir_porter = '/Users/throop/Dropbox/Data/NH_KEM_Hazard/Porter_Sep17/'
dir_synthetic = '/Users/throop/Dropbox/Data/NH_KEM_Hazard/synthetic/'

# Load the images into a table

images_raw = image_stack(dir_porter)
images_syn = image_stack(dir_synthetic)

data_raw = images_raw.data
data_syn = images_syn.data

t_raw = images_raw.t
t_syn = images_syn.t

num_images_raw = (np.shape(t_raw))[0]
num_images_syn = (np.shape(t_syn))[0]

# Create a bunch of possible image sets, based on various parameters

# Indices for 'raw' images
  
indices_sep17_raw = t_raw['et'] > sp.utc2et('15 sep 2017')  # The positon of MU69 has changed a few pixels.
                                                            # We can't blindly co-add between sep and pre-sep
indices_jan17_raw = t_raw['et'] < sp.utc2et('1 sep 2017')
                                                    
indices_rot0_raw  = t_raw['angle'] < 180   # One rotation angle
indices_rot90_raw = t_raw['angle'] > 180   # The other rotation angle
indices_10sec_raw = np.logical_and( t_raw['exptime'] < 10, t_raw['exptime'] > 5  )
indices_20sec_raw = np.logical_and( t_raw['exptime'] < 20, t_raw['exptime'] > 10 )

indices_30sec_raw = np.logical_and( t_raw['exptime'] < 30, t_raw['exptime'] > 20 )

indices_1x1_raw = t_raw['naxis1'] == 1024
indices_4x4_raw = t_raw['naxis1'] == 256

indices_30sec_4x4_raw = np.logical_and(indices_4x4_raw, indices_30sec_raw) # 94

# Indices for synthetic images

indices_ring_small_syn = t_syn['size_ring'] == 'small'
indices_ring_large_syn = t_syn['size_ring'] == 'large'

indices_iof_1em7_syn = t_syn['iof_ring'] == 1e-7
indices_iof_1em6_syn = t_syn['iof_ring'] == 1e-6
indices_iof_1em5_syn = t_syn['iof_ring'] == 1e-5
indices_iof_1em4_syn = t_syn['iof_ring'] == 1e-4

indices_small_1em7_syn = np.logical_and(indices_iof_1em7_syn, indices_ring_small_syn)
indices_small_1em6_syn = np.logical_and(indices_iof_1em7_syn, indices_ring_small_syn)
indices_small_1em5_syn = np.logical_and(indices_iof_1em7_syn, indices_ring_small_syn)
indices_small_1em4_syn = np.logical_and(indices_iof_1em7_syn, indices_ring_small_syn)
indices_large_1em7_syn = np.logical_and(indices_iof_1em7_syn, indices_ring_large_syn)
indices_large_1em6_syn = np.logical_and(indices_iof_1em6_syn, indices_ring_large_syn)
indices_large_1em5_syn = np.logical_and(indices_iof_1em5_syn, indices_ring_large_syn)
indices_large_1em4_syn = np.logical_and(indices_iof_1em4_syn, indices_ring_large_syn)

indices_raw = indices_30sec_4x4_raw.copy()   # 94 of 344
indices_syn = indices_small_1em5_syn.copy()  # 94 of 752

# Now take the first half of the synthetic indices, and the second half of the raw ones
# This is to assure that we are using different images for the two stacks! Otherwise, the results are trivial.

frames_max = int(np.sum(indices_raw) / 2)         # Total number of frames (94)

w = np.where(indices_raw)[0]
indices_raw[w[frames_max]:] = False          # De-activate all frames *below* frames_max

w = np.where(indices_syn)[0]
indices_syn[:w[frames_max]] = False          # De-activate all frames above frames_max

# Set the indices

images_raw.set_indices(indices_raw)
images_syn.set_indices(indices_syn)

# Do the flattening

arr_raw = images_raw.flatten()    
arr_syn = images_syn.flatten()

# Look up the I/F from one of the synthetic images (doesn't matter which one -- they will all be the same)

t_syn = images_syn.t
iof_ring = t_syn[indices_syn]['iof_ring'][0]

# The two flattened images need some offsetting. Do that.

shift = ird.translation(arr_raw, arr_syn)['tvec']
shift = np.round(shift).astype('int')
arr_syn_shift = np.roll(np.roll(arr_syn, shift[0], axis=0), shift[1], axis=1)

arr_diff = arr_syn_shift - arr_raw

pos = (images.y_pix_mean*4, images.x_pix_mean*4)

(dist_pix_1d, profile_1d_median) = get_radial_profile_circular(arr_diff, pos, method='median')
(dist_pix_1d, profile_1d_mean)   = get_radial_profile_circular(arr_diff, pos, method='mean')

# Make a plot of the radial profile

plt.plot(dist_pix_1d, profile_1d_median, label = 'Median')
plt.plot(dist_pix_1d, profile_1d_mean, label = 'Mean')

#    plt.plot(dist_pix_1d, profile_1d_mean,   label = 'Mean')
plt.xlabel('Distance [km]')
plt.ylabel('DN per pixel')
plt.title('I/F = {:.0e}'.format(iof_ring))
plt.legend()
plt.show()