# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:58:26 2016

READ_ALICE_RING_OCC
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
import spiceypy as sp # was cspice
import skimage
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils import daofind
import wcsaxes
import time
import sys  # For stdout.write, without newline
from scipy.interpolate import griddata

from mpl_toolkits.axes_grid1 import host_subplot # For adding a second axis to a plot
import mpl_toolkits.axisartist as AA             # For adding a second axis to a plot

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

# Imports for Tk

from astropy import units as u
from astropy.coordinates import SkyCoord

# HBT imports
import hbt

def read_alice_occ_data(file_list, xlim, ylim, verbose=True, short=False, DO_MASK_LYA=False):

#==============================================================================
# Read the Alice data
#==============================================================================

# The deal here is that we can read the data in several ways.
#  o COUNT_RATE.data .  This gives total rate, for each timestep, integrated across the detector.
#  o PRIMARY.data .     This gives an image, for each file. It is integrated in time (~7 sec per file)
#  o PIXEL_LIST_TABLE . This is the 'raw' data, with a time-tagged list of photons.
#                       Need to use this, if we want to get high-time-cadence photon counts for individual rows/columns.
#  o DO_MASK_LYA = if requested, mask out the LyA data (generally columns 166:186])
        
    met_all = []               # A long array with a list of all of the timestamps, one per 4 ms (i.e., at 250 hz)
    count_rate_fits_all = []   # The count rate as read from the COUNT_RATE extension directly
    count_rate_all = []        # Count rate computed from the PIXEL_LIST_TABLE. Should match COUNT_RATE extension
    count_rate_target_all = [] # Count rate for the target only, extracted by spatially filtering the PIXEL_LIST_TABLE
 
    dn_bg = 0                  # Background DN level. For filling in empty data.

    bin_lya = [536,556]        # Column range for LyA radiation, in bins.
    
    # O_RING_OC2, O_RING_OC3
    
    if (True):
        print("read_alice_occ_data: ylim = " + repr(ylim))
        print("DO_MASK_LYA = {}".format(DO_MASK_LYA))
        
        image_target_summed = np.zeros((ylim[1]-ylim[0],xlim[1]-xlim[0])) 
                                            # The 2D extracted region of the full spectral-spatial array.
                                            # This is just the region we care about, around the star.
        image_summed        = np.zeros((32, 1024))  # A full image

        # For debugging purposes: if requested, extract just a short subset of the data
        
        file_list_use = file_list
        if (short):
            file_list_use = file_list_use[0:100]

        for i, file in enumerate(file_list_use):
            if (verbose): sys.stdout.write('.')
            
            hdulist = fits.open(file)
            image = hdulist['PRIMARY'].data # Units of this are float, but I'm not sure what they are. 
                                        # I prefer raw counts.
                                        # d is 32 x 1024 
            
            if (DO_MASK_LYA):           # If requested, mask out the LyA right at the start
                                        # The columns are hardcoded -- so be careful if we ever change them.
                image[:,bin_lya[0]:bin_lya[1]] = dn_bg
                                                  
            # First go into the 'image', and sum photons at the right wavelength and spatial position
            
            image_target = image[ylim[0]:ylim[1], xlim[0]:xlim[1]] # Extract the right rows and columns, as per Steffl.
                                         # First range is for the position (just a few rows)
                                         # Second range is for wavelength  (broad, many columns)
                                         
            image_target_summed += image_target  # A 2D array of the extracted region
    
            image_summed        += image
             
            # Now go through and process the FITS file in a second way, by summing all the photons in the list.
            # This should in theory give us the same results as above, if we spatially integrate.
            # Also, it gives us the ability to *not* spatially integrate, and extract individual rows/columns.
            
            p = hdulist['PIXEL_LIST_TABLE'].data           # Array of length 29460 (for first file)
            count_rate_fits_i = hdulist['COUNT_RATE'].data # Array of length 3308  (for first file)
            
            num_samples = hdulist['COUNT_RATE'].header['NAXIS1'] # Number of samples in this file
            dt          = hdulist['COUNT_RATE'].header['SAMPLINT']  # Count rate sampling interval [sec]
            
            bins = hbt.frange(0, num_samples) # Get a list of all of the timestep bins, inclusive, for this file.
                                                # Use '+1' so we create the histogram upper size bin.
            
            # Now downselect the pixel list for just the photons in the proper X and Y position on the detector
            
            is_target = (p['Y_INDEX'] < ylim[1]) & (p['Y_INDEX'] >= ylim[0]) & (p['X_INDEX'] >= xlim[0]) & \
                        (p['X_INDEX'] < xlim[1])
        
            
            is_lya = (p['X_INDEX'] >= bin_lya[0]) & (p['X_INDEX'] <= bin_lya[1]) # This is an array of booleans.
            # Not a numpy array.
            
            # Flag each pixel as good or bad, based on its position, and based on the LyA mask flag
            
            if (DO_MASK_LYA):
                is_good = is_target & (is_lya == False)   # Strangely, have to do ==False, rather than not().
                timesteps_all  = p['TIMESTEP'][(is_lya == False)]

            else:
                is_good = is_target
                timesteps_all  = p['TIMESTEP']
           
            # Now we have a list of all of the good pixels. For each of these, now we want to grab its timestep.
        
            timesteps_good = p['TIMESTEP'][is_good]
        
        # Now count how many photons are in each timestep bin. I have defined those timestep bins up above.
        
            (count_rate_target_i, junk) = np.histogram(timesteps_good, bins)
            (count_rate_i, junk)        = np.histogram(timesteps_all,  bins) 
            met_i = hdulist['PRIMARY'].header['STARTMET'] + dt * np.array(range(num_samples))
                  
        #    print "File " + os.path.basename(file) + ' : ' + repr(np.sum(count_rate_target_i)) + ' / ' + 
        #       repr(np.sum(count_rate_i))
        
        # Now append these into the output lists
        
            count_rate_fits_all.append(count_rate_fits_i)
            count_rate_all.append(count_rate_i)
            count_rate_target_all.append(count_rate_target_i)
            met_all.append(met_i)

            hdulist.close()

        count_rate_fits  = np.array([item for sublist in count_rate_fits_all for item in sublist])  # Flatten count rate 
                                                                                                    # (from 2D, to 1D)
        count_rate_fits  = np.array(count_rate_fits, dtype=float)					           # Convert to float - 
                                                                                                    # Otherwise,
                                                                                                    # get wraparound.
        
        count_rate  = np.array([item for sublist in count_rate_all for item in sublist])
        count_rate  = np.array(count_rate, dtype=float)					  
        
        count_rate_target  = np.array([item for sublist in count_rate_target_all for item in sublist]) 
        count_rate_target  = np.array(count_rate_target, dtype=float)					  
        
        met         = np.array([item for sublist in met_all for item in sublist])   # Flatten the MET array 2D to 1D
        met         = np.array(met, dtype=float)
    
    return (met, count_rate_target, count_rate, image_target_summed, image_summed)

 
#==============================================================================
# Start of main program
#==============================================================================

##########
# Select which of the sequences we want to read in
#
# NB: The spelling / capitalization is inconsistent in SAPNAME vs. VISITNAM. I have standardized it here.
##########

#sequence 	= 'O_RING_OC3'
#sequence 	= 'O_RING_OC2'
sequence   = 'OCCSTAR1'

DO_ABBREVIATED = False       # For debugging. Just use a subset of the data?
DO_HICAD       = False
DO_MASK_LYA    = True        # Mask out the LyA lines in the Alice data?
#DO_LOAD_PICKLE = True        # Load pre-processed data from Pickle file?

dir_out = '/Users/throop/Data/NH_Alice_Ring/out/' # For saving plots, etc.

file_pickle = dir_out + sequence + '.pkl'

#if (DO_LOAD_PICKLE) and os.path.isfile(file_pickle):
    
binning      = 25000		# Smoothing. 25000 is too much (shows opposite trend!). 5000 and 1000 look roughly similar.
                            # To be meaningful, the binning timescale must be less than the deadband 
                            #  timescale (~20-30 sec RT).
                            # At binning=3000, it's 12 sec... so that is the longest we'd really want to go.
                            #
                            # Use binning = 25,000 for the data analysis.
                            # Use binning = 3000 to make a plot of the trend of DN vs. RA
                            
fs           	= 15		# Font size
plt.rc('image', interpolation='None')       # Turn of interpolation for imshow. Could also do hbt.set_plot_defaults()
hbt.set_fontsize(fs)        # Set global matplotlib font size
plt.set_cmap('plasma')

# Figure out the directory, stellar positions, etc. based on which sequence we are using

if (sequence == 'O_RING_OC2') or (sequence == 'O_RING_OC3'):
    xlim = np.array([370,910])  # Wavelength
    ylim = np.array([15,18])    # Spatial of the star. Was [13,19]

if (sequence == 'O_RING_OC2'): 
    DO_HICAD = False             # Used hicad since that's where the data was when I grabbed it
    
if (sequence == 'OCCSTAR1'):
    DO_HICAD = False
    xlim = np.array([370,910]) # Wavelength 
    ylim = np.array([10,13])    # Spatial, star #1, brighter, v = 4.9, and closer to lollipop
    ylim_2 = np.array([7,10])   # Spatial, star #2, fainter,  v = 5.5, and further toward bottom of stick.
                                # NB: [10,13] means three rows (10,11,12) -- not four.

    ylim_12 = np.array([6,14])  # Include both targets in one.
    
if not DO_HICAD:
    dir_images = '/Users/throop/Data/NH_Alice_Ring/' + sequence + '/data/pluto/level2/ali/all'
else:
    dir_images = '/Users/throop/Data/NH_Alice_Ring/' + sequence + '/data/pluto/level2/ali_hicad'
    
file_tm    = '/Users/throop/gv/dev/kernels/15sci_rhr_25aug16.tm' # Def need to use an encounter kernel here.

sp.furnsh(file_tm)

#########
# Read the Alice data
#########

file_list = glob.glob(dir_images + '/*fit')

print("Found {} Alice data files for sequence {} in {}".format(np.size(file_list), sequence, dir_images))
print("Reading...")

# Read the entire Alice sequence. Extract the total count rate (across frame),
# and count rate within a box defined by xlim, ylim.

(met, count_rate_target, count_rate, image_target_summed, image_summed) = \
  read_alice_occ_data(file_list, xlim, ylim, verbose=True, short=DO_ABBREVIATED,DO_MASK_LYA=DO_MASK_LYA)

# If there are two stars, then also read count rate of second one.
# Programmatically this is a bit inefficient since I read the entire dataset 
# a second time. Function could be rewritten to extract two stars at once, 
# but that is a lot of work.

if (sequence == 'OCCSTAR1'):
    print()
    print('Reading star 2...')
    (met, count_rate_target_2, count_rate_2, image_target_summed_2, image_summed) = \
        read_alice_occ_data(file_list, xlim, ylim_2, verbose=True, short=DO_ABBREVIATED, DO_MASK_LYA=DO_MASK_LYA)

# Look up the sub-obs latitude. This is for the ring tilt angle. 
# Do this for just the first timestep, since it doesn't change much during the scan.

vec,lt = sp.spkezr('New Horizons', sp.utc2et(hbt.met2utc(met[0])), 'IAU_PLUTO', 'LT', 'Pluto')
(junk,subobslon,subobslat) = sp.reclat(vec[0:3])    # Return sub-obs latitude (aka ring tilt angle), in radians

dt = int((met[1] - met[0])*1000)/1000.          # Interval between consecutive samples

# Make a quick plot of the full image and the target

stretch = astropy.visualization.PercentileInterval(99)
 
hbt.figsize((10,5))
plt.imshow(stretch(image_summed), aspect=10)
plt.title(sequence + ', full')
plt.show()

plt.imshow(stretch(image_target_summed), aspect=10)
plt.title(sequence + ', target')
plt.show()

if (sequence == 'OCCSTAR1'): # Plot second target as well, if we have it
    plt.imshow(stretch(image_target_summed_2), aspect=10)
    plt.title(sequence + ', target #2')
    plt.show()

# Now start the analysis

# Generate fake count rate data

print("Generating fake data...")
 
count_rate_target_fake = np.random.poisson(np.mean(count_rate_target), np.size(count_rate))
count_rate_fake        = np.random.poisson(np.mean(count_rate), np.size(count_rate))

# Compute UTC and ET for the initial timestep

utc_start= hbt.met2utc(np.min(met)) # Crashes kernel
et_start = sp.utc2et(utc_start)

# Compute MET and ET for all timesteps. Both ET and MET are in seconds, but their offset is different.

et       = met - met[0] + et_start    # ET is now fully populated and correct
t        = et - et[0]                 # Seconds since start of observation

num_dt   = np.size(et)

#==============================================================================
# Compute the Alice boresight RA/Dec position for every timestep
#==============================================================================

# NB: For FSS, I got the FSS-Sun angle directly from Gabe -- I didn't get it from SPICE.

# Get vector to star. 67 Ori = HR 2159 = HD 41753, a V=4.2 B3V. RA~91.89, Dec~14.77.

print('Computing boresight positions...')

if (sequence == 'O_RING_OC2') or (sequence == 'O_RING_OC3'):
    pos_star_str = "06 07 34.326 +14 46 06.51"  # Vizier coords in FK5 = J2K

if (sequence == 'OCCSTAR1'):
    pos_star_str = "06 12 03.27955 +16 07 49.4614"  # HD 42545, brighter, closer to lollipop.  
    pos_star2_str= "06 15 25.12824 +16 08 35.4219"  # HD 42153, fainter,  closer to tip of stick
    pos_star2 = SkyCoord(pos_star2_str, unit=(u.hourangle, u.deg))
    vec_star2_j2k = sp.radrec(1, pos_star2.ra.rad, pos_star2.dec.rad)
         
pos_star = SkyCoord(pos_star_str, unit=(u.hourangle, u.deg))
ra_star  = pos_star.ra.rad
dec_star = pos_star.dec.rad

ra       = np.zeros(num_dt)
dec      = np.zeros(num_dt)

vec_star_j2k = sp.radrec(1., ra_star, dec_star)

name_fov = 'NH_ALICE_AIRGLOW'

vec_bsight_alice = (-1, 0, 0)  # -X defines Alice SOC FOV

vsep = np.zeros(np.size(et))

for i,et_i in enumerate(et):
  mx = sp.pxform(name_fov, 'J2000', et[i])
#  vec_alice_j2k = cspice.mxvg(mx, vec_bsight_alice)		# ICY did not have mxvg(), but pdstools does.
  vec_alice_j2k = sp.mxvg(mx, vec_bsight_alice, 3, 3)		# ICY did not have mxvg(), but pdstools does.

  vsep[i] = sp.vsep(vec_star_j2k, vec_alice_j2k)   # Angular separation, in radians

  (junk, ra[i], dec[i]) = sp.recrad(vec_alice_j2k)

#==============================================================================
# Save results into a Pickle file
#==============================================================================

#lun = open(file_pickle, 'wb')
#pickle.dump((ra, dec, ra_star, dec_star, 
#             t, et, met, num_dt, utc_start, et_start, 
#             count_rate_fake, count_rate_data), lun)
#lun.close()
#print("Wrote: " + file_pickle)
        
#==============================================================================
# Use linear fit to compute correlation between count rate and RA / Dec
#==============================================================================

# I think this is not used any more??
  
coeffs_ra   = linregress(ra, count_rate_target)
coeffs_dec  = linregress(dec, count_rate_target)
count_rate_nonlinear = coeffs_ra.intercept + coeffs_ra.slope * ra - np.mean(count_rate_target)
count_rate_target_fixed     = count_rate_target - count_rate_nonlinear 

#==============================================================================
# For OCCSTAR1, do a polynomial fit to remove effect of motion within deadband
#==============================================================================

if (sequence in ['O_RING_OC2', 'O_RING_OC3']):     # No vignetting for these sequences
    count_rate_target_fixed = count_rate_target
    
if (sequence == 'OCCSTAR1'):

    count_rate_target_3000   = hbt.smooth_boxcar(count_rate_target,   3000)
    count_rate_target_2_3000 = hbt.smooth_boxcar(count_rate_target_2, 3000)
    
    # For star 1, look at samples 0 .. 500K (ie, before the distance pluto appulse starts). 
    # For these, do a polynomial fit between (dec) and (DN).
    
    y = count_rate_target_3000[0:500000] / 5.1 # y is defined as a normalization factor
    x = dec[0:500000]
    order = 5
    poly = np.polyfit(x, y, order) # Get the coefficients
    f_norm    = np.poly1d(poly) # Define a function of x. This is the normalization function,
                                # so that (corrected flux) = (measured flux) / f_norm(dec)

    count_rate_target_fixed = count_rate_target / f_norm(dec)
                              
    # For star 2, look at samples from 1.1M to 2.0M -- the samples after the Pluto occ.
    # Do fit between (dec) and (DN).

    y = count_rate_target_2_3000[1100000:2000000] / 1.7 # y is defined as a normalization factor
    x = dec[1100000:2000000]
    order = 5
    poly_2 = np.polyfit(x, y, order)
    f_2_norm    = np.poly1d(poly_2)
    
    x1 = 850000 # Approximate end of occultation. Use the 'fixed' data afer this,
                # and the raw data before it (since during occultation, there is no deadband effect)
    
    count_rate_target_2_fixed = count_rate_target_2.copy() / f_2_norm(dec)
    count_rate_target_2_fixed[0:x1] = count_rate_target_2[0:x1]          

#==============================================================================
# Smooth the data at several different binnings
#==============================================================================

print("Smoothing to various binnings...")

# Do statistics for first star in the aperture

# 'count_rate'                Is raw out of the Alice frame, including all pixels
# 'count_rate_target'         Extracts just the proper rows and columns
# 'count_rate_target_fixed'   Is same as above, but corrected for vignetting
 
#count_rate_fixed_30000 = hbt.smooth_boxcar(count_rate_fixed, 30000)
#count_rate_fixed_3000 = hbt.smooth_boxcar(count_rate_fixed, 3000)
#count_rate_fixed_300 = hbt.smooth_boxcar(count_rate_fixed, 300)
#count_rate_fixed_30 = hbt.smooth_boxcar(count_rate_fixed, 30)

count_rate_target_30000 = hbt.smooth_boxcar(count_rate_target, 30000)

count_rate_target_fixed_30000 = hbt.smooth_boxcar(count_rate_target_fixed, 30000)
count_rate_target_fixed_3000 = hbt.smooth_boxcar(count_rate_target_fixed, 3000)
count_rate_target_fixed_300 = hbt.smooth_boxcar(count_rate_target_fixed, 300)
count_rate_target_fixed_30 = hbt.smooth_boxcar(count_rate_target_fixed, 30)
count_rate_target_fixed_5 = hbt.smooth_boxcar(count_rate_target_fixed, 5)
count_rate_target_fixed_3 = hbt.smooth_boxcar(count_rate_target_fixed, 3)

count_rate_3000 = hbt.smooth_boxcar(count_rate, 3000)
count_rate_30000 = hbt.smooth_boxcar(count_rate, 30000)

count_rate_fake_3000 = hbt.smooth_boxcar(count_rate_fake, 3000)
count_rate_fake_30000 = hbt.smooth_boxcar(count_rate_fake, 30000)
count_rate_target_fake_30000 = hbt.smooth_boxcar(count_rate_target_fake, 30000)

if (sequence == 'OCCSTAR1'): # Do statistics for second star in the aperture, if there is one
    
    count_rate_target_2_30000 = hbt.smooth_boxcar(count_rate_target_2, 30000)
    count_rate_target_2_3000 = hbt.smooth_boxcar(count_rate_target_2, 3000)
    count_rate_target_2_300 = hbt.smooth_boxcar(count_rate_target_2, 300)
    count_rate_target_2_30 = hbt.smooth_boxcar(count_rate_target_2, 30)
    count_rate_target_2_5 = hbt.smooth_boxcar(count_rate_target_2, 5)
    count_rate_target_2_3 = hbt.smooth_boxcar(count_rate_target_2, 3)

#    count_rate_target_fixed_30000 = hbt.smooth_boxcar(count_rate_target_fixed, 30000)
#    count_rate_target_fixed_3000 = hbt.smooth_boxcar(count_rate_target_fixed, 3000)
#    count_rate_target_fixed_300 = hbt.smooth_boxcar(count_rate_target_fixed, 300)
#    count_rate_target_fixed_30 = hbt.smooth_boxcar(count_rate_target_fixed, 30)
#    count_rate_target_fixed_5 = hbt.smooth_boxcar(count_rate_target_fixed, 5)
#    count_rate_target_fixed_3 = hbt.smooth_boxcar(count_rate_target_fixed, 3)

    count_rate_target_2_fixed_30000 = hbt.smooth_boxcar(count_rate_target_2_fixed, 30000)
    count_rate_target_2_fixed_3000 = hbt.smooth_boxcar(count_rate_target_2_fixed, 3000)
    count_rate_target_2_fixed_300 = hbt.smooth_boxcar(count_rate_target_2_fixed, 300)
    count_rate_target_2_fixed_30 = hbt.smooth_boxcar(count_rate_target_2_fixed, 30)
    count_rate_target_2_fixed_5 = hbt.smooth_boxcar(count_rate_target_2_fixed, 5)
    count_rate_target_2_fixed_3 = hbt.smooth_boxcar(count_rate_target_2_fixed, 3)
    
##########
# Calculate statistics
##########

# Sigma_s = S / SNR = (mean) / (sqrt(n * mean) / (n * mean))
# This is the expected stdev for binned data (i.e., 63% of data should be within this amount of the mean, etc.)

sigma_s = np.mean(count_rate) * np.sqrt(np.mean(count_rate * binning)) / (np.mean(count_rate * binning))

#==============================================================================
# Do some geometry calculations: Get the projected radii from Pluto.
#==============================================================================

# Define a plane centered on Pluto. This is in Pluto coords (not J2K).
# The orientation of IAU_PLUTO will change from start to end of observation, and this will give us a 
# different answer in the end vs. start.
         
#et_start = np.min(et)
#et_end   = np.max(et) + duration[-1]

print("Getting distance to Pluto...")

plane_plu = sp.nvp2pl([0,0,1], [0,0,0])    # nvp2pl: Normal Vec + Point to Plane

# Get vector from Pluto to star, in IAU_PLUTO. This will define the *vector* portion of the ray.

vec_plu_star = sp.radrec(1., ra_star, dec_star) # Vector from Pluto to star, in J2K

# Set the values of ET for which we look up the projected radius. We could do every ET, but that's pretty slow, and the
# behavior is essentially linear. So, we do the opposite extreme, which is just the starting and ending values.

et_vals = hbt.mm(et)

radius_pluto = np.zeros(np.size(et_vals))
radius_bary  = np.zeros(np.size(et_vals)) # Distance from the barycenter
lon          = np.zeros(np.size(et_vals))
lat          = np.zeros(np.size(et_vals))

for i,et_i in enumerate(et_vals):

# Define a ray from NH, toward the star.
# Ray starts at NH position and points toward star (*not* boresight). It should be in IAU_PLUTO coordxmas.

# Get vector from Pluto to S/C, in IAU_PLUTO. This will define the *point* (the ray is a point plus a vector)

    (st_plu_sc_plu_i, lt_i) = sp.spkezr('NEW_HORIZONS', et_i, 'IAU_PLUTO', 'LT+S', 
                                        'PLUTO BARYCENTER') # Need to change this to barycenter
    pt_plu_sc_plu_i = st_plu_sc_plu_i[0:3]

    mx_i = sp.pxform('J2000', 'IAU_PLUTO', et_i)
    vec_plu_star_plu_i = sp.mxvg(mx_i, vec_plu_star, 3, 3)

# Now find the intersection. A ray is *not* a datatype. 
# Instead, we just pass a point and a vector, and that defines the ray.
# Point = vector from Pluto to NH in Pluto coords.
# Vector = vector from Pluto star in Pluto coords

    (npts, pt_intersect_plu_i) = sp.inrypl(pt_plu_sc_plu_i, vec_plu_star_plu_i, plane_plu) 
                                                       # intersect ray and plane. Jup coords.

    (radius_bary[i], lon[i], lat[i]) = sp.reclat(pt_intersect_plu_i)  # Convert to lon/lat and distance from Pluto

#==============================================================================
# Compute angular distance from star to Pluto, for each timestep
#==============================================================================

r_pluto_km = 1187

if (sequence == 'OCCSTAR1'):
        
    # Loop over every point in scan, and do some geometry.
    # But, SPICE calcs here are very slow, and the curves are smooth. It takes > 1 hour to do all
    # the calcs explicitly. Though I don't do it much, here we just sample it at low-res, and that's fine.
    
    skip = 1000 # Pick one element every 'skip' points to keep.
    t_skip = t[::skip]
    num_dt_skip = np.size(t_skip)
    
    print("Computing angular distance to Pluto for both occ stars...")
    ang_target_center   = np.zeros(num_dt_skip)
    ang_target_2_center = np.zeros(num_dt_skip)
    ang_target_limb     = np.zeros(num_dt_skip)
    ang_target_2_limb   = np.zeros(num_dt_skip)
    ang_target_center_radii = np.zeros(num_dt_skip) # Center-to-center angle, in Pluto radii
    ang_target_2_center_radii = np.zeros(num_dt_skip)
    
    for i in range(num_dt_skip): # This loop is very slow. Like 20 minutes.
        (st_pl, lt) = sp.spkezr('Pluto', et[i*skip], 'J2000', 'LT+S', 'New Horizons')
        vec_pl_j2k = st_pl[0:3]

        angle_radius_pluto =  (r_pluto_km / sp.vnorm(st_pl[0:3])) # Size in radians
        
        ang_target_center[i]   = sp.vsep(vec_pl_j2k, vec_star_j2k) # Angle in radians
        ang_target_2_center[i] = sp.vsep(vec_pl_j2k, vec_star2_j2k)
        ang_target_limb[i]     = ang_target_center[i] - angle_radius_pluto
        ang_target_2_limb[i]   = ang_target_2_center[i] - angle_radius_pluto
        ang_target_center_radii = ang_target_center / angle_radius_pluto
        ang_target_center_radii[i] = ang_target_center[i] / angle_radius_pluto
        ang_target_2_center_radii[i] = ang_target_2_center[i] / angle_radius_pluto
            
    
#==============================================================================
# Make a time-series plot of Counts vs. Time, for Target, Off-Target, Fake Data, etc.
#==============================================================================

# Plot of count rate vs. time

if (sequence == 'O_RING_OC3') or (sequence == 'O_RING_OC2'): # Complex plot -- do this only for the single stellar occs

    plt.rcParams['figure.figsize'] = 15,5
    
    if (sequence == 'O_RING_OC3'):
        offset_fake = 60
        offset_total = 765
        offset_diff = 2545
        
    if (sequence == 'O_RING_OC2'):
        offset_fake = 45
        offset_total = 765
        offset_diff = 2630
    
    binning = 30000
    
    do_fake = False

    # Calc stellar velocity. It's not constant, but this gives it in terms of (total motion) / (total time).
    # That is fine for O_Ring_Oc2 and O_Ring_Oc3
    
    vel = (radius_bary[1]-radius_bary[0]) / (np.amax(t) - np.amin(t)) # Velocity in km/sec
    
    # Jump through some hoops to place a second x-axis here: et vs. radius_pluto
    
    host = host_subplot(111, axes_class=AA.Axes) # Set up the host axis
    par = host.twiny()                           # Set up the parasite axis
    plt.subplots_adjust(bottom=0.2)              # Adjusts overall height of the whole plot in y direction 
                                                 #  (like figure.figsize)
    offset = 50                                  # How far away from the main plot the parasite axis is. 
                                                 #   Don't know the units.
    new_fixed_axis     = par.get_grid_helper().new_fixed_axis
    par.axis["bottom"] = new_fixed_axis(loc="bottom", axes=par,
                                            offset=(0,-offset))
    par.axis["bottom"].toggle(all=True)          # Make sure the bottom axis is displayed
    par.axis["top"].set_visible(False)           # Do not display the axis on *top* of the plot.
            
    p1, = host.plot(t, 1/dt * count_rate_target_30000, linewidth=0.5, ms=0.1, label='Alice, 67 Ori only')
    host.plot(t, 1/dt * count_rate_30000 - offset_total,               marker = '.', linewidth=0.5, 
              ms=0.1, label='Alice, Total [+offset]')
    host.plot(t, 1/dt * count_rate_30000 - count_rate_target_30000/dt + offset_diff, linewidth=0.5, 
              ms=0.1, label='Alice, Total - 67 Ori [+offset]')
    plt.plot(t, 1/dt * count_rate_target_fake_30000 - offset_fake, linewidth=0.5, ms=0.1, 
             label='Fake Poisson data [+offset]')
    
    plt.title(sequence + ', dt = ' + repr(dt) + ' sec, smoothed x ' + repr(binning) + ' = ' + 
                         repr(dt * binning) + ' sec', \
                         fontsize=fs)
    host.get_xaxis().get_major_formatter().set_useOffset(False)

    plt.legend(loc='lower right', framealpha=0.95)
    plt.ylabel('Counts/sec (smoothed)', fontsize=fs)
    plt.xlabel('Time since ' + utc_start + ' [sec]', fontsize=fs)
    
    plt.xlim((t[binning], t[-binning])) # X limits for main plot. Remove transients at edges. 
    par.set_xlim(radius_bary)               # X limits for parasite plot. 
                                            # *** This is slightly off, since I am ignoring edges above
    par.set_xlabel('Distance from Pluto barycenter [km]', fontsize=fs)
    
    if (sequence == 'O_RING_OC3'):
      plt.ylim((3100,3450))
    
    if (sequence == 'O_RING_OC2'):
      plt.ylim((3250,3500))

    file_out = (dir_out + '/' + sequence + '_vignetting_fixed' +
        ('_mask-lya'   if DO_MASK_LYA   else '') + '.png')
        
    plt.savefig(file_out)
    print("Wrote: " + file_out)
    plt.show()

if (sequence == 'OCCSTAR1'):
    
        binning = 3000 
        DO_PLOT_NORMALIZED = True
        DO_PLOT_RAW = True
          
        plt.rcParams['figure.figsize'] = 15,5
    
        fix, ax1 = plt.subplots()
        
#        ax2 = ax1.twinx()
#        
#        ax2.plot(t_skip, ang_target_center_radii, linestyle = 'dashed', linewidth=2,
#                 label = 'Sep from Pluto Center, Star 1', color='blue')
#        ax2.plot(t_skip, ang_target_2_center_radii, linestyle = 'dashed', linewidth=2,
#                 label = 'Sep from Pluto Center, Star 2', color='green')
    
        if (DO_PLOT_RAW):
            ax1.plot(t, 1/dt * count_rate_target_3000, label='Count Rate, Star 1', color='lightblue')
        
        if (DO_PLOT_NORMALIZED):
            ax1.plot(t, 1/dt * count_rate_target_3000 /f_norm(dec), label='Count Rate, Star 1', color='darkblue')
        
        if (DO_PLOT_RAW):
            ax1.plot(t, 1/dt * count_rate_target_2_3000, label='Count Rate Fixed, Star 2', color='lightgreen')     
            
        if (DO_PLOT_NORMALIZED):
            x1 = 850000
            ax1.plot(t[0:x1], 1/dt * count_rate_target_2_3000[0:x1], 
                     label='Count Rate Fixed, Star 2', color='darkgreen')
            ax1.plot(t[x1:], 1/dt * count_rate_target_2_3000[x1:] / f_2_norm(dec[x1:]), color='darkgreen')
            
        ax1.set_xlabel('Time since ' + utc_start + ' [sec]', fontsize=fs)
        ax1.set_ylabel('Counts/sec')
        ax1.set_title(sequence + ', dt = ' + repr(dt) + ' sec, smoothed x ' + repr(binning) + ' = ' + 
                  repr(int(dt * binning)) + ' sec', fontsize=fs)
    
#        ax2.set_ylabel('Pluto-Star Separation [$r_P$]')
        ax1.legend(framealpha=0.8, loc='center left', fontsize=fs*0.7)
#        ax2.legend(framealpha=0.8, loc='center right',fontsize=fs*0.7)
        ax1.set_xlim(np.array(hbt.mm(t)))
#        ax2.set_ylim([0.04, 5])
        ax1.set_ylim([-290, 1400])
        ax1.text(1000, 1050, 'HD 42545')
        ax1.text(7000, 100,  'HD 42153')
        
        file_out = (dir_out + '/' + sequence + '_vignetting_fixed' +
            ('_mask-lya'   if DO_MASK_LYA   else '') + '.png')
        
        plt.savefig(file_out)
        print("Wrote: " + file_out)
        plt.show()

        # 
        # Now make a second set of plots -- split into two, for clarity.
        #
        
        # Plot 1 (Star 1)
        # Calc velocity of star. It moves in & out, so use a ballpark total distance.

        vel = 5 * r_pluto_km / (np.amax(t) - np.amin(t)) # Velocity is not constant, but this is ballpark km/sec
        
        hbt.figsize((15,3))
        
        fix, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax2.plot(t_skip, ang_target_center_radii, linestyle = 'dashed', linewidth=2,
                 label = 'Sep from Pluto Center', color='blue') 
        
        ax1.plot(t, 1/dt * count_rate_target_3000 /f_norm(dec), label='Count Rate, Fixed', color='darkblue')
        
        ax1.set_ylabel('Counts/sec')
        ax1.set_ylabel('Counts/sec')
        ax1.set_title(sequence + ', dt = ' + repr(dt) + ' sec, smoothed x ' + repr(binning) + ' = ' + 
                  repr(int(dt * binning)) + ' sec', fontsize=fs)
    
        ax1.text(100, 1080, 'HD 42545, v = {:4.2f} km/sec'.format(vel))
        ax2.set_ylabel('Pluto-Star Separation [$r_P$]')
        ax1.legend(framealpha=0.8, loc='upper left', fontsize=fs*0.7)
        ax2.legend(framealpha=0.8, loc='lower right',fontsize=fs*0.7)
        ax1.set_xlim(np.array(hbt.mm(t)))
        ax2.set_ylim([0, 10])
        ax1.set_ylim([1000, 1400])
        ax1.axes.get_xaxis().set_ticks([])
        plt.show()

        # Plot 2 (Star 2)
        
        # Calc velocity of star. It moves in & out, so use a ballpark total distance.
        
        vel_2 = 6 * r_pluto_km / (np.amax(t) - np.amin(t)) # Velocity is not constant, but this is ballpark km/sec
        
        fix, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax2.plot(t_skip, ang_target_2_center_radii, linestyle = 'dashed', linewidth=2,
                 label = 'Sep from Pluto Center', color='green') 

        x1 = 850000 # Leftward of this x value, plot the non-fixed count rate (ie, in the occultation).
                    # Rightward of this x value, plot the fixed count rate
                    
        ax1.plot(t[0:x1], 1/dt * count_rate_target_2_3000[0:x1], 
                     label='Count Rate, Fixed', color='darkgreen')
        ax1.plot(t[x1:], 1/dt * count_rate_target_2_3000[x1:] / f_2_norm(dec[x1:]), color='darkgreen')     
                    
        ax1.set_ylabel('Counts/sec')
        ax1.set_xlabel('Time since ' + utc_start + ' [sec]', fontsize=fs)
        ax1.set_ylabel('Counts/sec')

        ax1.text(100, 200, 'HD 42153, {:4.2f} km/sec'.format(vel_2))
        ax2.set_ylabel('Pluto-Star Separation [$r_P$]')
        ax1.legend(framealpha=0.8, loc='upper left', fontsize=fs*0.7)
        ax2.legend(framealpha=0.8, loc='lower right',fontsize=fs*0.7)
        ax1.set_xlim(np.array(hbt.mm(t)))
        ax2.set_ylim([0.2, 6])
        ax1.set_ylim([000, 500])
        plt.show()
        
        
        
#==============================================================================
# Make a plot of (angle from star) vs (time). This is a line-plot of thruster firings.
#==============================================================================

plt.rcParams['figure.figsize'] = 15,5

plt.plot(et - et[0], (vsep % 0.02) * hbt.r2d)
plt.xlabel('Seconds since ' + utc_start, fontsize=fs)
plt.ylabel('Degrees from star', fontsize=fs)
plt.show()

#==============================================================================
# Make a line plot of motion thru the deadband -- following the FOV thru every thruster firing
#==============================================================================

hbt.set_fontsize(9)
plt.rcParams['figure.figsize'] = 5,5
plt.plot(ra*hbt.r2d, dec*hbt.r2d, linestyle='none', marker='.', ms=1)
plt.title(sequence,fontsize=12)
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
ax = plt.gca()
ax.ticklabel_format(useOffset=False)  # Turn off the 'offset' that matplotlib can use
plt.show()

hbt.set_fontsize(12)

#==============================================================================
# Make a plot of RA vs DN. 
# This is to see if there is a trend in brightness as we move from one side of slit to the other.
#==============================================================================

ra_s = ra[binning:-binning]    # Same as RA, but cropped, so edges missing. 
                               # Has same # elements as count_rate_s (smoothed)
dec_s = dec[binning:-binning]

crop = 3000  # Don't plot quite at the start and end of array, to avoid transients

plt.rcParams['figure.figsize'] = 16,8

plt.subplot(1,2,1)
plt.rcParams['figure.figsize'] = 10,10
plt.plot(ra[crop:-crop]*hbt.r2d, count_rate_3000[crop:-crop], linestyle='none', marker='.', ms=0.1)
#plt.plot(ra[crop:-crop]*hbt.r2d, count_rate_nonlinear[crop:-crop] + 
#         np.mean(count_rate), color='red') # Add linear fit on top
plt.xlabel('RA [deg]', fontsize=fs)
plt.title(sequence + ': DN vs. Position', fontsize=fs*1.5)
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.ylabel('DN (smoothed x ' + repr(crop) + ')', fontsize=fs)

plt.subplot(1,2,2)
plt.rcParams['figure.figsize'] = 10,10
plt.plot(dec[crop:-crop]*hbt.r2d, count_rate_3000[crop:-crop], linestyle='none', marker='.', ms=0.1)
plt.xlabel('Dec [deg]', fontsize=fs)
plt.title(sequence + ': DN vs. Position', fontsize=fs*1.5)
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.ylabel('DN (smoothed x ' + repr(crop) + ')', fontsize=fs)
plt.show()

hbt.set_fontsize(12)

if (sequence == 'OCCSTAR1'): # Made plot for first star above. Now do it for the second star
    plt.rcParams['figure.figsize'] = 16,8
    plt.subplot(1,2,1)
    plt.rcParams['figure.figsize'] = 10,10
    plt.plot(ra[crop:-crop]*hbt.r2d, count_rate_target_2_3000[crop:-crop], linestyle='none', marker='.', ms=0.1)
    #plt.plot(ra[crop:-crop]*hbt.r2d, count_rate_nonlinear[crop:-crop] + 
    #         np.mean(count_rate), color='red') # Add linear fit on top
    plt.xlabel('RA [deg]', fontsize=fs)
    plt.title(sequence + ': DN vs. Position, #2', fontsize=fs*1.5)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    plt.ylabel('DN (smoothed x ' + repr(crop) + ')', fontsize=fs)
    
    plt.subplot(1,2,2) # Plot DN vs. Position. Exclude the occultation, so only take a subset of the data.
    plt.rcParams['figure.figsize'] = 10,10
    plt.plot(dec[1200000:2000000]*hbt.r2d, count_rate_target_2_3000[1200000:2000000], 
             linestyle='none', marker='.', ms=0.1)
    plt.xlabel('Dec [deg]', fontsize=fs)
    plt.title(sequence + ': DN vs. Position, #2', fontsize=fs*1.5)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    plt.ylabel('DN (smoothed x ' + repr(crop) + ')', fontsize=fs)
    plt.show()

#==============================================================================
# Regrid the data, and make a plot of the actual spatial variation.
# This is similar to plots above, but in a color image.
# Must use binning = 3000 (or something small) for this.
#==============================================================================
# NB: The effect that this routine was plotting was caused by not subtracting
# the field star. Now that I am extracting only the flux of the target star, then
# this routine is not needed.

plt.rcParams['figure.figsize'] = 5,5
plt.subplot(1,1,1)

num_dx = 100
num_dy = num_dx

ra_arr  = np.linspace(np.min(ra),  np.max(ra),  num_dx)
dec_arr = np.linspace(np.min(dec), np.max(dec), num_dy)

count_rate_s_arr = griddata((ra[crop:-crop], dec[crop:-crop]), count_rate_3000[crop:-crop], 
                            (ra_arr[None,:], dec_arr[:,None]), method='cubic')

# Make the plot, scaled vertically as per the sequence read in

if (sequence == 'O_RING_OC3'):
    plt.imshow(count_rate_s_arr, interpolation='none', vmin=13,vmax=14)
if (sequence == 'O_RING_OC2'):
    plt.imshow(count_rate_s_arr, interpolation='none', vmin=16.4,vmax=16.8)
if (sequence == 'OCCSTAR1'):
    plt.imshow(count_rate_s_arr, interpolation='none', vmin=3,vmax=6) # Star 1. Not a really useful plot, though.
    plt.imshow(count_rate_s_arr, interpolation='none', vmin=0,vmax=3) # Star 2
    
plt.title(sequence + ', Mean DN vs position, binning=' + repr(binning), fontsize=fs)
plt.xlabel('RA bin', fontsize=fs)
plt.ylabel('Dec bin', fontsize=fs)
plt.colorbar()
plt.show()

####################
# Now do some testing to read in the housekeeping data. Randy says the same data from the main FITS header should be
# available here, but in 'analog' form rather than 'digital.'
####################

# Just test one file (or file list) for now

file = '/Users/throop/Data/NH_Alice_Ring/O_RING_OC3/data/pluto/level2/ali/all/ali_0299391413_0x4b5_sci_1.fit'

# Open the file directly
# Can use hdulist.info() to get a list of all the different fields available.

hdulist = fits.open(file)
        
spect         = hdulist['PRIMARY'].data            # 32           x 1024              
hk            = hdulist['HOUSEKEEPING_TABLE']
pix           = hdulist['PIXEL_LIST_TABLE']          # A time-tagged list of every photon that has come in. 
                                                     # 30883 x 5. columns X_INDEX, Y_INDEX, WAVELEN, TIMESTEP, MET.
                                                     # During this 7-second observation, 30883 photons were 
                                                     # received & counted.

timestep_list = pix.data['TIMESTEP']                 # Extract the list of 30883 timesteps.

# Print a few fields from the HK data. However, these seem to be mostly flags and settings -- 
# no real 'raw' data here useful for me.
# This 'COUNT_RATE' field is a bit of a misnomer. I don't know what it is. It seems to be integrated 
# over all bins or something.

print('Housekeeping: MET = ' + repr(int(hk.data['MET'])))
print('Housekeeping: COUNT_RATE = ' + repr(int(hk.data['COUNT_RATE'])))

# Now try to extract the 'raw' data from the time-tagged pixel list.

# First create bins -- one for each timestep, including start and end... so N+1 elements in this total. 
# We feed this to histogram() to define the bins.

bins = hbt.frange(0, np.max(timestep_list)+1)

# Now count how many photons are in each timestep bin. I have defined those timestep bins up above with range().

(count_rate_raw, junk) = np.histogram(timestep_list, bins)

# Now read the count rate as processed by the SOC directly, from the FITS 'COUNT_RATE' extension.

count_rate_soc = hbt.read_alice(file)['count_rate']

# Concl: the values in the FITS 'COUNT_RATE' extension, and in the time-tagged photon list, are 100% identical.

hk.data

#==============================================================================
# Now make a time-series plot, but with several binning widths. To look for narrow rings.
#==============================================================================

if ((sequence == 'O_RING_OC2') or (sequence == 'O_RING_OC3')):
    
    DO_PLOT_LEGEND_TIME_SEQUENCE = True  # Boolean: Do we draw the legend on the plot?
    DO_PLOT_TITLE_TIME_SEQUENCE = False  # Boolean: Do we draw the title on the plot?
    
    plt.rcParams['figure.figsize'] = 20,12
    #plt.rcParams['figure.figsize'] = 10,6
    
    #fs = 25
    # Jump through some hoops to place a second x-axis here: et vs. radius_pluto
    
    plt.rc('font', size=20)
    
    host = host_subplot(111, axes_class=AA.Axes) # Set up the host axis
    par = host.twiny()                           # Set up the parasite axis
    plt.subplots_adjust(bottom=0.2)              # Adjusts overall height of the whole plot in y direction 
    offset = 50                                  # How far away from the main plot the parasite axis is.
    new_fixed_axis     = par.get_grid_helper().new_fixed_axis
    par.axis["bottom"] = new_fixed_axis(loc="bottom", axes=par,
                                        offset=(0,-offset))
    par.axis["bottom"].toggle(all=True)          # Make sure the bottom axis is displayed
    par.axis["top"].set_visible(False)           # Do not display the axis on *top* of the plot.
    
    res_1 = 1 * dt*u.s * vel*u.km/u.s
    res_30 = 30 * dt*u.s * vel*u.km/u.s
    res_300 = 300 * dt*u.s * vel*u.km/u.s
    res_3000 = 3000 * dt*u.s * vel*u.km/u.s
    res_30000 = 30000 * dt*u.s * vel*u.km/u.s
    
    p1, = host.plot(t, count_rate_target_fixed/dt, ms=0.1, linestyle='none', marker='.', color='orange', 
             label = 'Raw, {:6.3f} sec = {:3.1f} m [orange]'.format(dt, res_1.to('meter').value))
    
    host.plot(t, count_rate_target_fixed_30/dt, linestyle='none', marker='.', linewidth=0.1, ms=1, color='green', 
             label = 'Binned x30, {:5.2f} sec = {:3.1f} m [green]'.format(30*dt, res_30.to('m').value))
    
    host.plot(t, count_rate_target_fixed_300/dt, linewidth=0.1, ms=0.1, color='blue', 
             label = 'Binned x300, {:4.1f} sec = {:3.1f} km [blue]'.format(300*dt, res_300.to('km').value))
    
    host.plot(t, count_rate_target_fixed_3000/dt, linewidth=0.3, ms=0.1, color='red', 
             label = 'Binned x3000, {:3.0f} sec = {:2.0f} km [red]'.format(3000*dt, res_3000.to('km').value))
    
    host.plot(t, count_rate_target_fixed_30000/dt, linewidth=0.2, ms=0.1, color='yellow', 
             label = 'Binned x30000, {:3.0f} sec = {:3.0f} km [yellow]'.format(30000*dt, res_30000.to('km').value))
    
    #host.set_ylim((2800,4000))
    host.set_ylim((2500,4200))
    plt.xlim(hbt.mm(t))
    #plt.xlim((1750,1850))
    par.set_xlim(hbt.mm(radius_bary))
    if (DO_PLOT_TITLE_TIME_SEQUENCE):
        plt.title(sequence, fontsize=fs)
        
    plt.ylabel('DN/sec', fontsize=fs) # Fontsize is ignored here, probably because of the twin-axis thing...
    plt.xlabel('Seconds since ' + utc_start.split('.')[0], fontsize=fs)
    par.set_xlabel('Distance from Pluto barycenter [km]', fontsize=fs)
    
    # Plot the legend, if requested
    
    if (DO_PLOT_LEGEND_TIME_SEQUENCE):
        plt.legend(loc = 'upper left', framealpha=0.90)
        leg = host.get_legend()
        leg.legendHandles[0].set_linestyle('solid')  # Change symbols from points, to solid lines
        leg.legendHandles[1].set_linestyle('solid')
        leg.legendHandles[2].set_linestyle('solid')
        leg.legendHandles[3].set_linestyle('solid')
        leg.legendHandles[4].set_linestyle('solid')
        leg.legendHandles[4].set_linewidth(4)
        leg.legendHandles[3].set_linewidth(4)
        leg.legendHandles[2].set_linewidth(4)
        leg.legendHandles[1].set_linewidth(4)
        leg.legendHandles[0].set_linewidth(4)

    file_out = (dir_out + sequence + '_all' +
      ('_legend'  if DO_PLOT_LEGEND_TIME_SEQUENCE else '') +
      ('_title'   if DO_PLOT_TITLE_TIME_SEQUENCE  else '') + 
      ('_mask-lya'if DO_MASK_LYA                  else '') +
      '.png')

    plt.savefig(file_out)
    print("Wrote: " + file_out)
    plt.show()

if (sequence == 'OCCSTAR1') and True:

    # Plot first subplot, for star #1

    DO_PLOT_LEGEND_TIME_SEQUENCE = True  # Boolean: Do we draw the legend on the plot?
    DO_PLOT_TITLE_TIME_SEQUENCE = False  # Boolean: Do we draw the title on the plot?
        
    fs = 18 # Fontsize
    
    plt.rcParams['figure.figsize'] = 17,8 # Plot size
    hbt.set_fontsize(size=fs)

    # Do a very rough calculation of shadow velocity.
    
    host = host_subplot(111, axes_class=AA.Axes) # Set up the host axis
    plt.subplots_adjust(bottom=0.2)              # Adjusts overall height of the whole plot in y direction 
    offset = 50                                  # How far away from the main plot the parasite axis is.
#    new_fixed_axis     = par.get_grid_helper().new_fixed_axis

    p1, = host.plot(t, count_rate_target/dt, ms=0.1, linestyle='none', marker='.', color='orange', 
             label = 'Raw, ' + repr(dt) + ' sec = {:.0f} m'.format(dt*vel*1000))
    
    host.plot(t, count_rate_target_fixed_30/dt, linestyle='none', marker='.', linewidth=0.1, ms=1, color='green', 
             label='Binned x30, ' + repr(30*dt) + ' sec = {:.0f} m'.format(dt*vel*30*1000) )
    
    host.plot(t, count_rate_target_fixed_300/dt, linewidth=0.1, ms=0.1, color='blue', 
             label='Binned x300, ' + repr(300*dt) + ' sec = {:.0f} km'.format(dt*vel*300))
    
    host.plot(t, count_rate_target_fixed_3000/dt, linewidth=0.3, ms=0.1, color='red', 
             label='Binned x3000, ' + repr(3000*dt) + ' sec = {:.0f} km'.format(dt*vel*3000))
    
    host.plot(t, count_rate_target_fixed_30000/dt, linewidth=0.2, ms=0.1, color='yellow', 
             label='Binned x30000, ' + repr(30000 * dt) + 
             ' sec = {:.0f} km'.format(dt*vel*30000) )
    
    host.set_ylim((500,1650))
    
    plt.xlim(hbt.mm(t))

    plt.text(5000, 1525, "HD 42545", fontsize=fs*1.4) 

    plt.ylabel('Counts/sec', fontsize=fs) # Fontsize is ignored here, probably because of the twin-axis thing...
#    plt.xlabel('Seconds since ' + utc_start, fontsize=fs)

    if (DO_PLOT_LEGEND_TIME_SEQUENCE):
    
        plt.legend(framealpha=0.8, loc='lower left', fontsize=fs*0.75)
        
    # Retrieve the individual lines for the legend(), and tweak them a bit. Make them lines, not points
    # This was hard to figure out, but got w help from http://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html
    
        leg = host.get_legend()
        leg.legendHandles[0].set_linestyle('solid')
        leg.legendHandles[1].set_linestyle('solid')
        leg.legendHandles[2].set_linestyle('solid')
        leg.legendHandles[3].set_linestyle('solid')
        leg.legendHandles[4].set_linestyle('solid')
        leg.legendHandles[4].set_linewidth(4)
        leg.legendHandles[3].set_linewidth(4)
        leg.legendHandles[2].set_linewidth(4)
        leg.legendHandles[1].set_linewidth(4)
        leg.legendHandles[0].set_linewidth(4)
    
    DO_SUPPRESS_X_AXIS = False

    if (DO_SUPPRESS_X_AXIS):
        host.get_xaxis().set_ticks([]) # Do not print the x axis at all!
    else:
        plt.xlabel('Seconds since ' + utc_start.split('.')[0])    

    file_out = (dir_out + sequence + '_star1_all' +
      ('_legend'   if DO_PLOT_LEGEND_TIME_SEQUENCE else '') +
      ('_title'    if DO_PLOT_TITLE_TIME_SEQUENCE  else '') + 
      ('_mask-lya' if DO_MASK_LYA                  else '') +
      '.png')

    plt.savefig(file_out)
    print("Wrote: " + file_out)

    plt.show()


########### Now plot second subplot, for Star #2

# NB: For some reason this plot is messed up the second time it gets plotted. Works the first time only.
#     Something changes with count_rate_target_2 apparently.

#    dist_2 = (5 * r_pluto_km) - (-1*r_pluto_km) # From -1 RP, to 5 RP
#    v_2 = dist_2 / (np.amax(t) - np.amin(t))    # km/sec

    host = host_subplot(111, axes_class=AA.Axes) # Set up the host axis
    plt.subplots_adjust(bottom=0.2)              # Adjusts overall height of the whole plot in y direction 
    offset = 50                                  # How far away from the main plot the parasite axis is.
    par = host.twiny()                           # Set up the parasite axis
    new_fixed_axis     = par.get_grid_helper().new_fixed_axis
    
    p1, = host.plot(t, count_rate_target_2/dt, ms=0.1, linestyle='none', marker='.', color='orange', 
             label = 'Raw, ' + repr(dt) + ' sec = {:.0f} m'.format(dt*vel_2*1000))
    
    host.plot(t, count_rate_target_2_fixed_30/dt, linestyle='none', marker='.', linewidth=0.2, ms=1, color='green', 
             label='Binned x30 bins, ' + repr(30*dt) + ' sec = {:.0f} m'.format(dt*vel_2*30*1000))
    
    host.plot(t, count_rate_target_2_fixed_300/dt, linewidth=0.2, ms=0.1, color='blue', 
             label='Binned x300 bins, ' + repr(300*dt) + ' sec = {:.0f} km'.format(dt*vel_2*300))
    
    host.plot(t, count_rate_target_2_fixed_3000/dt, linewidth=0.3, ms=0.1, color='red', 
             label='Binned x3000 bins, ' + repr(3000*dt) + ' sec = {:.0f} km'.format(dt*vel_2*3000))
    
    host.plot(t, count_rate_target_2_fixed_30000/dt, linewidth=0.3, ms=0.1, color='yellow', 
             label='Binned x30000 bins, ' + repr(int(30000 * dt)) + 
             ' sec = {:.0f} km'.format(dt*vel_2*30000) )
    
    host.set_ylim((000,800))
    
    plt.xlim(hbt.mm(t))

    plt.text(1000, 400, "HD 42153", fontsize=fs*1.4) 
    plt.ylabel('Counts/sec', fontsize=fs) # Fontsize is ignored here, probably because of the twin-axis thing...
    plt.xlabel('Seconds since ' + utc_start.split('.')[0])    

    plt.legend(framealpha=0.8, loc = 'upper left', fontsize=fs*0.75)
    
    leg = host.get_legend()
    leg.legendHandles[0].set_linestyle('solid')
    leg.legendHandles[1].set_linestyle('solid')
    leg.legendHandles[2].set_linestyle('solid')
    leg.legendHandles[3].set_linestyle('solid')
    leg.legendHandles[4].set_linestyle('solid')
    leg.legendHandles[4].set_linewidth(4)
    leg.legendHandles[3].set_linewidth(4)
    leg.legendHandles[2].set_linewidth(4)
    leg.legendHandles[1].set_linewidth(4)
    leg.legendHandles[0].set_linewidth(4)
    
    file_out = file_out.replace('star1', 'star2')
    plt.savefig(file_out)
    print("Wrote: " + file_out)
    
    plt.show()

#==============================================================================
# Zoom in on one area of interest in OC2
#==============================================================================

if (sequence == 'O_RING_OC2'):
    plt.rcParams['figure.figsize'] = 20,12
    
    #fs = 25
    # Jump through some hoops to place a second x-axis here: et vs. radius_pluto
    
    host = host_subplot(111, axes_class=AA.Axes) # Set up the host axis
    par = host.twiny()                           # Set up the parasite axis
    plt.subplots_adjust(bottom=0.2)              # Adjusts overall height of the whole plot in y direction 
                                                 #  (~ figure.figsize)
    offset = 50                                  # How far away from the main plot the parasite axis is. 
                                                 #  Don't know units.
                                                 
    new_fixed_axis     = par.get_grid_helper().new_fixed_axis
    par.axis["bottom"] = new_fixed_axis(loc="bottom", axes=par,
                                        offset=(0,-offset))
    par.axis["bottom"].toggle(all=True)          # Make sure the bottom axis is displayed
    par.axis["top"].set_visible(False)           # Do not display the axis on *top* of the plot.
    
    p1, = host.plot(t, count_rate_target_fixed/dt, ms=0.2, linestyle='none', marker='+', color='orange', 
             label = 'Alice, raw, ' + repr(dt) + ' sec = 1.5 m [orange]')

    host.plot(t, count_rate_target_fixed_5/dt, linestyle='none', marker='+', linewidth=0.1, ms=5, color='black', 
             label='Alice, smoothed 5 bins = ' + repr(5*dt) + ' sec = 8 m [black]' )
    
    host.plot(t, count_rate_target_fixed_30/dt, marker='.', linewidth=2, ms=1, color='green', 
             label='Alice, smoothed 30 bins = ' + repr(30*dt) + ' sec = 40 m [green]' )
    
    host.plot(t, count_rate_target_fixed_300/dt, linewidth=2, ms=0.1, color='blue', 
             label='Alice, smoothed 300 bins = ' + repr(300*dt) + ' sec = 0.4 km [blue]')
    
    host.plot(t, count_rate_target_fixed_3000/dt, linewidth=1, ms=0.1, color='red', 
             label='Alice, smoothed 3000 bins = ' + repr(3000*dt) + ' sec = 4 km [red]')
    
#    host.plot(t, count_rate_target_30000/dt, linewidth=2, ms=0.1, color='yellow', 
#             label='Alice, smoothed 30000 bins = ' + repr(30000 * dt) + 
#             ' sec = 40 km [yellow]' )
    
    host.get_xaxis().get_major_formatter().set_useOffset(False)

    plt.xlim(hbt.mm(t))

# Select the x range to zoom in on. Both of these ranges below are interesting areas to check out.
# I don't think either are statistically significant.
    
    plt.xlim((1792,1798)) # Use 1792 .. 1798
#    plt.xlim((1370,1380)) # Use 1370 .. 1380


    plt.ylim((1500,5000))
#    par.set_xlim(hbt.mm(radius_bary)) # Not correct for the zoomed plot!
    plt.title(sequence + ' zoom', fontsize=fs)
    plt.ylabel('Counts/sec', fontsize=fs) # Fontsize is ignored here, probably because of the twin-axis thing...
    plt.xlabel('Seconds since ' + utc_start, fontsize=fs)
    par.set_xlabel('Distance from Pluto barycenter [km]', fontsize=fs)
    
    plt.legend()
    plt.show()
 
# Calculate some statistics for that blip at met0 + 1796 sec

    std_5 = np.std(count_rate_target_fixed_5)/dt
    std_30 = np.std(count_rate_target_fixed_30)/dt
    
    depth_5 = 3400-2000
    depth_30 = 3400 - 2800
    
    print("At binning = 5, depth = " + repr(depth_5) + " = " + repr(depth_5 / std_5) + " sigma")
    print("At binning = 30, depth = " + repr(depth_30) + " = " + repr(depth_30 / std_30) + " sigma")
    
#==============================================================================
# Do a histogram of the count rate
#==============================================================================

bins = hbt.frange(0, int(np.max(count_rate))+1)/dt
(h,b) = np.histogram(count_rate_target/dt,bins=bins)
plt.rcParams['figure.figsize'] = 5,5
plt.plot(bins[0:-1],h, drawstyle='steps')
plt.xlabel('67 Ori, Counts/sec')
plt.ylabel('# of samples')
plt.title(sequence)
plt.xlim((0,7000))
#plt.yscale('log')
#plt.set_yscale('log')
plt.show()

#==============================================================================
# Look for any count=1 bins right next to each other
#==============================================================================

plt.rcParams['figure.figsize'] = 5,5
w = np.where(count_rate < 4)[0]
plt.plot(range(np.size(w)), w)
plt.show()

## This is kind of an auto-correlation function. Maybe I should look at that too.
## But honestly it looks unlikely.
#==============================================================================
# Calculate the Fresnel limit
#==============================================================================

(st, lt) = sp.spkezr('NEW_HORIZONS', et[0], 'IAU_PLUTO', 'LT+S', 'PLUTO')
dist_start = sp.vnorm(st[0:3]) * u.km # Distance in km
alam = 100 * u.nm
d_fresnel_start = np.sqrt(dist_start * alam/2).to('m').value

(st, lt) = sp.spkezr('NEW_HORIZONS', et[-1], 'IAU_PLUTO', 'LT+S', 'PLUTO')
dist_end = sp.vnorm(st[0:3]) * u.km # Distance in km

d_fresnel_end = np.sqrt(dist_end * alam/2).to('m').value

#==============================================================================
# Make a binned plot right at the fresnel limit
#==============================================================================

plt.rcParams['figure.figsize'] = 15,5
plt.plot(t, count_rate_target_fixed_5/dt,linestyle='none', ms=0.5, marker='.')
plt.title(sequence + ', binning = 5 = Fresnel limit')
plt.xlabel('Seconds')
plt.ylabel('Counts/sec')
plt.xlim(hbt.mm(t))
plt.show()

#==============================================================================
# Do some autocorrelation
#==============================================================================

plt.rcParams['figure.figsize'] = 5,5

offsets = np.array(hbt.frange(-100,100),dtype=int)
count_rate_target_2 = hbt.smooth_boxcar(count_rate_target_fixed,2)
count_rate_target_3 = hbt.smooth_boxcar(count_rate_target_fixed,3)
count_rate_target_10 = hbt.smooth_boxcar(count_rate_target_fixed,10)

corr_full = np.zeros(np.size(offsets))
for i in range(np.size(offsets)):
    corr_full[i] = np.correlate(count_rate_target_fixed_30, np.roll(count_rate_target_fixed_30,offsets[i]))

corr = corr_full[100:]    

plt.plot(hbt.ln01(corr),drawstyle='steps')
plt.xlabel('Offset [bins]')
plt.ylabel('ln(Correlation)')
plt.title(sequence + ", Binning = 30")
plt.show()
    
#==============================================================================
# Now make a data table for the data, for the final optical depth
#==============================================================================

if ((sequence == 'O_RING_OC2') or (sequence == 'O_RING_OC3')):

    x0 = 100000 # Starting sample number
    x1 = 400000 # Ending sample number

if (sequence == 'OCCSTAR1'):
    
    x0   = 400000  # For 'Star 1' of OCCSTAR1. Skip the Pluto occultation / appulse.
    x1   = 700000

    x0_2 = 1500000  # For 'Star 2' of OCCSTAR1. Skip the Pluto occultation / appulse.
    x1_2 = 2000000
    
# Calculate the fractional 3-sigma variance from the mean.
# F = F0 * exp(-tau / cos(theta))
#   Where F0 is the unocculted flux, and F is the occulted.
#   I assume here that F = F0 - 3*stdev(F0)
#
# Solve this and it gives tau = -cos * ln(F/F0)
#   Where F = F_3S = F0 - 3*stdev(F0)

# Star 1

nsig =  3  # 3 sigma? 5 sigma? Plug it in here (3, 4, 5, etc)

# f0 is the mean flux, in bins x0:x1
# f_3s = f0 - 3*stdev(f0)  -- see above equation

f0 = np.mean(count_rate_target[x0:x1])
f_3s       = f0 - nsig * np.std(count_rate_target_fixed[x0:x1]) # This one might be negative due to noise. Not a problem
f_3s_3     = f0 - nsig * np.std(count_rate_target_fixed_3[x0:x1])
f_3s_30    = f0 - nsig * np.std(count_rate_target_fixed_30[x0:x1])
f_3s_300   = f0 - nsig * np.std(count_rate_target_fixed_300[x0:x1])
f_3s_3000  = f0 - nsig * np.std(count_rate_target_fixed_3000[x0:x1])
f_3s_30000 = f0 - nsig * np.std(count_rate_target_fixed_30000[x0:x1])
    
tau_norm       = -np.cos(subobslat) * np.log(f_3s       / f0)
tau_norm_3     = -np.cos(subobslat) * np.log(f_3s_3     / f0)
tau_norm_30    = -np.cos(subobslat) * np.log(f_3s_30    / f0)
tau_norm_300   = -np.cos(subobslat) * np.log(f_3s_300   / f0)
tau_norm_3000  = -np.cos(subobslat) * np.log(f_3s_3000  / f0)
tau_norm_30000 = -np.cos(subobslat) * np.log(f_3s_30000 / f0)

# We want to calculate optical depth at the fresnel limit. We could just figure out the proper binning width
# for this, and make an array in the code for it. But instead, we just interpolate from other widths, and say
# what the optical depth should be. 
# We assume the signal is basically all shot noise, so tau ~ 1/sqrt(width), which it does.

# How many bins wide is the fresnel limit?

d_fresnel = (d_fresnel_start + d_fresnel_end)/2 # Fresnel limit, in meters

bins_fresnel = d_fresnel / (vel * dt *1000)

# Now interpolate logarithmically .

x = np.array([3,          30,          300,          3000])  # Binning width, in bins
y = np.array([tau_norm_3, tau_norm_30, tau_norm_300, tau_norm_3000]) # Optical depth

# Modify these to be something we can do a linear fit against

x = np.sqrt(x)
y = 1/y

# Do a linear fit

m,b = np.polyfit(x, y, 1)

# And calculate the derived optical depth for the given binning width

tau_fresnel = 1 / (m * np.sqrt(bins_fresnel) + b)

# Make a plot of binning width vs. tau. This is just a check for linearity. It goes exactly as it should.
#plt.plot(x,y, marker='+')
#plt.xlabel('sqrt(bining width, bins)')
#plt.ylabel('1/tau')

print()
print("Sequence = {}, LyA masking = {}".format(sequence, DO_MASK_LYA))
print("----------------------------")
print("Star 1, Binning 3     = {:.3f} m, tau <= {:.3f}".format(dt*vel * 3*1000, tau_norm_3))
print("Star 1, Binning 30    = {:.3f} km, tau <= {:.3f}".format(dt*vel * 30, tau_norm_30))
print("Star 1, Binning 300   = {:.1f} km, tau <= {:.3f}".format(dt*vel * 300, tau_norm_300))
print("Star 1, Binning 3000  = {:.1f} km, tau <= {:.3f}".format(dt*vel * 3000, tau_norm_3000))
print("Star 1, Binning 30000 = {:.1f} km, tau <= {:.3f}".format(dt*vel * 30000, tau_norm_30000))
print("Star 1, Binning {:.1f} = {:.1f} m, tau <= {:.3f}    **** Fresnel limit ***".format(
      bins_fresnel, bins_fresnel * vel * dt * 1000, tau_fresnel))

# If there are two stars, now do statistics for the second one

if (sequence == 'OCCSTAR1'): # If there are two stars, then now do the second star
    # Star 2

    
    f0_2         = np.mean(count_rate_target_2[x0_2:x1_2])
    f_3s_2       = f0_2 - nsig * np.std(count_rate_target_2_fixed[x0_2:x1_2])
    f_3s_2_3     = f0_2 - nsig * np.std(count_rate_target_2_fixed_3[x0_2:x1_2])
    f_3s_2_30    = f0_2 - nsig * np.std(count_rate_target_2_fixed_30[x0_2:x1_2])
    f_3s_2_300   = f0_2 - nsig * np.std(count_rate_target_2_fixed_300[x0_2:x1_2])
    f_3s_2_3000  = f0_2 - nsig * np.std(count_rate_target_2_fixed_3000[x0_2:x1_2])
    f_3s_2_30000 = f0_2 - nsig * np.std(count_rate_target_2_fixed_30000[x0_2:x1_2])
    
    tau_norm_2       = -np.cos(subobslat) * np.log(f_3s_2       / f0_2)
    tau_norm_2_3     = -np.cos(subobslat) * np.log(f_3s_2_3     / f0_2)
    tau_norm_2_30    = -np.cos(subobslat) * np.log(f_3s_2_30    / f0_2)
    tau_norm_2_300   = -np.cos(subobslat) * np.log(f_3s_2_300   / f0_2)
    tau_norm_2_3000  = -np.cos(subobslat) * np.log(f_3s_2_3000  / f0_2)
    tau_norm_2_30000 = -np.cos(subobslat) * np.log(f_3s_2_30000 / f0_2)

    bins_fresnel_2 = d_fresnel / (vel_2 * dt *1000)
    x = np.array([3, 30, 300, 3000])  # Binning width, in bins
    y = np.array([tau_norm_2_3, tau_norm_2_30, tau_norm_2_300, tau_norm_2_3000]) # Optical depth    
    x = np.sqrt(x)
    y = 1/y
    m,b = np.polyfit(x, y, 1)
    tau_fresnel_2 = 1 / (m * np.sqrt(bins_fresnel_2) + b)
    
    print()
    print("Star 2, Binning 30    = {:.3f} km, tau <= {:.3f}".format(dt*vel_2 * 30, tau_norm_2_30))
    print("Star 2, Binning 300   = {:.1f} km, tau <= {:.3f}".format(dt*vel_2 * 300, tau_norm_2_300))
    print("Star 2, Binning 3000  = {:.0f} km, tau <= {:.3f}".format(dt*vel_2 * 3000, tau_norm_2_3000))
    print("Star 2, Binning 30000 = {:.0f} km, tau <= {:.3f}".format(dt*vel_2 * 30000, tau_norm_2_30000))
    print("Star 2, Binning {:.1f} = {:.1f} m, tau <= {:.3f}    **** Fresnel limit ***".format(bins_fresnel_2, 
          bins_fresnel_2 * vel_2 * dt * 1000, tau_fresnel_2))
    print()

print("These are {}-sigma values".format(nsig))
        
#
#==============================================================================
# Print some statistics, for use in tables.
#==============================================================================

print("Start time = " + sp.et2utc(et[0], 'C', 2))
print("End time = " + sp.et2utc(et[-1], 'C', 2))
print("Duration = {:.2f} s".format(et[-1]-et[0]))
print("Fresnel scale START = {:.2f} m".format(d_fresnel_start))
print("Fresnel scale END   = {:.2f} m".format(d_fresnel_end))

if (sequence == 'OCCSTAR1'):
    print("Star 1 dist from Pluto center: {:.1f} .. {:.1f} km. HD 42545.".format(
            np.min(ang_target_center_radii)*r_pluto_km, np.max(ang_target_center_radii)*r_pluto_km))
    
    print("Star 2 dist from Pluto center: {:.1f} .. {:.1f} km. HD 42153.".format(
            np.min(ang_target_2_center_radii)*r_pluto_km, np.max(ang_target_2_center_radii)*r_pluto_km))

    print("Star 1 dist from Pluto center: {:.2f} .. {:.2f} R_P".format(
            np.min(ang_target_center_radii), np.max(ang_target_center_radii)))
    
    print("Star 2 dist from Pluto center: {:.2f} .. {:.2f} R_P.".format(
            np.min(ang_target_2_center_radii), np.max(ang_target_2_center_radii)))
    
    print("Star 1 velocity typical = {:.2f} km/s".format(vel))
    print("Star 2 velocity typical = {:.2f} km/s".format(vel_2))
    print("Star 1 Per-sample resolution typical = {:.2f} m".format(vel * dt * 1000))
    print("Star 2 Per-sample resolution typical = {:.2f} m".format(vel_2 * dt * 1000))
    
else:
    print("Distance from Pluto Barycenter: {:.1f} .. {:.1f} km".format(radius_bary[0], radius_bary[1]))
    print("Velocity: {:.2f} km/s".format(vel))
    print("Per-sample resolution typical = {:.2f} m".format( vel*dt *1000   )   )

#==============================================================================
# Testing for LyA masking
#==============================================================================

# This is a simple test to make sure that the DO_MASK_LYA flag is working as it is supposed to.
# I get strange results. For the fainter star in OCCSTAR1 (star 2), omitting the LyA flux reduces 
# flux in the image to about 85% of nominal. But when masking out the same bands in the pixel-list
# (aka count rate), the change is only to about 98%.
# Why the difference? Not sure. My only guess is that perhaps the image arrage is not counting # of photons.
# It is a floating-point value. So maybe it is wavelength-weighted??
# If I am wrong and not masking, it's not really a big deal. It will weaken my optical depth limits by maybe 10%
# for one star, and lesss for the others. Not a big deal at all. Not worth obsessing over.
#
# requires sequence = OCCSTAR1

DO_TEST_LYA_MASKING = False

if DO_TEST_LYA_MASKING:
    print('Reading star 2, mask False...')
    (met, count_rate_target_2, count_rate_2, image_target_summed_2, image_summed) = \
        read_alice_occ_data(file_list, xlim, ylim_2, verbose=True, short=DO_ABBREVIATED, DO_MASK_LYA=False)
    
    print()
    print('Reading star 2, mask True...')  # _m -> masked (ie, with LyA removed)
    (met_m, count_rate_target_2_m, count_rate_2_m, image_target_summed_2_m, image_summed_m) = \
        read_alice_occ_data(file_list, xlim, ylim_2, verbose=True, short=DO_ABBREVIATED, DO_MASK_LYA=True)
    
    stretch = astropy.visualization.PercentileInterval(99)
     
    hbt.figsize((10,5))
    
    plt.imshow(stretch(image_target_summed_2), aspect=10)
    plt.title('target_summed')
    plt.show()
    
    plt.imshow(stretch(image_target_summed_2_m), aspect=10)
    plt.title('target_summed')
    plt.show()
    
    print("Ratio of images = {}".format(np.sum(image_target_summed_2_m) / np.sum(image_target_summed_2)))
    print("Ratio of time series = {}".format(np.sum(count_rate_2_m) / np.sum(count_rate_2)))
    print("Ratio of time series target = {}".format(np.sum(count_rate_target_2_m) / np.sum(count_rate_target_2)))
       