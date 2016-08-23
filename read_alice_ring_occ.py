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
import cspice
import skimage
from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
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
from scipy.interpolate import griddata

from mpl_toolkits.axes_grid1 import host_subplot # For adding a second axis to a plot
import mpl_toolkits.axisartist as AA             # For adding a second axis to a plot

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

# Imports for Tk

import Tkinter
import ttk
import tkMessageBox
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure

# HBT imports
import hbt

def read_alice_occ_data(file_list):

#==============================================================================
# Read the Alice data
#==============================================================================
    
    met_all = []               # A long array with a list of all of the timestamps, one per 4 ms (i.e., at 250 hz)
    count_rate_fits_all = []   # The count rate as read from the COUNT_RATE extension directly
    count_rate_all = []        # Count rate computed from the PIXEL_LIST_TABLE. Should match that in COUNT_RATE extension
    count_rate_target_all = [] # Count rate for the target only, extracted by spatially filtering the PIXEL_LIST_TABLE
    
    d_target_summed = np.zeros((5,540))
    
    for i,file in enumerate(file_list):
        
        hdulist = fits.open(file)
        d = hdulist['PRIMARY'].data # Units of this are float, but I'm not sure what they are. I would prefer raw counts.
        d_target = d[13:18, 370:910]
        d_target_summed += d_target
        p = hdulist['PIXEL_LIST_TABLE'].data
        count_rate_fits_i = hdulist['COUNT_RATE'].data
        num_samples = hdulist['COUNT_RATE'].header['NAXIS1'] # Number of samples in this file
        dt          = hdulist['COUNT_RATE'].header['SAMPLINT']  # Count rate sampling interval [sec]
        
        bins = hbt.frange(0, num_samples) # Get a list of all of the timestep bins, inclusive, for this file.
                                            # Use '+1' so we create the histogram upper size bin.
        
        # Now downselect the pixel list for just the photons in the proper X and Y position on the detector
        
        is_good = (p['Y_INDEX'] < 19) & (p['Y_INDEX'] >= 13) & (p['X_INDEX'] > 370) & (p['X_INDEX'] < 910)
    
        # Now we have a list of all of the good pixels. For each of these, now we want to grab its timestep.
    
        timesteps_good = p['TIMESTEP'][is_good]
        timesteps_all  = p['TIMESTEP']
    
    # Now count how many photons are in each timestep bin. I have defined those timestep bins up above.
    
        (count_rate_target_i, junk) = np.histogram(timesteps_good, bins)
        (count_rate_i, junk)        = np.histogram(timesteps_all,  bins) 
        met_i = hdulist['PRIMARY'].header['STARTMET'] + dt * np.array(range(num_samples))
              
    #    print "File " + os.path.basename(file) + ' : ' + repr(np.sum(count_rate_target_i)) + ' / ' + repr(np.sum(count_rate_i))
    
    # Now append these into the output lists
    
        count_rate_fits_all.append(count_rate_fits_i)
        count_rate_all.append(count_rate_i)
        count_rate_target_all.append(count_rate_target_i)
        met_all.append(met_i)
        
        hdulist.close()
    
    count_rate_fits  = np.array([item for sublist in count_rate_fits_all for item in sublist])  # Flatten the count rate (from 2D, to 1D)
    count_rate_fits  = np.array(count_rate_fits, dtype=float)					          # Convert to float. Otherwise get wraparound.
    
    count_rate  = np.array([item for sublist in count_rate_all for item in sublist])
    count_rate  = np.array(count_rate, dtype=float)					  
    
    count_rate_target  = np.array([item for sublist in count_rate_target_all for item in sublist]) 
    count_rate_target  = np.array(count_rate_target, dtype=float)					  
    
    met         = np.array([item for sublist in met_all for item in sublist])         # Flatten the MET array  (from 2D, to 1D)
    met         = np.array(met, dtype=float)
    
    return (met, count_rate_target, count_rate)

#==============================================================================
# Start of main program
#==============================================================================

##########
# Select which of the sequences we want to read in
##########
# (NB: The spelling / capitalization is inconsistent in SAPNAME vs. VISITNAM. I have standardized it here.)

sequence 	= 'O_RING_OC3'
#sequence 	= 'O_RING_OC2'


binning      = 25000		# Smoothing. 25000 is too much (shows opposite trend!). 5000 and 1000 look roughly similar.
                            # To be meaningful, the binning timescale must be less than the deadband timescale (~20-30 sec RT).
                            # At binning=3000, it's 12 sec... so that is the longest we'd really want to go.
                            #
                            # Use binning = 25,000 for the data analysis.
                            # Use binning = 3000 to make a plot of the trend of DN vs. RA
                            
fs           	= 15		# Font size

dir_images = '/Users/throop/Data/NH_Alice_Ring/' + sequence + '/data/pluto/level2/ali/all'
file_tm = '/Users/throop/gv/dev/gv_kernels_new_horizons.txt'

cspice.furnsh(file_tm)

#########
# Read the Alice data
#########

file_list = glob.glob(dir_images + '/*fit')

(met, count_rate_target, count_rate) = read_alice_occ_data(file_list)

#file_list = file_list[0:10]
#met_all = []  # A long array with a list of all of the timestamps
#count_rate_all = []
#duration_all = []
#
#for file in file_list:
#    data = hbt.read_alice(file)
##    print "."
#    startmet = data['header_spect']['STARTMET']
#    spcutcal = data['header_spect']['SPCUTCAL']	# Mid-obs time
#
#    visitnam = data['header_spect']['VISITNAM']
#    sapname   = data['header_spect']['SAPNAME']
#
#    dt       = data['header_count_rate']['SAMPLINT']
#    count_rate_i = data['count_rate']              # Extract count rate from file
#
#    met_i = startmet + dt * hbt.frange(0, np.shape(count_rate_i)[0]-1)   # Create the MET by interpolation
#    duration_i = dt * len(count_rate_i)    
#
#    
#    print "Read " + os.path.basename(file) + ", MET " + repr(startmet) + ' = ' + hbt.met2utc(startmet) + \
#      ', N = ' + repr(len(count_rate_i)) + ' samples' + \
#      ', duration = ' + hbt.trunc(duration_i,3) + ' sec'
##       print "  " + visitnam + " " + spcutcal
#
#    met_all.append(met_i)
#    count_rate_all.append(count_rate_i)
#    duration_all.append(duration_i)
#
#count_rate  = np.array([item for sublist in count_rate_all for item in sublist])  # Flatten the count rate (from 2D, to 1D)
#count_rate  = np.array(count_rate, dtype=float)					  # Convert to float. Otherwise get wraparound.
#met         = np.array([item for sublist in met_all for item in sublist])         # Flatten the MET array  (from 2D, to 1D)
#duration    = np.array(duration_all)

count_rate_fake = np.random.poisson(np.mean(count_rate), np.size(count_rate))

# Compute UTC and ET for the initial timestep

utc_start= hbt.met2utc(np.min(met))
et_start = cspice.utc2et(utc_start)

# Compute MET and ET for all timesteps. Both ET and MET are in seconds, but their offset is different.

et       = met - met[0] + et_start    # ET is now fully populated and correct
t        = et - et[0]                 # Seconds since start of observation

num_dt   = np.size(et)

#==============================================================================
# Compute the Alice boresight RA/Dec position for every timestep
#==============================================================================

# NB: For FSS, I got the FSS-Sun angle directly from Gabe -- I didn't get it from SPICE.

# Get vector to star. 67 Ori = HR 2159 = HD 41753, a V=4.2 B3V. RA=91.89, Dec=14.77.

ra_star  = 91.89 * hbt.d2r
dec_star = 14.77 * hbt.d2r

ra       = np.zeros(num_dt)
dec      = np.zeros(num_dt)

vec_star_j2k = cspice.radrec(1., ra_star, dec_star)

name_fov = 'NH_ALICE_AIRGLOW'

vec_bsight_alice = (-1, 0, 0)  # -X defines Alice SOC FOV

vsep = np.zeros(np.size(et))

for i,et_i in enumerate(et):
  mx = cspice.pxform(name_fov, 'J2000', et[i])
  vec_alice_j2k = cspice.mxvg(mx, vec_bsight_alice)		# ICY did not have mxvg(), but pdstools does.
  vsep[i] = cspice.vsep(vec_star_j2k, vec_alice_j2k)   # Angular separation, in radians

  (junk, ra[i], dec[i]) = cspice.recrad(vec_alice_j2k)

#==============================================================================
# Use linear fit to compute correlation between count rate and RA / Dec
#==============================================================================
  
coeffs_ra   = linregress(ra, count_rate)
coeffs_dec  = linregress(dec, count_rate)
count_rate_nonlinear = coeffs_ra.intercept + coeffs_ra.slope * ra - np.mean(count_rate)
count_rate_fixed     = count_rate - count_rate_nonlinear 

#==============================================================================
# Smooth the data at several different binnings
#==============================================================================

count_rate_fixed_30000 = hbt.smooth_boxcar(count_rate_fixed, 30000)
count_rate_fixed_3000 = hbt.smooth_boxcar(count_rate_fixed, 3000)
count_rate_fixed_300 = hbt.smooth_boxcar(count_rate_fixed, 300)
count_rate_fixed_30 = hbt.smooth_boxcar(count_rate_fixed, 30)

count_rate_target_30000 = hbt.smooth_boxcar(count_rate_target, 30000)
count_rate_target_3000 = hbt.smooth_boxcar(count_rate_target, 3000)
count_rate_target_300 = hbt.smooth_boxcar(count_rate_target, 300)
count_rate_target_30 = hbt.smooth_boxcar(count_rate_target, 30)

count_rate_3000 = hbt.smooth_boxcar(count_rate, 3000)
count_rate_30000 = hbt.smooth_boxcar(count_rate, 30000)

count_rate_fake_3000 = hbt.smooth_boxcar(count_rate_fake, 3000)
count_rate_fake_30000 = hbt.smooth_boxcar(count_rate_fake, 30000)

# Compute truncated versions of the time arrays, just in case they are useful
# _s extension = 'smoothed'

#t_s          = t[binning:-binning]
#et_s         = et[binning:-binning]
#met_s        = met[binning:-binning]

##########
# Calculate statistics
##########

# Sigma_s = S / SNR = (mean) / (sqrt(n * mean) / (n * mean))
# This is the expected stdev for binned data (i.e., 63% of data should be within this amount of the mean, etc.)

sigma_s = np.mean(count_rate) * np.sqrt(np.mean(count_rate * binning)) / (np.mean(count_rate * binning))

#==============================================================================
# Make plots
#==============================================================================

plt.rcParams['figure.figsize'] = 15,5

# Plot of count rate vs. time

if (sequence == 'O_RING_OC3'):
    offset_fake = 0.15
if (sequence == 'O_RING_OC2'):
    offset_fake = 0.15
    
#binning = 30000
    
plt.plot(t, count_rate_30000, marker = '.', linewidth=0.5, ms=0.1, label='Alice, Raw')
plt.plot(t, count_rate_fixed_30000, linewidth=0.5, ms=0.1, label='Alice, Linear Row Trend Removed')
plt.plot(t, count_rate_fake_30000 - offset_fake, linewidth=0.5, ms=0.1, label='Fake Poisson Data')

plt.title(sequence + ', dt = ' + repr(dt) + ' sec, smoothed x ' + repr(binning) + ' = ' + repr(dt * binning) + ' sec' + \
                     ', 1$\sigma$ = ' + hbt.trunc(sigma_s,4), fontsize=fs)

plt.errorbar(100, np.mean(count_rate_3000) + 5 * sigma_s, xerr=binning*dt/2, yerr=None, label='Binning Width', linewidth=2) 
			# X 'error bar' -- show the bin width
plt.errorbar(300, np.mean(count_rate_3000) + 5 * sigma_s, xerr=None, yerr=sigma_s/2, label='1$\sigma$ shot noise', linewidth=2) 
			# Y 'error bar' -- show the binned shot noise error
plt.legend()
plt.ylabel('Count Rate (smoothed)', fontsize=fs)
plt.xlabel('Time since ' + utc_start + ' [sec]', fontsize=fs)
plt.xlim((t[binning], t[-binning])) # Remove transients at edges

if (sequence == 'O_RING_OC3'):
  plt.ylim((np.mean(count_rate_3000) -10*sigma_s, np.mean(count_rate_3000) +11*sigma_s))

if (sequence == 'O_RING_OC2'):
  plt.ylim((np.mean(count_rate_3000) -9*sigma_s, np.mean(count_rate_3000) +11*sigma_s))
    
plt.show()

#==============================================================================
# Make a plot of (angle from star) vs (time). This is a line-plot of thruster firings.
#==============================================================================

plt.plot(et - et[0], (vsep % 0.02) * hbt.r2d)
plt.xlabel('Seconds since ' + utc_start, fontsize=fs)
plt.ylabel('Degrees from star', fontsize=fs)
plt.show()

#==============================================================================
# Make a line plot of motion thru the deadband -- following the FOV thru every thruster firing
#==============================================================================

plt.rcParams['figure.figsize'] = 5,5
plt.plot(ra*hbt.r2d, dec*hbt.r2d, linestyle='none', marker='.', ms=1)
plt.title(sequence + ', start = ' + utc_start)
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
ax = plt.gca()
ax.ticklabel_format(useOffset=False)  # Turn off the 'offset' that matplotlib can use
plt.show()

## Now do some correlation between position and pointing, to see if there is anything I should unwrap there.
#
# Concl: there is quite a bit of trend across the detector. So, a lot of the variation I'm seeing
# in DN value is probably not due to statistics, but due to sensitivity across the detector.

#==============================================================================
# Make a plot of RA vs DN. 
# This is to see if there is a trend in brightness as we move from one side of slit to the other.
#==============================================================================

ra_s = ra[binning:-binning]    # Same as RA, but cropped, so the edges are missing. Has same # elements as count_rate_s (smoothed)
dec_s = dec[binning:-binning]

crop = 3000  # Don't plot quite at the start and end of array, to avoid transients

plt.rcParams['figure.figsize'] = 16,8

plt.subplot(1,2,1)
plt.rcParams['figure.figsize'] = 10,10
plt.plot(ra[crop:-crop]*hbt.r2d, count_rate_3000[crop:-crop], linestyle='none', marker='.', ms=0.1)
plt.plot(ra[crop:-crop]*hbt.r2d, count_rate_nonlinear[crop:-crop] + np.mean(count_rate), color='red') # Add linear fit on top
plt.xlabel('RA [deg]', fontsize=fs)
plt.title(sequence + ': DN vs. Position', fontsize=fs*1.5)
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.ylabel('DN (smoothed x ' + repr(binning) + ')', fontsize=fs)

plt.subplot(1,2,2)
plt.rcParams['figure.figsize'] = 10,10
plt.plot(dec[crop:-crop]*hbt.r2d, count_rate_3000[crop:-crop], linestyle='none', marker='.', ms=0.1)
plt.xlabel('Dec [deg]', fontsize=fs)
plt.title(sequence + ': DN vs. Position', fontsize=fs*1.5)
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.ylabel('DN (smoothed x ' + repr(crop) + ')', fontsize=fs)
plt.show()

#==============================================================================
# Regrid the data, and make a plot of the actual spatial variation.
# This is similar to plots above, but in a color image.
# Must use binning = 3000 (or something small) for this.
#==============================================================================

plt.rcParams['figure.figsize'] = 5,5
plt.subplot(1,1,1)

num_dx = 100
num_dy = num_dx

ra_arr  = np.linspace(np.min(ra),  np.max(ra),  num_dx)
dec_arr = np.linspace(np.min(dec), np.max(dec), num_dy)

count_rate_s_arr = griddata((ra[crop:-crop], dec[crop:-crop]), count_rate_3000[crop:-crop], (ra_arr[None,:], dec_arr[:,None]), method='cubic')

# Make the plot, scaled vertically as per the sequence read in

if (sequence == 'O_RING_OC3'):
    plt.imshow(count_rate_s_arr, interpolation='none', vmin=16.2,vmax=16.45)
if (sequence == 'O_RING_OC2'):
    plt.imshow(count_rate_s_arr, interpolation='none', vmin=16.4,vmax=16.8)

plt.title(sequence + ', Mean DN vs position, binning=' + repr(binning), fontsize=fs)
plt.xlabel('RA bin', fontsize=fs)
plt.ylabel('Dec bin', fontsize=fs)
plt.colorbar()
plt.show()

#==============================================================================
# Now do some geometry calculations: Get the Pluto radii for this particular sequence.
#==============================================================================

# Overall picture

# Define a plane centered on Pluto. This is in Pluto coords (not J2K).
# The orientation of IAU_PLUTO will change from start to end of observation, and this will give us a 
# different answer in the end vs. start.
         
#et_start = np.min(et)
#et_end   = np.max(et) + duration[-1]

plane_plu = cspice.nvp2pl([0,0,1], [0,0,0])    # nvp2pl: Normal Vec + Point to Plane

# Get vector from Pluto to star, in IAU_PLUTO. This will define the *vector* portion of the ray.

vec_plu_star = cspice.radrec(1., ra_star, dec_star) # Vector from Pluto to star, in J2K

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

    (st_plu_sc_plu_i, lt_i) = cspice.spkezr('NEW_HORIZONS', et_i, 'IAU_PLUTO', 'LT+S', 'PLUTO BARYCENTER') # Need to change this to barycenter
    pt_plu_sc_plu_i = st_plu_sc_plu_i[0:3]

    mx_i = cspice.pxform('J2000', 'IAU_PLUTO', et_i)
    vec_plu_star_plu_i = cspice.mxvg(mx_i, vec_plu_star)

# Now find the intersection. A ray is *not* a datatype. Instead, we just pass a point and a vector, and that defines the ray.
# Point = vector from Pluto to NH in Pluto coords.
# Vector = vector from Pluto star in Pluto coords

    (npts, pt_intersect_plu_i) = cspice.inrypl(pt_plu_sc_plu_i, vec_plu_star_plu_i, plane_plu) # intersect ray and plane. Jup coords.

    (radius_bary[i], lon[i], lat[i]) = cspice.reclat(pt_intersect_plu_i)  # Convert to lon/lat and distance from Pluto

####################
# Now do some testing to read in the housekeeping data. Randy says that the same data from the main FITS header should be
# available here, but in 'analog' form rather than 'digital.'
####################

# Just test one file for now

file = '/Users/throop/Data/NH_Alice_Ring/O_RING_OC3/data/pluto/level2/ali/all/ali_0299391413_0x4b5_sci_1.fit' # Or file_list[0]

# Open the file directly
# Can use hdulist.info() to get a list of all the different fields available.

hdulist = fits.open(file)
        
spect         = hdulist['PRIMARY'].data            # 32           x 1024              
hk            = hdulist['HOUSEKEEPING_TABLE']
pix           = hdulist['PIXEL_LIST_TABLE']          # A time-tagged list of every photon that has come in. 30883 x 5. columns X_INDEX, Y_INDEX, WAVELEN, TIMESTEP, MET.
                                                     # During this 7-second observation, 30883 photons were received & counted.

timestep_list = pix.data['TIMESTEP']                 # Extract the list of 30883 timesteps.

# Print a few fields from the HK data. However, these seem to be mostly flags and settings -- no real 'raw' data here useful for me.
# This 'COUNT_RATE' field is a bit of a misnomer. I don't know what it is. It seems to be integrated over all bins or something.

print 'Housekeeping: MET = ' + repr(int(hk.data['MET']))
print 'Housekeeping: COUNT_RATE = ' + repr(int(hk.data['COUNT_RATE']))

# Now try to extract the 'raw' data from the time-tagged pixel list.

# First create bins -- one for each timestep, including start and end... so N+1 elements in this total. We feed this to histogram() to define the bins.

bins = hbt.frange(0, np.max(timestep_list)+1)

# Now count how many photons are in each timestep bin. I have defined those timestep bins up above with range().

(count_rate_raw, junk) = np.histogram(timestep_list, bins)

# Now read the count rate as processed by the SOC directly, from the FITS 'COUNT_RATE' extension.

count_rate_soc = hbt.read_alice(file)['count_rate']

# Concl: the values in the FITS 'COUNT_RATE' extension, and in the time-tagged photon list, are 100% identical.

hk.data

#==============================================================================
# Now do a plot with differeng binning widths. To look for narrow rings.
#==============================================================================

# Put another axis here: et vs. radius_pluto


#plt.rcParams['figure.figsize'] = 20,12
plt.rcParams['figure.figsize'] = 10,6

host = host_subplot(111, axes_class=AA.Axes) # Set up the host axis
par = host.twiny()                           # Set up the parasite axis
plt.subplots_adjust(bottom=0.2)              # Adjusts overall height of the whole plot in y direction (like figure.figsize)
offset = 50                                  # How far away from the main plot the parasite axis is. Don't know the units.
new_fixed_axis     = par.get_grid_helper().new_fixed_axis
par.axis["bottom"] = new_fixed_axis(loc="bottom", axes=par,
                                    offset=(0,-offset))
par.axis["bottom"].toggle(all=True)          # Make sure the bottom axis is displayed
par.axis["top"].set_visible(False)           # Do not display the axis on *top* of the plot.

p1, = host.plot(t, count_rate, ms=0.1, linestyle='none', marker='.', color='orange', 
         label = 'Alice, Raw, 250/sec = 1.5 m [orange]')
host.plot(t, count_rate_fixed_30, linewidth=0.1, ms=0.1, color='green', 
         label='Alice, smoothed 30 bins = ' + repr(30*dt) + ' = 4 sec = 40 m[green]' )
host.plot(t, count_rate_fixed_300, linewidth=0.1, ms=0.1, color='blue', 
         label='Alice, smoothed 300 bins = ' + repr(300*dt) + ' sec = 0.4 km [blue]')
host.plot(t, count_rate_fixed_3000, linewidth=0.3, ms=0.1, color='red', 
         label='Alice, smoothed 3000 bins = ' + repr(3000*dt) + ' sec = 4 km [red]')
host.plot(t, count_rate_fixed_30000, linewidth=0.2, ms=0.1, color='yellow', 
         label='Alice, smoothed 30000 bins = ' + repr(30000 * dt) + 
         ' sec = 40 km [yellow]' )
host.set_ylim((14,20))
plt.xlim(hbt.mm(t))
par.set_xlim(hbt.mm(radius_bary))
plt.title(sequence, fontsize=fs)
plt.ylabel('DN', fontsize=fs)
plt.xlabel('Seconds since ' + utc_start, fontsize=fs)
par.set_xlabel('Distance from Pluto barycenter [km]')

plt.legend()
plt.show()

 
#==============================================================================
# Do a histogram of the count rate
#==============================================================================

bins = hbt.frange(0, int(np.max(count_rate))+1)
(h,b) = np.histogram(count_rate,bins=bins)
plt.rcParams['figure.figsize'] = 5,5
plt.plot(bins[0:-1],h)
plt.xlabel('COUNT RATE')
plt.ylabel('# of samples')
plt.title(os.path.basename(file))
plt.yscale('log')
#plt.set_yscale('log')
plt.show()

#==============================================================================
# Look for any count=1 bins right next to each other
#==============================================================================

plt.rcParams['figure.figsize'] = 15,15
w = np.where(count_rate < 4)[0]
plt.plot(range(np.size(w)), w)

## This is kind of an auto-correlation function. Maybe I should look at that too.
## But honestly it looks unlikely.