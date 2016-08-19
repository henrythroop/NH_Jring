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

##########
# Select which of the sequences we want to read in
##########
# (NB: The spellling / capitalization is inconsistent in SAPNAME vs. VISITNAM. I have standardized it here.)

#sequence 	= 'O_RING_OC3'
sequence 	= 'O_RING_OC2'


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

#file_list = file_list[0:10]
met_all = []  # A long array with a list of all of the timestamps
count_rate_all = []
duration_all = []

for file in file_list:
    data = hbt.read_alice(file)
#    print "."
    startmet = data['header_spect']['STARTMET']
    spcutcal = data['header_spect']['SPCUTCAL']	# Mid-obs time

    visitnam = data['header_spect']['VISITNAM']
    sapname   = data['header_spect']['SAPNAME']

    dt       = data['header_count_rate']['SAMPLINT']
    count_rate_i = data['count_rate']              # Extract count rate from file

    met_i = startmet + dt * hbt.frange(0, np.shape(count_rate_i)[0]-1)   # Create the MET by interpolation
    duration_i = dt * len(count_rate_i)    

    
    print "Read " + os.path.basename(file) + ", MET " + repr(startmet) + ' = ' + hbt.met2utc(startmet) + \
      ', N = ' + repr(len(count_rate_i)) + ' samples' + \
      ', duration = ' + hbt.trunc(duration_i,3) + ' sec'
#       print "  " + visitnam + " " + spcutcal

    met_all.append(met_i)
    count_rate_all.append(count_rate_i)
    duration_all.append(duration_i)

count_rate  = np.array([item for sublist in count_rate_all for item in sublist])  # Flatten the count rate (from 2D, to 1D)
count_rate  = np.array(count_rate, dtype=float)					  # Convert to float. Otherwise get wraparound.
met         = np.array([item for sublist in met_all for item in sublist])         # Flatten the MET array  (from 2D, to 1D)
duration    = np.array(duration_all)

count_rate_fake = np.random.poisson(np.mean(count_rate), np.size(count_rate))

# Compute UTC and ET for the initial timestep

utc_start= hbt.met2utc(np.min(met))
et_start = cspice.utc2et(utc_start)

# Compute MET and ET for all timesteps. Both ET and MET are in seconds, but their offset is different.

et       = met - met[0] + et_start    # ET is now fully populated and correct
t        = et - et[0]                 # Seconds since start of observation

num_dt   = np.size(et)

##########
# Compute the angle from the star to Alice boresight, for every timestep.
# Also compute the RA and Dec for each timestep
##########

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
  
# Do some linear regression on the count rate vs. RA, to remove nonlinearity across the detector

coeffs_ra = linregress(ra, count_rate)
coeffs_dec = linregress(dec, count_rate)
count_rate_nonlinear = coeffs_ra.intercept + coeffs_ra.slope * ra - np.mean(count_rate)
count_rate_fixed = count_rate - count_rate_nonlinear 

# Computed a smoothed count rate. _s indicates smoothed. Array is cropped at edges too.

                            
count_rate_s = hbt.smooth_boxcar(count_rate, binning)[binning:-binning]

count_rate_fake_s = hbt.smooth_boxcar(count_rate_fake, binning)[binning:-binning]

count_rate_fixed_s = hbt.smooth_boxcar(count_rate_fixed, binning)[binning:-binning]

# Compute truncated versions of the time arrays, just in case they are useful

t_s          = t[binning:-binning]
et_s         = et[binning:-binning]
met_s        = met[binning:-binning]

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
    
plt.plot(t_s, count_rate_s, marker = '.', linewidth=0.5, ms=0.1, label='Alice, Raw')
plt.plot(t_s, count_rate_fixed_s, linewidth=0.5, ms=0.1, label='Alice, Linear Row Trend Removed')
plt.plot(t_s, count_rate_fake_s - offset_fake, linewidth=0.5, ms=0.1, label='Fake Poisson Data')

plt.title(sequence + ', dt = ' + repr(dt) + ' sec, smoothed x ' + repr(binning) + ' = ' + repr(dt * binning) + ' sec' + \
                     ', 1$\sigma$ = ' + hbt.trunc(sigma_s,4), fontsize=fs)

plt.errorbar(100, np.mean(count_rate_s) + 5 * sigma_s, xerr=binning*dt/2, yerr=None, label='Binning Width', linewidth=2) 
			# X 'error bar' -- show the bin width
plt.errorbar(300, np.mean(count_rate_s) + 5 * sigma_s, xerr=None, yerr=sigma_s/2, label='1$\sigma$ shot noise', linewidth=2) 
			# Y 'error bar' -- show the binned shot noise error
plt.legend()
plt.ylabel('Count Rate (smoothed)', fontsize=fs)
plt.xlabel('Time since ' + utc_start + ' [sec]', fontsize=fs)
plt.xlim((0, np.max(t_s)+50))

if (sequence == 'O_RING_OC3'):
  plt.ylim((np.mean(count_rate_s) -10*sigma_s, np.mean(count_rate_s) +11*sigma_s))

if (sequence == 'O_RING_OC2'):
  plt.ylim((np.mean(count_rate_s) -9*sigma_s, np.mean(count_rate_s) +11*sigma_s))
    
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

ra_s = ra[binning:-binning]    # Same as RA, but cropped, so the edges are missing. Has same # elements as count_rate_s (smoothed)
dec_s = dec[binning:-binning]

#==============================================================================
# Make a plot of RA vs DN. 
# This is to see if there is a trend in brightness as we move from one side of slit to the other.
#==============================================================================

plt.rcParams['figure.figsize'] = 16,8

plt.subplot(1,2,1)
plt.rcParams['figure.figsize'] = 10,10
plt.plot(ra_s*hbt.r2d, count_rate_s, linestyle='none', marker='.', ms=0.1)
plt.plot(ra*hbt.r2d, count_rate_nonlinear + np.mean(count_rate), color='red')
plt.xlabel('RA [deg]', fontsize=fs)
plt.title(sequence + ': DN vs. Position', fontsize=fs*1.5)
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.ylabel('DN (smoothed x ' + repr(binning) + ')', fontsize=fs)

plt.subplot(1,2,2)

plt.rcParams['figure.figsize'] = 10,10
plt.plot(dec_s*hbt.r2d, count_rate_s, linestyle='none', marker='.', ms=0.1)
plt.xlabel('Dec [deg]', fontsize=fs)
plt.title(sequence + ': DN vs. Position', fontsize=fs*1.5)
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.ylabel('DN (smoothed x ' + repr(binning) + ')', fontsize=fs)
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

count_rate_s_arr = griddata((ra_s, dec_s), count_rate_s, (ra_arr[None,:], dec_arr[:,None]), method='cubic')

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
         
et_start = np.min(et)
et_end   = np.max(et) + duration[-1]

for et_i in (et_start, et_end):
    plane_plu = cspice.nvp2pl([0,0,1], [0,0,0])    # nvp2pl: Normal Vec + Point to Plane

# Define a ray from NH, toward the star.
# Ray starts at NH position and points toward star (*not* boresight). It should be in IAU_PLUTO coordxmas.

# Get vector from Pluto to S/C, in IAU_PLUTO. This will define the *point* (the ray is a point plus a vector)

    (st_plu_sc_plu, lt) = cspice.spkezr('NEW_HORIZONS', et_i, 'IAU_PLUTO', 'LT+S', 'PLUTO')
    pt_plu_sc_plu = st_plu_sc_plu[0:3]

# Get vector from Pluto to star, in IAU_PLUTO. This will define the *vector* portion of the ray.

    vec_plu_star = cspice.radrec(1., ra_star, dec_star) # Vector from Pluto to star, in J2K
    mx = cspice.pxform('J2000', 'IAU_PLUTO', et_i)
    vec_plu_star_plu = cspice.mxvg(mx, vec_plu_star)

# Now find the intersection. A ray is *not* a datatype. Instead, we just pass a point and a vector, and that defines the ray.
# Point = vector from Pluto to NH in Pluto coords.
# Vector = vector from Pluto star in Pluto coords

    (npts, pt_intersect_plu) = cspice.inrypl(pt_plu_sc_plu, vec_plu_star_plu, plane_plu) # intersect ray and plane. Jup coords.

    (radius, lon, lat) = cspice.reclat(pt_intersect_plu)  # Convert to lon/lat and distance from Pluto

    print cspice.et2utc(et_i, 'C', 1)
    print "Intercept distance = " + hbt.trunc(radius,1) + " km from Pluto"
#    print "NH distance = " + hbt.trunc(cspice.vnorm(pt_plu_sc_plu),1) + " km from Pluto"
    print
    
# Call CSPICE_INRYPL to get the intersection between ray and plane.

# Then get the distance from this point, to Pluto (or alternatively, the Pluto system barycenter)

quit

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

binning = [1, 30, 300, 3000, 30000]

plt.rcParams['figure.figsize'] = 20,12

count_rate_fixed_30000 = hbt.smooth_boxcar(count_rate_fixed, 30000)
count_rate_fixed_3000 = hbt.smooth_boxcar(count_rate_fixed, 3000)
count_rate_fixed_300 = hbt.smooth_boxcar(count_rate_fixed, 300)
count_rate_fixed_30 = hbt.smooth_boxcar(count_rate_fixed, 30)

plt.plot(t, count_rate, ms=0.1, linestyle='none', marker='.', color='orange', 
         label = 'Alice, Raw, 250/sec [orange]')
plt.plot(t, count_rate_fixed_30, linewidth=0.1, ms=0.1, color='green', 
         label='Alice, smoothed 30 bins = ' + repr(30*dt) + ' sec [green]' )
plt.plot(t, count_rate_fixed_300, linewidth=0.1, ms=0.1, color='blue', 
         label='Alice, smoothed 300 bins = ' + repr(300*dt) + ' sec [blue]')
plt.plot(t, count_rate_fixed_3000, linewidth=0.3, ms=0.1, color='red', 
         label='Alice, smoothed 3000 bins = ' + repr(3000*dt) + ' sec [red]')
plt.plot(t, count_rate_fixed_30000, linewidth=0.2, ms=0.1, color='yellow', 
         label='Alice, smoothed 30000 bins = ' + repr(30000 * dt) + ' sec [yellow]' )
plt.ylim((14,20))
plt.xlim(hbt.mm(t))
plt.title(sequence, fontsize=fs)
plt.ylabel('DN', fontsize=fs)
plt.xlabel('Seconds since ' + utc_start, fontsize=fs)
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