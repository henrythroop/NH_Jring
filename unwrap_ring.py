# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 09:01:30 2016

@author: throop
"""

##########
# Unwrap the ring image, based on coordinates provided by the backplane
##########

# This is just a test development routine... not the real function
# HBT 28-Jun-2016

import hbt
import pickle
from   astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os.path
import astropy
from scipy.interpolate import griddata
import math
import spiceypy as sp

#==============================================================================
# Function for DN to I/F conversion
#==============================================================================

def dn2iof(dn, exptime, pixfov, rsolar):
    
# Convert DN to I/F. 
# This is not a general routine -- it's specific to LORRI @ Jupiter.

    # Calculate LORRI pixel size, in sr

    sr_pix       = (pixfov*(1e-6))**2  # Angular size of each MVIC pixel, in sr

    # Look up the solar flux at 1 AU, using the solar spectral irradiance at 
    #    http://rredc.nrel.gov/solar/spectra/am1.5/astmg173/astmg173.html
    # or http://www.pas.rochester.edu/~emamajek/AST453/AST453_stellarprops.pdf

    f_solar_1au_si     = 1.77                 # W/m2/nm. At 600 nm. 
    f_solar_1au_cgs    = f_solar_1au_si * 100 # Convert from W/m2/nm to erg/sec/cm2/Angstrom

    f_solar_1au        = f_solar_1au_cgs

    # Calculate the solar flux at Pluto's distance. [Actual distance is 33.8 AU]

    f_solar_jup       = f_solar_1au / (5**2)  # Use 1/r^2, not 1/(4pi r^2)

    # Use constants (from Level-2 files) to convert from DN, to intensity at detector.

    # This is the equation in ICD @ 62.

    i_per_sr    = dn / exptime / rsolar # Convert into erg/cm2/s/sr/Angstrom.

    # Because the ring is spatially extended, it fills the pixel, so we mult by pixel size
    # to get the full irradiance on the pixel.
    # *** But, somehow, this is not working. I get the right answer only if I ignore this 'per sr' factor.

    DO_OVERRIDE = True  # If true, ignore the missing factor of 'per sr' that seems to be in the conversion

    i           = i_per_sr * sr_pix     # Convert into erg/cm2/s/Angstrom

    if (DO_OVERRIDE):
        i = i_per_sr

    # Calculate "I/F". This is not simply the ratio of I over F, because F in the eq is not actually Flux.
    # "pi F is the incident solar flux density" -- ie, "F = solar flux density / pi"
    # SC93 @ 125

    iof = i / (f_solar_jup / math.pi)
    
    return iof

#==============================================================================
# Define input quantities and constants
#==============================================================================
    
file_pickle       = 'nh_jring_read_params_571.pkl' # Filename to read to get filenames, etc.
file_tm           = "/Users/throop/git/NH_Rings/kernels_nh_jupiter.tm"  # SPICE metakernel
dir_images        = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'
dir_backplanes    = '/Users/throop/data/NH_Jring/out/'
dir_out           = dir_backplanes

rj                = 71492           # Jupiter radius, km

abcorr            = 'LT+S'
frame             = 'J2000'

stretch_percent   = 95

stretch = astropy.visualization.PercentileInterval(stretch_percent)  # PI(90) scales array to 5th .. 95th %ile. 

# Define the size of the output array

num_bins_azimuth = 1000 
num_bins_radius  = 400

plt.rcParams['figure.figsize'] = 10,10

lun = open(dir_out + file_pickle, 'rb')
t = pickle.load(lun)
lun.close()
groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

index_group = 5
index_image = 5 # Which frame from the group do we extract?
index_images_stray = range(1)

#index_group = 7
#index_image = 42 # Which frame from the group do we extract?
#index_images_stray = hbt.frange(1,2)

#index_group = 5
#index_image = 1
#index_images_stray = hbt.frange(1,6,6)

groupmask = (t['Desc'] == groups[index_group])
t_group = t[groupmask]

# Load stray light (or generate if this one has not been made yet)

image_stray = hbt.nh_get_straylight_median(index_group, index_images_stray, do_sfit=True, power1=5, power2=5)

f       = t_group['Filename'][index_image] # Look up filename
f_short = t_group['Shortname'][index_image]

file_backplane = dir_backplanes + t_group['Shortname'][index_image].replace('.fit', '_planes.pkl')

# Load backplane, which is required for navigation
# ** For the MVIC observations of the Pluto rings, I put the backplanes into the FITS files,
# which had been created by Tod.
# But for Jupiter, since we have hundreds of images, I am not going to modify them. Safer to 
# analyze the virgin FITS files from the SOC, and use unrelated backplane files. 
# Also, since the Jupiter FITS file have all the SOC header info in them, there is really less
# need to mess with them than with Tod's 
				
if (os.path.isfile(file_backplane)):
    print('load_backplane: loading ' + file_backplane) 
    lun = open(file_backplane, 'rb')
    planes = pickle.load(lun)
    planes = planes # This will load self.t
    lun.close()
    
# Load image
# Then do basic image processing on it to remove stray light, etc.
# To remove stray light, we first remove low-frq terms using an sfit. 
# Then we remove hi-frq terms with the stray light image.
# This is my standard cookbook routine
# One change though: I should make my image clamping (remove_brightest) 
# only apply to the plotted image -- not the internal one
    
image = hbt.read_lorri(t_group['Filename'][index_image])
image = hbt.lorri_destripe(image)
image_processed = hbt.remove_brightest(image, 0.97, symmetric=True)
image_processed = image_processed - hbt.sfit(image_processed, 5)
image_processed = hbt.remove_brightest(image_processed, 0.97, symmetric=True)

DO_REMOVE_STRAY = False

if (DO_REMOVE_STRAY):
    image_processed -= image_stray 

plt.imshow(stretch(image_processed))
plt.title(f_short)
plt.show()

# Load some fields from the FITS image header

hdulist = fits.open(t_group['Filename'][index_image])
exptime = hdulist['PRIMARY'].header['EXPTIME']
rsolar  = hdulist['PRIMARY'].header['RSOLAR']
pixfov = 0.3 * hbt.d2r * 1e6 / 1024  # This keyword is missing. "Plate scale in microrad/pix "
et      = hdulist['PRIMARY'].header['SPCSCET']

# Roll the image as per the saved navigation offset values
												
#dx_total =  -(t_group['dx_offset'][index_image] + t_group['dx_opnav'][index_image])
#dy_total =  -(t_group['dy_offset'][index_image] + t_group['dy_opnav'][index_image])

#image_roll = np.roll(np.roll(image_processed, dx_total, axis=1), dy_total, axis=0)

# Read in values from the backplane

radius  = planes['Radius_eq']    # Radius in km
azimuth = planes['Longitude_eq'] # Azimuth in radians
phase   = planes['Phase']        # Phase angle 

r_ring_inner = 1.6 * rj
r_ring_outer = 1.9 * rj

plt.rc('image', cmap='Greys_r')               # Default color table for imshow

#==============================================================================
# Start up SPICE and extract some geometrical quantities
#==============================================================================

sp.furnsh(file_tm)
utc     = sp.et2utc(et, 'C', 0)

# Look up the ring tilt angle (ie, observer latitude)

(vec, lt)        = sp.spkezr('New Horizons', et, 'IAU_JUPITER', abcorr, 'Jupiter')
(junk, lon, lat) = sp.reclat(vec[0:3])
elev             = np.abs(lat)          # Elevation angle (aka 'B')  
emis             = math.pi/2 - elev     # Emission angle (ie, angle down from normal) 
mu               = abs(math.cos(emis))  # mu. See definitions of all these Throop 2004 @ 63 

#==============================================================================
# Make a plot with the ring and boundary shown
#==============================================================================

plt.rcParams['figure.figsize'] = 15,15

fs = 15
scaling = 5  # How much to mulitply two image by to superimpose them
image_mask = ( np.array(radius > r_ring_inner) & np.array(radius < r_ring_outer))
image_ring_filtered = hbt.remove_brightest(image_roll,0.99,symmetric=True) / np.max(image_roll)

plt.subplot(1,3,1)
plt.imshow(                      image_ring_filtered)

plt.subplot(1,3,2)
plt.imshow( image_mask + 0      *image_ring_filtered)
plt.title(f_short, fontsize=fs)              

plt.subplot(1,3,3)
plt.imshow( image_mask + scaling*image_ring_filtered)

plt.show()

#==============================================================================
# Examine backplane to figure out azimuthal limits of the ring image
#==============================================================================

# Select the ring points -- that is, everything inside the mask

is_ring_all = ( np.array(radius > r_ring_inner) & np.array(radius < r_ring_outer))

radius_all  = planes['Radius_eq'][is_ring_all]     # Make a list of all of the radius values
azimuth_all = planes['Longitude_eq'][is_ring_all]  # Make a list of all of the azimuth points for all pixels
dn_all      = image_roll[is_ring_all]              # DN values, from the rolled image

# Now take these raw data, and rearrange them so that we can take the longest continuous segment
# We do this by appending the timeseries to itself, looking for the largest gap (of no az data), 
# and then the data will start immediately after that.

# _2 indicates a double-length array (ie, with [azimuth, azimuth + 2pi])
# _s indicates sorted
# _d indicates delta

azimuth_all_3 = np.concatenate((azimuth_all, azimuth_all + 2*math.pi, azimuth_all + 4*math.pi))
dn_all_3      = np.concatenate((dn_all, dn_all, dn_all))
radius_all_3  = np.concatenate((radius_all, radius_all, radius_all))

azimuth_all_3_s = np.sort(azimuth_all_3, kind = 'heapsort')
azimuth_all_3_s_d = azimuth_all_3_s - np.roll(azimuth_all_3_s, 1)

# Look for the indices where the largest gaps (in azimuth) start

index_seg_start_3_s = (np.where(azimuth_all_3_s_d > 0.999* np.max(azimuth_all_3_s_d)))[0][0]
index_seg_end_3_s = (np.where(azimuth_all_3_s_d > 0.999* np.max(azimuth_all_3_s_d)))[0][1]-1

# Get proper azimithal limits. We want them to be a single clump of monotonic points.
# Initial point is in [0, 2pi) and values increase from there.
                                                   
azimuth_seg_start = azimuth_all_3_s[index_seg_start_3_s] # Azimuth value at the segment start
azimuth_seg_end   = azimuth_all_3_s[index_seg_end_3_s]   # Azimuth value at the segment end

indices_3_good = (azimuth_all_3 >= azimuth_seg_start) & (azimuth_all_3 < azimuth_seg_end)

azimuth_all_good = azimuth_all_3[indices_3_good]
radius_all_good  = radius_all_3[indices_3_good]
dn_all_good      = dn_all_3[indices_3_good]

# Extract arrays with the proper pixel values, and proper azimuthal values

azimuth_all = azimuth_all_good
radius_all  = radius_all_good
dn_all      = dn_all_good

#stop

#==============================================================================
#  Now regrid the data from xy position, to an unrolled map in (azimuth, radius)
#==============================================================================

# Method #1: Construct the gridded image all at once 

az_arr  = np.linspace(azimuth_seg_start, azimuth_seg_end,     num_bins_azimuth)
rad_arr = np.linspace(np.min(radius_all), np.max(radius_all), num_bins_radius)

# In order for griddata to work, d_radius and d_azimuth values should be basically equal.
# In Jupiter data, they are not at all. Easy soln: tweak azimuth values by multiplying
# by a constant to make them large.

f = (np.max(radius_all) - np.min(radius_all)) / (np.max(azimuth_all) - np.min(azimuth_all))
aspect = 1/f

# NB: griddata() returns an array in opposite row,column order than I would naively think.
# I want results in order (radius, azimuth) = (y, x)

dn_grid = griddata((f*azimuth_all, radius_all), dn_all, 
                   (f*az_arr[None,:], rad_arr[:,None]), method='linear')



# Method #2: Construct the gridded image line-by-line

dn_grid_2 = np.zeros((num_bins_radius, num_bins_azimuth))  # Row, column
bins_azimuth    = hbt.frange(azimuth_seg_start, azimuth_seg_end, num_bins_azimuth)
bins_radius     = hbt.frange(r_ring_inner, r_ring_outer, num_bins_radius)

for i in range(num_bins_radius-1):  # Loop over radius -- inner to outer
    
    # Select only bins with right radius and azimuth
    is_ring_i = np.array(radius_all > bins_radius[i]) & np.array(radius_all < bins_radius[i+1]) & \
                np.array(azimuth_all > azimuth_seg_start) & np.array(azimuth_all < azimuth_seg_end) 
    
    if np.sum(is_ring_i) > 0:
        dn_i = dn_all[is_ring_i]  # Get the DN values from the image (adjusted by navigation position error)
        radius_i = radius_all[is_ring_i]
        azimuth_i = azimuth_all[is_ring_i]
        grid_lin_i   = griddata(azimuth_i, dn_i, bins_azimuth, method='linear')
        
        dn_grid_2[i,:] = grid_lin_i

# Set limits for where the extraction should happen

# For the radial profile, we take only a portion of the unwrapped frame. 
# e.g., the inner 30%. We exclude the region on the edges, since it is smeared, and has contamination.
# We are not starved for photons. Take the best portion of the signal and use it.

frac_profile_radial = 0.3  # Of the available azimuth range, what fraction do we use for extracting radial profile?

# For azimuthal profile, focus on the main ring. Don't focus on the diffuse inner region.
# It is harder to do that photometry, more artifacts, fainter, and probalby more spread out anyhow.

# Define distances of [outer box, outer ring, inner ring, inner box] in km

limits_profile_azimuth = np.array([131e3,130e3,127e3,126e3]) 

limits_profile_azimuth_bins = limits_profile_azimuth.astype(int).copy() * 0
for i,r in enumerate(limits_profile_azimuth):
    limits_profile_azimuth_bins[i] = int(hbt.wheremin(abs(bins_radius - r)))

    limits_profile_radial_bins = int(np.shape(dn_grid)[1]) * \
      np.array([0.5-frac_profile_radial/2, 0.5+frac_profile_radial/2])

limits_profile_radial_bins = limits_profile_radial_bins.astype(int)
      
#==============================================================================
# Extract radial and azimuthal profiles, using entire reprojected image
#==============================================================================

profile_azimuth = np.nanmean(dn_grid, axis=0)
profile_radius  = np.nanmean(dn_grid, axis=1)

profile_azimuth_2 = np.nanmean(dn_grid_2, axis=0)
profile_radius_2  = np.nanmean(dn_grid_2, axis=1)

#==============================================================================
# Extract radial and azimuthal profiles, using subsections
#==============================================================================

# I am using the *mean* here along each row and column. That means that the final value
# in the profiles is 
# Azimuthal

#profile_azimuth_bg_inner = np.nansum(dn_grid[limits_profile_azimuth_bins[1]:limits_profile_azimuth_bins[0],:],0)

plt.rcParams['figure.figsize'] = 10,5

profile_azimuth_bg_inner = np.nanmean(dn_grid[limits_profile_azimuth_bins[1]:limits_profile_azimuth_bins[0],:],axis=0)
profile_azimuth_core     = np.nanmean(dn_grid[limits_profile_azimuth_bins[2]:limits_profile_azimuth_bins[1],:],axis=0)
profile_azimuth_bg_outer = np.nanmean(dn_grid[limits_profile_azimuth_bins[3]:limits_profile_azimuth_bins[2],:],axis=0)

# Get profile in DN
profile_azimuth_subtracted = profile_azimuth_core - (profile_azimuth_bg_inner + profile_azimuth_bg_outer)/2
profile_radius_central   = np.nanmean(dn_grid[:,limits_profile_radial_bins[0]:limits_profile_radial_bins[1]],1)

# Convert profile to I/F

profile_azimuth_subtracted_iof = dn2iof(profile_azimuth_subtracted, exptime, pixfov, rsolar)
profile_radius_central_iof     = dn2iof(profile_radius_central,  exptime, pixfov, rsolar)

# Convert profile to normalized I/F (ie, seen from above)

profile_radius_central_iof_norm = profile_radius_central_iof * mu

plt.plot(bins_azimuth, profile_azimuth_bg_inner, label=\
  'Background: Rows {} .. {}.'.format(limits_profile_azimuth_bins[1], limits_profile_azimuth_bins[0]))
plt.plot(bins_azimuth, profile_azimuth_core,     label='Core')
plt.plot(bins_azimuth, profile_azimuth_bg_outer, label=\
  'Background: Rows {} .. {}.'.format(limits_profile_azimuth_bins[1], limits_profile_azimuth_bins[0]))
plt.plot(bins_azimuth, profile_azimuth_subtracted + 4, color='black', label = 'Subtracted')
plt.legend()
plt.show()

#==============================================================================
#  Plot the remapped 2D images, with extraction boxes
#==============================================================================

plt.rcParams['figure.figsize'] = 10,5

extent = [azimuth_seg_start, azimuth_seg_end, np.min(radius_all),np.max(radius_all)]

# Method #1 -- using griddata() over the whole array monolithically

plt.subplot(1,2,1)
plt.imshow(dn_grid, extent=extent, aspect=aspect, vmin=-15, vmax=20, origin='lower')
plt.title('griddata [2d]')
plt.xlabel('Azimuth [radians]')
plt.ylabel('Radius [km]')

plt.hlines(limits_profile_azimuth[0], -10, 10, color='purple')
plt.hlines(limits_profile_azimuth[1], -10, 10, color='purple')
plt.hlines(limits_profile_azimuth[2], -10, 10, color='purple')
plt.hlines(limits_profile_azimuth[3], -10, 10, color='purple')
plt.xlim(hbt.mm(bins_azimuth))

plt.vlines(bins_azimuth[limits_profile_radial_bins[0]],-1e10, 1e10)
plt.vlines(bins_azimuth[limits_profile_radial_bins[1]],-1e10, 1e10)
plt.ylim(hbt.mm(bins_radius))

# Method #2 -- using griddata() over one line at a time

plt.subplot(1,2,2)
plt.imshow(hbt.remove_brightest(dn_grid_2, 0.95, symmetric=True), aspect=aspect, vmin=-10, vmax=20, \
           extent=extent, origin='lower')
plt.title('griddata [line-by-line]')
plt.xlabel('Azimuth [radians]')
plt.ylabel('Radius [km]')

plt.hlines(limits_profile_azimuth[0], -10, 10, color='purple')
plt.hlines(limits_profile_azimuth[1], -10, 10, color='purple')
plt.hlines(limits_profile_azimuth[2], -10, 10, color='purple')
plt.hlines(limits_profile_azimuth[3], -10, 10, color='purple')
plt.xlim(hbt.mm(bins_azimuth))

plt.vlines(bins_azimuth[limits_profile_radial_bins[0]],-1e10, 1e10)
plt.vlines(bins_azimuth[limits_profile_radial_bins[1]],-1e10, 1e10)
plt.ylim(hbt.mm(bins_radius))

plt.show()

#==============================================================================
# Plot the radial and azimuthal profiles
#==============================================================================

plt.rcParams['figure.figsize'] = 15,5  # <- this is in dx, dy... which is opposite from array order!

plt.subplot(1,2,1)
plt.plot(bins_azimuth, profile_azimuth, label = '1D')
plt.plot(bins_azimuth, profile_azimuth_2+4, label='2D')
plt.plot(bins_azimuth, profile_azimuth_subtracted, color='yellow', label='bg-sub')
plt.plot(bins_azimuth, hbt.smooth_boxcar(profile_azimuth_subtracted,30), label='bg-sub smooth', color='purple')

plt.title('Azimuthal Profile')
plt.xlabel('Azimuth [radians]')
plt.legend()

plt.subplot(1,2,2)
plt.plot(bins_radius, profile_radius, label = '1D')
plt.plot(bins_radius, profile_radius_2 + 1, label='2D')
plt.plot(bins_radius, profile_radius_central + 2, label='Small box')
plt.xlabel('Radius [km]')
plt.title('Radial Profile')
plt.xlim(hbt.mm(bins_radius))
plt.legend()
plt.show()

#
plt.rcParams['figure.figsize'] = 15,5  # <- this is in dx, dy... which is opposite from array order!

scalefac = 1e6

plt.subplot(1,2,1)
plt.plot(bins_azimuth, scalefac * profile_azimuth_subtracted_iof, color='yellow', label='bg-sub')
plt.plot(bins_azimuth, scalefac * hbt.smooth_boxcar(profile_azimuth_subtracted_iof,30), \
         label='bg-sub smooth', color='purple')

plt.title('Azimuthal Profile')
plt.xlabel('Azimuth [radians]')
plt.ylabel('I/F * {:.0e}'.format(scalefac))
plt.legend()

plt.subplot(1,2,2)
plt.plot(bins_radius, scalefac * profile_radius_central_iof, label='Small box')
plt.xlabel('Radius [km]')
plt.title('Radial Profile')
plt.ylabel('I/F * {:.0e}'.format(scalefac))
plt.xlim(hbt.mm(bins_radius))
plt.legend()
plt.show()


#==============================================================================
# Now measure the I/F, EW, etc. from the radial profile
#==============================================================================

# First convert DN to I/F

# Then integrate over I/F dr to get EW

# Set up an array with the bin width, in km

# List of all of the quantities I want to get:
# DN as function of radius
# I/F as function of radius
# Normal I/F
# Radially averaged normal I/F
# tau omega0 P ('taupp')

# EW = \int( (I/F)_norm dr)   [Throop 2004 @ 70]
# Radial Ave EW = EW / Dr.  Dr = 6500 km  [Throop 2004 @ 70]
# (I/F)_norm = (I/F) * mu [ie, corrected for ring opening angle] -- [Throop 2004 @ 2]

# Calc ew

width_ring = 6500 # in km. This is a constant to normalize by -- not the same as actual boundaries.

ew_edge     = [120000,130000]  # Integrate over this range. This is wider than the official ring width.
ew_edge_bin = hbt.x2bin(ew_edge,bins_radius)
    
dr = np.roll(bins_radius,-1) - bins_radius # Bin width, in km
dr[-1] = 0

iof_norm = profile_radius_central_iof_norm  # Just a shorter alias
ew_norm  = np.sum((iof_norm * dr)[ew_edge_bin[0] : ew_edge_bin[1]]) # Normalized EW (ie, from top)
ew_mean  = ew_norm / 6500                                           # Mean normalized EW

taupp    = iof_norm * 4 * mu  # Radially averaged tau_omega0_P

plt.plot(bins_radius, iof_norm) # Compare to Throop 2004 @ 65
plt.title('Phase = {:.2f}'.format(np.mean(phase) * hbt.r2d))

print("Ring mean EW = {}".format())
#print "DN = {}; I/F = {:.2g} .".format(dn_ring, iof_ring)


ew = np.sum((iof_ring * dr)[ew_edge_bin[0] : ew_edge_bin[1]])

stop



