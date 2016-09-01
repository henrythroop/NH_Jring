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
import matplotlib.pyplot as plt
import numpy as np
import os.path
import astropy
from scipy.interpolate import griddata
import math

file_pickle = 'nh_jring_read_params_571.pkl' # Filename to read to get filenames, etc.
dir_images =         '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'
dir_backplanes =     '/Users/throop/data/NH_Jring/out/'
rj = 71492 # km

# Define the size of the output array

num_bins_azimuth = 1000 
num_bins_radius  = 400

plt.rcParams['figure.figsize'] = 10,10

lun = open(file_pickle, 'rb')
t = pickle.load(lun)
lun.close()
groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

index_group = 8
index_image = 0 # Which frame from the group do we extract?
index_images_stray = hbt.frange(0,48)

index_group = 7
index_image = 42 # Which frame from the group do we extract?
index_images_stray = hbt.frange(1,2)

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
				
if (os.path.isfile(file_backplane)):
    print 'load_backplane: loading ' + file_backplane 
    lun = open(file_backplane, 'rb')
    planes = pickle.load(lun)
    planes = planes # This will load self.t
    lun.close()
    
# Load image
# Then do basic image processing on it to remove stray light, etc.
# To remove stray light, we first remove low-frq terms using an sfit. Then we remove hi-frq terms with the stray light image.
# This is my standard cookbook routine
# One change though: I should make my image clamping (remove_brightest) only apply to the plotted image -- not the internal one
    
image = hbt.read_lorri(t_group['Filename'][index_image])
image_processed = hbt.remove_brightest(image, 0.97, symmetric=True)
image_processed = image_processed - hbt.sfit(image_processed, 5)
image_processed = hbt.remove_brightest(image_processed, 0.97, symmetric=True)
image_processed -= image_stray 

plt.imshow(image_processed)
												
dx_total =  -(t_group['dx_offset'][index_image] + t_group['dx_opnav'][index_image])
dy_total =  -(t_group['dy_offset'][index_image] + t_group['dy_opnav'][index_image])

# Roll the image as per the saved navigation offset values

image_roll = np.roll(np.roll(image_processed, dx_total, axis=1), dy_total, axis=0)

# Read in values from the backplane

radius  = planes['Radius_eq']    # Radius in km
azimuth = planes['Longitude_eq'] # Azimuth in radians

r_ring_inner = 1.7 * rj
r_ring_outer = 1.81 * rj

plt.rc('image', cmap='Greys_r')               # Default color table for imshow

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

#==============================================================================
# Extract radial and azimuthal profiles
#==============================================================================

profile_azimuth = np.nansum(dn_grid, 0)
profile_radius  = np.nansum(dn_grid, 1)

profile_azimuth_2 = np.nansum(dn_grid_2, 0)
profile_radius_2  = np.nansum(dn_grid_2, 1)

#==============================================================================
#  Plot the remapped 2D images
#==============================================================================

plt.rcParams['figure.figsize'] = 15,15

extent = [azimuth_seg_start, azimuth_seg_end, np.min(radius_all),np.max(radius_all)]

# Method #1 -- using griddata() over the whole array monolithically

plt.subplot(1,2,1)
plt.imshow(dn_grid, extent=extent, aspect=aspect, vmin=-15, vmax=20, origin='lower')
plt.title('griddata [2d]')
plt.xlabel('Azimuth [radians]')
plt.ylabel('Radius [km]')

# Method #2 -- using griddata() over one line at a time

plt.subplot(1,2,2)
plt.imshow(hbt.remove_brightest(dn_grid_2, 0.95, symmetric=True), aspect=aspect, vmin=-10, vmax=20, extent=extent, origin='lower')
plt.title('griddata [line-by-line]')
plt.xlabel('Azimuth [radians]')
plt.ylabel('Radius [km]')

plt.show()

#==============================================================================
# Plot the radial and azimuthal profiles
#==============================================================================

plt.rcParams['figure.figsize'] = 15,6  # <- this is in dx, dy... which is opposite from array order!

plt.subplot(1,2,1)
plt.plot(bins_azimuth, profile_azimuth, label = '1D')
plt.plot(bins_azimuth, profile_azimuth_2+400, label='2D')
plt.title('Azimuthal Profile')
plt.xlabel('Azimuth [radians]')
plt.legend()

plt.subplot(1,2,2)
plt.plot(bins_radius, profile_radius, label = '1D')
plt.plot(bins_radius, profile_radius_2 + 800, label='2D')
plt.xlabel('Radius [km]')
plt.title('Radial Profile')
plt.xlim(hbt.mm(bins_radius))
plt.legend()
plt.show()
