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

plt.rcParams['figure.figsize'] = 12, 12

lun = open(file_pickle, 'rb')
t = pickle.load(lun)
lun.close()
groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

index_group = 8
index_image = 10 # Which frame from the group do we extract?
index_image_stray = [49, 50, 51, 52, 53]

index_group = 7
index_image = 0 # Which frame from the group do we extract?
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

image_roll = np.roll(np.roll(image_processed, dx_total, axis=1), dy_total, axis=0)

radius  = planes['Radius_eq'] / rj   # Radius stored as km, but do the calc in R_J
azimuth = planes['Longitude_eq'] * hbt.r2d # Azimuth is stored as radians, but do the calc in degrees

print "radius[100,100] = " + repr(radius[100,100])

num_bins_azimuth = 10000 # Might make this a user parameter later
num_bins_radius  = 400

r_ring_inner_rj = 1.7
r_ring_outer_rj = 1.81

plt.rc('image', cmap='Greys_r')               # Default color table for imshow

fs = 15
plt.imshow( ( np.array(radius > r_ring_inner_rj) & np.array(radius < r_ring_outer_rj)) + 
              1* image_roll / np.max(image_roll))
plt.title(f_short, fontsize=fs)              
plt.show()

plt.imshow( ( np.array(radius > r_ring_inner_rj) & np.array(radius < r_ring_outer_rj)) + 
              20* image_roll / np.max(image_roll))
plt.title(f_short, fontsize=fs)              
plt.show()

# Extract and plot radial profile

bins_radius = hbt.frange(r_ring_inner_rj, r_ring_outer_rj, num_bins_radius)
bins_azimuth = hbt.frange(-180,180,num_bins_azimuth+1)
    
# Select the ring points

is_ring_all = ( np.array(radius > r_ring_inner_rj) & np.array(radius < r_ring_outer_rj))

#==============================================================================
# Create our final 'raw data' set: radius_all, azimuth_all, dn_all
#==============================================================================

radius_all  = planes['Radius_eq'][is_ring_all]
azimuth_all = planes['Longitude_eq'][is_ring_all]  # Make a list of all of the azimuth points for all pixels
dn_all      = image_roll[is_ring_all]

# Now take these raw data, and rearrange them so that we can take the longest continuous segment
# We do this by appending the timeseries to itself, looking for the largest gap (of no az data), 
# and then the data will start immediately after that.

# _2 indicates a double-length array (ie, with [azimuth, azimuth + 2pi])
# _s indicates sorted
# _d indicates delta

azimuth_all_2 = np.concatenate((azimuth_all, azimuth_all + 2*math.pi))
dn_all_2      = np.concatenate((dn_all, dn_all))
radius_all_2  = np.concatenate((radius_all, radius_all))

azimuth_all_2_s = np.sort(azimuth_all_2, kind = 'heapsort')
azimuth_all_2_s_d = azimuth_all_2_s - np.roll(azimuth_all_2_s, 1)
index_max_2_s_d = hbt.wheremax(azimuth_all_2_s_d)  # Search for the index of the largest gap.
                                                   # Just after this starts the data
                                              
azimuth_seg_start = azimuth_all_2_s[index_max_2_s_d] # Azimuth value at the segment start
azimuth_seg_end   = azimuth_seg_start + 2 * math.pi  # Azimuth value at the segment end XXX No that's not right!

indices_2_good = (azimuth_all_2 >= azimuth_seg_start) & (azimuth_all_2 < azimuth_seg_end)
azimuth_all_good = azimuth_all_2[indices_2_good]
radius_all_good  = radius_all_2[indices_2_good]
dn_all_good      = dn_all_2[indices_2_good]

azimuth_all = azimuth_all_good
radius_all  = radius_all_good
dn_all      = dn_all_good

plt.plot(azimuth_all_good,marker='.', linestyle='none')

#quit

grid_azimuth_1d = hbt.frange(0, 4. * math.pi, num_bins_azimuth)

# Create output array
# First we try using griddata() on the whole 2D array. 

#num_bins_azimuth = 1000

#az_arr  = np.linspace(np.min(azimuth_all),  np.max(azimuth_all),  num_bins_azimuth)
#rad_arr = np.linspace(np.min(radius_all), np.max(radius_all), num_bins_radius)
#
#dn_grid = griddata((radius_all, azimuth_all), dn_all, 
#                   (rad_arr[None,:], az_arr[:,None]), method='cubic')

# We should be able to use a single call to griddata() to grid the entire dataset. But it gives
# screwy results, obviously wrong. So instead, I am doing it line-by-line (in radial bins).

grid_lin_2d = np.zeros((num_bins_radius, num_bins_azimuth))

for i in range(num_bins_radius-1):  # Loop over radius -- inner to outer
    
    is_ring_i = np.array(radius > bins_radius[i]) & np.array(radius < bins_radius[i+1]) # Select only bins with right radius
    dn_i = image_roll[is_ring_i]
    radius_i = planes['Radius_eq'][is_ring_i]
    azimuth_i = planes['Longitude_eq'][is_ring_i]
    azimuth_i[azimuth_i < 0] += 2*math.pi
    grid_lin_i   = griddata(azimuth_i, dn_i, grid_azimuth_1d, method='linear')
    
    grid_lin_2d[i,:] = grid_lin_i

# Remove the NaN's

grid_lin_2d[np.isnan(grid_lin_2d)] = 0

# Figure out the proper xlimits

grid_lin_2d_extend = np.concatenate((grid_lin_2d, grid_lin_2d), axis=1)  # Stuff two copies together, in radial direction

# Clamp the result

grid_lin_2d = hbt.remove_brightest(grid_lin_2d, 0.97, symmetric=True)

profile_azimuth = np.sum(grid_lin_2d, 0)
profile_radius  = np.sum(grid_lin_2d, 1)

bin_az_min = np.where(profile_azimuth > 0)[0][0]
bin_az_max = np.where(profile_azimuth > 0)[0][-1]

az_min = grid_azimuth_1d[bin_az_min]
az_max = grid_azimuth_1d[bin_az_max]

# And make a plot of the unwrapped ring!
# Set the aspect ratio manually -- we do not really need square pixels here.

plt.imshow(grid_lin_2d, aspect=00.1)
plt.xlim((azimuth_seg_start, azimuth_seg_end))
plt.show()

plt.plot(grid_azimuth_1d, profile_azimuth, label='Azimuthal')
plt.xlim(az_min, az_max)
plt.plot(bins_radius, profile_radius, label='Radial')
plt.legend()
plt.show()



fs = 15

plt.rcParams['figure.figsize'] = 16, 10




stop


# Now I want to map these values of dn, radius, azimuth into my output array.
# Different ways I could do this:
#  1. For each data pixel, find the output image pixel it corresponds to.
#     Map it there. Then do interpolation on the output pixels with no input pix.
#
#  2. For each output pixel, find the input pixel it corresponds to.
#
# I need to keep track of how many pixels go where. For instance, on the far edges, I know
# one pixel might get stretched over 100. And in the core of the image, I want to make sure that we're not
# under-sampling the highest res.
#
# So, I want to create a backplane to the output image, which is basically 'how many input pixels 
# went into this output pixel.'
# 
# There must be some general term for this: re-gridding irregular data?

plt.plot(bins_radius, profile_radius)
fs = 20

#        self.ax2.set_xlim([0,100])  # This is an array and not a tuple. Beats me, like so many things with mpl.

plt.show()

# Extract and plot azimuthal profile

# Create the azimuthal bins. 


###    
# Grid the data in 2D. This causes no errors, but the result doesn't look right.
# Thus, I went to gridding data line-by-line in 1D instead, rather than a single call to 2D-grid it.
###

# Create the output grid. This is kind of funny but it is what sample code does
# This mpgrid function creates a pair of 2D grids: one with all the x values, and one with all the y.

#grid_radius, grid_azimuth = np.mgrid[0:1:(num_bins_radius)*1j, 0:1:(num_bins_azimuth)*1j]
#
#grid_radius *= 0.1*rj
#grid_radius += 1.7*rj
#grid_azimuth *= 4 * math.pi

#grid = griddata(np.transpose((azimuth_all, radius_all)), dn_all, (bins_azimuth, bins_radius), method='nearest')

#points = np.transpose((azimuth_all, radius_all))
#values = dn_all
#grid_x = grid_azimuth
#grid_y = grid_radius
#
#grid_lin_2d   = griddata(points, values, (grid_x, grid_y), method='linear')
