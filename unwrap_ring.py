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

groupmask = (t['Desc'] == 'Jupiter ring - search for embedded moons')
index_image = 49 # Which frame from the group do we extract?

t_group = t[groupmask]	

f = t_group['Filename'][index_image] # Look up filename

file_backplane = dir_backplanes + t_group['Shortname'][index_image].replace('.fit', '_planes.pkl')

# Load backplane
				
if (os.path.isfile(file_backplane)):
    print 'load_backplane: loading ' + file_backplane 
    lun = open(file_backplane, 'rb')
    planes = pickle.load(lun)
    planes = planes # This will load self.t
    lun.close()
    
# Load image

image = hbt.get_image_nh(t_group['Filename'][index_image])
image_processed = hbt.remove_brightest(image, 0.97, symmetric=True)
image_processed = image_processed - hbt.sfit(image_processed, 5)
image_processed = hbt.remove_brightest(image_processed, 0.97, symmetric=True)

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
r_ring_outer_rj = 1.8

plt.rc('image', cmap='Greys_r')               # Default color table for imshow

plt.imshow( ( np.array(radius > 1.7) & np.array(radius < 1.8)) + 1* image_roll / np.max(image_roll))
plt.show()

# Extract and plot radial profile

bins_radius = hbt.frange(1.7, 1.8, num_bins_radius)
bins_azimuth = hbt.frange(-180,180,num_bins_azimuth+1)
    
# Select the ring points

is_ring_all = ( np.array(radius > r_ring_inner_rj) & np.array(radius < r_ring_outer_rj))

radius_all  = planes['Radius_eq'][is_ring_all]
azimuth_all = planes['Longitude_eq'][is_ring_all]
dn_all      = image_roll[is_ring_all]

# Unwrap and fix the azimuth angles.

azimuth_all[azimuth_all < 0] += 2 * math.pi

# For diagnostic purposes, draw a line at 1.75 Rj

dr = 0.003
r_mid = 1.78

is_ring_mid =  ( np.array(radius > 1.78-dr) & np.array(radius < 1.78+dr))
image_roll[is_ring_mid] = np.max(image_roll)

#is_ring_azmid = 
# Create the output grid. This is kind of funny but it is what sample code does
# This mpgrid function creates a pair of 2D grids: one with all the x values, and one with all the y.

#grid_radius, grid_azimuth = np.mgrid[0:1:(num_bins_radius)*1j, 0:1:(num_bins_azimuth)*1j]
#
#grid_radius *= 0.1*rj
#grid_radius += 1.7*rj
#grid_azimuth *= 4 * math.pi

grid_azimuth_1d = hbt.frange(0, 4. * math.pi, num_bins_azimuth)

# Grid the data in 1D (in azimuth, single radius)

#grid_lin     = griddata(azimuth_i, dn_i, grid_azimuth_1d, method='linear')
#grid_lin_1   = griddata(azimuth_1, dn_1, grid_azimuth_1d, method='linear')
#grid_lin_2   = griddata(azimuth_2, dn_2, grid_azimuth_1d, method='linear')
#grid_lin_3   = griddata(azimuth_3, dn_3, grid_azimuth_1d, method='linear')

#grid_lin_123 = 
#plt.plot(azimuth_i, dn_i, label='ungridded')
#plt.plot(grid_azimuth_1d, grid_lin, label='gridded', color='red', marker='o' )
#plt.xlabel('Azimuth [rad]')
#plt.ylabel('DN')
#
#plt.legend()
#plt.show()

# Grid the data in 2D

#grid = griddata(np.transpose((azimuth_all, radius_all)), dn_all, (bins_azimuth, bins_radius), method='nearest')

#points = np.transpose((azimuth_all, radius_all))
#values = dn_all
#grid_x = grid_azimuth
#grid_y = grid_radius
#
#grid_lin_2d   = griddata(points, values, (grid_x, grid_y), method='linear')
#
#grid_lin_2d   = griddata(np.transpose((azimuth_all, radius_all)), dn_all, (grid_azimuth, grid_radius), method='linear')
#grid_near_2d  = griddata(np.transpose((azimuth_all, radius_all)), dn_all, (grid_azimuth, grid_radius), method='nearest')
#grid_cub_2d   = griddata(np.transpose((azimuth_all, radius_all)), dn_all, (grid_azimuth, grid_radius), method='cubic')

# Create output array

grid_lin_2d = np.zeros((num_bins_radius, num_bins_azimuth))

# We should be able to use a single call to griddata() to grid the entire dataset. But it gives
# screwy results, obviously wrong. So instead, I am doing it line-by-line (in radial bins).

for i in range(num_bins_radius-1):
    
    is_ring_i = np.array(radius > bins_radius[i]) & np.array(radius < bins_radius[i+1])
    dn_i = image_roll[is_ring_i]
    radius_i = planes['Radius_eq'][is_ring_i]
    azimuth_i = planes['Longitude_eq'][is_ring_i]
    azimuth_i[azimuth_i < 0] += 2*math.pi
    grid_lin_i   = griddata(azimuth_i, dn_i, grid_azimuth_1d, method='linear')
    
    grid_lin_2d[i,:] = grid_lin_i

# Remove the NaN's

grid_lin_2d[np.isnan(grid_lin_2d)] = 0

profile_azimuth = np.sum(grid_lin_2d, 0)
profile_radius  = np.sum(grid_lin_2d, 1)

bin_az_min = np.where(profile_azimuth > 0)[0][0]
bin_az_max = np.where(profile_azimuth > 0)[0][-1]

az_min = grid_azimuth_1d[bin_az_min]
az_max = grid_azimuth_1d[bin_az_max]

#plt.plot(grid_azimuth_1d, profile_azimuth)
#plt.xlim(az_min, az_max)
#plt.plot(bins_radius, profile_radius)
#plt.show()

# And make a plot of the unwrapped ring!

plt.rcParams['figure.figsize'] = 16, 10

plt.imshow(grid_lin_2d, 
           extent = [np.min(grid_azimuth_1d), np.max(grid_azimuth_1d), 
                     r_ring_outer_rj, r_ring_inner_rj])
           
plt.xlim((az_min, az_max))
plt.show()

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

        
#print
#for i in range(num_bins_azimuth-1):
#    
#    is_good = np.array(radius > 1.5)                  & np.array(radius < 2.2) & \
#              np.array(azimuth > self.bins_azimuth[i]) & np.array(azimuth < self.bins_azimuth[i+1])
#              
#    self.profile_azimuth[i] = np.mean(image_roll[is_good])
##            print "{:<3}. {:<7} {:<7}".format(i, self.bins_azimuth[i], self.profile_azimuth[i])
#
## Now we do some crazy logic to unwrap the azimuth, so that start will be in -180 .. 180, and end after that.
## First, copy the azimuth and intensity, so we have two full loops in the array (720 deg)
#
#profile_azimuth_2      = np.concatenate((self.profile_azimuth[:-2], self.profile_azimuth[:-2]))
#bins_azimuth_2         = np.concatenate((self.bins_azimuth[:-2], self.bins_azimuth[:-2]+360))     
#
## Now, search this for the first pattern of [nan, non-nan]. That will be the start of the valid data.
#
#i = np.array(range(np.size(profile_azimuth_2)-1)) # Just an index. Chop off last entry so we can use i+1
#
#is_start = np.isnan(profile_azimuth_2[i]) & np.logical_not(np.isnan(profile_azimuth_2[i+1]))
#bin_start = i[is_start][0]+1 # get first match
#
## And then look for the end pattern: [non-nan, nan], starting *after* bin_start
#
#is_end = np.logical_not(np.isnan(profile_azimuth_2[i])) & 
        # np.isnan(profile_azimuth_2[i+1]) & np.array(i > bin_start)
#bin_end = i[is_end][0]
#
## Now we have the proper indices for the start and end.
# 
#az_start = bins_azimuth_2[bin_start]
#az_end   = bins_azimuth_2[bin_end]
#
#print "az = " + repr(bins_azimuth_2[bin_start]) + ' .. ' + repr(bins_azimuth_2[bin_end])
#
#self.ax3.plot(bins_azimuth_2, profile_azimuth_2)
#fs = 20
#
##        self.ax3.set_xlim([az_start-10, az_end+10])  
        # This is an array and not a tuple. Beats me, like so many things with mpl.
##        self.ax3.set_ylim([-5, 5])           # Hard-code this.... not sure what is best.
#self.canvas3.show()
#
#plot4 = self.ax4.imshow(azimuth)
##        plt.title('Lon_eq [deg]', fontsize=fs)
##        plt.xlim((0,1000))
##        plt.ylim((0,1000))
#
##        self.fig4.colorbar(plot4)  # Need to get this to updated itself. Not sure how.
#self.canvas4.show()
    
