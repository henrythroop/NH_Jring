#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:51:59 2017

Program to analyze rings data

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
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling

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

from astropy.convolution import Box1DKernel, Gaussian1DKernel, convolve

from   scipy.stats import linregress

import re # Regexp
import pickle # For load/save

from   matplotlib.figure import Figure

# HBT imports

import hbt

# Local NH rings imports

from  nh_jring_mask_from_objectlist import nh_jring_mask_from_objectlist

from nh_jring_mask_from_objectlist             import nh_jring_mask_from_objectlist
from nh_jring_unwrap_ring_image                import nh_jring_unwrap_ring_image

# For fun, I'm going to try to do this whole thing as a class. That makes it easier to pass 
# params back and forth, I guess.

class ring_profile:
    
# =============================================================================
#     Init method: load the main table file
# =============================================================================
    
    def __init__(self):

        file_pickle = 'nh_jring_read_params_571.pkl' # Filename to read to get filenames, etc.
        dir_out     = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
            
        lun = open(dir_out + file_pickle, 'rb')
        self.t = pickle.load(lun)
        lun.close()
        
        # Process the group names. Some of this is duplicated logic -- depends on how we want to use it.
        
        self.groups = astropy.table.unique(self.t, keys=(['Desc']))['Desc']

        stretch_percent = 90    
        self.stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

# =============================================================================
# Lookup the analysis filename
# =============================================================================

    def get_export_analysis_filename(self, index_group = None, index_image = None):

        if (index_image is None):
            index_image = self.index_image  # Use the current image, unless one is passed

        if (index_group is None):
            index_group = self.index_group

#        else:
                           # Use the passed-in image name
        
        groupmask = self.t['Desc'] == self.groups[index_group]
        t_group = self.t[groupmask]  # 
        
        t = t_group[index_image]  # Grab this, read-only, since we use it a lot.

        dir_export = '/Users/throop/data/NH_Jring/out/'
        file_export = dir_export + t['Shortname'].replace('_opnav', '').replace('.fit', '_analysis.pkl')

        return(file_export)
        
# =============================================================================
#     Read the profiles from disk
# =============================================================================
    
    def load(self, index_group, index_images, 
                      key_radius = 'core',  # Of the radial profiles, which one to read?
                      key_azimuth = 'net',  # Of the azimuthal profiles, which one to read?
                      **kwargs):

        # Each file on disk has several different extractions:
        #
        #   Radial profile has 'core', 'full', 'half, etc.
        #   Az     profile has 'inner', 'outer', 'net', etc. 
        #
        # Since we usually need just one of these, this routine *only* reads one.
        # The individual one read can be controlled with key_radius and key_azimuth
        
        import humanize
        
        # Define t_group, since we'll reference it a lot
        
        self.groupmask = self.t['Desc'] == self.groups[index_group]
        self.t_group = self.t[self.groupmask]  # 
        
        self.num_images_group = np.size(self.t_group)
        
        self.index_group            = index_group

        self.profile_radius_dn_arr  = []
        self.profile_azimuth_dn_arr = []
        self.exptime_arr            = []
        self.ang_phase_arr          = []
        self.radius_arr             = []
        self.azimuth_arr            = []
        self.index_image_arr        = []
        self.index_group_arr        = []
        self.ang_phase_arr          = []
        self.dt_arr                 = []
        self.dt_str_arr             = []

        for index_image in index_images:
            file = self.get_export_analysis_filename(index_group, index_image)
            print("{}".format(file))
            
            lun = open(file, 'rb')
            vals = pickle.load(lun)
            lun.close()
        
            (image_unwrapped,     # Unwrapped image itself
                        mask_unwrapped,      # Boolean mask
                        radius,              # Axis values for the unwrapped image
                        azimuth,             # Axis values for the unwrapped image 
                        profile_radius_dn,   # Radial profile (several, in a dictionary)
                        profile_azimuth_dn,  # Az profile (several, in a dictionary)
                        range_of_azimuth,
                        range_of_radius,
                        exptime,             # Exposure time
                        et,                  # ET
                        ang_elev,            # Elevation angle above ring
                        ang_phase,           # Phase angle (mean to rings -- not to planet center)
                        bg_method,
                        bg_argument,
                        index_image,         # Index of image
                        index_group) = vals  # Index of image group

            # Figure out when the analysis file was written, and turn into human readable form ('10 minute ago')
            
            time_file = os.path.getmtime(file)
            time_now  = time.time()
            dt        = time_now - time_file
            dt_str    = humanize.naturaltime(dt) 
            
            file_short = file.split('/')[-1].replace('_analysis.pkl', '')
            
            # Save all these in the ring object. This makes it easy to export them or use in other functions.
            
            self.profile_radius_dn_arr.append(profile_radius_dn[key_radius])  # This is now (eg) 8 x 30 array
            self.profile_azimuth_dn_arr.append(profile_azimuth_dn[key_azimuth])
            self.exptime_arr.append(exptime)
            self.ang_phase_arr.append(ang_phase)
            self.radius_arr.append(radius)
            self.azimuth_arr.append(azimuth)
            self.index_image_arr.append(index_image)
            self.index_group_arr.append(index_group)
            self.dt_arr.append(dt)
            self.dt_str_arr.append(dt_str)

            print("Read image {}/{}, phase = {:0.1f}°, {}, {}, {}".format(
                    index_group, index_image, ang_phase*hbt.r2d, file_short, bg_method, bg_argument))
        
        #==============================================================================
        # Now put these into arrays (not lists). Ideally we'd put these into an astropy table (so we can sort, etc.)
        #==============================================================================
        
        self.ang_phase_arr = np.array(self.ang_phase_arr)
        self.azimuth_arr   = np.array(self.azimuth_arr)
                   # Convert some of these from lists, to NP arrays
            
        self.profile_radius_dn_arr = np.array(self.profile_radius_dn_arr)
        self.profile_azimuth_dn_arr = np.array(self.profile_azimuth_dn_arr) 
        
        return self

# =============================================================================
# Return number of profiles
# =============================================================================
        
    def num_profiles(self):
        
        return np.size(self.dt_arr)
    
# =============================================================================
# Convert from DN to I/F. Or at least create the I/F curves.
# =============================================================================

      
# =============================================================================
# Smooth the profiles
# This smooths profiles in place. It doesn't return a new object. It destroys the 
# current state, overwriting it with a new state.
# =============================================================================

    def smooth(self, width=1, kernel=None):

        if (kernel is None):
            kernel = Gaussian1DKernel(width)
                
        num_profiles = self.num_profiles()

        for i in range(num_profiles):
        
            # First do radial kernels

            self.profile_radius_dn_arr[i,:]  = convolve(self.profile_radius_dn_arr[i,:],  kernel)

            # Then do azimuthal kernels
            
            self.profile_azimuth_dn_arr[i,:] = convolve(self.profile_azimuth_dn_arr[i,:], kernel)

        return self
    
#            for key in self.profile_azimuth_dn_arr[i]:
#                profile = convolve(self.profile_azimuth_dn_arr[i][key], kernel)
#                self.profile_azimuth_dn_arr[i][key] = profile
                
# =============================================================================
# Sum the profiles
# =============================================================================

    def sum(self):
        
        # Add up all of the radial profiles, to get a merged radial profile
        # Because of how this is implemented, we can't just do an np.sum(). 
        # Instead, we have to loop. 
        # I guess next time we could use AstroPy table, rather than a list of dictionaries.
        
        num_profiles = np.size(self.profile_radius_dn_arr)

        profile_out = {}

        # First do radial profiles
                
        self.profile_radius_dn_arr[0] = np.sum(self.profile_radius_dn_arr, axis=0) # Save the summed profile

        # Then do the azimuthal profiles
        
        self.profile_azimuth_dn_arr[0] = np.sum(self.profile_azimuth_dn_arr, axis=0) # Save the summed profile
        
        # Then collapse the other fields (from N elements, to 1). This is very rough, and not a good way to do it.
        # We collapse down into a 1-element array. Alternatively, we could collapse to a scalar.
        
        self.profile_radius_dn_arr = np.array([self.profile_radius_dn_arr[0]])
        self.profile_azimuth_dn_arr = np.array([self.profile_azimuth_dn_arr[0]])
        
        self.exptime_arr     = np.array([np.sum(self.exptime_arr)])
        self.azimuth_arr     = np.array([self.azimuth_arr[0]])
        self.radius_arr      = np.array([self.radius_arr[0]])
        self.index_image_arr = np.array(["{}-{}".format(self.index_image_arr[0], self.index_image_arr[-1])])
        self.dt_arr          = np.array([self.dt_arr[0]])
        
        return self
 
# =============================================================================
# Copy the object
# =============================================================================

    def copy(self):
        import copy
        
        return copy.deepcopy(self)   # copy.copy() copies by reference, which is def not what we want

# =============================================================================
# Remove background slope (trend) from radial profiles
# =============================================================================

    def remove_background_radial(self, xranges, do_plot=False):

        """
    This routine takes a 1D radial profile, and a list of radii. It fits a linear trend to the brightness
    at those radii. It then subtracts that linear trend. The object then contains the same radial profile,
    but with the linear trend removed.
    
    This works on all radial profiles in the object.
    
    With do_plot=True, then make a one-off plot showing the removed trend.
    
    """
    
         # Remove background from all the radial profiles
      
        num_ranges = hbt.sizex(xranges)  # The range of radii to use for fitting the bg level. 
                                         # Total of (radius_bg) x (2) elements.

        radius  = self.radius_arr[0]    # Assume that all profiles share a common radius array 

        is_radius_good = np.zeros((num_ranges, np.size(radius))) # 

        # For each set of distance ranges, convert from distance, to bins. Assume all radial profiles use same bins.
        
        for i in range(num_ranges):  # e.g., if we have three ranges, loop over all of them
            range_i = xranges[i,:]
            is_radius_good[i,:] = np.logical_and( (radius > range_i[0]), (radius < range_i[1]) )

        is_radius_good = (np.sum(is_radius_good, axis=0) == 1)
        
        # Now we have made our radius flag array. We assume (for now) that all profiles have the same radius values.
        
        # Now loop over each profile, and find the right fit coefficients for it.
        
        for i in range(self.num_profiles()):  
            profile = self.profile_radius_dn_arr[i]

            r = linregress(radius[is_radius_good], profile[is_radius_good])
            m = r[0]
            b = r[1]
            profile_fit = radius*m + b
            
            if (do_plot):
                plt.plot(radius, profile_fit)
                plt.plot(radius,profile, label='Before')
                plt.plot(radius[is_radius_good], profile[is_radius_good], marker='+', label='Fit Vals')
                plt.plot(radius, profile-profile_fit, label='After')
                plt.plot(radius[is_radius_good], (profile-profile_fit)[is_radius_good], marker='+', label='Fit Vals')
                plt.legend()
                
                plt.show()
         
            # Save the new subtracted profile
            
            self.profile_radius_dn_arr[i] -= profile_fit

        return self 

# =============================================================================
# Plot the profiles
# =============================================================================
    
    def plot(self,    plot_radial=True, 
                      plot_azimuthal=False, 
                      dy_radial=1.5,  # Default vertical offset between plots 
                      dy_azimuthal=3,
                      smooth=None,
                      title=None,
                      **kwargs):

        #==============================================================================
        # Make a consolidated plot of radial profile
        #==============================================================================
        
        if (plot_radial):
            hbt.figsize((10,8))
                                                
            for i,index_image in enumerate(self.index_image_arr):

                radius = self.radius_arr[i]
                profile_radius = self.profile_radius_dn_arr[i]
                
                plt.plot(radius, 
                         (i * dy_radial) + profile_radius, 
                         label = '{}/{}, {:.2f}°, {}'.format(
                                       self.index_group_arr[i], 
                                       self.index_image_arr[i], 
                                       self.ang_phase_arr[i]*hbt.r2d,
                                       self.dt_str_arr[i]),
                         **kwargs)
                
            plt.xlabel('Radial Distance [km]')
            plt.ylabel('DN')
            plt.legend()
            if (title is not None):
                plt.title(title)    
            plt.show()
        
        #==============================================================================
        # Make a consolidated plot of azimuthal profile
        #==============================================================================
        
        if (plot_azimuthal):
            hbt.figsize((18,2))
            
            for i,profile_azimuth in enumerate(self.profile_azimuth_dn_arr):
                plt.plot(self.azimuth_arr[i,:]*hbt.r2d, 
                         (i * dy_azimuthal) + profile_azimuth, 
                         label = '{}/{}, {:.2f}°'.format(
                                 self.index_group_arr[i], 
                                 self.index_image_arr[i], 
                                 self.ang_phase_arr[i]*hbt.r2d),
                         **kwargs)
                
            plt.xlabel('Azimuth [deg]')
            plt.ylabel('DN')
            plt.legend()
            if (title is not None):
                plt.title(title)
            plt.show()

        return self

# =============================================================================
# Calculate the radial area under a curve
# =============================================================================

    def area_radial(self, limit):
 
        radius = self.radius_arr[0]     # Assume that all profiles use the same radius
        
        dradius = radius - np.roll(radius,1)  # Width of each bin, in km
        
        bin0   = np.where(radius > limit[0])[0][0]
        bin1   = np.where(radius > limit[1])[0][0]
        
        area = []
        
        for i in range(self.num_profiles()):
            area_i = np.sum((self.profile_radius_dn_arr[i] * dradius)[bin0:bin1])
            area.append(area_i)
            
        return area    
        
        
# =============================================================================
# Now do some tests!
# =============================================================================
    
# Invoke the class
        
# Q: Why can't I do ring.load().sum().plot()  ?
        
        
ring = ring_profile()

# Define the off-ring locations. Pixels in this region are used to set the background level.

radius_bg = np.array([[127500,127900], [129200,129700]])

# Make a plot of all of the data

limit_radial = (128000,130000)

# Set up output arrays

# Create an astropy table to dump results into
#phase_all = np.array([])
#area_all  = np.array([])

t = Table([[], [], [], [], [], [], []], 
          names=('index_group', 'index_image', 'phase', 'area_dn', 'exptime', 'label', 'sequence'),
          dtype = ('int', 'int', 'float64', 'float64', 'float64', 'U30', 'U30'))

# Make a copy of this, which we put the 'summed' profiles in -- that is, one brightness, per image sequence

#t_summed = Table([[], [], [], [],
#                  names = ('area_dn', 'exptime', 'label', 'area_dn'),
#                  dtype = ('float64', '))

params        = [(7, hbt.frange(0,7),   'core'),  # For each plot, we list a tuple: index_group, index_image, key
                 (7, hbt.frange(8,15),  'core'),
                 (7, hbt.frange(16,23), 'core'),
                 (7, hbt.frange(24,31), 'core'),
                 (7, hbt.frange(32,35), 'core'),
                 (7, hbt.frange(36,39), 'outer'), # Lots of 
                 (7, hbt.frange(40,42), 'core'),
                 (7, hbt.frange(52,54), 'core'),
                 (7, hbt.frange(61,63), 'core'),
                 (7, hbt.frange(91,93), 'core'),
                 (7, hbt.frange(91,93), 'outer-30')]

# Loop over each of the sets of images

for param_i in params:

    # Unwrap the tuple, describing this set of files (e.g., 7/12 .. 7/15)
    
    (index_group, index_images, key_radius) = param_i   # Unwrap the tuple

    # Create a string to describe this set of images ("7/0-7 core")
    
    sequence = ("{}/{}-{} {}".format(index_group, np.amin(index_images), np.amax(index_images), key_radius))

    # Load the profile from disk.
    
    ring = ring_profile()
    ring.load(index_group, index_images, key_radius = key_radius).smooth(1)
    
    # Copy the profile, and sum all the individual curves in it.
    
    ring_summed = ring.copy().sum()
    
    # Remove the background, and measure the area of each curve
    
    ring.remove_background_radial(radius_bg, do_plot=False)
    area = ring.area_radial(limit_radial)
    
    # Plot summed radial profile
    
#    ring_summed.plot(plot_azimuthal=False)
    
    # Make a plot of the phase curve from just this file
    
    phase = ring.ang_phase_arr
#    plt.plot(phase*hbt.r2d, area, marker='+', linestyle='none')
    plt.show()
    
    # Save the area and phase angle for each individual phase curve

    for i in range(np.size(index_images)):
        
      label = ("{}/{} {}".format(index_group, index_images[i], key_radius))
        
      t.add_row([index_group, index_images[i], phase[i]*u.radian, area[i], ring.exptime_arr[i], label, sequence])

# Aggregate the groups. This merges things and takes means of all numeric columns. Any non-numeric columns are dropped.
# The argument to groups.aggregate is a function to apply (e.g., numpy)
      
t_mean = t.group_by('sequence').groups.aggregate(np.mean)  # Mean
t_std  = t.group_by('sequence').groups.aggregate(np.std)   # Stdev
      
#    phase_all = np.concatenate((phase_all, phase))
#    area_all  = np.concatenate((area_all, area))

# Plot the phase curve, with one point per image

hbt.set_fontsize(size=15)
for i in range(len(t_mean['phase'])):
    plt.errorbar(t_mean['phase'][i] * hbt.r2d, t_mean['area_dn'][i], yerr=t_std['area_dn'][i],
             marker = 'o', ms=10, linestyle = 'none', label=t_mean['sequence'][i])
plt.xlabel('Angle [deg]')
plt.ylabel('Ring Area [DN]')
plt.legend()
plt.show()

# Plot the phase curve, with one point per sequence
      
plt.plot(t['phase'] * hbt.r2d, t['area_dn'], marker = 'o', linestyle = 'none', label=t['label'])
plt.legend()
plt.show()

# Now we want to do some astropy table manipulation. Take all thigns 
# 
# Plot several individual radial profiles

a = ring_profile()
a.load(7, hbt.frange(91,93),key_radius='outer-30').smooth(1).remove_background_radial(radius_bg,do_plot=True).plot()
plt.show()

a = ring_profile()
a.load(7, hbt.frange(91,93),key_radius='core').plot()
plt.show()
