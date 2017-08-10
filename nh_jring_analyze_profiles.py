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
from nh_jring_extract_profile_from_unwrapped   import nh_jring_extract_profile_from_unwrapped   

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
# Remove Background from radial profiles
# =============================================================================

    def remove_background_radial(self, xranges, do_plot=False):

         # Remove background from all the radial profiles

        num_ranges = hbt.sizex(xranges)  # Radius_bg has total of (radius_bg) x (2) elements

        radius  = self.radius_arr[0]    # Assume that all profiles share a common radius array 

        is_radius_good = np.zeros((num_ranges, np.size(radius))) # 

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
    
    def plot(self, plot_radial=True, plot_azimuthal=True, 
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
radius_bg = np.array([[127500,127900], [129200,129700]])

# Make a plot

limit = (128000,130000)

ring1 = ring_profile()
ring1.load(7, hbt.frange(0,7)).smooth(1)
ring11 = ring1.copy().sum()
ring1.remove_background_radial(radius_bg, do_plot=False)
area1 = ring1.area_radial(limit)
ring1.plot(plot_azimuthal=False)
phase1 = ring1.ang_phase_arr
plt.plot(phase1*hbt.r2d, area1, marker = '+', linestyle='none')
plt.plot(ring11.radius_arr[0], ring11.profile_radius_dn_arr[0])

ring2 = ring_profile()
ring2.load(7, hbt.frange(8,15)).smooth(1).remove_background_radial(radius_bg, do_plot=False)
ring21 = ring2.copy().sum()
area2 = ring2.area_radial(limit)
ring2.plot(plot_azimuthal=False)
phase2 = ring2.ang_phase_arr
plt.plot(phase2*hbt.r2d, area2, marker = '+', linestyle='none')
plt.show()
plt.plot(ring21.radius_arr[0], ring21.profile_radius_dn_arr[0])


ring3 = ring_profile()
ring3.load(7, hbt.frange(16,23)).smooth(1).remove_background_radial(radius_bg, do_plot=False)
ring31 = ring3.copy().sum()
area3 = ring3.area_radial(limit)
ring3.plot(plot_azimuthal=False)
phase3 = ring3.ang_phase_arr
plt.plot(phase3*hbt.r2d, area3, marker = '+', linestyle='none')
plt.show()
plt.plot(ring31.radius_arr[0], ring31.profile_radius_dn_arr[0])

ring4 = ring_profile()
ring4.load(7, hbt.frange(24,31)).smooth(1).remove_background_radial(radius_bg, do_plot=False)
ring41 = ring4.copy().sum()
area4 = ring4.area_radial(limit)
ring4.plot(plot_azimuthal=False)
phase4 = ring4.ang_phase_arr
plt.plot(phase4*hbt.r2d, area4, marker = '+', linestyle='none')
plt.show()
plt.plot(ring41.radius_arr[0], ring41.profile_radius_dn_arr[0])


ring5 = ring_profile()
ring5.load(7, hbt.frange(32,35)).smooth(1).remove_background_radial(radius_bg, do_plot=False)
ring51 = ring5.copy().sum()
area5 = ring5.area_radial(limit)
ring5.plot(plot_azimuthal=False)
phase5 = ring5.ang_phase_arr
plt.plot(phase5*hbt.r2d, area5, marker = '+', linestyle='none')
plt.show()
plt.plot(ring51.radius_arr[0], ring51.profile_radius_dn_arr[0])


ring6 = ring_profile()
ring6.load(7, hbt.frange(36,39)).smooth(1).remove_background_radial(radius_bg, do_plot=False)
ring61 = ring6.copy().sum()
area6 = ring6.area_radial(limit)
ring6.plot(plot_azimuthal=False)
phase6 = ring6.ang_phase_arr
plt.plot(phase6*hbt.r2d, area6, marker = '+', linestyle='none')
plt.show()
plt.plot(ring61.radius_arr[0], ring61.profile_radius_dn_arr[0])


ring7 = ring_profile()
ring7.load(7, hbt.frange(40,42)).smooth(1).remove_background_radial(radius_bg, do_plot=False)
ring71 = ring7.copy().sum()
area7 = ring7.area_radial(limit)
ring7.plot(plot_azimuthal=False)
phase7 = ring7.ang_phase_arr
plt.plot(phase7*hbt.r2d, area7, marker = '+', linestyle='none')
plt.show()
plt.plot(ring71.radius_arr[0], ring71.profile_radius_dn_arr[0])

ring8 = ring_profile()
ring8.load(7, hbt.frange(52,55)).smooth(1).remove_background_radial(radius_bg, do_plot=False)
ring81 = ring8.copy().sum()
area8 = ring8.area_radial(limit)
ring8.plot(plot_azimuthal=False)
phase8 = ring8.ang_phase_arr
plt.plot(phase8*hbt.r2d, area8, marker = '+', linestyle='none')
plt.show()
plt.plot(ring81.radius_arr[0], ring81.profile_radius_dn_arr[0])

ring9 = ring_profile()
ring9.load(7, hbt.frange(91,93)).smooth(1).remove_background_radial(radius_bg, do_plot=False)
ring91 = ring9.copy().sum()
area9 = ring9.area_radial(limit)
ring9.plot(plot_azimuthal=False)
phase9 = ring9.ang_phase_arr
plt.plot(phase9*hbt.r2d, area9, marker = '+', linestyle='none')
plt.show()
plt.plot(ring91.radius_arr[0], ring91.profile_radius_dn_arr[0])

plt.plot(ring11.radius_arr[0], ring11.profile_radius_dn_arr[0] / ring1.num_profiles())
plt.plot(ring21.radius_arr[0], ring21.profile_radius_dn_arr[0] / ring2.num_profiles())
plt.plot(ring31.radius_arr[0], ring31.profile_radius_dn_arr[0] / ring3.num_profiles())
plt.plot(ring41.radius_arr[0], ring41.profile_radius_dn_arr[0] / ring4.num_profiles())
plt.plot(ring51.radius_arr[0], ring51.profile_radius_dn_arr[0] / ring5.num_profiles())
plt.plot(ring61.radius_arr[0], ring61.profile_radius_dn_arr[0] / ring6.num_profiles())
plt.plot(ring71.radius_arr[0], ring71.profile_radius_dn_arr[0] / ring7.num_profiles())
plt.plot(ring81.radius_arr[0], ring81.profile_radius_dn_arr[0] / ring8.num_profiles())
plt.plot(ring91.radius_arr[0], ring91.profile_radius_dn_arr[0] / ring9.num_profiles())
plt.ylim((-10,20))
plt.xlabel('Radius [km]')
plt.ylabel('DN sum')
plt.xlim((127000,130000))
plt.title('NH Radial Profiles, Group 7')
plt.show()

phase_all = np.concatenate((phase1, phase2, phase3, phase4, phase5, phase6, phase7, phase8, phase9))
area_all  = np.concatenate((area1,  area2,  area3,  area4,  area5,  area6,  area7,  area8,  area9))

hbt.set_fontsize(15)
plt.plot(phase_all*hbt.r2d, area_all, marker = 'o', linestyle='none')
plt.ylabel('DN')
plt.xlabel('Phase angle [deg]')
plt.title('NH Phase Curve, Group 7')
plt.show()


ring.load(7, hbt.frange(16,23)).smooth(1).plot(plot_azimuthal=False)
ring.remove_background_radial(radius_bg, do_plot=False).plot()


area = ring.area_radial(limit)
phase = ring.ang_phase_arr
plt.plot(phase*hbt.r2d, area, marker = '+', linestyle='none')



ring.plot(plot_azimuthal=False)
ring.sum().plot()

ring.plot(plot_radial = True, plot_azimuthal=False, title = 'Core')
ring.smooth()
ring.sum()
ring.smooth()

ring_full = ring_profile()
ring_full.load(7, hbt.frange(16,23), key_radius='full')
ring_full.plot(plot_azimuthal=False, title='Full')

ring.sum()
ring.plot(plot_radial = True, plot_azimuthal=False)  # This makes a really nice radial profile

# Try some background subtraction


profile = ring.profile_radius_dn_arr[0]
radius  = ring.radius_arr[0] 
vals_radius = radius_bg
bins_x = radius
num_ranges = hbt.sizex(radius_bg)  # Radius_bg has total of (radius_bg) x (2) elements


# Now loop over all of the radius ranges (e.g., three ranges).
# In each range, flag the x values that are within that range.
# Then, we'll sum these, to get a list of xvals that match *any* of the specified ranges.
# And then we'll do the linfit, using just these xvals, and the corresponding yvals.



#for i in 
#r = linregress(arr1_filter.flatten(), arr2_filter.flatten())
#
#m = r[0] # Multiplier = slope
#b = r[1] # Offset = intercept
#
#arr2_fixed = arr2 * m + b
#return (m,b)
ring.plot()

ring_sum = ring.copy()
ring_sum.sum()
ring_sum.smooth()
ring_sum.plot(plot_radial = True, plot_azimuthal=False)

