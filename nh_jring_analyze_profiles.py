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

class App:
    
# =============================================================================
#     Init method: load the main table file
# =============================================================================
    
    def __init__(self, master):

        file_pickle = 'nh_jring_read_params_571.pkl' # Filename to read to get filenames, etc.
        dir_out     = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
            
        lun = open(dir_out + file_pickle, 'rb')
        self.t = pickle.load(lun)
        lun.close()
        
        # Process the group names. Some of this is duplicated logic -- depends on how we want to use it.
        
        self.groups = astropy.table.unique(self.t, keys=(['Desc']))['Desc']

        stretch_percent = 90    
        stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

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
    
    def load_profiles(self, index_group, index_images, plot_radial=True, plot_azimuthal=True, 
                      **kwargs):

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
            
            self.profile_radius_dn_arr.append(profile_radius_dn)
            self.profile_azimuth_dn_arr.append(profile_azimuth_dn)
            self.exptime_arr.append(exptime)
            self.ang_phase_arr.append(ang_phase)
            self.radius_arr.append(radius)
            self.azimuth_arr.append(azimuth)
            self.index_image_arr.append(index_image)
            self.dt_arr.append(dt)
            self.dt_str_arr.append(dt_str)
            
            print("Read image {}/{}, phase = {:0.1f}°, {}, {}, {}".format(
                    index_group, index_image, ang_phase*hbt.r2d, file_short, bg_method, bg_argument))
        
        #==============================================================================
        # Now put these into arrays (not lists). Ideally we'd put these into an astropy table (so we can sort, etc.)
        #==============================================================================
        
        self.ang_phase_arr = np.array(self.ang_phase_arr)
        self.azimuth_arr   = np.array(self.azimuth_arr)
        
        return

# =============================================================================
# Smooth the profiles
# This smooths profiles in place. It doesn't return a new object. It destroys the 
# current state, overwriting it with a new state.
# =============================================================================

    def smooth_profiles(self, width=1, kernel=None):

        if (kernel is None):
            kernel = Gaussian1DKernel(width)
                
        num_profiles = np.size(self.profile_radius_dn_arr)

        for i in range(num_profiles):
        
            # First do radial kernels

            for key in self.profile_radius_dn_arr[i]: # Loop over every key ('core', 'edge', etc)
                profile = convolve(self.profile_radius_dn_arr[i][key], kernel)
                self.profile_radius_dn_arr[i][key] = profile

            # Then do azimuthal kernels

            for key in self.profile_azimuth_dn_arr[i]:
                profile = convolve(self.profile_azimuth_dn_arr[i][key], kernel)
                self.profile_azimuth_dn_arr[i][key] = profile
                
# =============================================================================
# Sum the profiles
# =============================================================================

    def sum_profiles(self):
        
        num_profiles = np.size(self.profile_radius_dn_arr)

        profile_out = {}

        # First do radial profiles
        
        
        for key in self.profile_radius_dn_arr[0]: # Loop over every key ('core', 'edge', etc)
            
            profile_out[key] = np.array(self.profile_radius_dn_arr[0][key] * 0)
            for i in range(num_profiles):
                profile_out[key] += np.array(self.profile_radius_dn_arr[i][key])
                
            self.profile_radius_dn_arr[0][key] = profile_out[key] # Save the summed profile

        # Then do the azimuthal profiles
        
        for key in self.profile_azimuth_dn_arr[0]: # Loop over every key ('core', 'edge', etc)
            
            profile_out[key] = np.array(self.profile_azimuth_dn_arr[0][key] * 0)
            for i in range(num_profiles):
                profile_out[key] += np.array(self.profile_azimuth_dn_arr[i][key])
                
            self.profile_azimuth_dn_arr[0][key] = profile_out[key] # Save the summed profile
        
        # Then collapse the other fields (from N elements, to 1). This is very rough, and not a good way to do it.
        
        self.profile_radius_dn_arr = np.array([self.profile_radius_dn_arr[0]])
        self.profile_azimuth_dn_arr = np.array([self.profile_azimuth_dn_arr[0]])
        
        self.exptime_arr     = np.array([self.exptime_arr[0] * num_profiles])
        self.azimuth_arr     = np.array([self.azimuth_arr[0]])
        self.radius_arr      = np.array([self.radius_arr[0]])
        self.index_image_arr = np.array(["{}-{}".format(self.index_image_arr[0], self.index_image_arr[-1])])
 
# =============================================================================
# Copy the object
# =============================================================================

    def copy(self):
        import copy
        
        return copy.deepcopy(self)
    
# =============================================================================
# Plot the profiles
# =============================================================================
    
    def plot_profiles(self, plot_radial=True, plot_azimuthal=True, 
                      dy_radial=1.5,  # Default vertical offset between plots 
                      dy_azimuthal=3,
                      smooth=None,
                      **kwargs):

        #==============================================================================
        # Make a consolidated plot of radial profile
        #==============================================================================
        
        if (plot_radial):
            hbt.figsize((10,8))
            
            alpha_radial = 0.5
            
            # Set up an empty array to hold a summed profile
            
            profile_radius_sum = profile_radius_dn_arr[0]['core']*0
            
            for i,index_image in enumerate(self.index_image_arr):
#            for i,profile_radius in enumerate(self.profile_radius_dn_arr):
                radius = self.radius_arr[i]
                profile_radius = self.profile_radius_dn_arr[i]['core']
                
                plt.plot(radius, 
                         (i * dy_radial) + profile_radius, 
                         alpha = alpha_radial, 
                         label = '{}/{}, {:.2f}°, {}'.format(
                                       index_group, 
                                       self.index_image_arr[i], 
                                       self.ang_phase_arr[i]*hbt.r2d,
                                       self.dt_str_arr[i]),
                         **kwargs)
                profile_radius_sum += profile_radius
                
            plt.xlabel('Radial Distance [km]')
            plt.ylabel('DN')
            plt.legend()
            plt.show()
        
        #==============================================================================
        # Make a consolidated plot of azimuthal profile
        #==============================================================================
        
        if (plot_azimuthal):
            hbt.figsize((18,2))
            alpha_azimuthal = 0.5
            
            for i,profile_azimuth in enumerate(self.profile_azimuth_dn_arr):
                plt.plot(self.azimuth_arr[i,:]*hbt.r2d, 
                         (i * dy_azimuthal) + profile_azimuth['net'], 
                         alpha = alpha_azimuthal,
                         label = '{}/{}, {:.2f}°'.format(
                                 index_group, self.index_image_arr[i], self.ang_phase_arr[i]*hbt.r2d),
                         **kwargs)
                
            plt.xlabel('Azimuth [deg]')
            plt.ylabel('DN')
            plt.legend()
            plt.show()

        return

# Invoke the class
        
ring = App(0)

# Make a plot
ring.load_profiles(7, hbt.frange(16,23))
ring.plot_profiles(plot_radial = True, plot_azimuthal=False)
ring.smooth_profiles(1)
ring.plot_profiles(plot_radial = True, plot_azimuthal=False)

ring2 = ring.copy()
ring2.smooth_profiles(10)

ring.sum_profiles()
ring.plot_profiles(plot_radial = True, plot_azimuthal=False)
ring2.plot_profiles(plot_radial = True, plot_azimuthal=False)
ring2.sum_profiles()
ring2.plot_profiles(plot_radial = True, plot_azimuthal=False)



    
    

       