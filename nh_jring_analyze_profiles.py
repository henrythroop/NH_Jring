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
from   pymiecoated import Mie

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

from scatter_mie_ensemble                      import scatter_mie_ensemble
#from area_between_line_curve                    import area_between_line_curve

# For fun, I'm going to try to do this whole thing as a class. That makes it easier to pass 
# params back and forth, I guess. [A: Yes, this is a great place for a class -- very clean.]

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
                      verbose = False,      # Verbose: List the filename as loading?
                      **kwargs):

        # Each file on disk has several different extractions:
        #
        #   Radial profile has 'core', 'full', 'half, etc.
        #   Az     profile has 'inner', 'outer', 'net', etc. 
        #
        # Since we usually need just one of these, this routine *only* reads one.
        # The individual one read can be controlled with key_radius and key_azimuth
        #
        # Usually the values of the profiles will be read in units of DN. 
        # But, it is up to the user to keep straight if it is DN, IoF, normal IoF, EW, etc.
        
        import humanize
        
        # Define t_group, since we'll reference it a lot
        
        self.groupmask = self.t['Desc'] == self.groups[index_group]
        self.t_group = self.t[self.groupmask]  # 
        
        self.num_images_group = np.size(self.t_group)
        
        self.index_group            = index_group

        self.profile_radius_arr     = []  # This can be in DN, *or* IoF.
        self.profile_azimuth_arr    = []  # This can be in DN, *or* IoF.
        self.exptime_arr            = []
        self.ang_phase_arr          = []
        self.radius_arr             = []
        self.ang_elev_arr           = []
        self.azimuth_arr            = []
        self.index_image_arr        = []
        self.index_group_arr        = []
        self.dt_arr                 = [] # dt is the time since the datafile was written. 
        self.dt_str_arr             = []
        
        self.profile_azimuth_units  = 'DN per pixel'
        self.profile_radius_units   = 'DN per pixel'

        for index_image in index_images:
            file = self.get_export_analysis_filename(index_group, index_image)
            
            lun = open(file, 'rb')
            vals = pickle.load(lun)
            lun.close()
        
            (image_unwrapped,     # Unwrapped image itself
                        mask_objects_unwrapped,  # Boolean mask: True = good pixel
                        mask_stray_unwrapped,
                        radius,              # Axis values for the unwrapped image
                        azimuth,             # Axis values for the unwrapped image 
                        profile_radius,      # Radial profile (several, in a dictionary)
                        profile_azimuth,     # Az profile (several, in a dictionary)
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
            
            self.profile_radius_arr.append(profile_radius[key_radius])  # This is now (eg) 8 x 30 array
            self.profile_azimuth_arr.append(profile_azimuth[key_azimuth])
            self.exptime_arr.append(exptime)
            self.ang_phase_arr.append(ang_phase)
            self.ang_elev_arr.append(ang_elev)
            self.radius_arr.append(radius)
            self.azimuth_arr.append(azimuth)
            self.index_image_arr.append(index_image)
            self.index_group_arr.append(index_group)
            self.dt_arr.append(dt)
            self.dt_str_arr.append(dt_str)

            if verbose:
                print("Read image {}/{}, phase = {:0.1f}°, {}, {}, {}".format(
                    index_group, index_image, ang_phase*hbt.r2d, file_short, bg_method, bg_argument))
                print("{}".format(file))
        
        #==============================================================================
        # Now put these into arrays (not lists). Ideally we'd put these into an astropy table (so we can sort, etc.)
        #==============================================================================
             
        self.ang_phase_arr = np.array(self.ang_phase_arr)
        self.ang_elev_arr  = np.array(self.ang_elev_arr)
        self.azimuth_arr   = np.array(self.azimuth_arr)
                 
        self.profile_radius_arr  = np.array(self.profile_radius_arr)
        self.profile_azimuth_arr = np.array(self.profile_azimuth_arr) 
        
        return self

# =============================================================================
# Return number of radial bins
# =============================================================================
        
    def num_bins_radius(self):
        
        return np.shape(self.profile_radius_arr)[1]

 
# =============================================================================
# Return size of azimuthal array
# =============================================================================
        
    def num_bins_azimuth(self):
        
        return np.shape(self.profile_azimuth_arr)[1]

 
# =============================================================================
# Return number of profiles
# =============================================================================
        
    def num_profiles(self):
        
        return np.size(self.dt_arr)

    
# =============================================================================
# Convert from DN to I/F. Or at least create the I/F curves.
# =============================================================================

    def dn2iof(self):
        
# The math here follows directly from that in "NH Ring Calibration.pynb", which I used for Pluto rings.
# That in turn follows from Hal Weaver's writeup.
# The key constants are RSOLAR and F_solar, which together convert from DN to I/F.
# 
# How to use: Call ring.dn2iof() once. Then after that, any reference to the profiles is in I/F, not DN.
        
        RSOLAR =  221999.98  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)

        C = self.profile_radius_arr  # Get the DN values of the ring. Typical value is 1 DN.

        # Define the solar flux, from Hal's paper.
        
        FSOLAR_LORRI  = 176.	     	    # We want to be sure to use LORRI value, not MVIC value!
        F_solar = FSOLAR_LORRI # Flux from Hal's paper
        
        # Calculate the Jupiter-Sun distance, in AU (or look it up). 
        
        r = 5.35 # Distance in AU. For the Jupiter encounter, dist = 5.30 AU .. 5.39 AU for the 10d surrounding it.
        
        TEXP_2D = np.transpose(np.tile(self.exptime_arr, (self.num_bins_radius(),1)))  # Increase this from 1D to 2D
        
        I = C / TEXP_2D / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. All v similar, except for spectrum assumed.
        
        # Apply Hal's conversion formula from p. 7, to compute I/F and print it.
        
        IoF = math.pi * I * r**2 / F_solar # Equation from Hal's paper
        
        # Now convert to 'normal I/F'
        
        # Define mu = cos(e), where e = emission angle, and e=0 is face-on.
        
        e = math.pi/2 - self.ang_elev_arr # NB: I made an error in this for the Pluto ring analysis! 
        
        mu = np.cos(e)
        
        mu_2D = np.transpose(np.tile(mu, (self.num_bins_radius(),1)))
        
        # Calculate the normal I/F
        
        IoF_normal = 4 * mu_2D * IoF
        
        # Save this result into 'self,' replacing the DN values with the normal I/F values
        
        self.profile_radius_arr = IoF_normal
        
        # Change the units
        
        self.profile_radius_units = 'Normal I/F'
        
        # Do the same for azimuthal profile

        TEXP_2D = np.transpose(np.tile(self.exptime_arr, (self.num_bins_azimuth(),1)))
        mu_2D = np.transpose(np.tile(mu, (self.num_bins_azimuth(),1)))
        self.profile_azimuth_arr = self.profile_azimuth_arr / TEXP_2D / RSOLAR * math.pi * r**2 / F_solar * 4 * mu_2D
        
        self.profile_azimuth_units = 'Normal I/F'
        
        return self
        
#   From Hal Weaver email 7-Dec-2016
# These are for 1X1. 4X4 is different, but since they're saturated, I am ignoring.
# These are intrinsic to detector, and not a function of distance    
#Diffuse sensitivity keywords
#Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
#RSOLAR =  221999.98
#RPLUTO =  214583.33
#RPHOLUS =  270250.00
#RCHARON =  219166.66
#RJUPITER =  195583.33
      
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

            self.profile_radius_arr[i,:]  = convolve(self.profile_radius_arr[i,:],  kernel)

            # Then do azimuthal kernels
            
            self.profile_azimuth_arr[i,:] = convolve(self.profile_azimuth_arr[i,:], kernel)

        return self
    
#            for key in self.profile_azimuth_arr[i]:
#                profile = convolve(self.profile_azimuth_arr[i][key], kernel)
#                self.profile_azimuth_arr[i][key] = profile
                
# =============================================================================
# Flatten the profiles -- that is, if there are N profiles, flatten to one, with mean value
# =============================================================================

    def flatten(self):
        
        """ Flatten the profile, from N profile, into one mean profile.
            Applies to radial and azimuthal and all other quantities.
            Result is returned, and internal values are changed as well.
            After flattening, most outputs are as 1-element arrays (not scalars)."""
            
        # Use mean() to get a merged radial profile

        # First do radial profiles
                
        self.profile_radius_arr[0] = np.mean(self.profile_radius_arr, axis=0) # Save the flattened profile

        # Then do the azimuthal profiles
        
        self.profile_azimuth_arr[0] = np.mean(self.profile_azimuth_arr, axis=0) # Save the flattened profile

        self.profile_radius_arr  = np.array([self.profile_radius_arr[0]])
        self.profile_azimuth_arr = np.array([self.profile_azimuth_arr[0]])

        # For some quantities, like angles, take the mean and save that.
        # NB: If we add any new fields to the object, we need to add a corresponding new line to this method!
        
        self.ang_elev_arr    = np.array([np.mean(self.ang_elev_arr)])
        self.ang_phase_arr   = np.array([np.mean(self.ang_phase_arr)])
        self.exptime_arr     = np.array([np.mean(self.exptime_arr)])
        
        # For other fields, take one entry, and save that.
        # This is a bit crude, but probably OK.
                        
        self.azimuth_arr     = np.array([self.azimuth_arr[0]])
        self.radius_arr      = np.array([self.radius_arr[0]])
        self.index_image_arr = np.array(["{}-{}".format(self.index_image_arr[0], self.index_image_arr[-1])])
        self.index_group_arr = np.array([self.index_group_arr[0]]) 
        self.dt_arr          = np.array([self.dt_arr[0]])
        self.dt_str_arr      = np.array([self.dt_str_arr[0]])
        
        return self          # The value is returned, so that commands can be chained.
                             # Also, the value of all quantities (e.g., profile_azimuth_arr) is changed internally.
 
# =============================================================================
# Copy the object
# =============================================================================

    def copy(self):
        
        """ Copy an object to a new object (not a reference). """
        
        import copy
        
        return copy.deepcopy(self)   # copy.copy() copies by reference, which is def not what we want

# =============================================================================
# Remove background slope (trend) from radial profiles
# =============================================================================

    def remove_background_radial(self, xranges, do_plot=False, verbose=False):

        """
    This routine takes a 1D radial profile, and a list of radii. It fits a linear trend to the brightness
    at those radii. It then subtracts that linear trend. The object then contains the same radial profile,
    but with the linear trend removed.
    
    This works on all radial profiles in the object.
    
    With do_plot=True, then make a one-off plot showing the removed trend.
    
    Parameters
    -----
    
    xranges: 
        Range of radii to use for fitting the background level. Total of (radius_bg) x (2) elements. 
        Usually 2x2 = 4 elements.    
    
    XXX Warning: If the input radial range is not within the data's radial range, then array returned will be empty!
                 Needs better error handling.  
    """
    
         # Remove background from all the radial profiles
      
        num_ranges = hbt.sizex(xranges)  # The range of radii to use for fitting the bg level. 
                                         # Total of (radius_bg) x (2) elements. Usually 2x2 = 4 elements.

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
            profile = self.profile_radius_arr[i]

            r = linregress(radius[is_radius_good], profile[is_radius_good])
            m = r[0]
            b = r[1]
            profile_fit = radius*m + b
            
            if (do_plot):
                plt.plot(radius, profile_fit)
                plt.plot(radius,profile, label='Before')
                plt.plot(radius[is_radius_good], profile[is_radius_good], marker='+', ls='none', label='Fit Vals')
                plt.plot(radius, profile-profile_fit, label='After')
                plt.plot(radius[is_radius_good], (profile-profile_fit)[is_radius_good], ls='none', marker='+', 
                                                 label='Fit Vals')
                plt.legend()
                
                plt.show()
         
            # Save the new subtracted profile
            
            self.profile_radius_arr[i] -= profile_fit
            
            if verbose:
                print("Subtracted fit [y={} x + {}] from profile {}".format(m, b, i))

        return self 

# =============================================================================
# Plot the profiles
# =============================================================================
    
    def plot(self,    plot_radial=True, 
                      plot_azimuthal=False, 
                      dy_radial=0,  # Default vertical offset between plots 
                      dy_azimuthal=3,
                      plot_sats=True,  # Plot the location of Metis / Adrastea on radial profile?
                      smooth=None,
                      title=None,
                      **kwargs):

        #==============================================================================
        # Make plot of radial profile. Plot one line per profile (or one, if flattened)
        #==============================================================================
        
        if (plot_radial):
#            hbt.figsize((10,8))
                                                
            for i,index_image in enumerate(self.index_image_arr):

                radius = self.radius_arr[i]
                profile_radius = self.profile_radius_arr[i]
                
                p = plt.plot(radius, 
                         (i * dy_radial) + profile_radius, 
                         label = r'{}/{}, $\alpha$={:.2f}°, e={:.2f}°, {}'.format(
                                       self.index_group_arr[i], 
                                       self.index_image_arr[i], 
                                       self.ang_phase_arr[i]*hbt.r2d,
                                       self.ang_elev_arr[i]*hbt.r2d,
                                       self.dt_str_arr[i]),
                         **kwargs)
            
            # Plot Metis and Adrastea
            
            if (plot_sats):
                a_metis    = 128000
                a_adrastea = 129000
            
#                args = {'linestyle': 'dash', 'alpha':0.5, 'color':'black'}
                plt.axvline(x=a_metis, linestyle='dashed', alpha=0.2, color='black')
                plt.axvline(x=a_adrastea, linestyle='dashed', alpha=0.2, color='black')
                
            plt.xlabel('Radial Distance [km]')
            plt.ylabel(self.profile_radius_units)
 
            # Set scientific notation on the Y axis, if appropriate
            
            axes = plt.gca()
            
            if ('I/F' in self.profile_radius_units):       
                axes.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
            
            plt.legend()
            if (title is not None):
                plt.title(title)    
            plt.show()
        
        #==============================================================================
        # Make a consolidated plot of azimuthal profile
        #==============================================================================
        
        if (plot_azimuthal):
            hbt.figsize((18,2))
            
            for i,profile_azimuth in enumerate(self.profile_azimuth_arr):
                plt.plot(self.azimuth_arr[i,:]*hbt.r2d, 
                         (i * dy_azimuthal) + profile_azimuth, 
                         label = '{}/{}, {:.2f}°'.format(
                                 self.index_group_arr[i], 
                                 self.index_image_arr[i], 
                                 self.ang_phase_arr[i]*hbt.r2d),
                         **kwargs)
                
            plt.xlabel('Azimuth [deg]')
            plt.ylabel(self.profile_azimuth_units)
            plt.legend()
            if (title is not None):
                plt.title(title)
            plt.show()

        return self

# =============================================================================
# Calculate the radial area under a curve
# =============================================================================

    def area_radial(self, limit):
 
        """
        Calculate the radial area under a ring curve. 
        
        This routine is not very flexible -- recommend using area_between_line_curve() instead.
        """
        
        radius = self.radius_arr[0]     # Assume that all profiles use the same radius
        
        dradius = radius - np.roll(radius,1)  # Width of each bin, in km
        
        bin0   = np.where(radius > limit[0])[0][0]
        bin1   = np.where(radius > limit[1])[0][0]
        
        area = []
        
        for i in range(self.num_profiles()):
            area_i = np.sum((self.profile_radius_arr[i] * dradius)[bin0:bin1])
            area.append(area_i)
            
        return area    
        
        
# =============================================================================
# Now do some tests!
# =============================================================================
    
# Invoke the class
        
ring = ring_profile()

# Define the off-ring locations. Pixels in this region are used to set the background level.

radius_bg_core      = np.array([[127.5,127.9], [129.0,129.7]])*1000  # Core only
radius_bg_halo_core = np.array([[122,122.5],     [130,131]])*1000  # Cover full width of halo and core
radius_bg_117       = np.array([[117,120],     [130,131]])*1000  # This distance range is a bit better - no edge effects
radius_bg_127       = np.array([[125,126],     [130, 131]])*1000

# Make a plot of all of the data

# Define the limits used for the I/F 'area under the curve' total equivalent width

limit_radial = (122000,130000)  

# Set up output arrays

# Create an astropy table to dump results into. This is so we can make a phase curve in the end

# First create a table, one entry per image

t = Table([[], [], [], [], [], [], [], [], [], [], []],
    names=('index_group', 'index_image', 'phase', 'elev', 'area_dn', 'exptime', 'label', 'sequence', 'radius', 
           'profile_radius', 'profile_radius_units'),
    dtype = ('int', 'int', 'float64', 'float64', 'float64', 'float64', 'U30', 'U30', 'object', 'object', 'U30'))

# Then create a table, summed at one entry per phase angle. XXX No, we really don't need this one right now.

# List all of the data that we want to read in from disk

params        = [(7, hbt.frange(0,7),   'full'),  # For each plot, we list a tuple: (index_group, index_image, key)
                 (7, hbt.frange(8,15),  'full'),
                 (7, hbt.frange(16,23), 'full'),
                 (7, hbt.frange(24,31), 'full'),
                 (7, hbt.frange(32,35), 'full'),
                 (7, hbt.frange(36,39), 'full'), # Lots of 
                 (7, hbt.frange(40,42), 'full'),
                 (7, hbt.frange(52,54), 'full'),
                 (7, hbt.frange(61,63), 'full'),
                 (7, hbt.frange(91,93), 'full'),
                 (7, hbt.frange(94,96), 'full')]

# Loop over each of the sets of images
# For each image set, sum the individual profiles, and get one final profile

radius = []
profile_radius = []  # This is a list.

radius_bg = radius_bg_halo_core

for param_i in params:

    # Unwrap the tuple, describing this set of files -- e.g.,  (7, array([94, 95, 96]), 'full')
    
    (index_group, index_images, key_radius) = param_i   # Unwrap the tuple

    # Create a string to describe this set of images -- e.g.,    '7/94-96 full'
    
    sequence = ("{}/{}-{} {}".format(index_group, np.amin(index_images), np.amax(index_images), key_radius))

    # Load the profile from disk.
    
    ring = ring_profile()
    ring.load(index_group, index_images, key_radius = key_radius).smooth(1)
    
    # Copy the profile, and sum all the individual curves in it.
    
    ring_flattened = ring.copy().flatten()
    
    # Convert from DN to I/F
    
    ring.dn2iof()
    ring_flattened.dn2iof()
    
    # Remove the background from the I/F profile curve
    
    ring.remove_background_radial(radius_bg)
    ring_flattened.remove_background_radial(radius_bg)
    
    # Measure the area under the curve
    # XXX Also we have the function area_between_line_curve which is similar but not identical 
        
    area           = ring.area_radial(limit_radial)
    area_flattened = ring_flattened.area_radial(limit_radial)
    
    # Plot summed radial profile
    
#    ring_flattened.plot(plot_azimuthal=False)
#    plt.plot(ring_flattened.radius_arr, ring_flattened.profile_radius_arr)
    
    # Make a plot of the phase curve from just this file
    
    phase = ring.ang_phase_arr
    phase_flattened = ring_flattened.ang_phase_arr  # Since this is a flattened observation, phase is a 1-element array.
    
    # Save the area and phase angle for each individual phase curve

    for i in range(np.size(index_images)):  # Loop over each image
        
      label = ("{}/{} {}".format(ring.index_group_arr[0], ring.index_image_arr[i], key_radius))
        
      t.add_row([ring.index_group_arr[0],   # XXX this is slow, now that we are including arrays into it!
                 ring.index_image_arr[i], 
                 ring.ang_phase_arr[i],
                 ring.ang_elev_arr[i],
                 area[i],                  # This goes into 'area_dn'
                 ring.exptime_arr[i],
                 label,
                 sequence,
                 ring.radius_arr[i],
                 ring.profile_radius_arr[i],
                 ring.profile_radius_units])
    
# Aggregate the groups. This merges things and takes means of all numeric columns. 
# Any non-numeric columns are dropped -- even if they are all identical.
#    
# The argument to groups.aggregate is a function to apply (e.g., numpy)
# Q: What does .aggregate() will do with 2D NumPy arrays?
# A: I have tested it and it looks to actually properly do the job. 
#        For instance, it takes all of the radial profiles from that sequence, and averages them, resulting in one
#        mean radial profile.
# 
#        For np.std, it actually takes the stdev of the area, based on (e.g.) five individual curves that have been
#        summed for the final curve. So, it's a pretty decent way to get an error bar!      

t_mean = t.group_by('sequence').groups.aggregate(np.mean)  # Mean
t_std  = t.group_by('sequence').groups.aggregate(np.std)   # Stdev

num_sequences = len(t_mean)

# Add some columns to the table

# Set the core and halo positions. I have individually measured each of these.

radius_core_default = np.array([127.5, 129.5])*1000   # Radius, in units of 1000 km.
radius_halo_default = np.array([118.0, 130.0])*1000

t_mean.add_column(Table.Column(np.tile(radius_core_default, (num_sequences,1)),name='radius_core'))
t_mean.add_column(Table.Column(np.tile(radius_halo_default, (num_sequences,1)),name='radius_halo'))
t_mean.add_column(Table.Column(np.tile((0,0), (num_sequences,1)), name='bin_core'))
# Core radii, used
t_mean.add_column(Table.Column(np.zeros(num_sequences), name='area_core'))
t_mean.add_column(Table.Column(np.tile((0,0), (num_sequences,1)), name='bin_halo'))
t_mean.add_column(Table.Column(np.zeros(num_sequences), name='area_halo'))

# Copy the units over (which don't get aggregated properly, since they are a string)

t_mean.add_column(Table.Column(np.tile(ring.profile_radius_units, num_sequences),name='profile_radius_units'))

# Remove any columns that don't make any sense because they cannot be merged!

t_mean.remove_column('index_image')
t_std.remove_column('index_image')

# Sort by phase angle

t_mean.sort('phase')

# Set the radial position of the ring core. Inner and outer limits of it.
# Do this manually for each profile. This is just more accurate than trying to do it automated.

#%%%

# =============================================================================
# Extract the brightness of the core from the brightness of the halo
# =============================================================================

t_mean['radius_core'][t_mean['sequence'] == '7/8-15 full'] = np.array((127.85, 129.5))*1000
t_mean['radius_core'][t_mean['sequence'] == '7/0-7 full']  = np.array((127.65, 129.35))*1000
t_mean['radius_core'][t_mean['sequence'] == '7/24-31 full'] = np.array((127.8, 129.5))*1000
t_mean['radius_core'][t_mean['sequence'] == '7/16-23 full'] = np.array((127.65, 129.55))*1000 # Red
t_mean['radius_core'][t_mean['sequence'] == '7/36-39 full'] = np.array((127.75, 129.25))*1000 # Purple
t_mean['radius_core'][t_mean['sequence'] == '7/32-35 full'] = np.array((127.8, 129.6))*1000   # Brown
t_mean['radius_core'][t_mean['sequence'] == '7/40-42 full'] = np.array((128.0, 129.70))*1000 # pink
t_mean['radius_core'][t_mean['sequence'] == '7/61-63 full'] = np.array((127.8, 129.50))*1000 # grey
t_mean['radius_core'][t_mean['sequence'] == '7/52-54 full'] = np.array((127.5, 129.50))*1000 # olive
t_mean['radius_core'][t_mean['sequence'] == '7/91-93 full'] = np.array((128.0, 129.55))*1000 # lt blue
t_mean['radius_core'][t_mean['sequence'] == '7/94-96 full'] = np.array((127.5, 129.3))*1000 # dk blue top

t_mean['radius_halo'][t_mean['sequence'] == '7/8-15 full'] = np.array((118, 129.5))*1000
t_mean['radius_halo'][t_mean['sequence'] == '7/0-7 full']  = np.array((118, 129.5))*1000
t_mean['radius_halo'][t_mean['sequence'] == '7/24-31 full'] = np.array((118, 129.5))*1000
t_mean['radius_halo'][t_mean['sequence'] == '7/16-23 full'] = np.array((118, 129.5))*1000 # Red
t_mean['radius_halo'][t_mean['sequence'] == '7/36-39 full'] = np.array((118, 129.25))*1000 # Purple
t_mean['radius_halo'][t_mean['sequence'] == '7/32-35 full'] = np.array((118.7, 129.6))*1000   # Brown
t_mean['radius_halo'][t_mean['sequence'] == '7/40-42 full'] = np.array((118.5, 129.70))*1000 # pink
t_mean['radius_halo'][t_mean['sequence'] == '7/61-63 full'] = np.array((119.2, 129.50))*1000 # grey
t_mean['radius_halo'][t_mean['sequence'] == '7/52-54 full'] = np.array((118.5, 129.50))*1000 # olive
t_mean['radius_halo'][t_mean['sequence'] == '7/91-93 full'] = np.array((118.5, 129.55))*1000 # lt blue
t_mean['radius_halo'][t_mean['sequence'] == '7/94-96 full'] = np.array((120.5, 130.5))*1000 # dk blue top

# Now that we have the radii of the inner core calculated, look up the bin limits for the inner core

for i in range(num_sequences):
    radii       = t_mean['radius'][i]
    radius_core = t_mean['radius_core'][i]
    radius_halo = t_mean['radius_halo'][i]

    # Look up bin limits of the core and store them
    
    bins = [ np.argmin(np.abs(radii - t_mean['radius_core'][i][0])),
             np.argmin(np.abs(radii - t_mean['radius_core'][i][1])) ]
            
    t_mean['bin_core'][i] = bins

    # Look up bin limits of the halo, and store them
    
    bins = [ np.argmin(np.abs(radii - t_mean['radius_halo'][i][0])),
             np.argmin(np.abs(radii - t_mean['radius_halo'][i][1])) ]
            
    t_mean['bin_halo'][i] = bins

# Calculate the area under the curve for the core. 
# XXX We should probably write a method for this.
#     We should also improve things so that we can get an stdev from the component curves that make this up.
    
    (area_i_core, line, bins) = area_between_line_curve(t_mean['radius'][i],
                                     t_mean['profile_radius'][i],
                                     t_mean['radius_core'][i])
    t_mean['area_core'][i] = area_i_core
    
# Calculate the area under the curve for the halo.
# For the halo, we subtract off the area already measured for the core

    (area_i_halo, line, bins) = area_between_line_curve(t_mean['radius'][i],
                                     t_mean['profile_radius'][i],
                                     t_mean['radius_halo'][i])
    t_mean['area_halo'][i] = area_i_halo - area_i_core
    
    
# =============================================================================
# Make a plot: Radial Profile (Summed, I/F) vs. Sequence
# =============================================================================

hbt.figsize((12,20))

DO_PLOT_CORE_LIMITS = True
DO_PLOT_HALO_LIMITS = True

dy = 13e-7

#colors = iter(cm.rainbow(np.linspace(0, 1, 100)))

for i,s in enumerate(t_mean['sequence']):  # Loop over the text sequence name (e.g., '7/0-7 full')
    
    plt.plot(t_mean['radius'][i]/1000, t_mean['profile_radius'][i] + i * dy,
             linewidth=4,
             label = r'{}, {:.1f}$^{{\circ}}$, {:.1f}$^{{\circ}}$'.format(
            s, t_mean['phase'][i]*hbt.r2d, t_mean['elev'][i]*hbt.r2d) )
    
    if (DO_PLOT_CORE_LIMITS):
        bins = t_mean['bin_core'][i]
        plt.plot( np.array( (t_mean['radius'][i][bins[0]], t_mean['radius'][i][bins[1]]) )/1000 ,
                  np.array( (t_mean['profile_radius'][i][bins[0]], t_mean['profile_radius'][i][bins[1]]))  + i*dy,
                  color = 'red')
                         
    if (DO_PLOT_CORE_LIMITS):
        bins = t_mean['bin_halo'][i]
        plt.plot( np.array( (t_mean['radius'][i][bins[0]], t_mean['radius'][i][bins[1]]) )/1000 ,
                  np.array( (t_mean['profile_radius'][i][bins[0]], t_mean['profile_radius'][i][bins[1]]))  + i*dy,
                  color = 'grey', linestyle='dotted')
                         
plt.ylim((-1e-6,3e-6 + dy*i))
plt.xlim((115, 133))
plt.xlabel('Radius [1000 km]')
plt.title('J-ring radial profiles, summed per sequence')
plt.ylabel('I/F')
plt.legend()
plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
plt.show()

#%%%


#    phase_all = np.concatenate((phase_all, phase))
#    area_all  = np.concatenate((area_all, area))

# =============================================================================
# Do some Mie calculations
# =============================================================================

# Set up the size distribution

alam   = 500  * u.nm
rmin   = 0.01 * u.micron
rmax   = 50   * u.micron     
num_r  = 150

pi     = math.pi

# Define the exponent of the size distribution

q_1    = 2
q_2    = 5
r_break=15*u.micron  

# Define ratio of small:large bodies

ratio_lambert_mie = 0.20  # Ratio of Lambertian:Mie (or actually, Lambertian:Total)

r      = hbt.frange(rmin, rmax, num_r, log=True)*u.micron  # Astropy bug? When I run frange(), it drops the units.
n      = hbt.powerdist_broken(r, r_break, q_1, q_2)        # Power law  

# Set up the index of refraction
n_refract = 1.5
m_refract = 0.001 # Must be positive, not negative.
nm_refract = complex(n_refract,m_refract)

# Set up the angular distribution. We actually have two sets of angles we deal with here:
#   ang_model: A smooth set of phase angles, evenly spaced, for making a plot
#   ang_data:  The actual phase angles that the actual data were observed at.

num_ang_model = 300
ang_model = hbt.frange(0, pi, num_ang_model)*u.rad # Phase angle
ang_data  = t_mean['phase']*u.rad

# Compute the Mie and Lambert values at the 'data' and 'model' angles.

(p11_mie_data,  p11_mie_raw, qsca) = hbt.scatter_mie_ensemble(nm_refract, n, r, ang_data,  alam)  
              # Model at observed angles
(p11_mie_model, p11_mie_raw, qsca) = hbt.scatter_mie_ensemble(nm_refract, n, r, ang_model, alam)  
              # Model at all angles 

p11_lambert_data  = hbt.scatter_lambert(ang_data)
p11_lambert_model = hbt.scatter_lambert(ang_model)

# Merge Mie + Lambertian in the specified ratios

p11_merged_model = p11_mie_model * (1-ratio_lambert_mie) + p11_lambert_model * ratio_lambert_mie
p11_merged_data  = p11_mie_data  * (1-ratio_lambert_mie) + p11_lambert_data  * ratio_lambert_mie

# Do a linfit to compute proper magnitude of model, to fit data

DO_FIT_LOGSPACE = True
(scalefac_merged_log, _) = hbt.linfit_origin(p11_merged_data.value, t_mean['area_halo'], log=DO_FIT_LOGSPACE)
(scalefac_merged_lin, _) = hbt.linfit_origin(p11_merged_data.value, t_mean['area_halo'])
    
# =============================================================================
# Plot the phase curve, with one point per sequence
# =============================================================================

hbt.figsize((12,8))

hbt.set_fontsize(size=15)

for i in range(len(t_mean['phase'])):
    plt.errorbar(t_mean['phase'][i] * hbt.r2d, t_mean['area_halo'][i], #yerr=t_std['area_dn'][i],
             marker = 'o', ms=10, linestyle = 'none', alpha=0.5)

    plt.errorbar(t_mean['phase'][i] * hbt.r2d, t_mean['area_core'][i],
             marker = '+', ms=12, markeredgewidth=3, lw=3, linestyle = 'none')

plt.plot(ang_model.to('deg'), p11_mie_model        * scalefac_merged_log * (1-ratio_lambert_mie), 
         label = 'q = ({},{}), rb={}, n={}'.format(q_1, q_2, r_break, nm_refract), color = 'black', 
         linestyle = 'dotted')

plt.plot(ang_model.to('deg'), p11_lambert_model    * scalefac_merged_log * ratio_lambert_mie, label = 'Lambert', 
         color='black', linestyle='dashed')

plt.plot(ang_model.to('deg'), p11_merged_model     * scalefac_merged_log, 
         label = 'Merged, ratio = {}'.format(ratio_lambert_mie), 
         color='black', linestyle = 'dashdot')

plt.xlabel('Angle [deg]')
plt.ylabel(t_mean['profile_radius_units'][0][0])
plt.yscale('log')
plt.ylim((1e-6,1e-3))
plt.legend()
plt.show()



# =============================================================================
# Plot the phase curve, with one point per image
# =============================================================================
      
plt.plot(t['phase'] * hbt.r2d, t['area_dn'], marker = 'o', linestyle = 'none', label=t['label'])
#plt.yscale('log')
plt.title('Phase Curve, New Horizons J-Ring')
plt.ylim((1e-5, 1e-2))
plt.show()
 
# Plot several individual radial profiles

# =============================================================================
# Now make some customized individual profiles. For a 'best merged set' profile thing.
# =============================================================================

# This is a good profile for 0-7! Shows moonlet belt really well.
# Background subtration needed: None.

a = ring_profile()
a.load(7, hbt.frange(0,7),key_radius='full').flatten().plot()
plt.show()

# This is a good profile for 8-15. Shows moonlet belt really well. Looks like 0-7 one. Has an extra bump = stray.
# We could consider masking they stray. Right now I just use horizontal and vertical regions to remove it.
# But, that would be a big task.
# I guess I could just pass an additional argument, like 'c20(120,200)' for a 20-pixel radius circle.

a = ring_profile()
# a.load(7, hbt.frange(8,15),key_radius='full').\
#     remove_background_radial(radius_bg_halo_core,do_plot=False).flatten().plot()
a.load(7, hbt.frange(8,15),key_radius='full').flatten().plot()
plt.show()

# Decent profile for 16-23.

a = ring_profile()
a.load(7, hbt.frange(16,23),key_radius='full').remove_background_radial(radius_bg_117,do_plot=False).flatten().plot()
plt.show()


# Decent profile for 24-31.
# It is a very narrow moonlet belt between A+M, plus some inner stuff. 

a = ring_profile()
a.load(7, hbt.frange(91,93),key_radius='full').plot()

a.dn2iof()
a.plot()
a.remove_background_radial(radius_bg_halo_core).flatten().plot()

plt.show()

# Marginal profile for 32-35.
# There is more stray light here. It really needs masking.
# I did increase the polynomial power from 2 -> 5, which helped (like it should have).
# The inner region is going to be hard to get. I might do better if I mask it.

a = ring_profile()
a.load(7, hbt.frange(32,35),key_radius='core').remove_background_radial(radius_bg_127,do_plot=False).flatten().plot()

# OK profile for 36-39. Similar to previous. Thet moonlet belt is visible but the dust is probably subtracted out here.

a = ring_profile()
radius_bg_127 = np.array([[125,126], [130, 131]])*1000
a.load(7, hbt.frange(36,39),key_radius='full').remove_background_radial(radius_bg_127,do_plot=False).flatten().plot()

# Profile for 40-42. Still working on this one.

a = ring_profile()
radius_bg_127 = np.array([[125,126], [130, 131]])*1000
a.load(7, hbt.frange(40,42),key_radius='full').remove_background_radial(radius_bg_127,do_plot=False).flatten().plot()

plt.show()


a_dn = a.copy()
a_dn.dn2iof().plot()
plt.show()

a = ring_profile()
a.load(7, hbt.frange(91,93),key_radius='core').plot()
plt.show()

a_dn = a.copy()
a_dn.dn2iof().plot()
plt.show()
