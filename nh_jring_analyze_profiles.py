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

from   scipy import signal, fftpack

import spiceypy as sp
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs
from   astropy.wcs import WCS
from   astropy import units as u           # Units library
from   astropy import constants as c
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   scipy.stats import mode
from   scipy.stats import linregress
from   astropy.visualization import wcsaxes
import time
from   scipy.interpolate import griddata
from   importlib import reload            # So I can do reload(module)

from   astropy.convolution import Box1DKernel, Gaussian1DKernel, convolve
from   pymiecoated import Mie

from   scipy.stats import linregress

import re # Regexp
import pickle # For load/save

from   matplotlib.figure import Figure
from   scipy import interpolate

# HBT imports

import hbt

# Local NH rings imports

from  nh_jring_mask_from_objectlist            import nh_jring_mask_from_objectlist

from  nh_jring_mask_from_objectlist             import nh_jring_mask_from_objectlist
from  nh_jring_unwrap_ring_image                import nh_jring_unwrap_ring_image

from  nh_jring_extract_profile_from_unwrapped   import nh_jring_extract_profile_from_unwrapped

from  scatter_mie_ensemble                      import scatter_mie_ensemble
from  area_between_line_curve                   import area_between_line_curve

# For fun, I'm going to try to do this whole thing as a class. That makes it easier to pass 
# params back and forth, I guess. [A: Yes, this is a great place for a class -- very clean.]

class ring_profile:

    """
    This is a class to store ring profiles. They may be printed, aligned, loaded, averaged, etc.
    """    
    
# =============================================================================
#     Init method: load the main table file
# =============================================================================
    
    def __init__(self):

        file_pickle = 'nh_jring_read_params_571.pkl'     # Filename to read to get filenames, etc.
        dir_out     = '/Users/throop/data/NH_Jring/out/' # Directory for saving of parameters, backplanes, etc.
            
        self.dir_out = dir_out

        self.is_flattened = False
        
        lun = open(self.dir_out + file_pickle, 'rb')
        self.t = pickle.load(lun)                        # Self.t is the *entire* table for all J-ring obs, not subset.
        lun.close()
        
        
        # Process the group names. Some of this is duplicated logic -- depends on how we want to use it.
        
        self.groups = astropy.table.unique(self.t, keys=(['Desc']))['Desc']

        stretch_percent = 90    
        self.stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

        # Define a few constants within the method. We're going to keep using these in many places,
        # so best to set them up here. As per PEP-8, constants should be ALL_CAPS.
        
        self.A_METIS    = 127979.8*u.km        # Orbital distance, in km. From SCW07
        self.A_ADRASTEA = 128980.5*u.km        # Orbital distance, in km. From SCW07
        

# =============================================================================
# Define the 'string' for the class. This is the human-readable value returned when we do   print(ring)
# =============================================================================

    def __str__(self):
        
        if self.is_flattened:
           return("{}/{}".format(self.index_group_arr[0], 
               self.index_image_arr[0]))

        else:
            return("{}/{}-{}".format(self.index_group_arr[0], 
               self.index_image_arr[0],
               self.index_image_arr[-1])) 

# =============================================================================
# Define the 'string' for the class, which is what is shown when we type its name on the commandline.
#   This in theory should be a string to reconstruct the object... hopeless, here,
#   so we just return the string instead.
# =============================================================================

    def __repr__(self):
        return(self.__str__())
                
# =============================================================================
# Define an __iter__ function so we can iterate over our individual radial profile
# =============================================================================

    def __iter__(self):
        
        # NOT WORKING YET (obviously)
        return(True)

# =============================================================================
# Define a __getitem__ method. This should nominally take an object which has 7 profiles in it, and extract just nth.
# I think I need this if I want to define __iter__, __next__, etc. in order to iterate over class.
# However, maybe they are not needed for basic usage of the class.        
# =============================================================================

    def __getitem__(self, index):
        # NOT WORKING YET
        obj = self.copy()
        obj.exptime_arr = obj.exptime_arr[index]
        obj.groups
        return obj
        
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
                      key_radius = 'core',  # Which of the radial profiles to read? 'full', 'center', or 'core'
                                            # 'Core' means that only the inner 10% of the az range is used
                                            # for rad profile. Fewer data points, but less smeared out
                                            # if it is mis-aligned.
                      key_azimuth = 'net',  # Of the azimuthal profiles, which one to read?
                                            # 'net' means the central region, minus (inner+outer)/2.
                      verbose = False,      # Verbose: List the filename as loading?
                      **kwargs):
        """
        This loads a set of profiles from disk, into the object.
        It does not append in memory. It overwrites anything already loaded into that object.
        
        Parameters
        ----
        
        key_radius:
            String. Which of several pre-extracted radial profiles to use.
            - 'core'   : Uses the central 10% of the full azimuthal range.
            - 'center' : Uses the central 25% of the full azimuthal range.
            - 'full'   : Uses all of the data.
            
        verbose: 
            Boolean. List each file as it loads?
        
        key_azimuth:
            String. Which of several pre-extracted azimuthal profiles to use. Almost always 'net' is best.
            
        """
        
        # Each file on disk has several different extractions:
        #
        #   Radial profile has 'core',  'full',  'half', etc.
        #   Az     profile has 'inner', 'outer', 'net',  etc. 
        #
        # Since we usually need just one of these, this routine *only* reads one.
        # The individual one read can be controlled with key_radius and key_azimuth
        #
        # Usually the values of the profiles will be read in units of DN. 
        
        import humanize
        
        # Define t_group, since we'll reference it a lot
        
        self.groupmask = self.t['Desc'] == self.groups[index_group]
        self.t_group = self.t[self.groupmask]
        
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
        self.et_arr                 = []
        self.image_unwrapped_arr    = []
        self.mask_objects_unwrapped_arr = []
        self.mask_stray_unwrapped_arr = []
        
        self.profile_azimuth_units  = 'DN per pixel'
        self.profile_radius_units   = 'DN per pixel'

        for index_image in index_images:
            file = self.get_export_analysis_filename(index_group, index_image)
            
            lun = open(file, 'rb')
            vals = pickle.load(lun)
            lun.close()

            if verbose: 
                print("Loaded pickle file {}/{} {}".format(index_group, index_image, file))
            
            (image_unwrapped,                # Unwrapped image itself
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
                        index_image,         # Index of imagep
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
            self.et_arr.append(et)
            self.dt_str_arr.append(dt_str)
            self.image_unwrapped_arr.append(image_unwrapped)
            self.mask_stray_unwrapped_arr.append(mask_stray_unwrapped)
            self.mask_objects_unwrapped_arr.append(mask_objects_unwrapped)

            # Define a new property, which is a shift to be applied to each curve. Shift is in radial bins (ie, pixels)
            
            self.shift_bin              = np.zeros(self.num_images_group,dtype=int)

            if verbose:
                print("Read image {}/{}, phase = {:0.1f}°, {}, {}, {}".format(
                    index_group, index_image, ang_phase*hbt.r2d, file_short, bg_method, bg_argument))
                print("{}".format(file))
        
        print(f'Loaded {self.__str__()}: {len(index_images)} images')
        
        #==============================================================================
        # Now put these into arrays (not lists). Ideally we'd put these into an astropy table (so we can sort, etc.)
        #==============================================================================
        # NB: Doing np.array() on a list of 1D arrays will convert it to a 2D array, *iff* arrays are of same length!
        
        self.ang_phase_arr = np.array(self.ang_phase_arr)
        self.ang_elev_arr  = np.array(self.ang_elev_arr)
        self.azimuth_arr   = np.array(self.azimuth_arr)
                 
        self.profile_radius_arr  = np.array(self.profile_radius_arr)
        self.profile_azimuth_arr = np.array(self.profile_azimuth_arr) 
        
        return self

# =============================================================================
# Return number of radial bins
# =============================================================================
    
    @property         # Property -- allow it to be num_bins_radius, not num_bins_radius()    
    def num_bins_radius(self):
        
        return np.shape(self.profile_radius_arr)[1]

# =============================================================================
# Return size of azimuthal array
# =============================================================================
    
    @property    
    def num_bins_azimuth(self):
        
        return np.shape(self.profile_azimuth_arr)[1]

 
# =============================================================================
# Return number of profiles
# =============================================================================
    
    @property      # Putting the 'property' decorator here means that this func appears like 
                   # ring.num_profiles, not ring.num_profiles(). It looks like a property, not
                   # the function that it really is.
                   
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
        
        TEXP_2D = np.transpose(np.tile(self.exptime_arr, (self.num_bins_radius,1)))  # Increase this from 1D to 2D
        
        I = C / TEXP_2D / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. All v similar, except for spectrum assumed.
        
        # Apply Hal's conversion formula from p. 7, to compute I/F and print it.
        
        IoF = math.pi * I * r**2 / F_solar # Equation from Hal's paper
        
        # Now convert to 'normal I/F'
        
        # Define mu = cos(e), where e = emission angle, and e=0 is face-on.
        
        e = math.pi/2 - self.ang_elev_arr # NB: I made an error in this for the Pluto ring analysis! 
        
        mu = np.cos(e)
        
        mu_2D = np.transpose(np.tile(mu, (self.num_bins_radius,1)))
        
        # Calculate the normal I/F
        
        IoF_normal = 4 * mu_2D * IoF
        
        # Save this result into 'self,' replacing the DN values with the normal I/F values
        
        self.profile_radius_arr = IoF_normal
        
        # Change the units
        
        self.profile_radius_units = 'Normal I/F'
        
        # Do the same for azimuthal profile

        TEXP_2D = np.transpose(np.tile(self.exptime_arr, (self.num_bins_azimuth,1)))
        mu_2D = np.transpose(np.tile(mu, (self.num_bins_azimuth,1)))
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

    def smooth_azimuthal(self, width=1, kernel=None):

        if not width:     # If missing, None, 0, etc.
            return self
        
        if (kernel is None):
            kernel = Gaussian1DKernel(width)
                
        num_profiles = self.num_profiles

        for i in range(num_profiles):
        
            # Smooth over azimuth. Do this carefully: wrap the boundary conditions, but any NaNs in the input,
            # remain NaN in the output. This is so we don't make up any data in the wide unobserved gaps.
            
            self.profile_azimuth_arr[i,:] = convolve(self.profile_azimuth_arr[i,:], 
                                    kernel, 
                                    boundary='wrap', 
                                    preserve_nan=True)
#                                    p)

        return self
    
#            for key in self.profile_azimuth_arr[i]:
#                profile = convolve(self.profile_azimuth_arr[i][key], kernel)
#                self.profile_azimuth_arr[i][key] = profile
                
# =============================================================================
# Flatten the profiles -- that is, if there are N profiles, flatten to one, with mean value
# =============================================================================

    def flatten(self, a = 127_900*u.km, et=0):
        
        """ Flatten the profile, from N profile, into one mean (or median) profile.
            Applies to radial and azimuthal and all other quantities.
            Result is returned, and internal values are changed as well.
            After flattening, most outputs are as 1-element arrays (not scalars).
            
            NB: This always unwinds, as well.
            
            Parameters
            ---    
            
            a:
                Reference distance for doing the unwinding in the azimuthal direction. Mandatory.
                Astropy units. e.g., 127,500 km.

            et:
                Reference epoch to unwind to.
                
            Return value:
                self : Newly flattened object, with one profile, rather than N.
                
        """
        
        # Set up the reference distance
        
        a_ref = a
        
        et_ref = et
        
        # If is is already flattened (or a single profile), return without changing anything.
        # Flag as a warning, but do not cause error.
        
        if (self.num_profiles == 1):
            warnings.warn('Warning: Profile has already been flattened')
            return self
        
        # Do a check to verify that arrays all have the same shape. 

        try:
            _ = np.stack(self.profile_radius_arr)    # Attempt to convert list of 1D arrays, into a 2D array.
#            _ = np.stack(self.profile_azimuth_arr)
        except ValueError:
            raise ValueError("Cannot flatten because profiles are of different sizes.")
            
        # First flatten radial profiles. This is easy -- we just take the mean of all the profiles.
                
        self.profile_radius_arr[0] = np.nanmean(self.profile_radius_arr, axis=0) # Save the flattened profile

        # Now do the azimuthal profiles. This is much harder, since we need to unwind them to make sense of them.

        ### XXX We can now dramatically simplify this az stacking, now that every file uses the same
        #   gridding, at 0.001 radians.
        
        # Set up a uniformly spaced azimuthal grid, going 0 .. 2pi in 0.001 rad increments.
       
        rad_per_pix       = self.azimuth_arr[0][1] - self.azimuth_arr[0][0]
        num_pts_az        = int(2*math.pi / rad_per_pix)    
        az_arr            = np.array(range(num_pts_az)) * rad_per_pix
        
        # azimuth_unwind    = hbt.frange(0,math.pi *2, num_az_unwind+1)[0:-1]
        azimuth_unwind         = az_arr.copy()
 
        # Also, set up a 2D grid, where we make an 'image' of the azimuthal output
        
        unwind_2d = np.zeros( (self.num_profiles, num_pts_az) )
        unwind_2d[:] = np.nan
        
        # Loop over every azimuthal profile. For each one, unwind it properly to a common ET, then put it on a 
        # common grid. Then merge all of them.
            
        for i,profile_azimuth in enumerate(self.profile_azimuth_arr):
            
            # Unwind the longitudes. self.azimuth_arr is in radians, and covers only the observed longitudes
            #   (ie, not 0 .. 2pi)
            
            theta_i = unwind_orbital_longitude(self.azimuth_arr[i], self.et_arr[i], 'Jupiter', a_ref, et_ref)
            
            # Quantize unwound to 0.001 radian steps
            
            theta_i = (theta_i / rad_per_pix).astype(int) * rad_per_pix
            
            # Calculate the offset
            
            offset_pix = np.where(theta_i[0] > az_arr)[0][-1]
            
            # Load the unwound 1D array into the 2D array
    
            unwind_2d[i, 0:hbt.sizex(profile_azimuth)] = profile_azimuth
            unwind_2d[i,:] = np.roll(unwind_2d[i,:], offset_pix)
             
        # Save the 2D azimuth 'image' 
        
        self.profile_azimuth_arr_2d = unwind_2d

        # Flatten the 'image' of azimuth down into an array of size (1,3600), for instance.
        # This is just so the profile can be extracted consistently.

        warnings.simplefilter(action = "ignore", category = RuntimeWarning)  # "<" causes a warning with NaN's...
        
        self.profile_azimuth_arr = np.array([np.nanmean(unwind_2d,        axis=0)])    # Save the flattened profile
        self.profile_radius_arr  = np.array([self.profile_radius_arr[0]])

        # For statistics, count up how many datapoints went into each az bin.
        
        self.num_profile_azimuth_arr = np.array([np.sum(np.logical_not(np.isnan(unwind_2d)), axis=0)])
        
        warnings.simplefilter(action = "default", category = RuntimeWarning)  # "<" or 'less' causes warning w/ NaN's...
        
        # For some quantities, like angles, take the mean and save that.
        # NB: If we add any new fields to the object, we need to add a corresponding new line to this method!
        
        self.ang_elev_arr    = np.array([np.mean(self.ang_elev_arr)])
        self.ang_phase_arr   = np.array([np.mean(self.ang_phase_arr)])
        self.exptime_arr     = np.array([np.mean(self.exptime_arr)])
        
        # For other fields, take one entry, and save that.
        # This is a bit crude, but probably OK.
                        
        self.azimuth_arr     = np.array([azimuth_unwind])  # The newly unwound az grid values
        self.radius_arr      = np.array([self.radius_arr[0]])
        self.index_image_arr = np.array(["{}-{}".format(self.index_image_arr[0], self.index_image_arr[-1])])
        self.index_group_arr = np.array([self.index_group_arr[0]]) 
        self.dt_arr          = np.array([self.dt_arr[0]])
        self.dt_str_arr      = np.array([self.dt_str_arr[0]])
        
        # Flatten the images and mask arrays, by taking median (that is, reducing from N images, to 1.)
        # XXX Stacking images is hard and requires more subtlety to unwrap. Not doing this yet.
        
#        self.image_unwrapped_arr = np.array([np.nanmedian(self.image_unwrapped_arr, axis=0)])
#        self.mask_objects_unwrapped_arr = np.array([np.nanmedian(self.mask_objects_unwrapped_arr, axis=0)])
#        self.mask_stray_unwrapped_arr = np.array([np.nanmedian(self.mask_stray_unwrapped_arr, axis=0)])
        
        self.image_unwrapped_arr = np.array([self.image_unwrapped_arr[0]])
        self.mask_objects_unwrapped_arr = np.array([self.mask_objects_unwrapped_arr[0]])
        self.mask_stray_unwrapped_arr = np.array([self.mask_stray_unwrapped_arr[0]])
        
        
        self.is_flattened    = True
        
        return self          # The value is returned, so that commands can be chained.
                             # Also, the value of all quantities (e.g., profile_azimuth_arr) is changed internally.

# =============================================================================
# Copy the object
# =============================================================================

    def copy(self):
        
        """ 
        Copy an object to a new object (not a reference).
        This is useful do to before flattening, so that we can retain a copy of the original unflattened object.
        """
        
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
        
        for i in range(self.num_profiles):  
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
# Plot the radial profile
# =============================================================================

    def plot_radial(self, 
                    dy = 0, 
                    plot_sats = True,
                    smooth = None,
                    title = None,
                    xlim = (120000, 131000),
                    plot_legend = False,
                    **kwargs):
        """
        Just a helper routine to plot the radial profile.
        """
        
        for i,index_image in enumerate(self.index_image_arr):

            radius = self.radius_arr[i]
            profile_radius = self.profile_radius_arr[i]
            
            p = plt.plot(radius, 
                     np.roll( (i * dy) + profile_radius, self.shift_bin[i]),  
                     label = r'{}/{}, $\alpha$={:.2f}°, e={:.2f}°, {}'.format(
                                   self.index_group_arr[i], 
                                   self.index_image_arr[i], 
                                   self.ang_phase_arr[i]*hbt.r2d,
                                   self.ang_elev_arr[i]*hbt.r2d,
                                   self.dt_str_arr[i]),
                     **kwargs)

            # Set the xlimit explicitly. This is stupid that it cannot be passed as a kwarg...
            
            if (xlim):
                plt.xlim(xlim)
               
        # Plot Metis and Adrastea on the radial plot
        
        if (plot_sats):

#                args = {'linestyle': 'dash', 'alpha':0.5, 'color':'black'}
            plt.axvline(x=self.A_METIS.to('km').value/1000,    linestyle='dashed', alpha=0.2, color='black')
            plt.axvline(x=self.A_ADRASTEA.to('km').value/1000, linestyle='dashed', alpha=0.2, color='black')
            
        plt.xlabel('Radial Distance [km]')
        plt.ylabel(self.profile_radius_units)
 
        # Set scientific notation on the Y axis, if appropriate
        
        axes = plt.gca()
        
        if ('I/F' in self.profile_radius_units):       
            axes.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
        
        if plot_legend:
            plt.legend(loc='upper left', fontsize=8)
        if (title is not None):
            plt.title(title)    
        plt.show()

# =============================================================================
# Plot the azimuthal profile 
# =============================================================================

    def plot_azimuthal(self, 
                    dy = 3, 
                    plot_sats = True,
                    smooth = None,
                    title = None,
                    xlim = None,
                    plot_legend = False,
                    **kwargs):
        """
        Plot the azimuthal profile.
        
        Plot is in degrees. Internally, everything is stored in radians.
        
        Parameters
        -----
        
        smooth:
            Integer. If set, this is the width of a Gaussian used for smoothing the data.
        """
        
        #==============================================================================
        # Make a plot of azimuthal profile. Don't flatten -- just plot all the data
        #==============================================================================
        
        hbt.figsize((18,5))
        
        # Loop over all the az profiles we have, and plot each one.
        # If array is flattened, there will be just one. If raw, then there might be 30 to do.

        if (smooth):
            kernel = Gaussian1DKernel(smooth)
        
#        for i,profile_azimuth in enumerate(self.profile_azimuth_arr):
#            y = profile_azimuth
#         
#            if smooth:
#                y = convolve(self.profile_azimuth_arr[i,:], kernel)

        num_profiles = self.num_profiles

        for i in range(num_profiles):
        
            # Loop over azimuthal profiles
            
            y = self.profile_azimuth_arr[i]
         
            if smooth:
                y = convolve(y, kernel, boundary='wrap', preserve_nan=True)

#            plt.plot(self.azimuth_arr[i]*hbt.r2d, 
            plt.plot(self.azimuth_arr[i]*1000, 
                     (i * dy) + y,  
                     label = '{}/{}, {:.2f}°'.format(
                             self.index_group_arr[i], 
                             self.index_image_arr[i], 
                             self.ang_phase_arr[i]*hbt.r2d),
                     **kwargs)

        # Set the xlimit explicitly. This is stupid that it cannot be passed as a kwarg...
            
        if (xlim):
            plt.xlim(xlim)
        
        # Calculate and plot the Metis and Adrastea positions.
        
        if plot_sats:
            
            pass
        
            # Look up sub-satellite longitude at time of the 'current' observation
            # Unwind this longitude back to 2000, using Metis or Adrastea's orbital distance
            # Plot it
        
            # XXX Right now I'm not really sure of the point of this, so we do not do it.
            
        plt.xlabel('Azimuth [mrad]')
        plt.ylabel(self.profile_azimuth_units)
        if plot_legend:
            plt.legend(fontsize=8)
        if (title is not None):
            plt.title(title)
        plt.show()
            
        hbt.figsize_restore()

        return self

# =============================================================================
# Plot the unwrapped azimuthal profile image
# This is essentially the raw profiles that goes into the flattened az profile        
# =============================================================================

    def plot_azimuthal_2d(self):
    
        dx = 2 * math.pi
        dy = self.num_profiles
        aspect = dx / dy * 8
        
        plt.imshow(self.profile_azimuth_arr_2d, aspect = aspect, origin = 'lower')
        
        plt.title(self)
        plt.xlabel('Azimuth [mrad]')
        plt.ylabel('Image number')
        
        plt.show()
        
        return self
    
# =============================================================================
# Plot the radial and azimuthal profile together
# =============================================================================
        
    def plot(self,    plot_radial    = True, 
                      plot_azimuthal = True, 
                      dy_radial      = 0,  # Default vertical offset between plots 
                      plot_sats      = True,  # Plot the location of Metis / Adrastea on radial profile?
                      smooth         = None,  # Width of any smoothing to apply, in pixel. For Az profile only. 
                      title          = None,      # String to print above the plots
                      xlim_a         = (120000,131000),     # Limit for radial distance
                      xlim_az        = None, 
                      plot_legend    = False,
#                      a_unwrap       = None,      # Reference distance for unwrapping
                      **kwargs):   
                              
        if (plot_radial):
            self.plot_radial(plot_sats=plot_sats, 
                             title=title, 
                             xlim=xlim_a, 
                             plot_legend=plot_legend,
                             smooth=smooth)

        if (plot_azimuthal):
            self.plot_azimuthal(plot_sats=plot_sats, 
                                title=title, 
                                xlim=xlim_az, 
                                plot_legend=plot_legend,
                                smooth=smooth)

# =============================================================================
# Calculate the radial area under a curve
# =============================================================================

    def area_radial(self, limit):
 
        """
        Calculate the radial area under a ring curve. 
        
        Parameters
        -----
        
        limit:
            (r_inner, r_outer): Tuple defining the x values of the radial limits
            
        """
        
        area = []
        
        for i in range(self.num_profiles):  # I guess this is where I should use __iter__
        
            (area_out, _, _) = area_between_line_curve(self.radius_arr[i], self.profile_radius_arr[i], limit)
            area.append(area_out)                 
            
        return area


# =============================================================================
# Extract strips from a set of images, showing only the rings
# =============================================================================

    def make_strip_mosaic(self, a_orbit=127_900*u.km, gap_y_pix=1, 
                          y0_extract = 415, dy_extract=50, do_plot=False, do_plot_masked=False,
                          do_plot_profile = False,
                          do_unwind = True, dwidth_chop = 0, xlim=None, et_ref = 0):
        """
        Create an image of individual strips, by extracting the ring region from many individual images.
        
        Optional parameters
        -----
        
        gap_y_pixel:
            Gap betwen strips, in the y direction.
            
        y0_extract:
            Y value, in pixels, of center of extraction region. **I should rewrite to pass in km, not pix.**
            
        dy_extract:
            Extent, in y dir, of extraction region. Full height.
        
        do_unwind:
            Do we unwind in longitude to a common epoch, or print in terms of observed System III longitude?
            
        width_chop:
            The width (on each end) to chop az profile by. This is to remove vignetting, rather than using whole width.
            Pixels.
            
        a_orbit:
            Orbital distance used for unwinding.
        
        gap_y_pix:
            Gap in pixels, bewteen successive observations in the output mosaic.

        do_plot:
            Boolean. If set, make some plots.
            
        et_ref:
            Epoch to unwind longitude to.
             
            
        Return values:
        ----
        
        Tuple:   `((im_mosaic, mask_mosaic), (profile_im, profile_masked), az_arr)`


            
        """  
        # Initialize things
        
        rad_per_pix       = self.azimuth_arr[0][1] - self.azimuth_arr[0][0]
        width_mosaic_pix  = int(2*math.pi / rad_per_pix)    
        height_mosaic_pix = self.num_profiles * (dy_extract)
        az_arr            = np.array(range(width_mosaic_pix)) * rad_per_pix
        stretch           = astropy.visualization.PercentileInterval(95)
       
        # Create the output arrays
        
        im_mosaic        = np.zeros((height_mosaic_pix, width_mosaic_pix))  # Images are y, x
        im_mosaic[:,:]   = np.nan
        mask_o_mosaic    = im_mosaic.copy()
        mask_s_mosaic    = im_mosaic.copy()
        
        im_strip         = np.zeros((dy_extract, width_mosaic_pix))
        im_strip_arr     = np.zeros((self.num_profiles, dy_extract, width_mosaic_pix))
        
        mask_strip       = im_strip.copy()
        mask_strip_arr   = im_strip_arr.copy()
        
        lon0_unwind_rad = np.zeros(self.num_profiles)
        lon0_unwind_pix = np.zeros(self.num_profiles)
                 
        # Loop and calculate the appropriate shift, for each image
                                                     
        for j in range(self.num_profiles):
        
            lon0_unwind_rad[j] = unwind_orbital_longitude(self.azimuth_arr[j][0], self.et_arr[j], 
                            'Jupiter', a_orbit=a_orbit, et_ref = et_ref)
        
            lon0_unwind_pix[j] = (lon0_unwind_rad[j] - self.azimuth_arr[j][0]) / rad_per_pix
        
        # Loop and add each strip to the output extraction array
            
        for j in range(self.num_profiles):
            
            # Calculate where the start should be
            
            x0 = np.where( np.mod(self.azimuth_arr[j][0], 2*math.pi) <= az_arr)[0][0] - 1
        
            # If asked to unwind the longitude, apply the appropriate additional correction 
            
            if do_unwind:
                x0 += int(lon0_unwind_pix[j])
                
            im   = self.image_unwrapped_arr[j][y0_extract-int(dy_extract/2):y0_extract+int(dy_extract/2)-gap_y_pix]
            
            # If requested, chop off the ends of the radial profiles, because many of them are vignetted.
            
            if (dwidth_chop):
                im = im[:, dwidth_chop : -dwidth_chop]  # Make it less wide
                x0 += dwidth_chop                           # Advance position of it    
            
            # Place the individual profile into the output array, in right vertical position
        
            im_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, 0:hbt.sizey(im)] = im
        
            # Roll it into place horizontally, by shifting by an appropriate amount. 
            # It will properly roll across the edge at 2pi and back to the beginning.
            
            im_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, :] = \
              np.roll(im_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, :], x0, axis=1)
        
            # Roll it into place horizontally, by shifting by an appropriate amount. 
            # It will properly roll across the edge at 2pi and back to the beginning.
        
            # Now do the same again, for the masks.
            # mask_o = object mask (from SPICE)
            # mask_s = stray light mask (from Photoshop), *and* object mask, combined
            
            mask_o = \
              self.mask_objects_unwrapped_arr[j][y0_extract-int(dy_extract/2):y0_extract+int(dy_extract/2)-gap_y_pix]
            mask_s = \
              self.mask_stray_unwrapped_arr[j][y0_extract-int(dy_extract/2):y0_extract+int(dy_extract/2)-gap_y_pix]
            
            if (dwidth_chop):
                mask_o = mask_o[:, dwidth_chop : -dwidth_chop]  # Make it less wide
                mask_s = mask_s[:, dwidth_chop : -dwidth_chop]  # Make it less wide
        
            mask_o_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, 0:hbt.sizey(mask_o)] = mask_o
            mask_s_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, 0:hbt.sizey(mask_s)] = mask_s
        
            mask_o_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, :] = \
              np.roll(mask_o_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, :], x0, axis=1)
        
            mask_s_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, :] = \
              np.roll(mask_s_mosaic[(j*dy_extract):((j+1)*dy_extract)-gap_y_pix, :], x0, axis=1)
         
            # Now grab these rolled images, and save them in a place where will stack them later. A single long strip.
            
            mask_strip_arr[j,:,:] = mask_s_mosaic[(j*dy_extract):((j+1)*dy_extract), :]
            im_strip_arr[j,:,:]   = im_mosaic[(j*dy_extract):((j+1)*dy_extract), :]
            
        mask_mosaic = mask_s_mosaic.copy()
        
        # Convert the mask from a Boolean mask, to a NaN mask. 1 = Good, NaN = bad, and nothing in between
        
        mask_mosaic[mask_mosaic != True] = np.nan   
        
        # Generate azimuthal profiles
        
        profile_im     = np.nanmedian(im_mosaic,               axis=0)
        profile_masked = np.nanmedian(im_mosaic * mask_mosaic, axis=0)
        
        # Construct a single long strip of the image
        
        mask_strip_arr[mask_strip_arr == 0] = np.nan   # Change 0 → nan
        im_strip     = np.nanmedian(im_strip_arr, axis=0)
        masked_strip = np.nanmedian(mask_strip_arr * im_strip_arr, axis=0)
        
        # Make plots, if requested
        
        aspect = hbt.sizey(im_mosaic) / hbt.sizex(im_mosaic)/5
        
        if do_plot:
            
            plt.imshow(stretch(im_mosaic), origin='lower', aspect = aspect)
            plt.title(self)
            plt.show()
        
        if do_plot_masked:
            plt.imshow(stretch(im_mosaic * mask_mosaic), origin='lower', aspect = aspect)
            plt.title(f'{self}, Masked')
            plt.show()
        
        if do_plot_profile:
            alpha = 0.5
            plt.plot(profile_im,     alpha=alpha, label = 'Raw')
            plt.plot(profile_masked, alpha=alpha, label = 'Masked')
            plt.legend()
            plt.title(f'{self}, Az profile')
            plt.xlim(xlim)
            plt.show()
            
        # Return results
        
        return((im_mosaic, mask_mosaic), (profile_im, profile_masked), (im_strip, masked_strip), az_arr)


# =============================================================================
# Unwind an orbital longitude back to a common time frame
# =============================================================================
#%%%
        
def unwind_orbital_longitude(lon_in, et_in, name_body, a_orbit = 127_900*u.km, et_ref=0):
    """
    Takes a set of longitudes, and unwind them back to a given frame, incorporating
    both body rotation, and keplerian orbital motion. Both are assumed to be constant.
    
    That is, it takes the preset longitudes lon_in at time et_in of a satellite 
    in orbit above a rotating body, and it maps them backward to determine their longitudes at time et_out.
    
    "The object is currently in orbit and is at this longitude above the body. 
      If we go back to ET_OUT, what longitude will it be at?"
      
    Or, as Showalter 2007 said,  
      "Longitudes are rotated to a common epoch, assuming the mean motion of a body at a = 129,000 km"
    
    Typically et_out will be some common reference, like '1 Jan 2000 12:00:00'. But, it can be anything.
    
    The calculation requires the central body mass and rotation rate. These are
    determined from SPICE.
    
    The orbital period is determined from the orbital distance `a_orbit`.
    
    SPICE is assumed to be loaded, and have all of the necessary kernels.
    
    This is a standalone function, not part of any class.
    
    Parameters
    -----
    
    lon_in:
        Longitude, in radians. The longitude of the observations.
    
    et_in:
        ET, in seconds. This is the time at which the observations were taken.
    
    name_body:
        String. Name of the central body, for which keplerian rotation will be calculated.
    
    Optional Parameters
    -----
    
    a_orbit: 
        Orbital distance, from which to calculate Keplerian motin.
        
    et_ref:
        ET, in seconds. This is the ET for which we want to put everything into the time base of (ie, the epoch).
        Default value is 0 (ie, roughly 1 Jan 2000 12:00:00, but about 64 sec from that, due to leap seconds etc.)
      
    """
    
    name_body = 'jupiter'
    r_jup = 69_911*u.km
#    a_orbit = 127_900*u.km
    
    masses = {'JUPITER': c.M_jup,   # Agrees w/ wikipedia. 1.89818e27 kg
              'EARTH':   c.M_earth,
              'PLUTO':   1.303e22*u.kg,
              'SUN':     c.M_sun}
    
    m_body = masses[name_body.upper()]
    
    dt = et_ref - et_in                # Time shift to apply. Typically will be negative, 
                                     # since we usually go back in time.
                                     
  
#    print(f'dt = {dt} sec')
    
    (_, pm) = sp.bodvrd(name_body.upper(), 'PM', 3)
    
#    v_orbit = (c.G * m_body / a_orbit)    
#    p_orbit = (2 * math.pi * a_orbit) / v_orbit
    
    
    p_orbit =  2*math.pi * (np.sqrt(a_orbit**3 / (c.G * m_body))).to('s')  # Orbital period at Metis: approx 7 hr.
    
    dtheta_dt_rot = pm[1] * hbt.d2r / (u.day.to(u.s))  # Convert rotation rate from deg/day (spice) to rad/sec (hbt)
                                                       # We ignore the pm[0] component since that's offset, and 
                                                       # it is already programmed into our original longitudes.

    dtheta_dt_orbit = (2*math.pi / p_orbit.to('s')).value
    
    # Get Jupiter J2, J4, and the equation for dtheta_dt, from 
    # http://astro.cornell.edu/academics/courses/astro6570/Rotation_precession_gravity_fields.pdf
    # Also see 
    # https://space.stackexchange.com/questions/25868/equation-for-orbital-period-around-oblate-bodies-based-on-j2
    
    j2 = 0.01473
    j4 = -587e-6
    mean_motion_orbit = dtheta_dt_orbit * hbt.r2d * 86400  # Get orbital mean motion, in deg/day
 
    # Apply the orbital perturbation. Note that this is kind of made up. The basic idea is from eq's in PDF below.
    # But I had to change exponents to fit right. I think in the end my results are right, and I was just unable
    # to find the proper expression. No one gives a closed form for mean motion as func of J2 + J4.
    # I bet I could ask DPH for it, but I haven't.
    
    if (name_body.upper() == 'JUPITER'):    
#        print(f'Mean motion at {a_orbit} is {mean_motion_orbit} deg/day (pre-J2)')
        dtheta_dt_orbit *= (1 + 3/2 * j2 * (r_jup / a_orbit)**(6/2))
         
                # From the Final-2 page at Rotaton_precession_gravity_fields PDF above.
                # The exponent I use is different than used in paper, so I am prob wrong.
                    
    mean_motion_orbit = dtheta_dt_orbit * hbt.r2d * 86400  # Get orbital mean motion, in deg/day
 
#    print(f'Mean motion at {a_orbit} is {mean_motion_orbit} deg/day')
#    print(f'Orbital period = {p_orbit/86400} d')
    
    lon_out = lon_in + (dtheta_dt_orbit - dtheta_dt_rot)*dt

    # Put in range 0 .. 2pi
    
    lon_out = np.mod(lon_out, 2*math.pi)
    
    # Remove any units, which might be left. Output should be just radians.
    
    lon_out = lon_out.value
    
    return lon_out

#%%%
    
def recut_strip(strip, dtheta = 1):
    
    """ Take a long unrolled strip plot, and cut it up into pieces for easier plotting.
    
    Assumed to be 1 mrad/pixel
    """
    
    num_strips = int(np.ceil(2*math.pi / dtheta))  # Round up
    
    dy_strip = hbt.sizey_im(strip)
    
    dx_out = int( np.ceil( hbt.sizex_im(strip) * (dtheta / (2*math.pi)) ) )
    dy_out = num_strips * dy_strip
    
    arr = np.zeros((dy_out, dx_out))
    
    # Now copy the strips in. 
    
    # First, do all but the final strip
    
    for i in range(num_strips-1):
        arr[i*dy_strip:(i+1)*dy_strip, 0:dx_out] = strip[:,i*dx_out:(i+1)*dx_out]
    
    i+=1
    
    strip_last = strip[:,(i*dx_out):-1]
    
    arr[i*dy_strip:(i+1)*dy_strip, 0:hbt.sizex_im(strip_last)] = strip_last
    
    return arr
    
    
    
#%%%
    
def merge_mosaic_strip(mosaic, strip, gap_y_pix = 1):
    """
    Make a meger of the mosaic, and the strip.
    The 1D strip is repeated over and over. The mosaic is placed on top of that.
    
    Typical sizes: strip  = (100,6300),   since it is one azimuthal map of the ring
                   mosaic = (5000, 6300)  since it is 50 individual (and incomplete) azimuthal maps of the ring
                   
    This is a regular function -- not part of the class.

    Parameters
    -----

    mosaic: Typical (n x 100) x 6300

    strip: Typically 100 x 6300     
           
    """
           
    ###
    
    stack = mosaic.copy()
    dy = hbt.sizey_im(strip)
    
    num_strips = int(hbt.sizey_im(mosaic) / hbt.sizey_im(strip))

    for i in range(num_strips):
        stack[i*dy:(i+1)*dy,:] = strip
    
    # Make an edge mask
    
    out = mosaic.copy()
    is_nan = np.isnan(out)
    is_val = np.logical_not(np.isnan(out))
    out[is_nan] = 0
    out[is_val] = 1  # out is now 0 for off-mosaic, and 1 for on mosaic
    kernel = astropy.convolution.Gaussian2DKernel(1)
    out_c = convolve(out, kernel) # 10 sec to execute
    mask = np.logical_and((out_c > 0), np.logical_not(out))  # True for good pixels
    is_edge = (mask - out_c) > 0
    
    # What I want is a mask of pixels which are off the mosaic, but have out_c > 0
    # I'll just set these values to np.nan.
    
    is_good = np.logical_not(np.isnan(mosaic))
    
    # Copy the good pixels over
    
    stack[is_good] = mosaic[is_good]
    
    # Set the edge pixels (ie, within a few pix of the good ones) to nan, to mark a border
    
    stack[is_edge] = np.nan
    
    return stack
    
    ###  Jupiter Prime Meridian:      W       =  284.95  + 870.5366420 d"  pm[0] = 284; pm[1] = 870.
#%%%    
# =============================================================================
# # ===========================================================================
# # END OF CLASS AND FUNCTION DEFINITION ###################################################
# # ===========================================================================
# =============================================================================


# =============================================================================
# Now read in some rings data and plot it
# =============================================================================

pass

stretch = astropy.visualization.PercentileInterval(95)

# Do a test with Metis in the frame
    
a = ring_profile()
self = a
plot_azimuthal=True
plot_legend=True
a.load(8,hbt.frange(0,18), key_radius='full')
a.plot_azimuthal(plot_legend=True)
a.plot_radial()
a.plot(plot_legend=False, plot_azimuthal=True, title=a,a_unwind=129700*u.km)
a.load(8,hbt.frange(24,48), key_radius='full').plot(plot_legend=False, plot_azimuthal=True, plot_radial=False, title=a, 
      a_unwind=129700*u.km)



# Load SPICE, if it isn't started already
        
file_tm = 'kernels_nh_jupiter.tm'  # SPICE metakernel
sp.unload(file_tm)                 # Pre-emptively unload .tm file, to prevent polluting the kernel pool 
sp.furnsh(file_tm)
    
# Invoke the class
        
ring = ring_profile()

# Define the off-ring locations. Pixels in this region are used to set the background level.
# XXX THESE VALUE ARE NOT VERY USEFUL. 
# We do much better when we define them individually for each profile.

radius_bg_core      = np.array([[127.5,127.9], [129.0,129.7]])*1000  # Core only
radius_bg_main_core = np.array([[122,122.5],   [130,131]])*1000  # Cover full width of main ring and core
radius_bg_117       = np.array([[117,120],     [130,131]])*1000  # This range is a bit better - no edge effects
radius_bg_127       = np.array([[125,126],     [130, 131]])*1000

BINS_SMOOTH = 0

# Make a plot of all of the data

# Define the limits used for the I/F 'area under the curve' total equivalent width

limit_radial = (118000,1295000)  

# Set up output arrays

# Create an astropy table to dump results into. This is so we can make a phase curve in the end

# First create a table, one entry per image

t = Table([[], [], [], [], [], [], [], [], [], [], []],
    names=('index_group', 'index_image', 'phase', 'elev', 'area_dn', 'exptime', 'label', 'sequence', 'radius', 
           'profile_radius', 'profile_radius_units'),
    dtype = ('int', 'int', 'float64', 'float64', 'float64', 'float64', 'U30', 'U30', 'object', 'object', 'U30'))

# List all of the data that we want to read in from disk

params        = [(7, hbt.frange(0,7),   'full'),  # For each plot, we list a tuple: (index_group, index_image, key)
                 (7, hbt.frange(8,15),  'full'),
                 (7, hbt.frange(16,23), 'full'),
                 (7, hbt.frange(24,31), 'full'),
                 (7, hbt.frange(32,35), 'full'),
                 (7, hbt.frange(36,39), 'full'), # Lots of 
                 (7, hbt.frange(40,42), 'full'),
                 (7, hbt.frange(52,54), 'full'),
                 (7, hbt.frange(55,57), 'full'),                 
                 (7, hbt.frange(58,60), 'full'),                 
                 (7, hbt.frange(61,63), 'full'),
                 (7, hbt.frange(91,93), 'full'),
                 (7, hbt.frange(94,96), 'full')]

# Loop over each of the sets of images
# For each image set, sum the individual profiles, and get one final profile

radius = []
profile_radius = []  # This is a list.

radius_bg = radius_bg_main_core

for param_i in params:

    # Unwrap the tuple, describing this set of files -- e.g.,  (7, array([94, 95, 96]), 'full')
    
    (index_group, index_images, key_radius) = param_i   # Unwrap the tuple

    # Create a string to describe this set of images -- e.g.,    '7/94-96 full'
    
    sequence = ("{}/{}-{} {}".format(index_group, np.amin(index_images), np.amax(index_images), key_radius))

    # Load the profile from disk.
    
    ring = ring_profile()
    ring.load(index_group, index_images, key_radius = key_radius).smooth_azimuthal(BINS_SMOOTH)
    
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

# Set the core and main ring positions. I have individually measured each of these.

radius_core_default = np.array([127.5, 129.5])*1000   # Radius, in units of 1000 km.
radius_main_default = np.array([118.0, 130.0])*1000

t_mean.add_column(Table.Column(np.tile(radius_core_default, (num_sequences,1)),name='radius_core'))
t_mean.add_column(Table.Column(np.tile(radius_main_default, (num_sequences,1)),name='radius_main'))
t_mean.add_column(Table.Column(np.tile((0,0), (num_sequences,1)), name='bin_core'))
# Core radii, used
t_mean.add_column(Table.Column(np.zeros(num_sequences), name='area_core'))
t_mean.add_column(Table.Column(np.tile((0,0), (num_sequences,1)), name='bin_main'))
t_mean.add_column(Table.Column(np.zeros(num_sequences), name='area_main'))

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
# Extract separately the brightness of the moonlet core vs the brightness of the main ring
# =============================================================================
# These radial limits have been set individually for each profile. I think they are very good.
t_mean['radius_core'][t_mean['sequence'] == '7/8-15 full']  = np.array((127.85, 129.5))*1000
t_mean['radius_core'][t_mean['sequence'] == '7/0-7 full']   = np.array((127.65, 129.35))*1000
t_mean['radius_core'][t_mean['sequence'] == '7/24-31 full'] = np.array((127.8,  129.5))*1000
t_mean['radius_core'][t_mean['sequence'] == '7/16-23 full'] = np.array((127.65, 129.55))*1000 # Red
t_mean['radius_core'][t_mean['sequence'] == '7/36-39 full'] = np.array((127.75, 129.25))*1000 # Purple
t_mean['radius_core'][t_mean['sequence'] == '7/32-35 full'] = np.array((127.8,  129.6))*1000   # Brown
t_mean['radius_core'][t_mean['sequence'] == '7/40-42 full'] = np.array((128.0,  129.70))*1000 # pink
t_mean['radius_core'][t_mean['sequence'] == '7/61-63 full'] = np.array((127.8,  129.50))*1000 # grey
t_mean['radius_core'][t_mean['sequence'] == '7/52-54 full'] = np.array((127.5,  129.50))*1000 # olive
t_mean['radius_core'][t_mean['sequence'] == '7/91-93 full'] = np.array((128.0,  129.55))*1000 # lt blue
t_mean['radius_core'][t_mean['sequence'] == '7/94-96 full'] = np.array((127.5,  129.3))*1000 # dk blue top

t_mean['radius_main'][t_mean['sequence'] == '7/8-15 full']  = np.array((118,    129.5))*1000
t_mean['radius_main'][t_mean['sequence'] == '7/0-7 full']   = np.array((118,    129.5))*1000
t_mean['radius_main'][t_mean['sequence'] == '7/24-31 full'] = np.array((118,    129.5))*1000
t_mean['radius_main'][t_mean['sequence'] == '7/16-23 full'] = np.array((118,    129.5))*1000 # Red
t_mean['radius_main'][t_mean['sequence'] == '7/36-39 full'] = np.array((118,    129.25))*1000 # Purple
t_mean['radius_main'][t_mean['sequence'] == '7/32-35 full'] = np.array((118.7,  129.6))*1000   # Brown
t_mean['radius_main'][t_mean['sequence'] == '7/40-42 full'] = np.array((118.5,  129.70))*1000 # pink
t_mean['radius_main'][t_mean['sequence'] == '7/61-63 full'] = np.array((119.2,  129.50))*1000 # grey
t_mean['radius_main'][t_mean['sequence'] == '7/52-54 full'] = np.array((118.5,  129.50))*1000 # olive
t_mean['radius_main'][t_mean['sequence'] == '7/91-93 full'] = np.array((118.5,  129.55))*1000 # lt blue
t_mean['radius_main'][t_mean['sequence'] == '7/94-96 full'] = np.array((120.5,  130.5))*1000 # dk blue top

# Now that we have the radii of the inner core calculated, look up the bin limits for the inner core

for i in range(num_sequences):
    radii       = t_mean['radius'][i]
    radius_core = t_mean['radius_core'][i]
    radius_main = t_mean['radius_main'][i]

    # Look up bin limits of the core and store them
    
    bins = [ np.argmin(np.abs(radii - t_mean['radius_core'][i][0])),
             np.argmin(np.abs(radii - t_mean['radius_core'][i][1])) ]
            
    t_mean['bin_core'][i] = bins

    # Look up bin limits of the main ring, and store them
    
    bins = [ np.argmin(np.abs(radii - t_mean['radius_main'][i][0])),
             np.argmin(np.abs(radii - t_mean['radius_main'][i][1])) ]
            
    t_mean['bin_main'][i] = bins

# Calculate the area under the curve for the core. 
# XXX We should probably write a method for this.
#     We should also improve things so that we can get an stdev from the component curves that make this up.
    
    (area_i_core, line, bins) = area_between_line_curve(t_mean['radius'][i],
                                     t_mean['profile_radius'][i],
                                     t_mean['radius_core'][i])
    t_mean['area_core'][i] = area_i_core
    
# Calculate the area under the curve for the main ring.
# For the main ring, we subtract off the area already measured for the core

    (area_i_main, line, bins) = area_between_line_curve(t_mean['radius'][i],
                                     t_mean['profile_radius'][i],
                                     t_mean['radius_main'][i])
    t_mean['area_main'][i] = area_i_main - area_i_core
    
# =============================================================================
# Do some Mie calculations
# =============================================================================

# Set up the size distribution

alam   = 500  * u.nm
rmin   = 0.01 * u.micron
rmax   = 50   * u.micron     
num_r  = 20

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
(scalefac_merged_log, _) = hbt.linfit_origin(p11_merged_data.value, t_mean['area_main'], log=DO_FIT_LOGSPACE)
(scalefac_merged_lin, _) = hbt.linfit_origin(p11_merged_data.value, t_mean['area_main'])
    
# =============================================================================
# Make a tall plot of Radial Profile (Summed, I/F) vs. Sequence
# =============================================================================

hbt.figsize((12,20))

DO_PLOT_CORE_LIMITS = True
DO_PLOT_MAIN_LIMITS = True

dy = 13e-7 # Set the y separation between adjacent lines

#hbt.figsize((10,10))
#plt.axes([0.2,0.4,0.7,0.9])

fig, ax1 = plt.subplots() # This returns a FIGURE object, plus 1 or more AXES objects.
ax2      = ax1.twiny()    # Create a second AXES, which is for the top (so we can plot both R_J and km)

ax2.locator_params(axis='x', numticks=3)
ax1.locator_params(axis='x', numticks=20)

for i,s in enumerate(t_mean['sequence']):  # Loop over the text sequence name (e.g., '7/0-7 full')
    
    # Plot the radial profile itself
    
    ax1.plot(t_mean['radius'][i]/1000, t_mean['profile_radius'][i] + i * dy,
             linewidth=3,
             label = r'{}, {:.1f}°, {:.1f}°'.format(
                  s, t_mean['phase'][i]*hbt.r2d, t_mean['elev'][i]*hbt.r2d) )
    
    # Plot the line denoting limits for main ring core
    
    if (DO_PLOT_CORE_LIMITS):
        bins = t_mean['bin_core'][i]
        ax1.plot( np.array( (t_mean['radius'][i][bins[0]], t_mean['radius'][i][bins[1]]) )/1000 ,
                  np.array( (t_mean['profile_radius'][i][bins[0]], t_mean['profile_radius'][i][bins[1]]))  + i*dy,
                  color = 'red')
    
    # Plot the line denoting limits for main ring full
                     
    if (DO_PLOT_CORE_LIMITS):
        bins = t_mean['bin_main'][i]
        ax1.plot( np.array( (t_mean['radius'][i][bins[0]], t_mean['radius'][i][bins[1]]) )/1000 ,
                  np.array( (t_mean['profile_radius'][i][bins[0]], t_mean['profile_radius'][i][bins[1]]))  + i*dy,
                  color = 'grey', linestyle='dotted')
                         
ax1.set_ylim((-1e-6,3e-6 + dy*i))

# Set the x axis, in 1000 km.

xlim_kkm = np.array((115,133))
ax1.set_xlim(xlim_kkm)
ax1.set_xlabel('Radius [1000 km]')
ax1.set_title('J-ring radial profiles, summed per sequence, binning={}'.format(BINS_SMOOTH), y=1.04) 
                                                            # Bump up the y axis
ax1.set_ylabel('I/F')
ax1.legend()
plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))

# Mark position of Metis + Adrastea

ax1.axvline(ring.A_METIS.to('km').value/1000, linestyle='dashed', color='grey', alpha=0.3)
ax1.axvline(ring.A_ADRASTEA.to('km').value/1000, linestyle='dashed', color='grey', alpha=0.3)

ax1.text(ring.A_METIS.to('km').value/1000, -5e-7, 'M')
ax1.text(ring.A_ADRASTEA.to('km').value/1000, -5e-7, 'A')

# Put the second axis on the top

ax2.plot([1,1])
ax2.set_xlabel('Radius [R_J]')
ax2.set_xlim(xlim_kkm * 1000 / 71492) # Create these automatically from values on the other X axis.
plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))

plt.show()

#%%%
    
# =============================================================================
# Plot the main phase curve, separated into contributions of Main ring + Core
# =============================================================================

hbt.figsize((12,8))

hbt.set_fontsize(size=15)

for i in range(len(t_mean['phase'])):
    plt.errorbar(t_mean['phase'][i] * hbt.r2d, t_mean['area_main'][i], #yerr=t_std['area_dn'][i],
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
import matplotlib.lines as mlines

#symbol_plus = mlines.Line2D([], [], color='black', marker='+', ms=12, mew=3, linestyle='none')
#symbol_circle = mlines.Line2D([], [], color = 'black', marker = 'o', ms=10, ls='none')

plt.plot([], [], marker = 'o', ms = 10,              linestyle='none', color='black', label='Main Ring, everything')
plt.plot([], [], marker = '+', ms = 10, lw=3, mew=3, linestyle='none', color='black', label='Main Ring, core only')

plt.title('Jupiter ring, by sequence')
plt.xlabel('Angle [deg]')
plt.ylabel(t_mean['profile_radius_units'][0][0])
plt.yscale('log')
plt.ylim((1e-6,1e-3))
plt.legend()
plt.show()

# =============================================================================
# Plot the phase curve, with one point per image XXX DO NOT PLOT
# =============================================================================

# This is duplicated by much better radial profiles above. Skip it.

if False:
    plt.plot(t['phase'] * hbt.r2d, t['area_dn'], marker = 'o', linestyle = 'none', label=t['label'])
    #plt.yscale('log')
    plt.title('Phase Curve, New Horizons J-Ring, one point per image')
    plt.ylim((1e-7, 3e-4))
    plt.show()
 
# =============================================================================
# Now make some customized individual profiles. For a 'best merged set' profile thing.
# =============================================================================

# Keep these as an example of how to use the class

if False:
    a = ring_profile()
    radius_bg_127 = np.array([[125,126], [130, 131]])*1000
    a.load(7, hbt.frange(40,42),key_radius='full').remove_background_radial(radius_bg_127,do_plot=False).flatten()
    a.dn2iof().plot()
    plt.show()
    
    a = ring_profile()
    a.load(7, hbt.frange(40,42),key_radius='full').flatten().plot()

# =============================================================================
# Now construct a 'best possible radial profile'
# =============================================================================

a = ring_profile()

hbt.set_fontsize(20)
lw                  = 3.5          # Line weight 
dy                  = 3           # Offset between adjacent lines 
num_bins_az_central = 15            # How many azimuthal bins do we use for the profile?


a.load(7, hbt.frange(40,42),key_radius='full')
shift               = [1, 0, 1]     # This is the # of bins to shift it *outward* by. Negative is inward.

a.load(7, hbt.frange(24,31),key_radius='full')
shift               = [0, 0, 0, 2, 0, 0, 0, 1]     # This is the # of bins to shift it *outward* by. Negative is inward.

images = np.delete(hbt.frange(0,47),16)  # 0-47, but exclude #16 (bad background level, maybe due to Metis)
images = np.delete(hbt.frange(0,47),(16,18,19))  # 0-47, but exclude #16 (bad background level, maybe due to Metis)
images = np.delete(hbt.frange(0,47),(16,18,19,29))  # 0-47, but exclude a few specific bad profiles (bg, or satellites)

a.load(8, images, key_radius='full')
shift = np.zeros(48)

# From the image, extract just the central few azimuthal bins and make a profile. These are the best and most reliable2
shift = np.array(shift,dtype=int)  # Must be integers, for np.roll() to work in np 1.11. Floats allowed in np 1.13.

profile = []

bin0  = int(a.num_bins_azimuth/2 - num_bins_az_central/2)
bin1  = int(a.num_bins_azimuth/2 + num_bins_az_central/2)

for i in range(a.num_profiles):        
    image_i = a.image_unwrapped_arr[i] 
    mask = hbt.nanlogical_and(a.mask_stray_unwrapped_arr[i][:,:], 
                              a.mask_objects_unwrapped_arr[i][:,:])  # XXX This mask mergering isn't working right.

    # Now actually sum to get the radial profile. These two should give identical results.
    profile_i = np.nanmean( image_i[:,bin0:bin1], axis=1)
    profile_i2 = nh_jring_extract_profile_from_unwrapped(a.image_unwrapped_arr[i],
                                                        a.radius_arr[i][:],
                                                        a.azimuth_arr[i][:],
                                                        0.047,
                                                        'radial',
                                                        mask_unwrapped=mask)
    profile.append(profile_i)

# Compute the the mean of all the profiles.

profile_sum = 0. * profile[0]

for i,p in enumerate(profile):
    profile_i = np.roll(p, shift[i]) # Grab the profile, and roll it.
    profile_sum += profile_i
i += 1
profile_mean = profile_sum / i

# Make a plot of the summed radial profile

plt.plot(a.radius_arr[0]/1000, profile_mean + dy*i, label = 'Mean', lw=lw)
plt.xlim((127.5, 130))
plt.title("{}, shifted by preset, central {} az bins".format(a, num_bins_az_central))
plt.xlabel('Radius [1000 km]')
plt.show()

# Make the plot of all of the individual radial profiles

for i,p in enumerate(profile):
    profile_i = np.roll(p, shift[i]) # Grab the profile, and roll it.
    plt.plot(a.radius_arr[0]/1000, profile_i + dy*i, label="{}, $\Delta$ = {}".format(a.index_image_arr[i], 
             shift[i]), lw=lw)
i+=1

plt.plot(a.radius_arr[0]/1000, profile_mean + dy*i, label = 'Mean', lw=lw)

#plt.legend(loc='upper left')
plt.xlim((127.5,130))
plt.ylabel(a.profile_radius_units)
plt.xlabel('Radius [1000 km]')
plt.title("{}, shifted by preset, central {} az bins".format(a, num_bins_az_central))
plt.ylim((0, np.amax(profile_mean + dy*(i+3))))
plt.show()

# Now that we have a sum, we can do the correlation
# We want to get the correlation between the sum, and each individual curve.
# If this value is zero, then we are in good shape.
# In theory it is best if we correlate not to the sum, but to the sum of everything *except this one*.
# But if we have enough curves, this should be effectively the same.

lw = 1

shift_computed = []

for i,p in enumerate(profile):
    profile_i = np.roll(p, shift[i]) # Grab the profile, and roll it.
    r = signal.correlate(profile_mean, profile_i)
    shift_computed_i = np.argmax(r) - (len(profile_i)-1)
    plt.plot(r, label="{}, $\Delta_c$ = {}".format(i, shift_computed_i), lw=lw)
    shift_computed.append( shift_computed_i )
plt.title("{}, central {} az bins, signal.correlate(mean, i)".format(a, num_bins_az_central))
plt.xlim((480,520)) 
#plt.ylim((np.amin(r), dy * (i+3) + 2*np.amax(r)))
plt.legend(loc = 'upper left')  
plt.show()

shift_computed = np.array(shift_computed, dtype=int)

## Need to validate sign of shift_computed. Negative, or positive?

for i,p in enumerate(profile):
    profile_i = np.roll(p, shift[i]) # Grab the profile, and roll it.
    r = signal.correlate(profile_mean, profile_i)
    plt.plot(np.roll(r, -shift_computed[i]), label="{}, $\Delta_c$ = {}".format(i, shift_computed[i]), lw=lw)
plt.title("{}, central {} az bins, signal.correlate(mean, i)".format(a, num_bins_az_central))
plt.xlim((480,520)) 
#plt.ylim((np.amin(r), dy * (i+3) + 2*np.amax(r)))
plt.legend(loc = 'upper left')  
plt.show()

### NOW REDO THE STUFF AGAIN, USING THE COMPUTED SHIFTS

# Compute the the mean of all the profiles.

profile_sum = 0. * profile[0]

for i,p in enumerate(profile):
    profile_i = np.roll(p, -shift_computed[i]) # Grab the profile, and roll it.
    profile_sum += profile_i
i += 1
profile_mean = profile_sum / i

# Make a plot of the summed radial profile

lw=3.5

plt.plot(a.radius_arr[0]/1000, profile_mean + dy*i, label = 'Mean', lw=lw)
plt.xlim((127.5, 130))
plt.title("{}, central {} az bins, aligned".format(a, num_bins_az_central))
plt.xlabel('Radius [1000 km]')
plt.show()

# Make the plot of all of the individual radial profiles

lw = 1

for i,p in enumerate(profile):
    profile_i = np.roll(p, shift_computed[i]) # Grab the profile, and roll it.
    plt.plot(a.radius_arr[0]/1000, profile_i + dy*i, label="{}, $\Delta$ = {}".format(a.index_image_arr[i], 
             shift_computed[i]), lw=lw)
i+=1

plt.plot(a.radius_arr[0]/1000, profile_mean + dy*i, label = 'Mean', lw=lw)

#plt.legend(loc='upper left')
plt.xlim((127.5,130))
plt.ylabel(a.profile_radius_units)
plt.xlabel('Radius [1000 km]')
plt.title("{}, shifted by centroid, central {} az bins".format(a, num_bins_az_central))
plt.ylim((0, np.amax(profile_mean + dy*(i+3))))
plt.show()


# =============================================================================
# Plot radial profile vs. Extraction Region
# =============================================================================

plt.set_cmap('plasma')
a_ref = 127_900*u.km

images0 = hbt.frange(0, 47)  # Entire range 
keys = ['center', 'full', 'core']        # Interate over the 'key_radius' field. This affects radial extraction.
do_flatten = [True, False]

for do_flatten_i in do_flatten:
    for key in keys:
        a0 = ring_profile()
        a0.load(8, images0, key_radius=key, verbose=False)
    #    a0.plot_radial(dy=1)
        if do_flatten_i:
            a0.flatten()
        plt.plot(a0.radius_arr[0], a0.profile_radius_arr[0], label = key)
    plt.title(f'Radial Profile vs Extraction Regidon, {a0.__str__()}')
    plt.xlabel('Radius [km]')
    plt.ylabel('DN')    
    plt.axvline(a0.A_METIS.to('km').value, alpha=0.2)
    plt.axvline(a0.A_ADRASTEA.to('km').value, alpha=0.2)
    plt.xlim((122_000, 132_000))
    plt.legend()    
    plt.show()    


#%%%
    
# =============================================================================
# Now make some azimuthal profile plots, with all the data
# =============================================================================

plt.set_cmap('plasma')
a_ref = 129_000*u.km
et_ref = sp.utc2et('2007 24 Feb 12:00:00')  # Epoch to unwind into. For best results, use a time close to the obs

# Look at sequence 8.
# Load it in three segments: Entire range, beginning, and end. Make plots of all of these, separately.
# Then make a nice plot with all of the data, so we can see the correlations

images0 = hbt.frange(0, 47)  # Entire range 
images1 = hbt.frange(0, 20)  # First half
images2 = hbt.frange(35,47)  # Second half

images = [images0, images1, images2]

for images_i in images:
    a0 = ring_profile()
    a0.load(8,images_i,key_radius='full', verbose=False)  #.plot_azimuthal(smooth=3)
    a0.flatten(a=a_ref, et=et_ref)
    plt.imshow(a0.profile_azimuth_arr_2d, aspect=20, origin='lower')
    plt.title(a0)
    plt.show()
    a0.plot_azimuthal(xlim={0,360}, title=a0,smooth=3)

# Make a plot superimposing the data (three curves, group 8)

images = [images0, images1, images2]
smoothing = 5
group = 8

for i,images_i in enumerate(images):
    a0 = ring_profile()
    a0.load(group,images_i,key_radius='full', verbose=False)  #.plot_azimuthal(smooth=3)
    a0.flatten(a=a_ref, et=et_ref)
    if (smoothing):
        a0.smooth_azimuthal(smoothing)
    plt.plot(a0.azimuth_arr[0], a0.profile_azimuth_arr[0], label=a0.__str__(), alpha=0.7)
    plt.title(f'Az Profile, unwind, {group}/, {a_ref}, smoothing {smoothing}')
plt.legend(fontsize=8)
plt.xlabel('Az [rad]')
plt.ylabel('DN')
plt.show()

#%%%    
# Now make a plot showing just a zoom of the range where the data actually overlap (two curves, group 8)
    
images1 = hbt.frange(0, 20)  # First half
images2 = hbt.frange(35,47)  # Second half

images = [images1, images2]

a_ref = 127_900*u.km # Increasing by 100 km shifts orange right by 0.0001 radians

smoothing = 5
group = 8
y = np.zeros((len(images), 3600))   # Save the data here so we can calc the Correlation Coeff afterwards
xlim = (3.5,5.7)                        # If desired, constrain the xlim to just a small region

for i,images_i in enumerate(images):
    a0 = ring_profile()
    a0.load(group,images_i,key_radius='full', verbose=False)  #.plot_azimuthal(smooth=3)
    a0.flatten(a=a_ref, et=et_ref)
    if (smoothing):
        a0.smooth_azimuthal(smoothing)
    plt.plot(a0.azimuth_arr[0], a0.profile_azimuth_arr[0], label=a0.__str__(), alpha=0.5)
        # Plot the DN values themselves
    plt.plot(a0.azimuth_arr[0], a0.num_profile_azimuth_arr[0], label=f'Number of pts for {a0.__str__()}')

    plt.title(f'Az Profile, unwind, {group}/, {a_ref}, smoothing {smoothing}')
plt.legend(fontsize=8)
plt.xlim(xlim)
plt.xlabel('Az [rad]')
plt.ylabel('DN')
plt.show()
    
#%%%

# =============================================================================
# Now see if we can find correlations between adjacent sequences (ie, not in the same sequence).
# Compare 8/0-48, and 8/54-107. These are separated by ~33 hours.
# =============================================================================
images = [
          hbt.frange(0, 48),   # One entire orbit
          hbt.frange(54,107)]  # A second entire orbit, 33 hours later

et_ref = sp.utc2et('2007 24 Feb 12:00:00')  # Epoch to unwind into. For best results, use a time close to the obs

group = 8
smoothing = 2
a_ref = a0.A_METIS
a_ref = a0.A_ADRASTEA

for i,images_i in enumerate(images):
    a0 = ring_profile()
    a0.load(group,images_i,key_radius='full', verbose=False)  #.plot_azimuthal(smooth=3)
    a0.flatten(a=a_ref, et=et_ref)
    if (smoothing):
        a0.smooth_azimuthal(smoothing)
    plt.plot(a0.azimuth_arr[0], a0.profile_azimuth_arr[0], label=a0.__str__(), alpha=0.7)
    plt.title(f'Az Profile, unwind, {group}/, {a_ref}, smoothing {smoothing}')
    x = a0.azimuth_arr[0]
    y[i] = a0.profile_azimuth_arr[0]
plt.legend(fontsize=8)
#plt.xlim(xlim)
plt.xlabel('Az [rad]')
plt.ylabel('DN')
plt.show()

#%%%

# =============================================================================
# Plot sets of nearby subsequent images, unwrapped in az. 
# This is a check to see that the same things are in next-door images!
# =============================================================================

images = [hbt.frange(0,10),
          hbt.frange(11,21)]
et_ref = sp.utc2et('2007 24 Feb 12:00:00')  # Epoch to unwind into. For best results, use a time close to the obs

hbt.figsize((18,6))
hbt.fontsize(15)
a_ref = 127_900*u.km
xlim  = (3.9, 5.2)
smoothing = None
group = 8

for i,images_i in enumerate(images):
    a0 = ring_profile()
    a0.load(group,images_i,key_radius='full', verbose=False)  #.plot_azimuthal(smooth=3)
    a0.flatten(a=a_ref, et=et_ref)
    if (smoothing):
        a0.smooth_azimuthal(smoothing)
    
    # Plot the # of data points that were summed for each DN value 
    plt.plot(a0.azimuth_arr[0], a0.profile_azimuth_arr[0], label=a0.__str__(), alpha=0.7)
    
    # Plot the DN values themselves
    plt.plot(a0.azimuth_arr[0], a0.num_profile_azimuth_arr[0], label=f'Number of pts for {a0.__str__()}')

plt.title(f'Az Profile, unwind, {group}/, {a_ref}, smoothing {smoothing}')
    
plt.legend(fontsize=8)
plt.xlim(xlim)
plt.xlabel('Az [rad]')
plt.ylabel('DN; # of Data Points')
plt.legend(fontsize=10)
plt.show()

#%%%

# =============================================================================
# Now make plot where we do show actual DN of Metis and Adrastea.
# This is as a dummy check to make sure we are summing on them properly, esp. over different sequences.
# =============================================================================

images = [hbt.frange(0,48),
          hbt.frange(54,107)]

hbt.figsize((18,6))
hbt.fontsize(15)
xlim  = (3.9, 5.2)
smoothing = None
group = 8

a_ref = ring_profile().A_METIS
body_ref = 'Metis'

a_ref = a0.A_ADRASTEA
body_ref = 'Adrastea'

et_ref = sp.utc2et('2007 24 Feb 12:00:00')  # Epoch to unwind into. For best results, use a time close to the obs

for i,images_i in enumerate(images):
    a0 = ring_profile()
    a0.load(group,images_i,key_radius='full', verbose=False)  #.plot_azimuthal(smooth=3)

    # Make a strip lot of all the data
    
    im_extract = a0.make_strip_mosaic(a=3,do_plot=True, do_unwind=True)

    # Now that we have read in all the profiles (and their pre-unwrapped images),
    # loop over them, and extract the profile from the images. This will not apply stray light corrections,
    # and is very crude. But it should show the big sats.
    
    profiles = a0.profile_azimuth_arr

    for i in range(len(profiles)):
        
        profile_i = nh_jring_extract_profile_from_unwrapped(a0.image_unwrapped_arr[i], 
                                                          a0.radius_arr[i],
                                                          a0.azimuth_arr[i],
                                                          1,
                                                          'azimuthal',
                                                          a0.mask_stray_unwrapped_arr[i])
        profiles[i] = profile_i
    a0.profile_azimuth_arr = profiles
    do_plot_raw_phased = False
    
    if do_plot_raw_phased:
        plt.imshow(profiles)
        plt.show()
    
    # Unroll the images, applying keplerian motion
    
    a0.flatten(a=a_ref, et=et_ref)

    # Map Metis, Adrasta back to this time
    # What we want to do is use SPICE to calc Metis' longitude at the time of observation, and then use
    # func to wrap this back to J2000 common time.
    
    (pt, et_pt, vec) = sp.subpnt('Intercept: Ellipsoid', 'Jupiter', a0.et_arr[0], 'IAU_JUPITER', 'LT', body_ref.upper())
    (rad, lon, lat)  = sp.reclat(vec)  # Returns longitude in radians
    lon_unwind = unwind_orbital_longitude(lon, a0.et_arr[0], 'Jupiter', a_ref, et_ref)
    print(f'Unwind lon + pi = {lon_unwind + math.pi}')
    
    # Plot the map
    
    plt.imshow(stretch(a0.profile_azimuth_arr_2d),aspect=.030, extent=(0,2*math.pi,0,len(profiles)))
    
    # Plot moon
    
    plt.axvline(np.mod(lon_unwind + math.pi, 2*math.pi), alpha=0.1, lw=10)
    
    plt.xlabel('Azimuth [rad]')
    plt.ylabel('Image #')
    plt.title(f'{a0.__str__()}, Unwinding wrt {body_ref}')
    plt.show()


#%%%
    
# Just a quick test routine to unwind orbital longitudes
# This function below should give the same result for any ET at all! 
# Basically, it moves positions forward with SPICE, and then backwards with my routine, and checks if they end up
# where they started
# At first it did not. But now, after I added the J2 correction to my unwinding routine, it works properly.
# (Not exactly -- there is still an error in my J2 calc -- but it is a lot closer than it was.)
    
    
ets = [0,10000, 20000, 30000]
#
a_ref = 127_979.8*u.km  # Metis    = larger, inner.  From table in SCW07
#a_ref = 128_980.5*u.km  # Adrastea = smaller, outer. From table in SCW07. 

p_metis = 0.294780*86400  # From wiki. Orbital period

for et in ets:
    (pt, et_pt, vec) = sp.subpnt('Intercept: Ellipsoid', 'Jupiter', et, 'IAU_JUPITER', 'LT', 'METIS')
    (rad, lon, lat)  = sp.reclat(vec)  # Returns longitude in radians
    lon_unwind = unwind_orbital_longitude(lon, et, 'Jupiter', a_ref)
    print(f'Unwind lon + pi = {lon_unwind + math.pi}')


#%%%
    
# =============================================================================
# Do some tests on high-resolution extractions. I have increased resolution
# (300,500) → (600,1000) for a few images, just to compare.
# =============================================================================

hbt.figsize((8,5))
hbt.fontsize(12)
images = [hbt.frange(0,90),
          hbt.frange(91,107)]

plt.set_cmap('Greys_r')
plt.set_cmap('plasma')

a_ref = 127_900*u.km
xlim  = (5, 6.5)
smoothing = None
group = 8

for i,images_i in enumerate(images):
    a0 = ring_profile()
    a0.load(group,images_i,key_radius='core', verbose=False)  #.plot_azimuthal(smooth=3)
    a0_flat = a0.copy()
    a0_flat.flatten(a=a_ref)
    if (smoothing):
        a0_flat.smooth_azimuthal(smoothing)
    
    # Make an initial plot of az profile
    
    plt.plot(a0_flat.azimuth_arr[0], a0_flat.profile_azimuth_arr[0], label=a0_flat.__str__(), alpha=0.7)
    plt.xlim((0,6.28))
    plt.title(f'Az Profile, unwind, {a0_flat}, {a_ref}, smoothing {smoothing}')
    plt.show()

#hbt.figsize(25,25)
for i,images_i in enumerate(images):
    a0 = ring_profile()
    a0.load(group,images_i,key_radius='full', verbose=False)  #.plot_azimuthal(smooth=3)
    plt.imshow(stretch(a0.image_unwrapped_arr[0]), aspect=0.5)
    plt.show()

#    
#    a0_flat = a0.copy()
#    a0_flat.flatten(a=a_ref)
#    if (smoothing):
#        a0_flat.smooth_azimuthal(smoothing)


    
#%%%
# =============================================================================
# Now look at the sequence where MRS extracted the clumps from. See if I can reproduce his result.
# His images were 8/97 and 8/100.
# =============================================================================

group = 8
images = [hbt.frange(97,100)]

plt.set_cmap('Greys_r')
plt.set_cmap('plasma')

a_ref = 127_900*u.km
xlim  = (5, 6.5)
smoothing = None

for i,images_i in enumerate(images):  # Load each image set
    a0 = ring_profile()
    a0.load(group,images_i,key_radius='full', verbose=False)  #.plot_azimuthal(smooth=3)
    a0_flat = a0.copy()
    a0_flat.flatten(a=a_ref)
    if (smoothing):
        a0_flat.smooth_azimuthal(smoothing)
    
    # Make an initial plot of az profile
      
    plt.plot(a0_flat.azimuth_arr[0], a0_flat.profile_azimuth_arr[0], label=a0_flat.__str__(), alpha=0.7)
    plt.title(f'Az Profile, unwind, {group}/, {a_ref}, smoothing {smoothing}')
    plt.show()
    
    # Extract all the strips into an image, and show it.
    
    ((im_mosaic, mask_mosaic), (profile_im, profile_masked), strip, az_arr)\
    = a0.make_strip_mosaic(do_plot=False, do_plot_masked=True, do_unwind=True, a_orbit=127_900*u.km)
    
    # Make a series of TV plots
    
    hbt.figsize((10,3))
    for j in range(a0.num_profiles):
        plt.subplot(5,1,j+1)
        plt.imshow(stretch(a0.image_unwrapped_arr[j][380:450,:]), aspect=0.4)
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
    
    plt.tight_layout()    
    plt.show()
    
    # Make a series of plots, summing the data manually, from the 'raw' image
    
    for j in range(a0.num_profiles):
        plt.plot(j*100 + np.sum( a0.image_unwrapped_arr[j][400:430,:], axis=0))        
    plt.show()

    # Re-extract the data in a formal / complete way
    
    # Set the limits for extraction
    
    profile_azimuth_masked  = {'inner' : np.array([0]), 'core'   : np.array([0]), 'outer' : np.array([0])}
#                               'net'   : np.array([0]) }
    profile_azimuth_raw     = {'inner' : np.array([0]), 'core'   : np.array([0]), 'outer' : np.array([0])}
#                               'net'   : np.array([0]) }

    profile_radius          = {'full'  : np.array([0]), 'center' : np.array([0]), 'core'  : np.array([0])}
#                                                     'outer-30'  : np.array([0]), 'outer-50' : np.array([0]) }
    
    range_of_radius  = {'inner' : (126500,127500), 'core' : (127500,129500), 'outer' : (130000,131000)} # for az
    range_of_azimuth = {'full'  : 1,               'center' : 0.25,          'core' : 0.1}
#                                                       'outer-30' : (1,-0.3), 'outer-50' : (1,-0.5)}    

    # Grab the masks (already done)

# Make azimuthal profiles
# Confirmed. I can extract these quite well. When I do this formula, the extracted profile is identical to 
# that which I already have. It is a masked profile.

    
    
    # Loop over all the images in the sequence
    
    for j in range(a0.num_profiles):

        # Loop over and extract profiles from inner, outer, core, etc. regions
        
        for key in profile_azimuth_masked:
            
            # Extract the masked profile
            
            profile_azimuth_masked[key] = nh_jring_extract_profile_from_unwrapped(a0.image_unwrapped_arr[j], 
                                                  self.radius_arr[j],   # Array defining bins of radius
                                                  self.azimuth_arr[j],  # Array defining bins of azimuth
                                                  range_of_radius[key],        # range of rad used for az profile
                                                  'azimuth',
                                                  mask_unwrapped = a0.mask_objects_unwrapped_arr[j])   

            # Extract the non-masked, raw profile
            
            profile_azimuth_raw[key] = nh_jring_extract_profile_from_unwrapped(a0.image_unwrapped_arr[j], 
                                                  self.radius_arr[j],   # Array defining bins of radius
                                                  self.azimuth_arr[j],  # Array defining bins of azimuth
                                                  range_of_radius[key],        # range of rad used for az profile
                                                  'azimuth',
                                                  mask_unwrapped = a0.mask_stray_unwrapped_arr[j])      # Remove the BG
            
        profile_azimuth_masked_net = profile_azimuth_masked['core'] - (profile_azimuth_masked['inner'] + 
                              profile_azimuth_masked['outer'])/2

        profile_azimuth_raw_net = profile_azimuth_raw['core'] - (profile_azimuth_raw['inner'] + 
                              profile_azimuth_raw['outer'])/2
        
#        plt.plot(a0.azimuth_arr[0], profile_azimuth_raw_net, label = 'Raw')
#        plt.plot(a0.azimuth_arr[0], profile_azimuth_masked_net, label = 'Masked')
#        plt.legend()
#        plt.ylim((-3,10))
#        plt.title(f'{a0.__str__()}, {a0.index_image_arr[j]}')
#        plt.show()
#        
                                             
# Now, extract the image strips. I don't have a very good az image of the ring -- just a 1D profile.
# Take the raw strip, and subtract mean of inner and outer profiles, converted to a strip
        
        strip_raw = a0.image_unwrapped_arr[j][340:360,:]
        strip_inner = np.add.outer(np.zeros(20), profile_azimuth_raw['inner'])
        strip_outer = np.add.outer(np.zeros(20), profile_azimuth_raw['outer'])
        
        strip_cleaned = strip_raw - (strip_inner + strip_outer) / 2
        
        plt.imshow(stretch(strip_cleaned))
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.show()
        
# Also, extract non-masked images
            
plt.legend(fontsize=8)
plt.xlim(xlim)
plt.xlabel('Az [rad]')
plt.ylabel('DN')
plt.legend(fontsize=10)
plt.show()

#%%%

# =============================================================================
# Make a plot superimposing the unwrapped images, with the az profiles
# This tests the unwrapping for consistency, by four different methods: 
#   flatten() on profiles
#   flatten() on a 2D image
#   make_strip_mosaic()
#   unwind_orbital_longitude() on Metis / Adrastea itself
# Concl: Unwrapping is working well.
# =============================================================================

group  = [8, 8]

images = [hbt.frange(0,48),  # No masking here. It was not done?
          hbt.frange(54,107)]


do_short = True
if do_short:
    group =  [group[0]]
    images = [images[0]]

profiles_good = []
profiles_bad  = []
im_mosaics    = []
labels        = []

et_ref = sp.utc2et('24 Feb 2007 00:00:00') # Reference ET to unwrap to.

for i in range(len(images)):
    group_i = group[i]
    images_i = images[i]
    
    ring = ring_profile()
    ring.load(group_i,images_i,key_radius='core')
    
    # Define the longitude to unwrap to
    # This can only be one body. If we unwrap to Metis, then Adrastea will be (slightly) mis-aligned, etc.
    
    a_ref = ring.A_ADRASTEA
    
    dwidth_chop = 200
    ((im_mosaic, mask_mosaic),(profile_im, profile_masked), strip, az_arr) = \
      ring.make_strip_mosaic(do_plot=False, do_plot_profile=False, dwidth_chop=dwidth_chop, 
                             a_orbit = a_ref, et_ref = et_ref)

    im_mosaics.append(im_mosaic)
    labels.append(ring.__str__())
    
    profiles_bad.append(profile_masked)  # Save the profile made by simply summing the image

    ring.flatten(a=a_ref, et=et_ref)
    profiles_good.append(ring.profile_azimuth_arr[0])  # Save the profile made by (data-(inner+outer)/2)
    
# Now make a plot superimposing the radial profiles, with the stacked mosaics

width     = 10
linewidth = 4

hbt.figsize(20,10)
hbt.fontsize(18)

for i in range(len(images)):
    plt.plot(hbt.smooth(profiles_good[i], width, boundary='wrap'), alpha=0.5, 
             label=labels[i] + ' (Extracted, hi-qual)',     linewidth=linewidth)
    plt.plot(hbt.smooth(profiles_bad[i],  width, boundary='wrap'), alpha=0.5, 
             label=labels[i] + ' (Summed Images, lo-qual)', linewidth=linewidth/2, ls='--')
    
plt.ylabel('DN/pixel')
plt.xlabel('Azimuth [mrad]')
plt.title(f'J-Ring Azimuthal Brightness Map, with Unwrapped Data, {ring}')

# Unwrap Metis and Adrastea

#bodies = ['Metis', 'Adrastea']
bodies = ['Adrastea']

for body in bodies:
    plt.plot([], [])  # Advance to next color, since axvline() doesn't do so.

    # Look up longitude of body using SPICE

    (pt, _, srfvec)  = sp.subpnt('intercept: ellipsoid', 'Jupiter', ring.et_arr[0], 'IAU_JUPITER', 'LT', body)
    (radius,lon,lat) = sp.reclat(pt)
    dist = sp.vnorm(srfvec - pt)
    
    # Unwind this longitude to the one in the past
    
    lon_unwind = unwind_orbital_longitude(lon, ring.et_arr[0], 'Jupiter', a_orbit=a_ref, et_ref = et_ref)
    
    # And make a plot
    
    plt.axvline(x=lon_unwind*1000, linestyle = '--', label=body)

plt.legend(loc = 'lower right')    

# Now, overlay the image mosaic stack on this plot

plt.imshow(stretch(mask_mosaic*im_mosaic), origin='lower', extent=(0,6300,1,5), aspect=600)
plt.show()    

# For comparison, plot the map generated not from images, but from extracted profiles

do_plot_map_extracted = True
if do_plot_map_extracted:    
    plt.imshow(ring.profile_azimuth_arr_2d, aspect=45, origin='lower')
    plt.title(f'J-Ring Azimuthal Brightness Map, from Individual Az Profiles, {ring}')
    plt.ylabel(f'Image Number, {ring}')
    plt.show()

#%%%

# =============================================================================
# Make a plot superimposing results from highest quality regions, 7 and 8
# =============================================================================

group  = [7, 8, 8]

images = [hbt.frange(0,42),
          hbt.frange(0,48),  # No masking here. It was not done?
          hbt.frange(54,107)]

profiles_good = []
profiles_bad  = []
im_mosaics    = []
labels        = []

et_ref = sp.utc2et('24 Feb 2007 00:00:00') # Reference ET to unwrap to.

for i in range(len(images)):
    group_i = group[i]
    images_i = images[i]
    
    ring = ring_profile()
    ring.load(group_i,images_i,key_radius='core')
    
    # Define the longitude to unwrap to
    # This can only be one body. If we unwrap to Metis, then Adrastea will be (slightly) mis-aligned, etc.
    
    a_ref = ring.A_ADRASTEA
    
    dwidth_chop = 200
    ((im_mosaic, mask_mosaic),(profile_im, profile_masked), strip, az_arr) = \
      ring.make_strip_mosaic(do_plot=False, do_plot_profile=False, dwidth_chop=dwidth_chop, 
                             a_orbit = a_ref, et_ref = et_ref)

    im_mosaics.append(im_mosaic)
    labels.append(ring.__str__())
    
    profiles_bad.append(profile_masked)  # Save the profile made by simply summing the image

    ring.flatten(a=a_ref, et=et_ref)
    profiles_good.append(ring.profile_azimuth_arr[0])  # Save the profile made by (data-(inner+outer)/2)
    
# Now make a plot superimposing the radial profiles, with the stacked mosaics

width     = 10
linewidth = 4

hbt.figsize(20,10)
hbt.fontsize(18)

for i in range(len(images)):
    plt.plot(hbt.smooth(profiles_good[i], width, boundary='wrap'), alpha=0.5, 
             label=labels[i] + ' (Extracted, hi-qual)',     linewidth=linewidth)
    plt.plot(hbt.smooth(profiles_bad[i],  width, boundary='wrap'), alpha=0.5, 
             label=labels[i] + ' (Summed Images, lo-qual)', linewidth=linewidth/2, ls='--')
    
plt.ylabel('DN/pixel')
plt.xlabel('Azimuth [mrad]')
plt.title(f'J-Ring Azimuthal Brightness Map, with Unwrapped Data')

# Unwrap Metis and Adrastea

#bodies = ['Metis', 'Adrastea']
bodies = ['Metis']

for body in bodies:
    plt.plot([], [])  # Advance to next color, since axvline() doesn't do so.

    # Look up longitude of body using SPICE

    (pt, _, srfvec)  = sp.subpnt('intercept: ellipsoid', 'Jupiter', ring.et_arr[0], 'IAU_JUPITER', 'LT', body)
    (radius,lon,lat) = sp.reclat(pt)
    dist = sp.vnorm(srfvec - pt)
    
    # Unwind this longitude to the one in the past
    
    lon_unwind = unwind_orbital_longitude(lon, ring.et_arr[0], 'Jupiter', a_orbit=a_ref, et_ref = et_ref)
    
    # And make a plot
    
    plt.axvline(x=lon_unwind*1000, linestyle = '--', label=body)

plt.legend(loc = 'lower right')    

# Now, overlay the image mosaic stack on this plot

#plt.imshow(stretch(mask_mosaic*im_mosaic), origin='lower', extent=(0,6300,1,5), aspect=600)
#plt.show()    
#
## For comparison, plot the map generated not from images, but from extracted profiles
#
#do_plot_map_extracted = True
#if do_plot_map_extracted:    
#    plt.imshow(ring.profile_azimuth_arr_2d, aspect=45, origin='lower')
#    plt.title(f'J-Ring Azimuthal Brightness Map, from Individual Az Profiles, {ring}')
#    plt.ylabel(f'Image Number, {ring}')
#    plt.show()

plt.show()



#%%%

# =============================================================================
# Make a high-quality large stacked image output of the best two sequences, 
# which we can take into Photoshop or DS9 for analysis.
# =============================================================================

group  = [\
#        5, 
          7, 8, 8]

images = [\
#        np.array([1,2,3,4,5,6,7,8,9,  11,12,13,14]),
          hbt.frange(0,42),
          hbt.frange(0,48),  
          hbt.frange(54,107)]

profiles_good = []
profiles_bad  = []
im_mosaics    = []
labels        = []

et_ref = sp.utc2et('24 Feb 2007 00:00:00') # Reference ET to unwrap to.

for i in range(len(images)):
    
    group_i = group[i]
    images_i = images[i]
    num_images = len(images_i)

    figsize_base = 10
    hbt.figsize((figsize_base,figsize_base))
    
    ring = ring_profile()
    ring.load(group_i,images_i,key_radius='core')
    
    num_az = hbt.sizex_im(im_mosaic)
    num_profiles = len(images_i)
    
    # Define the longitude to unwrap to
    # This can only be one body. If we unwrap to Metis, then Adrastea will be (slightly) mis-aligned, etc.
    
    a_ref = ring.A_ADRASTEA
    
    # Loop over dwidth. This is the width to *remove* from each end of each azimuthal scan. 
    # Larger dwidth → fewer pixels used in the output array, and closer to the ansa.
    
    dwidth_chop = [500, 200]  # Can set to 200 or 500
    
    for dw_i in dwidth_chop:
        ((im_mosaic, mask_mosaic),(profile_im, profile_masked), (im_strip, masked_strip), az_arr) = \
          ring.make_strip_mosaic(do_plot=False, do_plot_profile=False, dwidth_chop=dw_i, 
                                 a_orbit = a_ref, et_ref = et_ref)
    
        im_mosaics.append(im_mosaic)
        labels.append(ring.__str__())
    
        aspect = num_az / num_profiles
    
        plt.imshow(stretch(            im_mosaic), origin='lower', extent=(0,6300,0,num_profiles-1), aspect=aspect, alpha=0.5)
        plt.imshow(stretch(mask_mosaic*im_mosaic), origin='lower', extent=(0,6300,0,num_profiles-1), aspect=aspect)
        plt.xlabel('Azimuth [mrad]')
        plt.ylabel(f'Image Number, {ring}')
        
        # Generate the base name
        
        file_base = (f'mosaic_{ring}_dw{dw_i}.png').replace('/', '_')
        file_out = (os.path.join(ring.dir_out, file_base))
    
        plt.show()    
        
        # Save the mosaics to disk as PNGs
        
        plt.imsave(file_out, im_mosaic, cmap = plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out}')
        
        file_out_tmp = file_out.replace('.png', '_stretch.png')
        plt.imsave(file_out_tmp, stretch(im_mosaic), cmap=plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out_tmp}')
        
        file_out_tmp = file_out.replace('.png', '_mask_stretch.png')
        plt.imsave(file_out_tmp, stretch(im_mosaic * mask_mosaic), cmap=plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out_tmp}')
        
        # Save the strips to disk as PNGs
        
        file_out = file_out.replace('mosaic', 'strip')
        
        plt.imsave(file_out, im_strip, cmap = plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out}')
    
        file_out_tmp = file_out.replace('.png', '_mask.png')
        plt.imsave(file_out_tmp, masked_strip, cmap = plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out}')
    
        file_out_tmp = file_out.replace('.png', '_stretch.png')
        plt.imsave(file_out_tmp, stretch(im_strip), cmap=plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out_tmp}')
        
        file_out_tmp = file_out.replace('.png', '_mask_stretch.png')
        plt.imsave(file_out_tmp, stretch(masked_strip), cmap=plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out_tmp}')
        
        # Make a 2D map of the ring
        
        file_out_tmp = file_out.replace('strip', 'strip-recut')

        strip_recut = recut_strip(im_strip, dtheta = 1)
        plt.imsave(file_out, stretch(strip_recut), cmap=plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out_tmp}')
        
        # Make a merged plot, of strips and mosaic together
        
        merged = merge_mosaic_strip(im_mosaic * mask_mosaic, masked_strip)
        
        file_out = file_out.replace('strip-recut', 'merged')
        plt.imsave(file_out, stretch(merged), cmap=plt.cm.plasma, origin = 'lower')
        print(f'Wrote: {file_out}')
    
        
    
#%%%
    
# Now make a plot superimposing the radial profiles, with the stacked mosaics


hbt.figsize(20,10)
hbt.fontsize(18)


plt.xlabel('Azimuth [mrad]')
plt.title(f'J-Ring Azimuthal Brightness Map, with Unwrapped Data, {ring}')


plt.legend(loc = 'lower right')    

# Plot the image mosaic stack


# For comparison, plot the map generated not from images, but from extracted profiles

do_plot_map_extracted = True
if do_plot_map_extracted:    
    plt.imshow(ring.profile_azimuth_arr_2d, aspect=45, origin='lower')
    plt.title(f'J-Ring Azimuthal Brightness Map, from Individual Az Profiles, {ring}')
    plt.ylabel(f'Image Number, {ring}')
    plt.show()



















a2 = ring_profile()
a2.load(8,hbt.frange(0,48),key_radius='core').plot(plot_legend=False)
#a2.plot()
a2.flatten()
a2.plot_azimuthal(smooth=1)

# Test different smoothing algorithms

a3 = ring_profile()
a3.load(8,hbt.frange(0,28),key_radius='full').plot(plot_legend=False)
a3.flatten()
a3.plot_azimuthal(smooth=1)
a3.plot_azimuthal(smooth=10)
a3.plot_azimuthal(smooth=100)


a3 = ring_profile()
a3.load(8,[13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29],key_radius='core')
a3.flatten()
a3.plot()
a3.plot_azimuthal(smooth=10,xlim=(0,360))

#.flatten().plot_radial()

a1 = ring_profile()
a1.load(8,hbt.frange(0,48),key_radius='core').flatten(a=a1.A_ADRASTEA)

a2 = ring_profile()
a2.load(8,hbt.frange(54,107),key_radius='core').flatten(a=a1.A_ADRASTEA).plot_azimuthal(xlim=(0,360),smooth=3)

plt.plot(a1.azimuth_arr[0], a1.profile_azimuth_arr[0])
plt.plot(a2.azimuth_arr[0], a2.profile_azimuth_arr[0])
plt.plot(a2.azimuth_arr[0], np.roll(a2.profile_azimuth_arr[0],2150))
plt.show()

plt.xlim((0,360)*hbt.r2d)


flatten(a=127_900*u.km).plot_azimuthal(xlim=(0,360))

a2_flat1 = a2.copy()
a2_flat1.flatten(a=127_900*u.km).plot_azimuthal(xlim=(0,360))

a2 = ring_profile()
a2.load(8,hbt.frange(30,48),key_radius='core').plot(plot_legend=False)
a2_flat2 = a2.copy()
a2_flat2.flatten(a=127_900*u.km).plot_azimuthal()
a2_flat2.plot_azimuthal(smooth=5)


a3 = ring_profile()
a3.load(8,hbt.frange(49,53),key_radius='full', verbose=True)

a4 = ring_profile()
a4.load(8,hbt.frange(54,107),key_radius='full', verbose=False)
a4.flatten(a=110_000*u.km).plot_azimuthal()


a2.plot_azimuthal()
a2.flatten(a=127_900*u.km).plot_azimuthal()

a3 = ring_profile()
a3.load(8,hbt.frange(36,49),key_radius='full').flatten(a=127_900*u.km).plot_azimuthal()

a.plot(plot_legend=False)
a_flat = a.copy().flatten()
a_flat.plot()

a3 = ring_profile()
a3.load(8,hbt.frange(13,23)).flatten().plot_azimuthal()

# Do a test with Metis in the frame

a = ring_profile()
self = a
plot_azimuthal=True
plot_legend=True
a.load(8,hbt.frange(8,28), key_radius='full').plot(plot_legend=False, plot_azimuthal=True, plot_radial=True, title=a, 
      a_unwind=129700*u.km)
a.flatten()
plt.plot(a.azimuth_arr[0], a.profile_azimuth_arr[0])
plt.show()
plt.imshow(a.profile_azimuth_arr_2d, aspect=30)
plt.show()




ring2 = ring_profile()
ring2.load(7,hbt.frange(24,31), key_radius='full').flatten().plot()
ring2.plot_azimuthal(xlim=(0,340))

flatten().plot_azimuthal()

plot(plot_legend=False, plot_azimuthal=False, plot_radial=False, 
          title=a, a_unwind=129700*u.km)
ring2.flatten()
plt.plot(ring2.azimuth_arr[0], ring2.profile_azimuth_arr[0])
plt.show()
plt.imshow(ring2.profile_azimuth_arr_2d, aspect=30)
plt.show()


ring3 = ring_profile()
ring3.load(8,hbt.frange(54,107), key_radius='full').plot(plot_legend=False, plot_azimuthal=False, plot_radial=False, 
          title=a, a_unwind=129700*u.km)
ring3.flatten()
plt.plot(ring3.azimuth_arr[0], ring3.profile_azimuth_arr[0])
plt.show()
plt.imshow(ring3.profile_azimuth_arr_2d, aspect=30)
plt.show()

ring3.plot(plot_azimuthal=True)


ring2.load(7,hbt.frange(24,31), key_radius='full')
ring2.plot(plot_legend=False, plot_azimuthal=True, title=a, a_unwind=129700*u.km)


a.load(8,hbt.frange(24,48), key_radius='full')
a.plot(plot_legend=False, plot_azimuthal=True, title=a, a_unwind=129700*u.km)

a_8_1 = ring_profile()
a_8_1.load(8,hbt.frange(0,48), key_radius='full').flatten().plot_azimuthal()


a_8_2 = ring_profile()
a_8_2.load(8,hbt.frange(54,107), key_radius='full').flatten().plot_azimuthal()

a.load(7,hbt.frange(24,31), key_radius='full')
a.plot(plot_legend=False, plot_azimuthal=True, title=a, a_unwind=129700*u.km)

a.load(5,hbt.frange(1,7), key_radius='full')
a.plot(plot_legend=False, plot_azimuthal=True, title=a, a_unwind=129700*u.km)


# Plot some images and masks to explore what is happening here.

hbt.figsize((10,5))
plt.subplot(2,2,1)
plt.imshow(stretch(a.image_unwrapped_arr[8]), aspect=0.3)
plt.subplot(2,2,2)
plt.imshow(a.mask_objects_unwrapped_arr[8], aspect=0.3)
plt.subplot(2,2,3)
plt.imshow(stretch(a.image_unwrapped_arr[9]), aspect=0.3)
plt.subplot(2,2,4)
plt.imshow(a.mask_objects_unwrapped_arr[9], aspect=0.3)
plt.show()

((im,mask),(profile, profile_masked), strip, az_arr)=\
    ring.make_strip_mosaic(do_plot_masked=True, do_unwind=True, do_plot=True, do_plot_profiles=True, dwidth_chop=0)

profile8_1 = np.nanmedian(im, axis=0)
profile_mask8_1 = np.nanmedian(mask * im, axis=0)
#plt.plot(profile, alpha=0.3)
#plt.plot(profile_mask, alpha=0.3)
#plt.xlim((0,6300))
#plt.show()

((im,mask),(profile,profile_mask), strip, az_arr)=ring8_2.make_strip_mosaic(do_unwind=True, do_plot_image=True, 
                                                             do_plot_profile=True)

profile = np.nanmedian(im, axis=0)
profile_mask = np.nanmedian(mask * im, axis=0)
plt.plot(profile, alpha=0.3)
plt.plot(profile8_1, alpha=0.3)
plt.plot(profile_mask, alpha=0.3)
plt.plot(profile_mask8_1, alpha=0.3)
plt.xlim((0,6300))
plt.show()

profile = np.nanmedian(im, axis=0)
profile_mask = np.nanmedian(mask * im, axis=0)
plt.plot(profile, alpha=0.3)
plt.plot(profile_mask, alpha=0.3)
plt.show()


hbt.figsize(15,5)
plt.subplot(1,2,1)
indices = [0,1,5,8,10]
for index in indices:
    delta_az = ring3.azimuth_arr[index][1] - ring3.azimuth_arr[index][0]
    plt.plot(ring3.azimuth_arr[index], label=f'{index}, delta_az = {delta_az:8.5f} rad')
plt.ylabel('Radians [deg]')
plt.xlabel('Bin #')
plt.title(ring3)
plt.legend()

plt.subplot(1,2,2)
for index in indices:
    delta_radius = ring3.radius_arr[index][1] - ring3.radius_arr[index][0]
    plt.plot(ring3.radius_arr[index], label=f'{index}, delta_radius = {delta_radius:6.3f} km')
plt.ylabel('Radius [km]')
plt.xlabel('Bin #')
plt.title(ring3)
plt.legend()
plt.show()

(im,mask)=ring3.make_strip_mosaic(do_plot=False,do_mask=True, do_unwind=False)
(im_u,mask_u)=ring3.make_strip_mosaic(do_plot=False,do_mask=True, do_unwind=True)

hbt.figsize((20,20))
plt.imshow(stretch(im),origin='bottom')
plt.show()

hbt.figsize((20,20))
plt.imshow(stretch(im_u),origin='bottom')
plt.show()
pl

plt.plot(ring3.azimuth_arr[10])
plt.plot(ring3.azimuth_arr[20])
plt.plot(ring3.azimuth_arr[30])
plt.plot(ring3.azimuth_arr[50])
plt.show()


### On-off to test unwrapping and making a stacked mosaic

self = ring3
do_plot = False
do_mask = True
do_unwind = True
dy_extract = 50
y0_extract = 415
gap_y_pix = 1

# Calc the size of the output array

rad_per_pix       = self.azimuth_arr[0][1] - self.azimuth_arr[0][0]
width_mosaic_pix  = int(2*math.pi / rad_per_pix)    
height_mosaic_pix = self.num_profiles * (dy_extract)
az_arr            = np.array(range(width_mosaic_pix)) * rad_per_pix


# Create the output arrays

im_mosaic      = np.zeros((dy_extract*self.num_profiles, width_mosaic_pix))  # Images are y, x
im_mosaic[:,:] = np.nan
mask_o_mosaic    = im_mosaic.copy()
mask_s_mosaic    = im_mosaic.copy()

lon0_unwind_rad = np.zeros(self.num_profiles)
lon0_unwind_pix = np.zeros(self.num_profiles)
 
do_unwind = True

hbt.figsize((20,20)) 
plt.imshow(stretch(im_mosaic))
plt.title(f'Unwind: {do_unwind}')
plt.show()

im_mosaic_masked = im_mosaic.copy()
im_mosaic_masked[mask_s_mosaic == False] = np.nan
plt.imshow(im_mosaic_masked)

profile=np.nanmedian(im_mosaic, axis=0)
profile_masked = np.nanmedian(im_mosaic_masked, axis=0)

hbt.figsize((20,5))
plt.plot(profile, label = 'Masked, no stray or objects', alpha=0.5)
plt.plot(profile_masked, label='Raw data', alpha=0.5)
plt.title(self)
plt.legend()
plt.show()
        

## One-off to test flattening, which used to work but broke
group = 8
images_i = hbt.frange(97,100)
a0 = ring_profile()
a0.load(group,images_i,key_radius='full', verbose=False)
((im_mosaic, mask_mosaic), (profile_im, profile_masked), az_arr) = a0.make_strip_mosaic(do_plot=True, xlim = (0,6300))
a0.flatten()
a0.plot_azimuthal_2d()
a0.plot_azimuthal(xlim=(0,360))
