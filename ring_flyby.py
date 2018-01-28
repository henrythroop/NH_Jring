#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 10:53:57 2018

@author: throop
"""

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
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
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
import imreg_dft as ird                    # Image translation

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt

class ring_flyby:

    """
    Simulate passage of New Horizons through a ring system surrounding MU69.
    
    For Hazard team. Jan 2019. New Horizons KEM.
    
    Mark Showalter has his own code to do this. But I will use mine.
    
    List of all parameters:
        
    Albedo
    
    """

# =============================================================================
# Set up the grid, and basic initialization.
# =============================================================================
    
    def __init__(self, num_pts_3d, gridspacing, frame, name_target):

        """
        Initialize the grid.
        
        Parameters
        ----
            
            num_pts_3d:
                Tuple (n_x, n_y, n_z) describing the size of the grid. Each entry is the number of points in that
                dimension.
                
                The axes (XYZ) are the SPICE axes, in MRS's Sunflower frame. This is *not* the same as 
                normal ring axes. For instance, to sum vertically, we sum in the __ direction, not Z.
                
            gridspacing:
                Tuple (dx, dy, dz) describing the binsize in each dirction. Each entry is the width, in km,
                of a single bin.
            
            frame:
                The name of the SPICE frame (e.g., '2014_MU69_SUNFLOWER_ROT')
                
            name_target:
                The name of the body (e.g., 'MU69')
                
            q_dust:
                Power law exponent of the dust size distribution: n(r) dr = r^(-q) dr
            
        """
        
        # We want to set up this grid in the MU69 Sunflower frame.
    
        if frame:
            self.frame = frame
        if name_target:
            self.name_target = name_target
        
        self.abcorr = 'LT'
        
        # Define the grids
        
        # NB: "_arr" always refers to a 3D array.
        
        self.numberdensity_arr = np.array(num_pts_3d) # Number of dust grains in this bin. Normalized somehow, not sure how.
        self.azimuth_arr = np.array(num_pts_3d) # Azmuth angle. Not sure we need this. In body frame.
        self.radius_arr  = np.array(num_pts_3d) # Radius. This will be useful. In body frame.
        
        self.x_arr       = np.array(num_pts_3d)
        self.y_arr       = np.array(num_pts_3d)
        self.z_arr       = np.array(num_pts_3d)
    
        n = num_pts_3d[0] # Number of grid points on each side
        
        # Define 1D arrays for xyz position
        
        x_1d = hbt.frange(-n/2, n/2, n) * gridspacing[0]
        y_1d = hbt.frange(-n/2, n/2, n) * gridspacing[1]
        z_1d = hbt.frange(-n/2, n/2, n) * gridspacing[2]
        
        # Save the 1D arrays for future use
        
        self.x_1d = x_1d
        self.y_1d = y_1d
        self.z_1d = z_1d
        
        # Define 3D arrays for xyz position
        # Note that the LHS xyz order is opposite the RHS order. I don't know exactly why, but it has something
        # to do with python array indexing order. In any event, the order below gives expected results
        # for accessing an array with values    arr[x_index, y_index, z_index]
        
#        (y_arr, x_arr, z_arr) = np.meshgrid(x, y, z)
        
        (self.y_arr, self.x_arr, self.z_arr) = np.meshgrid(x_1d, y_1d, z_1d)

        # Define the radius (ie, distance from the center, 3D)
        # ** Do not confuse this with the ring radius!
        
        radius_arr = np.sqrt(self.x_arr**2 + self.y_arr**2 + self.z_arr**2)
        self.radius_arr = radius_arr
        
        # Define the 'ring radius' -- ie, distance from the center, in the XZ plane only!

        radius_ring_arr = np.sqrt(self.x_arr**2 + self.z_arr**2)
        self.radius_ring_arr = radius_ring_arr

                
# =============================================================================
# Set up the ring itself
# =============================================================================

    def set_ring_parameters(self, file_profile_ring, albedo, q_dust):
        
        """
        Define the ring radial profile, size distribution, albedo, size dist, etc.
        """
        
        # Read the radial profile
        
        t = Table.read(file_profile_ring, format = 'ascii')
        radius_km  = t['RadiusKM']
        IoF        = t['I/F']
        radius_pix = t['RadiusPixels']
        
        # Extrapolate the radial profile to the right bins, and set it.
        
        bin_profile = np.digitize(radius_km, self.x_1d)
        
        # Now, we need to loop over the value of radius_ring, in 1-pixel steps.
        # For each step, fill everything at that radius and outward with value from IoF.
        
        # numberdensity = # of particles per cm3.
        # It is the total number of particles, from min(self.r_dust) .. max(self.r_dust).
        # It does not include a size distribution. 
        
        # Q: I could alternatively describe numberdensity as the # of particles in the smallest bin.
        #    This seems mathematically easier.
        #    *** I will choose this. 
        
        self.numberdensity_arr = self.radius_arr.copy()*0 
                
        # Get the positive x values. Use these as the radial bins to loop over
        
        radii = self.x_1d[self.x_1d > 0]  # Loop over the radius_ring in the output grid
        
        for i,radius_i in enumerate(radii): # Radius is in km
            is_good = self.radius_ring_arr > radius_i
            self.numberdensity_arr[is_good] = IoF[ int(np.digitize(radius_i, self.x_1d)) ]
                        
        # Set the size distribution
        # We set one size distribution for the entire ring. It is uniform and does not change spatially or temporally.
        
        is_good = np.logical_and(self.radius_ring_arr > 5000, self.radius_ring_arr < 10000)

        # Save the albedo
        
        self.albedo = albedo
                      
        self.n_dust = []
        self.r_dust = []
        
        # Define the size distribution
        # We'll normalize it later if we need to. For now, just make it a power, law, spanning the bins.

        # We assume this size distribution to apply to the entire ring equally.
        
        r_dust = np.array([0.046, 0.1, 0.215, 0.46, 1, 2.15, 4.6])*u.mm.to('cm') # From MRS 28-Jan-2017
        
#        r_dust_min = 0.1*u.micron.to('cm')
#        r_dust_max = 5*u.mm.to('cm')         
#        nbins_dust = 10
#        r_dust = hbt.frange(0.1*u.micron.to('cm'), 5*u.mm.to('cm'), nbins_dust, log=True)
        
        n_dust = r_dust**-q_dust
        
        self.r_dust = r_dust
        self.n_dust = n_dust
        
        
# =============================================================================
# Plot the ring radial profile
# =============================================================================

    def plot_radial_profile(self):
        """
        Plot the ring's radial profile. This is mostly used for testing, just to 
        make sure that it has been loaded properly.
        
        Also might plot vertical profile, or other slices through it.
        """
        
        # What we want to plot is the radial profile, in terms of I/F.
        
        # We now want to take a radial slice. This will let us plot the I/F.
        
        # We use that to then normalize s.t. 'density' is correct, and I/F matches given.
        
        # Then, we output values of 'density' along the path.
        
        num_radius = len(self.x_1d) / 2
        
        index_center = int(num_radius)
        
        radius_ring = self.radius_ring_arr[index_center:,0,index_center]
        
        IoF = self.albedo * np.sum(self.n_dust * self.r_dust**2) * self.numberdensity_arr[index_center:,0,index_center]
        
        # Flatten in the vertical direction.
        
        IoF_2d = np.sum(IoF, axis=2)
        
        # Plot the size distribution, just for fun
        
        plt.plot(self.r_dust, self.n_dust, marker = 'o', linestyle='none')
        plt.title('Size Distribution')
        plt.xlabel('Size [cm]')
        plt.ylabel('Number')
        plt.yscale('log')
        plt.xscale('log')
        plt.show()
        
# =============================================================================
# Fly a trajectory through the ring and sample it
# =============================================================================
        
    def fly_trajectory(self, name_observer, et_start, et_end, dt):
        """
        Now that all parameters are set, sample the ring along a flight path.
        """
    
        # Save the passed-in parameters in case we need them 
        
        self.et_start = et_start
        self.et_end   = et_end
        self.dt       = dt
        
        self.name_observer = name_observer
        
        # Set up the output time array
        
        num = math.ceil( (et_end - et_start) / dt )
        
        et_t = hbt.frange(int(et_start), int(et_end), num)
        
        # Loop over et
        
        radius_t = []
        x_t      = []
        y_t      = []
        z_t      = []
        lon_t    = []
        lat_t    = []
        density_t= []
        bin_x_t  = []
        bin_y_t  = []
        bin_z_t  = []
        
        
        for i,et_i in enumerate(et_t):
#            (st, lt) = sp.spkezr(self.name_target, et_i, 'J2000', self.abcorr, self.name_observer)  
                                                                    # Gives RA/Dec of MU69 from NH
            (st, lt) = sp.spkezr(self.name_observer, et_i, self.frame, self.abcorr, self.name_target)
                                                                    # Get position of s/c in MU69 frame!
            (radius, lon, lat) = sp.reclat(st[0:3])
            
            # Get the lon/lat wrt time. Note that for MU69 flyby, the lat is always positive.
            # This is non-intuituve, but it is because the MU69 Sunflower frame is defined s.t. the ring
            # is in the XZ plane. 
            # In this plane, indeed NH stays 'above' MU69 the whole flyby, with lat always positive.
            
            radius_t.append(radius)
            lon_t.append(lon)
            lat_t.append(lat)

            # Get the XYZ positions wrt time.
        
            x_t.append(st[0])
            y_t.append(st[1])
            z_t.append(st[2])
        
        # Convert into bin values. This is vectorized.
        # ** If values exceed the largest bin, then the index returned will be too large for the density lookup!
        
        bin_x_t = np.digitize(x_t, self.x_1d)
        bin_y_t = np.digitize(y_t, self.y_1d)
        bin_z_t = np.digitize(z_t, self.z_1d)
        
        # If indices are too large, drop them by one. This just handles the edge cases so the code doesn't crash.
        
        bin_x_t[bin_x_t >= len(self.x_1d)] = len(self.x_1d)-1
        bin_y_t[bin_y_t >= len(self.y_1d)] = len(self.y_1d)-1
        bin_z_t[bin_z_t >= len(self.z_1d)] = len(self.z_1d)-1
        
        # Now that we have all the XYZ bins that the s/c travels in, get the density for each one.
        # We should be able to do this vectorized -- but can't figure it out, so using a loop.
        
        numberdensity_t = 0. * et_t
        
        for i in range(len(et_t)):
            numberdensity_t[i] = self.numberdensity_arr[bin_x_t[i], bin_y_t[i], bin_z_t[i]]

        # ** This gives us number density at the smallest binsize. We then need to apply n(r), to get 
        # density at all other bin sizes.
        # Make a cumulative sum of the density. We don't need this -- Doug Mehoke will make his own -- but for testing.
        
        numberdensity_cum_t       = np.cumsum(numberdensity_t)
                            
        # Save all variables so they can be retrieved
        
        self.radius_t      = radius_t
        self.et_t          = et_t        
        self.bin_x_t       = bin_x_t
        self.bin_y_t       = bin_y_t
        self.bin_z_t       = bin_z_t
        self.numberdensity_t = numberdensity_t
        self.numberdensity_cum_t = numberdensity_cum_t
        self.lon_t         = lon_t
        self.lat_t         = lat_t
        self.x_t           = x_t
        self.y_t           = y_t
        self.z_t           = z_t
        
        self.radius_t = radius_t  # This is distance from the center, in 3D coords. Not 'ring radius.'
        
        return()
        
# =============================================================================
# END OF METHOD DEFINITION
# =============================================================================
    
# =============================================================================
# The main function to call to run the simulation
# =============================================================================
    
def do_ring_flyby():
    
    albedo              = [0.05]
    q_dust              = [2]
    inclination_ring    = [0., 0.1]
    
    file_profile_ring = '/Users/throop/Dropbox/Data/ORT1/throop/backplaned/K1LR_HAZ04/stack_n40_z4_profiles.txt'
    
 # This defines I/F as a function of orbital distance.
                                       # n(r) and q are then derived from I/F and albedo.
    
    utc_ca = '2019 1 Jan 05:33:00'
    dt_before = 1*u.hour
    dt_after  = 1*u.hour
    
    frame = '2014_MU69_SUNFLOWER_ROT'
    
    name_target = 'MU69'
    
    name_observer = 'New Horizons'
    
    dt = 10         # Sampling time through the flyby. Assumed to be seconds.

    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem.tm')
    
    # Define the number of grid locations on each side.
    
    n_dx = n_dy = n_dz = 201   # This is 'multiple assignment' in Python. Unlike IDL, it is official and legal.
    
    # Define the size of each grid box
    
    dx_bin = dy_bin = dz_bin = 250  # Distances are km. Not worth tracking the units. Consistent with SPICE.
                                    # XXX Wiki said 25 km, but should be 250.
    
    # Initialize the grid 
    
    ring = ring_flyby((n_dx, n_dy, n_dz), (dx_bin, dy_bin, dz_bin), frame, name_target)

    # Load the trajectory

    et_ca = sp.utc2et(utc_ca)
    et_start = et_ca - dt_before.to('s').value
    et_end   = et_ca + dt_after.to('s').value
    
    # Loop over the input parameters
    
    for albedo_i in albedo:
        for q_dust_i in q_dust:
            for inclination_i in inclination_ring:
    
                albedo_i = albedo[0]
                q_dust_i = q_dust[0]
                inclination_i = inclination_ring[0]
                
                # Set up the ring itself
                
                ring.set_ring_parameters(file_profile_ring, albedo_i, q_dust_i)

                # And fly through it
                
                ring.fly_trajectory(name_observer, et_start, et_end, dt)    
                
                # Write the trajectory to a file, plot it, etc.
                
                plt.subplot(2,2,1)
                plt.plot(ring.radius_t)
                plt.title('Radius')
                
                plt.subplot(2,2,2)
                plt.plot(ring.lat_t)
                plt.title('Lat')
                
                plt.subplot(2,2,3)
                plt.plot(ring.lon_t)
                plt.title('Lon')
                
                plt.show()

                plt.subplot(2,1,1)                
                plt.plot(ring.et_t - np.amin(ring.et_t), ring.numberdensity_cum_t)
                plt.title('Number Density (cumulative)')
                plt.xlabel('ET')
                plt.ylabel('Density Cum')

                plt.subplot(2,1,2)
                plt.plot(ring.et_t - np.amin(ring.et_t), ring.numberdensity_t)
                plt.title('Density')
                plt.xlabel('ET')
                plt.ylabel('Density Cum')
                plt.show()
                
    self = ring
                
                
# =============================================================================
# Run the function
# =============================================================================
                
if (__name__ == '__main__'):
    do_ring_flyby()
    