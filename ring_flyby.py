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
    
    def __init__(self, num_pts_3d, binsize_km_3d, frame, name_target):

        """
        Initialize the grid.
        
        Parameters
        ----
            
            num_pts_3d:
                Tuple (n_x, n_y, n_z) describing the size of the grid. Each entry is the number of points in that
                dimension.
                
                The axes (XYZ) are the SPICE axes, in MRS's Sunflower frame. This is *not* the same as 
                normal ring axes. For instance, to sum vertically, we sum in the __ direction, not Z.
                
            binsize_km:
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
        
        # NB: _arr refers to a 3D array.
        #     _1d  refers to a 1D array.
        
        self.numberdensity_arr = np.array(num_pts_3d) # Number of dust grains in this bin.
        self.azimuth_arr = np.array(num_pts_3d) # Azmuth angle. Not sure we need this. In body frame.
        self.radius_arr  = np.array(num_pts_3d) # Radius. This will be useful. In body frame.
        
        self.x_arr       = np.array(num_pts_3d)
        self.y_arr       = np.array(num_pts_3d)
        self.z_arr       = np.array(num_pts_3d)
    
        n = num_pts_3d[0] # Number of grid points on each side
        
        # Define 1D arrays for xyz position
        
        x_1d = hbt.frange(-n/2, n/2, n) * binsize_km_3d[0]
        y_1d = hbt.frange(-n/2, n/2, n) * binsize_km_3d[1]
        z_1d = hbt.frange(-n/2, n/2, n) * binsize_km_3d[2]
        
        # Save the 1D arrays for future use
        
        self.x_1d = x_1d
        self.y_1d = y_1d
        self.z_1d = z_1d
        
#        x_1d = self.x_1d + 0
#        y_1d = self.y_1d * 100
#        z_1d = self.z_1d * -1
#        
        self.binsize_km_3d = binsize_km_3d
        
        # Define 3D arrays for xyz position. 
        # Note that when using meshgrid, need to pass 'indexing' keyword -- otherwise first two axes are swapped.
        # I don't quite understand it, but this works and documentation confirms it.
        # This gives the expected results for accessing an array with values    arr[x_index, y_index, z_index]
                
        (self.x_arr, self.y_arr, self.z_arr) = np.meshgrid(x_1d, y_1d, z_1d, indexing = 'ij')  # Nominal

        # Define the radius (ie, distance from the center, 3D)
        # ** Do not confuse this with the ring radius!
        
        radius_arr = np.sqrt(self.x_arr**2 + self.y_arr**2 + self.z_arr**2)
        self.radius_arr = radius_arr
        
        # Define the 'ring radius' -- ie, distance from the center, in the XZ plane only!

        radius_ring_arr = np.sqrt(self.x_arr**2 + self.z_arr**2)
        self.radius_ring_arr = radius_ring_arr

        # Define the 'ring vertical' -- ie, distance from the midplane, relative to XZ plane.
        
        vertical_ring_arr = self.y_arr - np.mean(self.y_1d)
        self.vertical_ring_arr = vertical_ring_arr
        
        ### Just some testing
        
        x_1d = hbt.frange(-10,10,21)
        y_1d = hbt.frange(-100,100,21)
        z_1d = hbt.frange(-1000,1000,21)
        
        (x_arr, y_arr, z_arr) = np.meshgrid(x_1d, y_1d, z_1d, indexing = 'xy')  # Testing
               
        return
                
# =============================================================================
# Set up the ring itself
# =============================================================================

    def set_ring_parameters(self, file_profile_ring, albedo, q_dust, inclination=None):
        
        """
        Define the ring radial profile, size distribution, albedo, size dist, etc.
        
        Parameters
        ----
        
        file_profile_ring:
            Name of a text file which defines the ring's radial profile.
            
        albedo:
            Float. Ring particle albedo.
            
        q_dust:
            Float. Dust size distribution exponent. Positive. n(r) dr = r^(-q) dr
            
        inclination:
            Float. Fraction (e.g., 0.1 = 10%). If None, then the ring takes the full vertical height of the box 
            at all radii.
            
        """
        
        # Read the radial profile from file. We will load this into the 3D numberdensity_arr array.
        
        t               = Table.read(file_profile_ring, format = 'ascii')
        
        radius_km_file  = t['RadiusKM']
        IoF_file        = t['I/F']
        radius_pix_file = t['RadiusPixels']
        
        IoF_file        = np.clip(IoF_file, 0, None)   # Clip values to be positive.

        # Save the I/F values, in case we want to do a reality check later, or use them for normalization
        
        self.IoF_file       = IoF_file
        self.radius_km_file = radius_km_file
        
        # Now, we need to loop over the value of radius_ring, in 1-pixel steps.
        # For each step, fill everything at that radius and outward with value from IoF.
        
        # numberdensity = # of particles per km3.
        # This is taken to be the # of particles for the smallest size in n(r).
        # The entire size dist can be easily calculated from that.
        
        self.numberdensity_arr = self.radius_arr.copy()*0 
                
        # Get the positive x values. Use these as the radial bins to loop over
        
        radii = self.x_1d[self.x_1d > 0]  # Loop over the radius_ring in the output grid
        
        for radius_i in radii:                         # Loop over the radius_ring bin in the output grid 
            is_good = self.radius_ring_arr >= radius_i # Flag all bins at or beyond current radial bin
            self.numberdensity_arr[is_good] = IoF_file[ int(np.digitize(radius_i, radius_km_file)) ]
        
        # We will normalize this later!
        
        # Now bring in the vertical profile.
        # For each grid location, we want to calculate the value (vertical from midplane) / (radius_ring from center).
        # If this is < 'inclination', then we set it to 1. Otherwise, 0.
        # If no inclination is passed, then fill everything in the box.

        if (inclination):
            
            is_good_vertical_arr = ((np.abs(self.vertical_ring_arr) / self.radius_ring_arr) < inclination).astype(float)

            self.numberdensity_arr *= is_good_vertical_arr
            
        self.inclination_ring = inclination
        
        # Set the albedo
        
        self.albedo = albedo
                      
        # Define the size distribution. This size dist applies to the entire ring equally.

        self.r_dust = []
        self.n_dust = []
        self.q_dust = q_dust
        
        # Define the bins as per MRS 28-Jan-2017. These are the fixed sizes we output for Doug Mehoke.
        # Lower limit: 0.046 mm = 46 micron.
        
        r_dust = np.array([0.046, 0.1, 0.215, 0.46, 1, 2.15, 4.6])*u.mm.to('cm') # From MRS 28-Jan-2017
        
        # Make it a power law, spanning the bins.

        n_dust = r_dust**-q_dust
        
        # Normalize the size dist, s.t. smallest bin has exactly one particle.
        # The ring itself is defined s.t. numberdensity_arr is the number of smallest-bin particles per km3.
        
        n_dust = n_dust / n_dust[0]
        
        # Save the n(r) where we can get it later
        
        self.r_dust = r_dust
        self.n_dust = n_dust
        
    
        return

# =============================================================================
# Normalize the ring population
# =============================================================================

    def normalize(self, IoF_max=2e-7):
        """
        This changes numberdensity_arr s.t. I/F is the proper values.
        """
        
        num_radius = len(self.x_1d) / 2
        
        index_center = int(num_radius)
        
        radius_ring = self.radius_ring_arr[index_center:,0,index_center]
        
        # Calculate the I/F formally. See TPW04 eq 1.
        
        mu  = 1.  # Assume a face-on ring, for sunflower

        # Calculate the I/F along a vertical-radial slice in the Z-Y plane, from center outward 
        
        IoF = self.albedo * \
              np.sum(math.pi * self.n_dust * self.r_dust**2) * \
              self.numberdensity_arr[index_center:,:,index_center] / (4 * mu)

        # Sum the I/F vertically
        
        IoF = np.sum(IoF, axis=1)
        
        # And normalize number density based on this
        
        self.numberdensity_arr *= (IoF_max / np.amax(IoF))
        
        # Convert this to number per km3, from number 
        
        self.numberdensity_arr *= ( ((1*u.km) / (1*u.cm)).to('1').value ) **3
        
        # XXX NB: This is not yet finished! We still need to normalize based on the particles that 
        # Kaufmann / Hamilton see that actually exist. However, that is a second step. For now, 
        # get the proper population here.
        
        return
    
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
        
#        IoF_2d = np.sum(IoF, axis=2)
        
        # Make 2D plot of radial profile
        
        plt.plot(radius_ring, IoF)
        plt.xlabel('Radius [km]')
        plt.ylabel('I/F [arbitrary]')
        plt.show()
        
        # Plot an image of the face-on ring
        
        plt.imshow(self.numberdensity_arr[:,index_center,:], extent =[np.amin(self.x_1d), np.amax(self.x_1d), 
                                                     np.amin(self.z_1d), np.amax(self.z_1d)])
        plt.title('Ring, Face on')
        plt.xlabel('Z [km]')
        plt.ylabel('X [km]')
        plt.show()
        
        # Make a plot of the edge profile
        
        is_populated_arr = (self.numberdensity_arr > 0).astype(int)
        
        plt.imshow(is_populated_arr[index_center,:,:], extent =[np.amin(self.x_1d), np.amax(self.x_1d), 
                                                     np.amin(self.z_1d), np.amax(self.z_1d)])
        plt.xlabel('Radius [km]')
        plt.ylabel('Vertical [km]')
    
        plt.title('Edge Slice Profile, i = {}'.format(self.inclination_ring))
#        plt.xlim((-10,10))
        plt.show()
        
        # Plot the size distribution, just for fun
        
        plt.plot(self.r_dust, self.n_dust, marker = 'o', linestyle='none')
        plt.title('Size Distribution')
        plt.xlabel('Size [cm]')
        plt.ylabel('Number')
        plt.yscale('log')
        plt.xscale('log')
        plt.show()

# =============================================================================
# Output the trajectory and particle intercept info to a file
# =============================================================================
    def output_trajectory(self, do_positions=True):

        """
        Write an output table. The filename is auto-generated based on the parameters already supplied.
        
        Parameters
        -----
        do_positions:
            If True, output three columns which describe the XYZ position of the spacecraft as function of time.
        """
        
        # Create the table
        
        # We really don't care about fractional values on any of these. So, truncate everything
        # to be an integer. There is no easy way to force all columns to print as int, so just do it manually.
        
        t = Table([np.array(self.et_t).astype(int), 
                   np.array(self.delta_et_t).astype(int)], 
                      names = ['ET', 'ET from CA'])

        if do_positions:
            t['X [km]'] = np.array(self.x_t).astype(int)
            t['Y [km]'] = np.array(self.y_t).astype(int)
            t['Z [km]'] = np.array(self.z_t).astype(int)
            
        for i in range(len(self.n_dust)):
                t['n_{}, {:.3f} mm, # per km3'.format(i, self.r_dust[i])] = \
                                       np.array(self.numberdensity_t * self.n_dust[i]).astype(int)
        

        # Create the output filename
        
        file_out = 'mu69_hazard_q{}_i{}_a{}.txt'.format(self.q_dust, self.inclination_ring, self.albedo)
        if do_positions:
            file_out = file_out.replace('.txt', '_pos.txt')
            
        # Write the file

        t.write(file_out, format = 'ascii', overwrite=True)
            
        print("Wrote: {}".format(file_out))
         
        return
        
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
        
        num = math.ceil( (et_end - et_start) / dt ) + 1
        
        et_t = hbt.frange(int(et_start), int(et_end), num)
        
        # Calc offset in DT from C/A time
            
        delta_et_t = et_t - np.mean(et_t)
        
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
        self.delta_et_t    = delta_et_t
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
        
#        return()


            
# =============================================================================
# END OF METHOD DEFINITION
# =============================================================================
    
# =============================================================================
# The main function to call to run the simulation
# =============================================================================
    
def do_ring_flyby():
    
    albedo              = [0.05]
    q_dust              = [2]
    inclination_ring    = [None]
    
    file_profile_ring = '/Users/throop/Dropbox/Data/ORT1/throop/backplaned/K1LR_HAZ04/stack_n40_z4_profiles.txt'
    
# This defines I/F as a function of orbital distance.
# n(r) and q are then derived from I/F and albedo.
    
    utc_ca = '2019 1 Jan 05:33:00'
    dt_before = 1*u.hour
    dt_after  = 1*u.hour
    
    frame = '2014_MU69_SUNFLOWER_ROT'
    
    name_target = 'MU69'
    
    name_observer = 'New Horizons'
    
    dt = 1         # Sampling time through the flyby. Assumed to be seconds.

    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem.tm')
    
    # Define the number of grid locations on each side.
    
    n_dx = n_dy = n_dz = 201   # This is 'multiple assignment' in Python. Unlike IDL, it is official and legal.
    
    # Define the size of each grid box
    
    dx_bin = dy_bin = dz_bin = 250  # Distances are km. Not worth tracking the units. Consistent with SPICE.
                                    # XXX Wiki said 25 km, but should be 250.
    
    # Initialize the grid 
    
    num_pts_3d    = (n_dx,   n_dy,   n_dz)
    binsize_km_3d = (dx_bin, dy_bin, dz_bin)
    
    ring = ring_flyby(num_pts_3d, binsize_km_3d, frame, name_target)

    # Load the trajectory

    et_ca = int( sp.utc2et(utc_ca) )  # Force this to be an integer, just to make output cleaner.
    
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
                
                ring.set_ring_parameters(file_profile_ring, albedo_i, q_dust_i, inclination_i)

                ring.normalize()
                
                # Plot the ring profile
                
                ring.plot_radial_profile()
                
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
                
                ring.output_trajectory()
                
    self = ring
                
                
# =============================================================================
# Run the function
# =============================================================================
                
if (__name__ == '__main__'):
    do_ring_flyby()
    