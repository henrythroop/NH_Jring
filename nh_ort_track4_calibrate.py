#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 14:41:44 2018

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
import struct

import re # Regexp
import pickle # For load/save

from   datetime import datetime

import scipy

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular

# HBT imports

import hbt

from nh_ort_track4_grid            import nh_ort_track4_grid    # Includes .read, .write, .plot, .flythru, etc.

from nh_ort_track3_plot_trajectory import nh_ort_track3_plot_trajectory
from nh_ort_track3_read            import nh_ort_track3_read
from nh_ort_track3_read            import stretch_hbt, stretch_hbt_invert  # Some stretch routines good for trajectories
from nh_ort_track3_read            import plot_flattened_grids_table

# =============================================================================
# 
# This is the main code to read Doug Hamilton's "Track 3" dust trajectories, and merge them to 
# create the 'Track 4' model rings. For NH MU69 Hazards Apr-2018.
#
# To run Track 4:
#
#   - Execute nh_track4_calibrate.py . This reads in all of DPH/DK's individual dust trajectories,
#     and merges them into '4D' dust grids, which are properly calibrated to match a given I/F.
#     Typically this reads in 108 files, and outputs 64 files, named *.grids4d.gz. These grids
#     are essentially just matrices (7, 200, 200, 200) with the dust density as a func of XYZ and grain size.
#
#   - Then execute nh_ort_track4_flyby.py. This reads all of the 64 grids files, and 
#     outputs a list of dust densities vs. time, for each one. 
#     Output is a table, essentially showing dust density (in # km-3) as a func of grain size, and time.
#     Typically 64 files, *.dust .

# =============================================================================
# Run the file
# =============================================================================
    
def nh_ort_track4_calibrate():
    
    """
    This file does the 'calibration' for NH MU69 ORT 'Track 4'.
    Calibration refers to merging all of the Track-3 trajectories (from DPH), in various 
    combinations, and saving them as grid files ('.grid4d.gz').
    
    I then use another program to fly thru these grids, and output as .dust file for Doug Mehoke.
    
    HBT 28-Mar-2018
    
    """
    
    plt.set_cmap('plasma')

    pi = math.pi
    
    # Define the limit for the I/F we want to match.
    
    iof_limit_ring = 2e-8  # The actual observed ring in ORT2 Track1 was 5e-8.
                           # So, for the ORT2 'test' case, since we are using DPH's ring which is oriented
                           # differently and was not detected, I'll take a bit less than this -- say, I/F = 2e-8.
    
    # Search for the input files to use
    
# For Hamilton, sample is ~/Data/ORT2/hamilton/deliveries/sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/
#                          ort2-0003/y2.2/beta2.2e-01/subset04/grid.array2'
    
    dir_base='/Users/throop/data/ORT2/hamilton/deliveries'    
    runs_full = glob.glob(os.path.join(dir_base, '*', '*', 'ort2-0003', '*', '*', '*')) # Hamilton - ORT actual
    runs_full = glob.glob(os.path.join(dir_base, '*', '*', 'ort2-0003_6apr18', '*', '*', '*')) # Hamilton - ORT actual

    stretch_percent = 98
    
    stretch = astropy.visualization.PercentileInterval(stretch_percent)
    
    # Define the axis to sum over: 0 → X axis, 1 → Y axis, 2 → Z axis.
    # Usually this will be Y dir (axis=1), but for ORT2, DPH's axes are incorrect, so we use X dir instead.
#    
#    if ('hamilton' in dir_base):
#        axis_sum = 0

    axis_sum = 1                                            # In general sum in the y direction. Do not correct for DPH.
    
    do_short = False
    
    if do_short:
        runs_full = runs_full[0:20]
        
    num_files             = len(runs_full)                  # Scalar, number of files
    density_flattened_arr = np.zeros((num_files, 200, 200)) # Actual data, vertically flattened
    beta_arr              = np.zeros(num_files)             # Beta used for each run. Scalar.
    time_step_arr         = np.zeros(num_files)             # Timestep. Each run uses a constant dt (!). Seconds.
    duration_arr          = np.zeros(num_files)             # Total time, in yr.
    speed_arr             = np.zeros(num_files)             # Exponent for speed ('y'). -3 or -2.2 
    body_id_arr           = np.zeros(num_files, dtype='U30')# This is a string, like "ort2-0003".
                                                            # All bodies are in all sims. But only body_id emits dust.  
    subset_arr            = np.zeros(num_files)             # Subset, ie where in orbit does it start. Int, 0 ..7.
    grains_arr            = np.zeros(num_files)             # Number of grains in the simulatioin. Typically 100. 
    obj_file_arr          = np.zeros(num_files, dtype='U90')# File listing the objects themselves - photometry, Track1
    state_file_arr        = np.zeros(num_files, dtype='U90')# State file -- orbit solution -- Track2
    
    weighting = np.zeros((len(runs_full), 1, 1)) # Weight each individual run. We make this into a funny shape
                                                 # so we can broadcast it with the full array easily.
    
#    run_full = runs_full[0].replace(dir_base, '')[1:]  # Populate one, just for cut & paste ease

    run_full = runs_full[0]                       # Populate one, just for cut & paste ease
    i = 0
    rings = []                                    # Make a list with all of the ring objects in it.
    
    hbt.figsize((10,10))
    
    # Loop over all of the files, and read each one in.
    # We read these into a list called 'ring' and rings'. Maybe 'grid' would be a better name, but we
    # use that below already.
    
    for i,run_full in enumerate(runs_full):

        run  = run_full.replace(dir_base, '')[1:]  # Remove the base pathname from this, and initial '/'
        ring = nh_ort_track3_read(run)             # Read the data array itself.
        ring.print_info()

        # Save this entire ring (including 3D density array) to memory. 
        # We will use it later when we reconstruct the densities to fly thru.
        
        rings.append(ring)
        
        # Read various values from the ring
        
        halfwidth_km = ring.km_per_cell_x * hbt.sizex(ring.density) / 2
        extent = [-halfwidth_km, halfwidth_km, -halfwidth_km, halfwidth_km]  # Make calibrated labels for X and Y axes
        
        beta_arr[i]                  = ring.beta  # Take the params from individual runs, and stuff into an array   
        time_step_arr[i]             = ring.time_step.value
        duration_arr[i]              = ring.duration.value
        speed_arr[i]                 = ring.speed
        body_id_arr[i]               = ring.body_id
        subset_arr[i]                = ring.subset
        grains_arr[i]                = ring.grains

        # Flatten the arrays, along the axis of the viewer (ie, the sun).
        # The output to this is bascially D2D of MRS slide 6.4 -- called 'density_flattened_arr' here.
        
        density_flattened_arr[i,:,:] = np.sum(ring.density, axis=axis_sum)

        # Take the 3D arrays and flatten them, just so we can visualize more easily.
        # To get a full picture, we sum along all three axes (X, Y, Z) and show each.
        # XXX This code should be deprecated. Use the code in nh_ort_track3_plot_trajectory instead.

        do_plot_xyz_views = False
        
#        if do_plot_xyz_views:
#            plt.subplot(1,3,1)
#            plt.imshow(stretch_hbt(np.sum(ring.density, axis=0)), extent=extent)
#            plt.title('Summed along X')
#            plt.subplot(1,3,2)
#            plt.imshow(stretch_hbt(np.sum(ring.density, axis=1)), extent=extent)
#            plt.title('Summed along Y')
#            plt.subplot(1,3,3)
#            plt.imshow(stretch_hbt(np.sum(ring.density, axis=2)), extent=extent)
#            plt.title('Summed along Z')
#            plt.show()
        
        print('-----')

    print(f"Read {num_files} files.")
    
    num_subsets = len(np.unique(subset_arr))
    
    # Now that we have read all of the input files, we combine these in various ways to make 
    # the output files. 
    # We have 4 parameters to iterate over → 4 x 2 x 4 x 2 = 64 output cases.
    # Each output case yields an image in I/F. 
    # I'll calibrate based on that I/F, and then create paths for Doug Mehoke to fly.

    # Units of everything are:
    #  rho:        Density, g/cm3       [Slide 6.4]
    #  delta_{xy}: Grid sizes, km.      [Slide 6.4]
    #  s:          Particle radius. mm. [Slide 6.4]
    #  density_flattened_arr (and DPH data): No units. Just particles. This is what I calibrate to.
    #  albedo:     No units
    #  ds:         Bin width. mm
    #  R:          Moon radius. km.     [Slide 6.4 inferred from E0.] [Also, 4.3 says that area is delivered in km2.]
    
    albedo = [0.05, 0.1, 0.3, 0.7]   # Q: Are albedo and beta independent? A: Yes, completely indep -- see eq @ 3.7. 
    q      = [-2, -3.5]              # Should be negative, since exponent in eq @ 6.4 is positive
    rho    = [1, 0.4641588, 0.2154434, 0.10000]   # 10**(-1/3), 10**(-2/3)
    speed  = [-2.2, -3]
    orb_sol= [1]                     # Just one case here for 14-Mar-2018 ORT2 case, but might have more in future.
        
    b      = hbt.frange(-12,0)  # Exponent for beta, and radius, and bin width. 
                                #  b=[-12,0] → beta=[1e-4,1], which is the range that DPH uses.
                                # We iterate over 'b', rather than over beta or size.
                                
#    s      = 10.**b/3            # Particle size, some units?
    
    epsilon = 1e-5               # A small value for testing equality of floats.    
    
    # We create an astropy table for the output.
    
    t = Table(names=['albedo', 'q', 'rho', 'speed', 'val_img_med', 'val_img_typical', 'val_img_max', 
                     'img_2d', 'E_0', 'profile'],
              dtype = [float, float, float, float,    float,        float,             float, 
                      'object', float, 'object'])
    
    sarea = pi * (5)**2  # Moon surface area. First term in MRS eq @ Slide 6.4 . This is in km2. 
                         # The actual value here should be taken from Track1, and I don't have it here. 
    
    # Q: Don't we need another albedo term on sarea? Since if moons are dark, they may provide a lot more suface area.
    #    And, if dust is small, there can be a lot more dust. So, albedo might well enter as a square. I only have it
    #    as a single power.

    (albedo_i, q_i, rho_i, speed_i) = (albedo[0], q[0], rho[0], speed[0]) # Set up values, for testing only. Can ignore.
    
    # Now, do the loop in output space
    
    for albedo_i in albedo:
        for q_i in q:
            for rho_i in rho:
                for speed_i in speed:
                
                    img = np.zeros((200, 200))  # Create the output array
                    
                    for b_i in b:  # This is the loop over particle size. b is index used for beta, particle size, etc.
                        
                        beta_i = 10**(b_i/3)
                      
                        s_i      = 5.7e-4 * (1 + 1.36 * albedo_i) / (rho_i * beta_i)  # Particle size. MRS slide 6.4
                        ds_i     = 0.7865*s_i

                        # Find all runs that match this set of parameters. Typically this will match 8 'subsets.'
                        
                        is_good = ( np.logical_and( (speed_arr == speed_i),
                                      np.logical_and( ((np.abs(beta_i - beta_arr) / beta_i) < epsilon),
                                                      (body_id_arr == 'ort2-0003') ) ) )
                        
                        print(f'Found {np.sum(is_good)} matching runs for speed={speed_i:4.1f},' + 
                              f' beta={beta_i:5.2e}, s={s_i:7.3f} mm; applying q={q_i:3.2f},' + 
                              f' rho={rho_i:4.2f}, albedo={albedo_i} ')
                        
                        # Now apply MRS eq @ slide 6.4. Divide by num_subsets so we don't sum.
                        # The np.sum(axis=0) here does the sum over all the matching subsets.
                        # They have already been flattened -- we are not summing in XYZ space again
                        # This time we want to pick ut the density_flattened_array for the proper size.
                        
                        # Axis=0 means to sum along all of the 'is_good' arrays -- that is, ones that
                        # match the proper subset.
                        
                        # Units of 'img' are s.t. img * E_0 = unitless.
                        
                        img_i = (sarea * albedo_i * pi * s_i**2 * ds_i * s_i**q_i *
                                 np.sum(density_flattened_arr[is_good], axis=0) /
                                 (ring.km_per_cell_x * ring.km_per_cell_y) * 1e-12) / num_subsets
                        
                        # Add this image (with one beta) to the existing image (with a dfft beta)
                        
                        img += img_i
#                        print(f'Adding bin for size {s_i} mm')
                        
                        # If requested, make a plot of the running total for this array
                
                        do_plot_running_total = False
                        
                        if do_plot_running_total:
                            plt.subplot(1,2,1)
                            plt.imshow(img_i)
                            plt.title(f'img_i, b_i = {b_i}, max={np.amax(img_i):5.2e}')
                            plt.subplot(1,2,2)
                            plt.imshow(img)
                            plt.title('Running total')
                            plt.show()

                    # Do some statistics on the image -- brightest pixel, etc.
                    # We want to get the brightest pixel, etc. 
                    
                    val_img_max     = np.amax(img)               # Get brightest pixel level. Too high.
                    val_img_typical = np.percentile(img, 99)     # Get 'brightest region' level. Good to use.
                    val_img_med     = np.median(img[img != 0])   # Get median of non-zero pixels. Usually too low.

                    # Now calculate the actual calibration coefficient. 
                    # E_0 is the dust production rate, particles/km2/sec. 
                    
                    E_0_i = iof_limit_ring / val_img_typical
                    
                    # XXX Validated! The quantity np.percentile(E_0 * img,99) now is exactly the target I/F.
                    # E0 is a constant for the ring. It is not size-dependent (not a func of b, aka s).
                    
                    # Take a radial profile
                    
                    (radius_pix, profile) = get_radial_profile_circular(img)
                    
                    # Plot the radial profile
                    # XXX For now, disable this, because DPH's inputs are not even in the right plane.
                    # It is pointless to take a radial profile.
                    
                    do_plot_radial_profile = False
                    
                    if do_plot_radial_profile:                    
                        plt.imshow(stretch_hbt(img))
                        plt.show()

                    # We are now finished with summing this array over particle size (ie, beta). Save it!
                    
                    t.add_row((albedo_i, q_i, rho_i, speed_i,
                               val_img_med, val_img_typical, val_img_max, 
                               img, E_0_i, profile))
                
# =============================================================================
#  Now that all combinations have been done, make some plots.
#  All of the output quantities are in the table 't' -- one entry per combination of paramters.
# =============================================================================

# Make a plot of radial profile for each run, all stacked

    do_plot_profile_radial = False
    
    if do_plot_profile_radial:
        hbt.figsize((8,6))
        for i in range(len(t)):        
            plt.plot(radius_pix * ring.km_per_cell_x, t['profile'][i])
        plt.yscale('log')    
        plt.xlabel('Radius [km]')
        plt.ylabel('I/F (non-normalized)')
        plt.ylim((1e-12, 1e-6))
        plt.show()
    
# Make a gridded layout of various quantities
    
    columns = ['albedo', 'q', 'rho', 'speed']
    for i,column in enumerate(columns):
        plt.subplot(2,2,i+1)
        plt.plot(t[column], t['val_img_max'], marker = 'o', markersize=3, ls = 'none')
        plt.yscale('log')
        plt.xlabel(column)
        plt.ylabel('Max I/F')
    plt.tight_layout()    
    plt.show()
    
    # Make a plot of E_0 vs typical brightness
    
    hbt.figsize((8,6))    
    plt.plot(t['val_img_typical'], t['E_0'], marker = 'o', ls='none', color = 'red', label = 'Typical')
    plt.plot(t['val_img_med'], t['E_0'], marker = 'o', ls='none', color = 'blue', label = 'Median', ms=5)
    plt.plot(t['val_img_max'], t['E_0'], marker = 'o', ls='none', color = 'green', label = 'Max', ms=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('E_0')
    plt.xlabel('val_imgl')
    plt.legend()
    plt.show()
    
#    plt.plot(t['val_img_max'], t['val_img_med'], marker = 'o', ls = 'none')
#    plt.yscale('log')
#    plt.xscale('log')
#    
    
# Make a plot showing all of the flattened disks. This plots the t['img_2d'] field from the table.
# Each disk is stretched individually.
    
    do_plot_flattened_grids = False
    
    if (do_plot_flattened_grids):
        hbt.figsize((30,30)) 
        plot_flattened_grids_table(t,stretch_percent=98)
     
# =============================================================================
# Now loop over the table of E0, and create the output files for DM
# This is a loop in output space. We could loop over (speed, albedo, q, rho).
# Or equivalently, we could loop over values in table 't', which has all combinations of these.
# ** We will do the latter! **    
# =============================================================================

    # Loop over all paramters in output space.
    # At each of these combinations, we will want to find the *input* parameters that match this --
    # specifically, all of the subsets.
    
    do_short = False

    if do_short:
        t = t[52:55]   # Index 4/64 is a good one to try - a classic 'question mark' to check proper orientation.
        
    for k,t_i in enumerate(t):   # Loop over every element in the combination of output parameters,
                                 # which have already been tabulated in table 't'.  
        
        print(f'Starting output {k}/{len(t)}')
        
        # Calculate the indices of the smallest and largest grains, based on table.
        # This is so that we use a consistent, limited size range for all obs -- rather than
        # a different set of sizes depending on albedo, etc.
        
        rho_i    = t_i['rho']
        albedo_i = t_i['albedo']
        speed_i  = t_i['speed']
        q_i      = t_i['q']
        E_0_i    = t_i['E_0'] 
        
        arr = np.array( [[-6, -5, -4, -3],     # Transcription of table @ 4.6
                         [-6, -5, -4, -3], 
                         [-5, -4, -3, -2], 
                         [-5, -4, -3, -2]] )
        
        # Given values for albedo and rho, look up b_max
        
        index_column = hbt.wheremin( np.abs(np.array(albedo) - albedo_i)) 
        index_row    = hbt.wheremin( np.abs(np.array(rho) - rho_i))
        
        b_max = arr[index_column, index_row]
        b_min = b_max - 6
        beta_min = 10**(b_min/3)
        beta_max = 10**(b_max/3)
        s_min = 5.7e-4 * (1 + 1.36 * albedo_i) / (rho_i * beta_max)
        s_max = 5.7e-4 * (1 + 1.36 * albedo_i) / (rho_i * beta_min)
        
        print(f'Using b = {b_min} .. {b_max} → s = {s_min:.3f} .. {s_max:.3f} mm')

        # Get the volume of each cell. It is fixed -- no need to calc in every loop iteration.
        
        dxdydz_km3 = rings[0].km_per_cell_x * rings[0].km_per_cell_y * rings[0].km_per_cell_z
        
        num_b = b_max - b_min + 1
        
        # Create a 4D output array, which will hold the 3D distribution, for each grain size.
        # They are different, because each grain size has a different trajectory.
        
        D4D   = np.zeros((num_b, 200, 200, 200))
        
        hbt.figsize((6,6))
        
        # Now that we have found the values of b (=size), loop over them
        
        i = 0  # i stores the particle size index, that we loop over to sum.
        
        for b_i in np.arange(b_min, b_max+1):  # Sum over particle size
            
            D3D     = np.zeros((200, 200, 200))   # Create the 3D distribution, for a given particle size
            
            beta_i  = 10**(b_i/3)   # Calc beta so we can get size. MRS slide 5.5
            s_i     = 5.7e-4 * (1 + 1.36 * albedo_i) / (rho_i * beta_i)  # Particle size. MRS slide 6.4. Millimeters.
            ds_i    = 0.7865*s_i   # MRS @ slide 6.4

            # Find all runs in the input that match this parameter set
            # Typically this will be 8 -- since this is the number of subsets.
            # The indices calculated here are in the list of input arrays from DPH's Track 3.
            
            is_good = ( np.logical_and( (speed_arr == speed_i),
                          np.logical_and( ((np.abs(beta_i - beta_arr) / beta_i) < epsilon),
                                          (body_id_arr == 'ort2-0003') ) ) )

            # Finally, sum up all of the appropriate subsets, weighted by q, E_0, sarea, etc., into an output array
            # This array is in units of # per km3.  MRS slide 6.6.
            
            for j in np.where(is_good)[0]:
                D3D += E_0_i * sarea * (s_i**q_i) * ds_i * rings[j].density / dxdydz_km3

# ** For reference, the equation for calculating E0 is:
#                        img_i = (sarea * albedo_i * pi * s_i**2 * ds_i * s_i**q_i *
#                                 np.sum(density_flattened_arr[is_good], axis=0) /
#                                 (ring.km_per_cell_x * ring.km_per_cell_y) * 1e-12) / num_subsets
#                        E_0_i = iof_limit_ring / max(img)

            # And save the 3D array for this particle size, into the 4D array
            
            D4D[i] = D3D
            i += 1
        
        # Plot slices thru these cubes, to the screen
        
        # Convert the input parameters from index 'b', to radiation param beta, and size s. Make all arrays
        
        b    = np.arange(b_min, b_max+1)                        # Particle index
        beta = 10**(b/3)                                        # QR parameter. Larger, for smaller grains
        s    = 5.7e-4 * (1 + 1.36 * albedo_i) / (rho_i * beta)  # Particle size, in mm
        
        # Now created a 'track4_grid' object to store this new 4D array  

        grids_i = nh_ort_track4_grid(D4D)
        grids_i.set_parameters(albedo = albedo_i, speed=speed_i, q=q_i, rho=rho_i, 
                               b = list(b), s = list(s), beta = list(beta),
                               resolution_km = (ring.km_per_cell_x, ring.km_per_cell_y, ring.km_per_cell_z))  
        
        # Make a plot of this array, in various slices
        
        do_plot_xyz_slices = True

        if do_plot_xyz_slices:
            hbt.figsize((20,20))
            hbt.set_fontsize(7)
            grids_i.plot(axis_sum=0)
            grids_i.plot(axis_sum=1)
            grids_i.plot(axis_sum=2)

        # Now plot the max optical depth summed for all sizes, as a reality check

        do_plot_tau = True
        
        if (do_plot_tau):
            hbt.set_fontsize(10)
            hbt.figsize((6,6))
            grids_i.plot_tau()
        
        # Save the array to disk. Files are written as .grid4d.gz
        # These .grids4d.gz files are read by nh_org_track4_flyby.py, and create output for Doug Mehoke.
        
        grids_i.write()
        
        print('---')
        
           
# =============================================================================
# Run the function
# =============================================================================
    
if __name__ == '__main__':

    nh_ort_track4_calibrate()
    