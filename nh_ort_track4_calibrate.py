#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 14:41:44 2018

@author: throop
"""

import glob
import math       # We use this to get pi. Documentation says math is 'always available' 
                  # but apparently it still must be imported.

import os.path
import os

import astropy
from   astropy.table import Table
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
import matplotlib.pyplot as plt # pyplot
import numpy as np
import astropy.modeling
                       # Pylab defines the 'plot' command
#from   itertools import izip    # To loop over groups in a table -- see astropy tables docs

import pickle # For load/save
import scipy
from   get_radial_profile_circular import get_radial_profile_circular, get_profile_linear
import zlib  # For pickle decompression.

# HBT imports

import hbt

from nh_ort_track4_grid            import nh_ort_track4_grid    # Includes .read, .write, .plot, .flythru, etc.

from nh_ort_track3_read            import nh_ort_track3_read
from nh_ort_track3_read            import stretch_hbt, stretch_hbt_invert  # Some stretch routines good for trajectories
from nh_ort_track3_read            import plot_flattened_grids_table

#%%%
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
    
def nh_ort_track4_calibrate(dir_in, dir_out, runs, do_force=False):
    
    """
    This file does the 'calibration' for NH MU69 ORT 'Track 4'.
    Calibration refers to merging all of the Track-3 trajectories (from DPH/DK), in various 
    combinations, and saving them as grid files ('.grid4d.gz', or 
    
    I then use another program to fly thru these grids, and output as .dust file for Doug Mehoke.
    
    HBT 28-Mar-2018
    
    Parameters
    -----
    
    dir:
        Directory 
        
    runs:
        A list of all of the subdirectories that have the input grid files
        
    """
#%%%    
        
    name_run = dir_in.split('/')[-1]  # Extract the DPH/DK run name. Have to do hoops to avoid final '/' in path.
    if len(name_run) == 0:
        name_run = dir_in.split('/')[-2]
        
    plt.set_cmap('plasma')

    pi = math.pi
    
    do_compress = False             # Do we do .gzip compression on the output grids?
                                    # This is good for flexibility, but slows things down a lot.
    
    # Define the limit for the I/F we want to match.
    
    iof_limit_ring = 2e-7  # The actual observed ring in ORT2 Track1 was 5e-8.
                           # So, for the ORT2 'test' case, since we are using DPH's ring which is oriented
                           # differently and was not detected, I'll take a bit less than this -- say, I/F = 2e-8.
                           
                           # For ORT4 (the shepherd ring), the peak I/F is 2e-7.
                           # For ORT5 (ie, CHORE3), I am assuming 2e-7?
                     
    origin = 'lower'
                      
    stretch_percent = 98
    
    stretch = astropy.visualization.PercentileInterval(stretch_percent)
    
    # Define the axis to sum over: 0 → X axis, 1 → Y axis, 2 → Z axis.
    # Usually this will be Y dir (axis=1), but for ORT2, DPH's axes are incorrect, so we use X dir instead.
#    
#    if ('hamilton' in dir_in):
#        axis_sum = 0

    axis_sum = 1                                            # In general sum in the y direction. Do not correct for DPH.
    
    do_short = False
    
    if do_short:
        runs = runs[0:2]
        
    num_files             = len(runs)                  # Scalar, number of files
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
    
    weighting = np.zeros((len(runs), 1, 1)) # Weight each individual run. We make this into a funny shape
                                                 # so we can broadcast it with the full array easily.
    
#    run = runs[0].replace(dir_in, '')[1:]  # Populate one, just for cut & paste ease

    run = runs[0]                       # Populate one, just for cut & paste ease
    i = 0
    rings = []                                    # Make a list with all of the ring objects in it.
    
    hbt.figsize((10,10))
    run_name_base = runs[0].split('/')[-6]       # Extract the 'ort4_bc3_10cbr2_dph' portion if we need it
    
    # Parse the run name and guess if we are doing a tunacan here. Not guaranteed, but probably a good guess for tuna.
    
    is_tuna      = ('TUNA'      in name_run.upper())
    is_sunflower = ('SUNFLOWER' in name_run.upper())
    
    # Create the pickle save name.
    
    file_pickle = os.path.join(dir_out, f'grids_n{len(runs)}.pkl')
    
# =============================================================================
#     Restore from pickle file, if it exists
#     In reality, this routine doesn't save much time -- uncompressing the file
#     is nearly as slow as reading all the grids from disk.
#     Mac python cannot deal with the 3 GB uncompressed file sizes.    
# =============================================================================
    
    if (os.path.isfile(file_pickle) and (not do_force)):
        lun = open(file_pickle, 'rb')
        print(f'Reading pickle file {file_pickle}')
        
        (run, ring, rings_Z, beta_arr, time_step_arr, duration_arr, speed_arr, body_id_arr, subset_arr, grains_arr,
                     density_flattened_arr, halfwidth_km) = pickle.load(lun)

        rings = pickle.loads(zlib.decompress(rings_Z))
        lun.close()
        print(f'Read pickle file {file_pickle}')

# =============================================================================
# Otherwise, read the raw grids from disk         
# =============================================================================

    else:    
        # Loop over all of the files, and read each one in.
        # We read these into a list called 'ring' and rings'. Maybe 'grid' would be a better name, but we
        # use that below already.
        
        for i,run in enumerate(runs):
    
            run_shortname  = run.replace(dir_in, '')  # Remove the base pathname from this, and initial '/'
            ring = nh_ort_track3_read(run_shortname, dir_base=dir_in)            # Read the data array itself.
            
            do_print_ring_info = False
            do_plot_ring       = False
            
            if do_print_ring_info:
                ring.print_info()
            else:
                print(f'Read {i}/{len(runs)}: {run}')
    
            if do_plot_ring:
                ring.plot()
                
            # Save this entire ring (including 3D density array) to memory. 
            # We will use it later when we reconstruct the densities to fly thru.
            
            rings.append(ring)
            
            # Read various values from the ring
            
            halfwidth_km = ring.km_per_cell_x * hbt.sizex(ring.density) / 2
            extent = [-halfwidth_km, halfwidth_km, -halfwidth_km, halfwidth_km]  # Make calibrated labels for X Y axes
            
            beta_arr[i]                  = ring.beta  # Take the params from individual runs, and stuff into an array   
            time_step_arr[i]             = ring.time_step.value
            duration_arr[i]              = ring.duration.value
            speed_arr[i]                 = ring.speed
            body_id_arr[i]               = ring.body_id
            subset_arr[i]                = ring.subset
            grains_arr[i]                = ring.grains
    
            # Flatten the arrays, along the axis of the viewer (ie, the sun).
            # The output to this is bascially D2D of MRS slide 6.4 -- called 'density_flattened_arr' here.
            # This is a 2D image. When displayed with imshow(origin='lower'), oriented approx as seen from Sun.
            
            density_flattened_arr[i,:,:] = np.transpose(np.sum(ring.density, axis=axis_sum))
            
            if do_print_ring_info:
                print('-----')
    
        print(f'Read {num_files} files.')
  
        # Write pickle file to disk
        # Note that due to 'bug' in Python + OSX, file size is limited to ~4 GB. This means that saving the full
        # 'rings' array is not possible. Solution: compress the object in memory before writing pickle file.
        # This slows things down a lot -- not is not clear that in the end it's worth it.
        
        print(f'Writing: {file_pickle}')
        
        lun = open(file_pickle, 'wb')
        rings_Z = zlib.compress(pickle.dumps(rings))  #  Compress 'rings' in memory.
        pickle.dump((run, ring, rings_Z, beta_arr, time_step_arr, duration_arr, speed_arr, \
                     body_id_arr, subset_arr, grains_arr,\
                     density_flattened_arr, halfwidth_km), lun) 
        print(f'Wrote: {file_pickle}')
        lun.close()
    
#%%    
# =============================================================================
#  Now that we have read all of the input files, we combine these to make the output files. 
# =============================================================================
    
    # Usually this means summing sets of 8 'subsets', which are identical inputs, but starting a dfft azimuth.
    # We also choose the correct set of beta values to use.
    #
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
    q      = [-2.5, -3.5]            # Should be negative, since exponent in eq @ 6.4 is positive
    rho    = [1, 0.4641588, 0.2154434, 0.10000]   # 10**(-1/3), 10**(-2/3)
    speed  = [-2.2, -3]
    orb_sol= [1]                     # Just one case here for 14-Mar-2018 ORT2 case, but might have more in future.
    
    # Set the range of 'b', which is basically the range of beta. 
    # This was a 13 bins for ORT1 .. ORT3, and extended to 19 bins for ORT4.
    # b is the exponent used for beta, and radius, and bin width. We iterate over b, not beta.
    #
    # b = [-12,0] → beta = [1e-4, 1]
    # b = [-18,0] → beta = [1e-6, 1]

    if 'ORT3' in file_pickle:  # For ORT3 (and before, but we don't care about rerunning those)
        b = hbt.frange(-12,0)
        
    else:                       # For ORT4 and beyond
        b = hbt.frange(-18,0)

    r_moon_km = 1            # Radius of moon. THis will scale out in the end, but will change value of E_0.
    
    sarea = 4 * pi * (r_moon_km)**2  # Moon surface area. First term in MRS eq @ Slide 6.4 . This is in km2. 
                                     # The actual value here should be taken from Track1, and I don't have it here. 
                                     # 4 pi r^2 because this is production, which is assumed to go equally from 
                                     # all surfaces of the moon. Not cross-sectional sweepup.
    
    epsilon = 1e-5               # A small value for testing equality of floats.    
 
    #%%%
    
    # We create an astropy table for the output.
    
    t = Table(names=['albedo', 'q', 'rho', 'speed', 'val_img_med', 'val_img_typical', 'val_img_max', 
                     'img_2d', 'E_0', 'profile'],
              dtype = [float, float, float, float,    float,        float,             float, 
                      'object', float, 'object'])
    
    # Q: Don't we need another albedo term on sarea? Since if moons are dark, they may provide a lot more suface area.
    #    And, if dust is small, there can be a lot more dust. So, albedo might well enter as a square. I only have it
    #    as a single power.

    # Create an array for 'extent', which is for plotting. Base this on the size of the grid, which we now know.
    
    halfwidth_km = ring.km_per_cell_x * hbt.sizex(ring.density) / 2
    extent = [-halfwidth_km, halfwidth_km, -halfwidth_km, halfwidth_km]  # Make calibrated labels for X and Y axes



    (albedo_i, q_i, rho_i, speed_i) = (albedo[0], q[0], rho[0], speed[0]) # Set up values, for testing only. Can ignore.
    
    # =============================================================================
    # Now combine the runs of different particle sizes into output grids, and normalize to a target I/F value.
    # =============================================================================

    # Do the loop in output space. For each of the combinations of output requested, find the input runs 
    # that should go into that grid. This just takes a few seconds to run.
    
    hbt.figsize((10,5))
    
    # Just for reference, calculate the number of output files we will make. We use this just for plotting.
    
    num_combos = len(albedo) * len(q) * len(rho) * len(speed)
    
    for albedo_i in albedo:
        for q_i in q:
            for rho_i in rho:
                for speed_i in speed:
                
                    print() 
                    img = np.zeros((200, 200))  # Create the output array. Zero it for a new albedo/q/rho/speed combo.
                    
                    for b_i in b:  # This is the loop over particle size. b is index used for beta, particle size, etc.
                        
                        beta_i = 10**(b_i/3)
                      
                        s_i      = 5.7e-4 * (1 + 1.36 * albedo_i) / (rho_i * beta_i)  # Particle size. MRS slide 6.4
                        ds_i     = 0.7865*s_i

                        # Find all runs that match this set of parameters. 
                        # Typically this will match 8 'subsets' -- or 16 if we have two moons, etc.

                        is_subset = ( (speed_arr == speed_i) & ((np.abs(beta_i - beta_arr) / beta_i) < epsilon) )   
                        
                        num_subsets = np.sum(is_subset)  # Number of grids that are selected, to sum. 
                                                    # (ie, # of subsets)
                        
                        print(f'Found {num_good} matching runs for speed={speed_i:4.1f},' + 
                              f' beta={beta_i:5.2e}, s={s_i:7.3f} mm; applying q={q_i:3.2f},' + 
                              f' rho={rho_i:4.2f}, albedo={albedo_i} ')
                        print(f'  Indices = {np.where(is_good)[0]}')
                        
                        # Now apply MRS eq @ slide 6.4. Divide by num_subsets so we don't sum.
                        #   [Huh? What does this comment above mean? It is the key issue I think.]
                        # The np.sum(axis=0) here does the sum over all the matching subsets.
                        # They have already been flattened -- we are not summing in XYZ space again
                        # This time we want to pick out the density_flattened_arr for the proper size.
                        
                        # density_flattened_arr must be summed along the line-of-sight. I use it to make the I/F
                        # images, which we do here, and calc E_0. Then we go back to the grids, read those,
                        # and use them and E_0 to calc the final output.
                        
                        # Axis=0 means to sum along all of the 'is_good' arrays -- that is, ones that
                        # match the proper subset.
                        
                        # I am summing over subsets here. I am not sure about why dividing by num_good. That 
                        # seems to make it *averaging* over subsets, rather than *summing*.
                        
                        # Units of 'img' are s.t. img * E_0 = unitless.
                        
                        # NB: density_flattened_arr = [304, 200, 200]
                        
                        # XX Confirmed. This summation here looks OK.
                        
                        # XX We need to divide by num_subsets, because that is the number of subsets.
                        # By adding all these subsets, we increase the dust production. But we really only
                        # have one moon, not 8. So, must divide out by 8.
                        # This factor of 8 is *not* in the writeups. It should be.
                        # MRS didn't have it, but as long as you omit the same term in two places (when computing E_0,
                        #  and when making the grids), then those things cancel out.
                        
                        img_i = (sarea * albedo_i * pi * s_i**2 * ds_i * s_i**q_i *
                                 np.sum(density_flattened_arr[is_good], axis=0) /
                                 (ring.km_per_cell_x * ring.km_per_cell_y) * 1e-12) / num_good 
                        
                        # Add this image (with one beta) to the existing image (with a dfft beta)
                        
                        img += img_i
                        
                        # If requested, make a plot of the running total for this array
                        # This will make one plot per beta bin. Each plot will have all 
                        # subsets already summed.
                
                        do_plot_running_total = False
                        
                        if do_plot_running_total:
                            plt.subplot(1,2,1)
                            plt.imshow(stretch(img_i))
                            plt.title(f'img_i, b_i = {b_i}, max={np.amax(img_i):5.2e}')
                            plt.xlabel('X')
                            plt.ylabel('Z')
                            plt.subplot(1,2,2)
                            plt.imshow(stretch(img))
                            plt.title('Running total, summed along Y, front view')
                            plt.xlabel('X')
                            plt.ylabel('Z')
                            plt.show()
                            
                    # Do some statistics on the image -- brightest pixel, etc.
                    # We want to get the brightest pixel, etc. 
                    
                    val_img_max     = np.amax(img)               # Get brightest pixel level. Too high.
                    val_img_typical = np.percentile(img, 99)     # Get 'brightest region' level. Good to use.
                    val_img_med     = np.median(img[img != 0])   # Get median of non-zero pixels. Usually too low.
                    val_img_max_blur= np.amax(scipy.ndimage.filters.gaussian_filter(img,4))
                                                                 # As per MRS, blur it and take max

                    # Now calculate the actual calibration coefficient. 
                    # E_0 is the dust production rate, particles/km2/sec. 
                    
                    E_0_i = iof_limit_ring / val_img_max_blur
                    
                    # XXX Validated! The quantity np.percentile(E_0 * img,99) now is exactly the target I/F.
                    # E0 is a constant for the ring. It is not size-dependent (not a func of b, aka s).
                    
                    # Take a profile. This a radial profile in the case of a sunflower, 
                    # but a vertical profile in the case of a tunacan. This is because we want to get the 
                    # density at the point of flyby, which will be very nearly directly below the center, 
                    # in the Z direction (as labeled on the plot).
                    
                    if is_tuna:
                        (radius_pix, profile) = get_profile_linear(img, method='median', axis_move=1)
                    else:    
                        (radius_pix, profile) = get_radial_profile_circular(img, method='median')
                    
                    # Plot the radial profile, next to an image of that run, as seen from sun. 
                    
                    do_plot_profile_individual = False
                    
                    vals_radial_fiducial = [3500, 10000]

                    if do_plot_profile_individual:
                        plt.subplot(1,2,1)
                        plt.imshow(stretch_hbt(img), extent=extent)  
                        plt.title(f'{name_run}, {i}/{num_combos}')
                        plt.xlabel('X')
                        plt.ylabel('Z')
                        plt.subplot(1,2,2)
                        plt.plot(radius_pix * ring.km_per_cell_x, profile * E_0_i)
                        
                        # Plot grid lines at fixed distances, for references
                        
                        for val in vals_radial_fiducial:                        
                            plt.axvline(val, color='blue', alpha=0.3)
                            if is_tuna:
                                plt.axvline(-val, color='blue', alpha=0.3)
                        if is_tuna:        
                            plt.xlabel('Z [km]')
                        else:            
                            plt.xlabel('Radius [km]')
                        plt.ylabel('I/F Median in annulus')
        
                        plt.tight_layout()
                        plt.show()

                    # We are now finished with summing this array over particle size (ie, beta). Save it!
                    
                    t.add_row((albedo_i, q_i, rho_i, speed_i,
                               val_img_med, val_img_typical, val_img_max, 
                               img, E_0_i, profile))

    # For debugging, print value for y2.2_q2.5_pv0.05_rho1

    print( t['speed', 'q', 'albedo', 'rho', 'E_0'][0] )  # Check. I am *not* getting same E_0 value. Low by 8x.
    print( f'Total combinations made: {len(t)} combins. Each one has a 13 beta, and 8 subsets.')
    
#%%%
                    
# =============================================================================
#  Now that all combinations have been done, make some plots.
#  All of the output quantities are in the table 't' -- one entry per combination of parameters.
# =============================================================================

# Make a plot of radial (or linear) profile for each run, all stacked

    do_plot_profile_merged = True
    
    # Find the maximum value of all of the profiles. This is cool: use an iterator and a lambda function!
    # profile * E_0 is the final I/F value (I think). We first find which index has the max, and then look up the max.
    
    index_profile_max = max(range(len(t)), key = lambda i:np.amax(t['profile'][i] * t['E_0'][i]))
    val_max = np.amax(t['profile'][index_profile_max] * t['E_0'][index_profile_max])
    
    num_plots = 2
    yscales = ['log', 'linear']
    
    # Define some fiducual markers
    
    vals_radius_fiducial = [0, 1000, 2000, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    
    if do_plot_profile_merged:
        hbt.figsize((14,5))
        hbt.fontsize(14)
        for num_plot in range(num_plots):
            plt.subplot(2,1,num_plot+1)
            for i in range(len(t)):                 # Loop over run number
                plt.plot(radius_pix * ring.km_per_cell_x, 
                         t['profile'][i] * t['E_0'][i],
                         color='black', alpha=0.3)
            plt.yscale(yscales[num_plot]) # ** This is the difference: one log, one linear
            plt.ylabel('I/F median')
            plt.ylim((1e-9, val_max*1.1))
            for val in vals_radius_fiducial:
                plt.axvline(val, color='blue', alpha=0.3)
                if is_tuna:
                    plt.axvline(-val, color='blue', alpha=0.3)            

            # Make an appropriate label, depending on if tunacan or not

            if is_tuna:  
                plt.title(f'Profiles, {len(t)} runs, {name_run}' )
                plt.xlabel('Z [km]')
                
            else: 
                plt.title(f'Radial profiles, {len(t)} runs, {name_run}' )
                plt.xlabel('Radius [km]')

            plt.tight_layout()
        plt.show()

        hbt.fontsize()
        hbt.figsize()
        
#%%%
        
# Make a some diagnostic plots, and a beautiful family portrait, of all the combined 64 disks
    
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
    
    do_plot_flattened_grids = True
    
    print(f'Making plot of all N={len(t)} flattened grids that will be output to Doug Mehoke')
    print(f'These are plots in the XZ plane, as seen from Y dir')
    
    if (do_plot_flattened_grids):
        hbt.figsize((30,30))
        percent = 99                             # with new normalization, 98 is too hard a stretch. 95 OK.
        file_out = os.path.join(dir_out, f'plot_flattened_grids_n{len(runs)}.png')
        plot_flattened_grids_table(t,stretch_percent=percent, file_out = file_out) 
        hbt.figsize()

#%%%
            
# =============================================================================
# Now loop over the table of E0, and create the output grids which we will fly through for Doug Mehoke.
# This is a loop in output space. We could loop over (speed, albedo, q, rho).
# Or equivalently, we could loop over values in table 't', which has all combinations of these.
# ** We will do the latter! **    
# =============================================================================

    # Loop over all parameters in output space.
    # At each of these combinations, we will want to find the *input* parameters that match this --
    # specifically, all of the subsets.
    
    do_short = False

    indices_short = [52,5]
    
    k = 3
    t_i = t[0]
    
    if do_short:
        print (f'For output to Doug Mehoke, reducing output number from {len(t)} → {len(indices_short)}')
        t = t[indices_short]   # Index 4/64 is a good one to try - a classic 'question mark' to check proper orientation.

#%%%
        
    for k,t_i in enumerate(t):   # Loop over every element in the combination of output parameters,
                                 # which have already been tabulated in table 't'.  

#%%%        
        print(f'Starting output {k}/{len(t)}')
        
        # Calculate the indices of the smallest and largest grains, based on table.
        # This is so that we use a consistent, limited size range for all obs -- rather than
        # a different set of sizes depending on albedo, etc.
        
        rho_i    = t_i['rho']
        albedo_i = t_i['albedo']
        speed_i  = t_i['speed']
        q_i      = t_i['q']
        E_0_i    = t_i['E_0'] 

        # Now take the array for b_max. This is a transcription of MRS table @ 4.6

        arr = np.array( [[-6, -5, -4, -3],     # Transcription of table @ 4.6
                         [-6, -5, -4, -3],     # XXX I will have to redo this for ORT4!
                         [-5, -4, -3, -2], 
                         [-5, -4, -3, -2]] )
        
        # Given values for albedo and rho, look up b_max
        
        index_column = hbt.wheremin( np.abs(np.array(albedo) - albedo_i)) 
        index_row    = hbt.wheremin( np.abs(np.array(rho) - rho_i))
        
        b_max = arr[index_column, index_row]

        # For ORT4, we are using a larger range of sizes on the ouptut.
        # So, instead of always outputting 7 bins wide, we output 13 bins wide.

        if 'ORT3' in file_pickle:
            b_min = b_max - 6

        else:                               # For ORT4 and beyond
            b_min = b_max - 12

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
        
        # Now that we have found the values of b (= index to beta = size) to be used for this 
        # particular grid (ie, combo of q, velocity, etc), loop over them
        
        i = 0  # i stores the particle size index, that we loop over to sum.
        
        for b_i in np.arange(b_min, b_max+1):  # Sum over particle size
            
            D3D     = np.zeros((200, 200, 200))   # Create the 3D distribution, for a given particle size
            
            beta_i  = 10**(b_i/3)   # Calc beta so we can get size. MRS slide 5.5
            s_i     = 5.7e-4 * (1 + 1.36 * albedo_i) / (rho_i * beta_i)  # Particle size. MRS slide 6.4. Millimeters.
            ds_i    = 0.7865*s_i   # MRS @ slide 6.4

            # Find all runs in the input that match this parameter set
            # Typically this will be 8 -- since this is the number of subsets.
            # The indices calculated here are in the list of input arrays from DPH's Track 3.

            is_subset = ( (speed_arr == speed_i) & ((np.abs(beta_i - beta_arr) / beta_i) < epsilon) )   
            
            print(f'Found {np.sum(is_subset)} matching subsets to combine, with E0 = {E_0_i:8.3},' +
                           f'v = {speed_i}, q={q_i}, beta={beta_i:.3e}, rho={rho_i:.3}')
            print(f'  Indices = {np.where(is_subset)[0]} are the matching subsets')
            
            num_good = np.sum(is_good)  # Number of grids that are selected, to sum XXX NOT USED SO COMMENTED OUT
            
            # Finally, sum up all of the appropriate subsets, weighted by q, E_0, sarea, etc., into an output array
            # This array is in units of # per km3.  MRS slide 6.6.
            
            # D3D will usually be the sum of 13 values of b (=size) * 8 subsets (locations)
            
            # XX In the end we divide by num_subsets. This is *not* in MRS eq! But it should be.
            # If the num_subsets is omitted above when we compute E_0, then it needs to be omitted here too.
            
            for j in np.where(is_subset)[0]:
                D3D += E_0_i * sarea * (s_i**q_i) * ds_i * rings[j].density / dxdydz_km3 / num_subsets
                print(f'  Adding subset with index {j}')

# ** For reference, the equation for calculating E0 is:
#                        img_i = (sarea * albedo_i * pi * s_i**2 * ds_i * s_i**q_i *
#                                 np.sum(density_flattened_arr[is_subset], axis=0) /
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
        # This will plot the final merged data arrays
        
        do_plot_xyz_slices_merged = False

        if do_plot_xyz_slices_merged:
            hbt.figsize((20,20))
            hbt.fontsize(7)
            grids_i.plot(axis_sum=0)
            grids_i.plot(axis_sum=1)
            grids_i.plot(axis_sum=2)

        # Now plot the max optical depth summed for all sizes, as a reality check
        # This tau should be 'almost' the same as the target tau. Usually it will be slightly less.
        # due to the fact that the target tau is normalized using all bins of beta, and this output grid
        # excludes the smallest grains (which don't stick around for long anyhow, so only minor difference).
        
        do_plot_tau_merged = True
        
        if (do_plot_tau_merged):
            hbt.fontsize(10)
            hbt.figsize((6,6))
            grids_i.plot_tau()
        
        # Save the array to disk. Files are written as .grid4d.gz
        # These .grids4d.gz files are read by nh_org_track4_flyby.py, and create output for Doug Mehoke.
        
        grids_i.write(do_compress=do_compress, dir = dir_out, style = 'pickle')
        
        # Save the array to disk, in a format that MeshLab can read

        grids_i.write(style = 'xyz', dir = dir_out)
        
        # Save the 2D array as an image. This is so I can merge them in NH_ORT_MAKE_SUPERSTACK.PY.
        
        grids_i.write(style = 'image', dir=dir_out)
        
        print('---')
        
#%%%
           
# =============================================================================
# Run the function if file is called directly
# =============================================================================
    
if __name__ == '__main__':

#%%%
        # Search for the input files to use
    
# For Hamilton, sample is ~/Data/ORT2/hamilton/deliveries/sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/
#                          ort2-0003/y2.2/beta2.2e-01/subset04/grid.array2'
        

    do_force = False   # If set, read the raw N-body grid files from disk, rather than from a pre-saved pickle file

    do_short = False  # If set, just read a subset of the data ← Do not use
    # num_runs_max = 3  # ← This does not make sense to use here. Do not use it. This limits the number of input
                      #   files, but # of output files remains the same, and is what takes the time.

    # Set which of DPH or DEK runs we are doing here. Each run has ~304 grids in it. I will usually do four runs.
    
#    dir_in = '/Users/throop/data/ORT4/hamilton/ort4_bc3_10cbr2_dph/'   # Doug retrograde 
    # dir_in = '/Users/throop/data/ORT4/kaufmann/ort4_bc3_10cbr2_dek/'   # David retrograde 
    

    # Kaufmann ORT5 [Nov 2018]
    
    ## These directories must end with '/'
    
    # dir_in = '/Users/throop/data/ORT5/kaufmann/deliveries/chr3-sunflower3.5k/' # done
    # dir_in = '/Users/throop/data/ORT5/kaufmann/deliveries/chr3-sunflower10k/'
    # dir_in = '/Users/throop/data/ORT5/kaufmann/deliveries/chr3-tunacan10k/'
    # dir_in = '/Users/throop/data/ORT5/kaufmann/deliveries/chr3-sunflower10k-subsets16-dek'
    
    # Hamilton ORT5 [Nov 2018]
    
    # dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/sun10k_a/'  # Not sure what the diff btwn sun10k_{ab} is
    # dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/sun10k_b/'
    # dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/tuna9k/'
    dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/sun10k-DPH/'
    # dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/sun10kfast-DPH/'

    # dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/dph-tunacan10k/'
    # dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/dph-tunacan3.5k/'
    # dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/dph-sunflower10k/'
    # dir_in = '/Users/throop/data/ORT5/hamilton/deliveries/dph-sunflower3.5k/'

    # Make sure directory is properly terminated. This is what os.path is supposed to do, but does not!
    
    if dir_in[-1] != '/':
        dir_in += '/'

    # Get a list of all of the individual runs in the input dir.
    # DPH supplies a directory for each moon, while DK is one step thinner since he doesn't.
    # [fixed 6-Dec-2018; they are both the same now, with a moon directory]
            
    if ('kauf') in dir_in:
        runs = glob.glob(os.path.join(dir_in, '*/*/*/subset*/'))
    if ('hamilton' in dir_in):
        runs = glob.glob(os.path.join(dir_in, '*/*/*/subset*/'))
        
    if runs:
    
        # Create the output directory, and make sure it exists on disk
        #   (Do this by changing hamilton→throop, and kaufmann→throop)
        
        dir_out = dir_in.replace('hamilton', 'throop').replace('kaufmann', 'throop')
        
        if not os.path.exists(dir_out):
            os.makedirs(dir_out)
    
        if do_short:
            runs = runs[0:num_runs_max]
        
 
#%%%
        
# =============================================================================
#         Now do the run!
# =============================================================================

        nh_ort_track4_calibrate(dir_in, dir_out, runs, do_force = do_force)
        
        # NB: Takes about 70 seconds to decompress 304 files from pickle.
        #     Takes 68 seconds to read raw files from disk. ie, my pickling and compressing saves zero time.

    else:
        print(f'No files found in {dir_in}')
        