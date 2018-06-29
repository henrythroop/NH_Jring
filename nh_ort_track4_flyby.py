#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 23:02:00 2018

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

#from nh_ort_track3_plot_trajectory import nh_ort_track3_plot_trajectory
#from nh_ort_track3_read            import nh_ort_track3_read
#from nh_ort_track3_read            import stretch_hbt, stretch_hbt_invert  # Some stretch routines good for traj's
#from nh_ort_track3_read            import plot_flattened_grids_table

# =============================================================================
# Main function to fly s/c through the grids. All of the 4D grids must already 
# be generated and saved to disk. This routine loops over them, and creates the 
# output files for Doug Mehoke, of particle density vs. time.   
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
#
# This is a regular function, but it calls the class method nh_ort_track4_grid.fly_trajectory().
# =============================================================================

def nh_ort_track4_flyby(dir_in=None, dir_out=None, name_trajectory = 'prime'):

    #%%%
#    name_trajectory = 'alternate'  # Can be 'prime' or 'alternate'
    

#    dir_in = '/Users/throop/
#    dir_in = '/Users/throop/data/ORT4/throop/ort4_bc3_10cbr2_dph/'
    
    stretch_percent = 99
    stretch = astropy.visualization.PercentileInterval(stretch_percent)

#    dir_data = os.path.expanduser('~/Data/')
    
    dir_in 
    do_compress = False   # Do we use .gzip compression on the Track-4 input grids?
                          # If we used compression on the track4_calibrate routine, we must use it here too.
    
#    dir_track4 = os.path.join(dir_data, name_ort, 'throop', 'track4')
    
    if do_compress:
        files = glob.glob(os.path.join(dir_in, '*.grid4d.gz'))
    
    else:
        files = glob.glob(os.path.join(dir_in, '*.grid4d'))
        files = glob.glob(os.path.join(dir_in, '*.dust.pkl'))

    # Alphabetize file list
    
    files = sorted(files)
    
    plt.set_cmap('plasma')

    utc_ca = '2019 1 Jan 05:33:00'
    dt_before = 1*u.hour
    dt_after  = 1*u.hour
    
#    area_sc = (1*u.m)**2

    frame = '2014_MU69_SUNFLOWER_ROT'
    
    name_target = 'MU69'
    origin = 'lower'   # Required plotting order for imshow
    
    name_observer = 'New Horizons'
    
    hbt.figsize((8,6))
    hbt.set_fontsize(12)

    dt = 1*u.s         # Sampling time through the flyby. Astropy units.

    # Create an output table, Astropy format
    
    t = Table(names = ['trajectory', 'speed', 'q_dust', 'albedo', 'rho',
                       'tau_max', 'tau_typical', 'iof_max', 'iof_typical'],
              dtype = ['U30', float, float, float, float, 
                       float, float, float, float]  )
    
    # Start up SPICE if needed. Unload old kernels just as a safety precaution.
    
    sp.unload('kernels_kem_prime.tm')
    sp.unload('kernels_kem_alternate.tm')
    
    sp.furnsh(f'kernels_kem_{name_trajectory}.tm')
    
    do_short = False
    
    if do_short:
        files = files[0]
    
    file = files[27]  # Just for testing w cut+paste. Can ignore these.
    
    i=0
#%%%    
    for i,file in enumerate(files):
        
#%%%        
        print(f'Starting file {i}/{len(files)}')
              
        grid = nh_ort_track4_grid(file)    # Load the grid from disk. Uses gzip, so it is quite slow (10 sec/file)
         
    # Load the trajectory parameters

        et_ca = int( sp.utc2et(utc_ca) )  # Force this to be an integer, just to make output cleaner.
        
        et_start = et_ca - dt_before.to('s').value
        et_end   = et_ca + dt_after.to('s').value
                        
        grid.frame       = frame
        grid.name_target = name_target
        grid.name_trajectory = name_trajectory
        
        # And call the method to fly through it!
        # The returned density values etc are available within the instance variables, not returned explicitly.

        grid.fly_trajectory(name_observer, et_start, et_end, dt)
 
        # If the first time thru loop, make plot of our path through the system
     
        do_plots_geometry = True
        
        if (do_plots_geometry and (i==0)):
            
            grid.plot_trajectory_geometry()

        # Make slice plots thru the grid
        
        do_plot_slices_xyz = False

        if do_plot_slices_xyz:
            hbt.fontsize(8)
            hbt.figsize((20,5))
            grid.plot(axis_sum=0)
            grid.plot(axis_sum=1)
            grid.plot(axis_sum=2)
            
            hbt.fontsize(10)
        
        # Make a plot of optical depth
        
        do_plot_tau = True
        if do_plot_tau:
            grid.plot_tau()

# =============================================================================
# Make some plots of count rate vs. time!
# =============================================================================
                    
        # Make a plot of the instantaneous count rate

        hbt.figsize((10,15))

        # Make a plot of the actual density that we give to Doug Mehoke
        
        # Define a list of colors. This is so we can use colors= argument to set
        # a marker to show grain size, rather than let plot() auto-assign.
        
#        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # This is the default color iterator.
        colors = ['antiquewhite',
                  'tomato', 
                  'blueviolet',
                  'skyblue',
                  'gold',
                  'darkcyan',
                  'thistle',
                  'olive',
                  'red',
                  'sienna',
                  'deepskyblue',
                  'lightsalmon',
                  'pink',
                  'aqua']
                  
#                  'antiquewhite4', 'aqua', 'aquamarine4', 'black', 'blue', 'blueviolet', 
#                  'brown1', 'chartreuse1', 'darkgreen', 'darkorange1', 'dodgerblue1', 'lightpink', 'magenta']
        plt.subplot(3,1,1)
        for j,s in enumerate(grid.s):
            plt.plot(grid.delta_et_t, grid.number_t[j],
                     label = 's={:.2f} mm'.format(s))
        plt.legend()
        plt.title('Dust number density'.format(grid.area_sc))
        plt.xlabel('ET from C/A')
        plt.yscale('log')
        plt.ylabel(r'Dust,, # km$^{-3}$')

        plt.subplot(3,1,2)
        for j,s in enumerate(grid.s):   # 's' is dust size
            plt.plot(grid.delta_et_t, grid.number_sc_t[j],
                     label = f's={s:.2f} mm')
        plt.title('Impact rate, A={}'.format(grid.area_sc))
        plt.yscale('log')    
        plt.xlabel('ET from C/A')
        plt.legend()
        plt.ylabel(r'Dust, # Impacts sec$^{{-1}}$')

        # Make a plot of the cumulative count rate. Mark grain sizes here too.
        
        plt.subplot(3,1,3)
        for j,s in enumerate(grid.s):                                             # Loop over size
            plt.plot(grid.delta_et_t, grid.number_sc_cum_t[j],                    # Main plot line
                     label = 's={:.2f} mm'.format(s), color=colors[j])
            plt.plot([grid.delta_et_t[-1]], [grid.number_sc_cum_t[j,-1].value],   # Circle to indicate grain size
                     markersize=(7-j)*2, marker = 'o',                            # Use same color as prev line! 
                     color=colors[j])


        hbt.figsize(5,5)
        plt.legend()
        plt.title('Number of impacts (cumulative), A={}'.format(grid.area_sc))
        plt.xlabel('ET from C/A')
        plt.yscale('log')
        plt.ylabel('# of Impacts')
        plt.axhline(y = 1, linestyle = '--', alpha = 0.1)    

        plt.tight_layout()
        
        plt.show()
        
        # Make a plot of size distibution. 
        # Make two curves: one for n(r) for the entire grid, and one for n(r) that hits s/c
        
        # Now add an entry to the table. This is a table that lists all of the results --
        #     e.g., max_tau, count rate etc
        # One line per grid.

        t.add_row(vals=[grid.name_trajectory, grid.speed, grid.q, grid.albedo, grid.rho, 
                        grid.tau_max, grid.tau_typ,
                        grid.iof_max, grid.iof_typ])
                                        
        # Get size dist along path
         
        number_path = grid.number_sc_cum_t[:,-1].value
        
        # Take the full particle grid, and sum along all spatial axes, leaving just the size axis left.
        
        number_grid = np.sum(np.sum(np.sum(grid.density, axis=1), axis=1), axis=1)
        
        # Normalize the size dists both
        number_grid = hbt.normalize(number_grid)
        number_path = hbt.normalize(number_path)
        
        plt.plot(grid.s, number_path, label = 'Along s/c path')
        plt.plot(grid.s, number_grid, label = 'In grid, total')
        plt.yscale('log')
        plt.xscale('log')
        plt.ylim( (hbt.minval(np.array([number_grid, number_path]))/2, 1) )
        plt.xlabel('Radius [mm]')
        plt.ylabel('Particle number [arbitrary]')
        plt.legend(loc = 'lower right')
        plt.show()

        # Output the dust population for this run to a file. This is the file that Doug Mehoke will read.
        
        grid.output_trajectory(suffix=f'{name_trajectory}', do_positions=False, dir_out=dir_out)
    
        print('---')
                        

#%%%        
        
    # Print the table
    
    t.pprint(max_width=-1)
    
    # Save the table as output
    
    file_out = os.path.join(dir_out, f'nh_{name_trajectory}_track4_table.pkl')
    
    lun = open(file_out, 'wb')
    pickle.dump(t,lun)
    lun.close()
    print(f'Wrote: {file_out}')
    
#%%%    

# =============================================================================
# Plot some tabulated results.
# One-off function.
# =============================================================================

def plot_table():
    
    file_in = '/Users/throop/Data/ORT2/throop/track4/nh_ort_track4_table.pkl'
    
    lun = open(file_in, 'rb')
    t = ge.load(lun)
    lun.close()
    
    hbt.figsize((10,8))
    
    xvals = ['albedo', 'speed', 'rho', 'q_dust']
    i = 0
    for xval in xvals:
        plt.subplot(2,2,i+1)
#        plt.plot(t[xval], t['iof_max'], marker = 'o', linestyle='none', label = 'I/F Max')
        plt.plot(t[xval], t['iof_typical'], marker = '+', linestyle='none', label = 'I/F Typical')
        plt.xlabel(xval)
        plt.yscale('log')
        plt.legend()
        i +=1
    plt.tight_layout()    
    plt.show()

    t.sort(['speed', 'q_dust', 'albedo', 'rho'])
    t.pprint(max_width=-1, max_lines=-1)

# =============================================================================
# Output table indices for MRS
# =============================================================================

def make_table_grid_positions():
    
    """
    This is a one-off utility function for MRS. 
    In it, I just do a flyby of MU69, and output the X Y Z grid indices (as well as positions and timestamps).
    I don't output density at all -- just the s/c positions.
    
    I do this for both prime and alternate trajectories.
    
    This is just because he hasn't integrated SPICE into his grid code.
    
    This function is stand-alone. It doesn't rely on the grid class.
    It is included in this file because it directly relates to the grids.
    
    """
     
    name_trajectory = 'alternate'  # ← Set this to 'prime' or 'alternate'
   
    sp.unload('kernels_kem_prime.tm')
    sp.unload('kernels_kem_alternate.tm')
   
    frame = '2014_MU69_SUNFLOWER_ROT'
    
    name_observer = 'New Horizons'

    name_target = 'MU69'
 
    sp.furnsh(f'kernels_kem_{name_trajectory}.tm')
    
    file_in = '/Users/throop/Data/ORT2/throop/track4/ort2-ring_v2.2_q2.0_pv0.10_rho0.22.grid4d.gz'
    
    grid = nh_ort_track4_grid(file_in)    # Load the grid from disk. Uses gzip, so it is quite slow (10 sec/file)

    utc_ca = '2019 1 Jan 05:33:00'
    dt_before = 1*u.hour
    dt_after  = 1*u.hour
    dt = 1*u.s         # Sampling time through the flyby. Astropy units.
    
    et_ca = int( sp.utc2et(utc_ca) )  # Force this to be an integer, just to make output cleaner.
        
    et_start = et_ca - dt_before.to('s').value
    et_end   = et_ca + dt_after.to('s').value
                    
    grid.frame       = frame
    grid.name_target = name_target
    grid.name_trajectory = name_trajectory
    
    # And call the method to fly through it!
    # The returned density values etc are available within the instance variables, not returned explicitly.

    grid.fly_trajectory(name_observer, et_start, et_end, dt)

    # Make plots
    
    hbt.figsize((9,9))
    hbt.fontsize(12)

    plt.subplot(3,2,1)
    plt.plot(grid.bin_x_t)
    plt.ylabel('X Bin #')
    plt.title(f'MU69, Trajectory = {name_trajectory}, frame = {frame}')
    
    plt.subplot(3,2,3)
    plt.plot(grid.bin_y_t)
    plt.ylabel('Y Bin #')

    plt.subplot(3,2,5)
    plt.plot(grid.bin_z_t)
    plt.ylabel('Z Bin #')
    plt.xlabel('Timestep #')

    t_t = grid.et_t - np.mean(grid.et_t)
    bin_t = range(len(t_t))
    plt.subplot(3,2,2)
    plt.axhline(0, color='pink')
    plt.axvline(0, color='pink')
    plt.plot(t_t, grid.x_t)

    plt.ylabel('X [km]')

    plt.xlabel('t [sec]')
    
    plt.subplot(3,2,4)
    plt.axhline(0, color='pink')
    plt.axvline(0, color='pink')
    plt.plot(t_t, grid.y_t)
    plt.ylabel('Y [km]')

    plt.subplot(3,2,6)
    plt.axhline(0, color='pink')
    plt.axvline(0, color='pink')
    plt.plot(t_t, grid.z_t)
    plt.ylabel('Z [km]')
    plt.xlabel('Time from C/A [sec]')
    plt.tight_layout()

    # Save the plot to a file
    
    file_out = f'positions_trajectory_{name_trajectory}.png'
    path_out = os.path.join(dir_out, file_out)
     
    plt.savefig(path_out)
    print(f'Wrote: {path_out}')
    plt.show()
    
    # Make a table
    
    arr = {'bin' : bin_t, 
           'delta_et' : t_t,
           'X_km' : grid.x_t,
           'Y_km' : grid.y_t,
           'Z_km' : grid.z_t,
           'Bin_X' : grid.bin_x_t,
           'Bin_Y' : grid.bin_y_t,
           'Bin_Z' : grid.bin_z_t}
   
    t = Table(arr, names=['bin', 'delta_et', 'X_km', 'Y_km', 'Z_km', 'Bin_X', 'Bin_Y', 'Bin_Z'],
              dtype=['int', 'int', 'float', 'float', 'float', 'int', 'int', 'int'])
  
    # Save the table to a file
    
    file_out = f'positions_trajectory_{name_trajectory}.txt'
    path_out = os.path.join(dir_out, file_out)
    
    t.write(path_out, format = 'ascii.csv', overwrite=True)
    print(f'Wrote: {path_out}')
         
# =============================================================================
# Call the main function when this program as run
# =============================================================================
    
if (__name__ == '__main__'):
    
    name_trajectory = 'prime'  # Can be 'prime' or 'alternate'

    dir_in  = '/Users/throop/data/ORT4/throop/ort4_bc3_10cbr2_dph/'
#    dir_in  = '/Users/throop/data/ORT4/throop/ort4_bc3_10cbr2_dek/'
    
    dir_out = os.path.join(dir_in, 'for_mehoke')
    
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
        
    nh_ort_track4_flyby(dir_in=dir_in, dir_out=dir_out, name_trajectory = name_trajectory)
    