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


    
def nh_ort_track4_flyby():

#    name_trajectory = 'alternate'  # Can be 'prime' or 'alternate'
    
    name_trajectory = 'prime'  # Can be 'prime' or 'alternate'

    stretch_percent = 99
    stretch = astropy.visualization.PercentileInterval(stretch_percent)

    dir = '~/Data/ORT2/throop/track4'
    
    files = glob.glob(os.path.join(os.path.expanduser(dir), '*.grid4d.gz'))
    
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
    
    # Start up SPICE if needed. If we change the kernel file, we need to restart python.
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh(f'kernels_kem_{name_trajectory}.tm')
    
    # Define the number of grid locations on each side.
    # In the 'simulated' version, we used 200, not 201.

    do_short = False
    
    if do_short:
        files = files[0:5]

    iof_ring = 2e-8
    
    file = files[54]  # Just for testing. Can ignore these.
    i=0
    
    for i,file in enumerate(files):
        print(f'Starting file {i}/{len(files)}')
              
        grid = nh_ort_track4_grid(file)    # Load the grid from disk. Uses gzip, so it is quite slow (10 sec/file)
     
        do_adjust_iof = False
        
        if do_adjust_iof:
            grid.calc_tau()
            print(f'XXX WARNING: GRID I/F_TYP = {grid.iof_typ:.2e} !!!')       
            factor = iof_ring / grid.iof_typ   
            grid.density *= factor               # Change the grid density to match the intended density
            grid.calc_tau()
            
            print(f'XXX WARNING: ADJUSTED GRID DENSITY TO MATCH I/F = {iof_ring} !!!')
    
    # Load the trajectory parameters

        et_ca = int( sp.utc2et(utc_ca) )  # Force this to be an integer, just to make output cleaner.
        
        et_start = et_ca - dt_before.to('s').value
        et_end   = et_ca + dt_after.to('s').value
                        
        grid.frame       = frame
        grid.name_target = name_target
        
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

        hbt.figsize((18,5))

        # Make a plot of the actual density that we give to Doug Mehoke
        
        # Define a list of colors. This is so we can use colors= argument to set
        # a marker to show grain size, rather than let plot() auto-assign.
        
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # This is the default color iterator.
        
        plt.subplot(1,3,1)
        for j,s in enumerate(grid.s):
            plt.plot(grid.delta_et_t, grid.number_t[j],
                     label = 's={:.2f} mm'.format(s))
        plt.legend()
        plt.title('Dust number density'.format(grid.area_sc))
        plt.xlabel('ET from C/A')
        plt.yscale('log')
        plt.ylabel(r'Dust,, # km$^{-3}$')

        plt.subplot(1,3,2)
        for j,s in enumerate(grid.s):   # 's' is dust size
            plt.plot(grid.delta_et_t, grid.number_sc_t[j],
                     label = f's={s:.2f} mm')
        plt.title('Impact rate, A={}'.format(grid.area_sc))
        plt.yscale('log')    
        plt.xlabel('ET from C/A')
        plt.legend()
        plt.ylabel(r'Dust, # Impacts sec$^{{-1}}$')

        # Make a plot of the cumulative count rate. Mark grain sizes here too.
        
        plt.subplot(1,3,3)
        for j,s in enumerate(grid.s):
            plt.plot(grid.delta_et_t, grid.number_sc_cum_t[j],                    # Main plot line
                     label = 's={:.2f} mm'.format(s), color=colors[j])
            plt.plot([grid.delta_et_t[-1]], [grid.number_sc_cum_t[j,-1].value],   # Circle to indicate grain size
                     markersize=(7-j)*2, marker = 'o',                            # Use same color as prev line! 
                     color=colors[j])

        plt.legend()
        plt.title('Number of impacts (cumulative), A={}'.format(grid.area_sc))
        plt.xlabel('ET from C/A')
        plt.yscale('log')
        plt.ylabel('# of Impacts')
        plt.axhline(y = 1, linestyle = '--', alpha = 0.1)    

        plt.tight_layout()
        
        plt.show()        
        
        # Now add an entry to the table. This is a table that lists all of the results -- e.g., max_tau, count rate etc
        # one line per grid

        t.add_row(vals=[grid.name_trajectory, grid.speed, grid.q, grid.albedo, grid.rho, 
                        grid.tau_max, grid.tau_typ,
                        grid.iof_max, grid.iof_typ])
                        

        # Output the dust population for this run to a file. This is the file that Doug Mehoke will read.
        
        grid.output_trajectory(suffix=f'{name_trajectory}', do_positions=False)
    
        print('---')
            
    # Print the table
    
    t.pprint(max_width=-1)
    
    # Save the table as output
    
    file_out = os.path.join(os.path.expanduser(dir), 'nh_ort_track4_table.pkl')
    
    lun = open(file_out, 'wb')
    pickle.dump(t,lun)
    lun.close()
    print(f'Wrote: {file_out}')

# =============================================================================
# Plot some tabulated results.
# One-off function.
# =============================================================================

def plot_table():
    
    file_in = '/Users/throop/Data/ORT2/throop/track4/nh_ort_track4_table.pkl'
    
    lun = open(file_in, 'rb')
    t = pickle.load(lun)
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
# Call the main function when this program as run
# =============================================================================
    
if (__name__ == '__main__'):
    nh_ort_track4_flyby()
    