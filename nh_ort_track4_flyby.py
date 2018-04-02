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
    
    t = Table(names = ['trajectory', 'q_dust', 'inclination', 'albedo', 'r_c [mm]', 
                       'IoF > r_c', 'IoF full', 'tau > r_c', 'n_hits > r_c', 'n_hits full', f'PrLOM > r_c'],
              dtype = ['U30', float, float, float, float, float, float, float, float, float, float]  )
    
    # Start up SPICE if needed. If we change the kernel file, we need to restart python.
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh(f'kernels_kem_{name_trajectory}.tm')
    
    # Define the number of grid locations on each side.
    # In the 'simulated' version, we used 200, not 201.

    for i,file in enumerate(files):
        print(f'Starting file {i}/{len(files)}')
              
        grid = nh_ort_track4_grid(file)    # Load the grid from disk. Uses gzip, so it is quite slow (10 sec/file)
     
        n_dx = hbt.sizex(grid.density[0])
        n_dy = hbt.sizey(grid.density[0]) 
        n_dz = hbt.sizez(grid.density[0]) 
    
    # Define the size of each grid box
    # Wiki said 25 km, but should be 250.
    # XYZ sizes are tracked separately, but really should be identical.
    
        dx_bin = grid.resolution_km[0]
        dy_bin = grid.resolution_km[1]
        dz_bin = grid.resolution_km[2]
        
    # Load the trajectory parameters

        et_ca = int( sp.utc2et(utc_ca) )  # Force this to be an integer, just to make output cleaner.
        
        et_start = et_ca - dt_before.to('s').value
        et_end   = et_ca + dt_after.to('s').value
                
        do_test = False
        
        grid.frame       = frame
        grid.name_target = name_target
        
        # And call the method to fly through it!
        # The returned density values etc are available within the instance variables, not returned explicitly.

        grid.fly_trajectory(name_observer, et_start, et_end, dt)
 
        # Make a few diagnostics plots of our path through the system
        # Only do this on the first time thru the loop.
        
        do_plots_geometry = True
        
        if (do_plots_geometry and (i==0)):
            
            hbt.figsize((8,5))
            plt.subplot(2,3,1)
            plt.plot(grid.delta_et_t, np.array(grid.radius_t)/1000)
            plt.ylim((0,np.amax(grid.radius_t)/1000))
            plt.title('Radius')
            plt.xlabel('dt from CA [sec]')
            plt.ylabel('Radius [Mm]')
            
            plt.subplot(2,3,2)
            plt.plot(grid.delta_et_t, np.array(grid.lat_t) * hbt.r2d)
            plt.title('Lat [deg]')
            
            plt.subplot(2,3,3)
            plt.plot(grid.delta_et_t, np.array(grid.lon_t) * hbt.r2d)
            plt.title('Lon [deg]')
            
            plt.subplot(2,3,4)
            plt.plot(grid.delta_et_t, np.array(grid.x_t)/1000)
            plt.title('X')
            plt.xlabel('dt from CA [sec]')
            plt.ylabel('Radius [Mm]')
            
            plt.subplot(2,3,5)
            plt.plot(grid.delta_et_t, np.array(grid.y_t)/1000)
            plt.title('Y')
            
            plt.subplot(2,3,6)
            plt.plot(grid.delta_et_t, np.array(grid.z_t)/1000)
            plt.title('Z')
            
            plt.tight_layout()
                        
            plt.show()

        # Make slice plots thru the grid
        
        hbt.fontsize(8)
        hbt.figsize((20,5))
        grid.plot(axis_sum=0)
        grid.plot(axis_sum=1)
        grid.plot(axis_sum=2)
        
        hbt.fontsize(10)
        
        # Make a plot of optical depth
        
        grid.plot_tau()
                
        # Make a plot of the instantaneous count rate

        hbt.figsize((12,5))

        plt.subplot(1,2,1)
        for j,s in enumerate(grid.s):   # 's' is dust size
            plt.plot(grid.delta_et_t, grid.number_sc_t[j],
                     label = f's={s:.2f} mm')
        plt.title('Number of Impacts per sec, A={}'.format(grid.area_sc))
        plt.yscale('log')    
        plt.xlabel('ET')
        plt.legend()
        plt.ylabel('# of Impacts per sec')

        # Make a plot of the cumulative count rate
        
        plt.subplot(1,2,2)
        for j,s in enumerate(grid.s):
            plt.plot(grid.delta_et_t, grid.number_sc_cum_t[j],
                     label = 's={:.2f} mm'.format(s))
        plt.legend()
        plt.title('Number of Impacts (cumulative), A={}'.format(grid.area_sc))
        plt.xlabel('ET')
        plt.yscale('log')
        plt.ylabel('# of Impacts')
        plt.axhline(y = 1, linestyle = '--', alpha = 0.1)    
        plt.show()
        
        # Now add an entry to the table. This is a table that lists all of the results -- e.g., 
        # one line per grid

#        t.add_row(vals=[name_trajectory, q_dust_i, inclination_i, albedo_i, r_crit,
#                        IoF_rc[bin_a_flyby], IoF[bin_a_flyby], 
#                        tau_rc[bin_a_flyby], 
#                        n_hits_rc[bin_a_flyby], n_hits[bin_a_flyby],
#                        pr_lom_rc[bin_a_flyby]])

        # Output the dust population for this run to a file. This is the file that Doug Mehoke will read.
        
        grid.output_trajectory(suffix=f'{name_trajectory}', do_positions=False)
    
        print('---')
            
    # Print the table
    
    t.pprint(max_width=-1)


if (__name__ == '__main__'):
    nh_ort_track4_flyby()
    