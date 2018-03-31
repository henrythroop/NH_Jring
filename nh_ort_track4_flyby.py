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
#from nh_ort_track3_read            import stretch_hbt, stretch_hbt_invert  # Some stretch routines good for trajectories
#from nh_ort_track3_read            import plot_flattened_grids_table

# =============================================================================
# Main function to fly s/c through the grids. All of the 4D grids must already 
# be generated and saved to disk. This routine loops over them, and creates the 
# output files for Doug Mehoke, of particle density vs. time.   
#
# This is a regular function, but it calls the class method nh_ort_track4_grid.fly_trajectory().
#      
# =============================================================================
    
def nh_ort_track4_flyby():

#    name_trajectory = 'alternate'  # Can be 'prime' or 'alternate'
    name_trajectory = 'prime'  # Can be 'prime' or 'alternate'

    dir = '~/Data/ORT2/throop/track4'
    
    files = glob.glob(os.path.join(os.path.expanduser(dir), '*.grid4d.gz'))
    
    utc_ca = '2019 1 Jan 05:33:00'
    dt_before = 1*u.hour
    dt_after  = 1*u.hour
    
#    area_sc = (1*u.m)**2

    frame = '2014_MU69_SUNFLOWER_ROT'
    
    name_target = 'MU69'
    
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

    for file in files:
        
        grid_i = nh_ort_track4_grid(file)    # Load the grid from disk. Uses gzip, so it is quite slow (20 sec/file)
     
        n_dx = hbt.sizex(grid_i.density[0])  # This is 'multiple assignment' in Python. Unlike IDL, it is official and legal.
        n_dy = hbt.sizey(grid_i.density[0]) 
        n_dz = hbt.sizez(grid_i.density[0]) 
    
    # Define the size of each grid box
    # Wiki said 25 km, but should be 250.
    # XYZ sizes are tracked separately, but really should be identical.
    
        dx_bin = grid_i.resolution_km[0]
        dy_bin = grid_i.resolution_km[1]
        dz_bin = grid_i.resolution_km[2]
        
    # Load the trajectory

        et_ca = int( sp.utc2et(utc_ca) )  # Force this to be an integer, just to make output cleaner.
        
        et_start = et_ca - dt_before.to('s').value
        et_end   = et_ca + dt_after.to('s').value
                
        do_test = False
        i       = 0            # Counter of index in which to store results
        
        # And fly through it
        
        grid_i.fly_trajectory(name_observer, et_start, et_end, dt)
        
        # Make a few diagnostics plots of our path through the system
        
        do_plots_geometry = True
        
        if do_plots_geometry:
            plt.subplot(1,3,1)
            plt.plot(ring.delta_et_t, ring.radius_t)
            plt.ylim((0,np.amax(ring.radius_t)))
            plt.title('Radius')
            plt.xlabel('dt from CA [sec]')
            
            plt.subplot(1,3,2)
            plt.plot(ring.delta_et_t, ring.lat_t)
            plt.title('Lat')
            
            plt.subplot(1,3,3)
            plt.plot(ring.delta_et_t, ring.lon_t)
            plt.title('Lon')
            
            plt.tight_layout()
            
            plt.show()

        # Make a plot of the instantaneous count rate

        plt.plot(ring.delta_et_t, ring.number_sc_t)
        plt.title('Number of Impacts per sec, A={}, i={}'.format(area_sc, inclination_i))
        for i,r_dust_i in enumerate(self.r_dust):
            plt.plot(ring.delta_et_t, ring.number_sc_t * ring.n_dust[i],
                     label = 'r={}'.format(r_dust_i))
        plt.yscale('log')    
        plt.xlabel('ET')
        plt.legend()
        plt.ylabel('# of Impacts per sec')
        plt.show()

        # Make a plot of the cumulative count rate
        
        plt.plot(ring.delta_et_t, ring.number_sc_cum_t)
        for i,r_dust_i in enumerate(self.r_dust):
            plt.plot(ring.delta_et_t, ring.number_sc_cum_t * ring.n_dust[i],
                     label = 'r={}'.format(r_dust_i))
        plt.legend()    
        plt.title('Number of Impacts (cumulative), A={}, i={}'.format(area_sc, inclination_i))
        plt.xlabel('ET')
        plt.yscale('log')
        plt.ylabel('# of Impacts')
        plt.axhline(y = 1, linestyle = '--', alpha = 0.1)    
        plt.show()
        
        # Calculate the radial profiles explicitly, so I can extract the I/F, etc at the flyby distance,
        # and save in a table as output.
                            
        (radius_ring_rc, IoF_rc, tau_rc, n_hits_rc, pr_lom_rc) = self.get_radial_profile(r_min=r_crit)
        (radius_ring,    IoF,    tau,    n_hits,    pr_lom)    = self.get_radial_profile()

        # Calculate the proper radial bin to extract
        # We want this to be the geometrically closest bin. 
        
        bin_a_flyby = np.digitize(ring.a_flyby[name_trajectory], radius_ring)
        
        # Now add an entry to the table

        t.add_row(vals=[name_trajectory, q_dust_i, inclination_i, albedo_i, r_crit,
                        IoF_rc[bin_a_flyby], IoF[bin_a_flyby], 
                        tau_rc[bin_a_flyby], 
                        n_hits_rc[bin_a_flyby], n_hits[bin_a_flyby],
                        pr_lom_rc[bin_a_flyby]])

        # Output the dust population for this run to a file. This is the file that Doug Mehoke will read.
        
        ring.output_trajectory(suffix=f'{name_trajectory}', do_positions=False)
                
    # Print the table
    
    t.pprint(max_width=-1)
    