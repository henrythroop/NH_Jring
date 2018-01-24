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
    
    def __init__(self, num_pts_3d, gridspacing, frame=None, name_target=None):

        """
        Initialize the grid.
        """
        
        # We want to set up this grid in the MU69 Sunflower frame.
    
        if frame:
            self.frame = frame
        if name_target:
            self.name_target = name_target
        
        # Define the grids
        
        self.density_arr = np.array(num_pts_3d) # Number of dust grains in this bin. Normalized somehow, not sure how.
        self.azimuth_arr = np.array(num_pts_3d) # Azmuth angle. Not sure we need this. In body frame.
        self.radius_arr  = np.array(num_pts_3d) # Radius. This will be useful. In body frame.
        
        self.x_arr       = np.array(num_pts_3d)
        self.y_arr       = np.array(num_pts_3d)
        self.z_arr       = np.array(num_pts_3d)
    
        n = num_pts_3d[0] # Number of grid points on each side
        
        # Define 1D arrays for xyz position
        
        x = hbt.frange(-n/2, n/2, n) * gridspacing[0]
        y = hbt.frange(-n/2, n/2, n) * gridspacing[1]
        z = hbt.frange(-n/2, n/2, n) * gridspacing[2]
    
        # Define 3D arrays for xyz position
        
        (x_arr, y_arr, z_arr) = np.meshgrid(x, y, z)


    def set_ring_parameters(self, file, albedo, q):
        
        """
        Define the ring radial profile, size distribution, albedo, size dist, etc.
        """
        
        # Read the radial profile
        
        t = Table.read(file)
        self.radius_km  = t['RadiusKM']
        self.IoF        = t['I/F']
        self.radius_pix = t['RadiusPixels']
        
        # Extrapolate the radial profile to the right bins, and set it.
        
        # Set the size distribution
        # We set one size distribution for the entire ring. It is uniform and does not change spatially or temporally.
        
        self.n_dust = {}
        self.r_dust = {}
        return(0)
        
        
    def fly_trajectory(self, name_observer, et_start, et_end, dt):
        """
        Now that all parameters are set, sample the ring along a flight path.
        """
        
        et_arr = hbt.frange(int(et_start), int(et_end))
        
        return(0)
        
# =============================================================================
# END OF METHOD DEFINITION
# =============================================================================
    
# =============================================================================
# The main function to call to run the simulation
# =============================================================================
    
def do_ring_flyby():
    
    albedo              = {0.05, 0.1, 0.3, 0.7}
    q_dust              = {2, 3.}
    inclination_ring    = {0., 0.1}
    
    file_profile_ring = '/Users/throop/Dropbox/Data/ORT1/throop/backplaned/K1LR_HAZ04/stack_n40_z4_profiles.txt'
    
 # This defines I/F as a function of orbital distance.
                                       # n(r) and q are then derived from I/F and albedo.
    
    utc_ca = '2019 1 Jan 05:33:00'
    dt_before = 1*u.hour
    dt_after  = 1*u.hour
    
    frame = '2014_MU69_SUNFLOWER_ROT'
    
    name_target = 'MU69'
    
    name_observer = 'New Horizons'
    
    dt_sample = 1         # Sampling time through the flyby. Assumed to be seconds.

    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem.tm')
    
    # Define the number of grid locations on each side.
    
    n_dx = n_dy = n_dz = 201   # This is 'multiple assignment' in Python. Unlike IDL, it is official and legal.
    
    # Define the size of each grid box
    
    dx_bin = dy_bin = dz_bin = 25  # Distances are km. Not worth tracking the units. Consistent with SPICE.
    
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
                
                # Set up the ring itself
                
                ring.set_ring_parameters(file=file_profile_ring, albedo=albedo_i, q = q_dust_i)

                # And fly through it
                
                out = ring.fly_trajectory(name_observer, et_start, et_end, dt_sample)    

                # Write the trajectory to a file, plot it, etc.
                    
# =============================================================================
# Run the function
# =============================================================================
                
if (__name__ == '__main__'):
    do_ring_flyby()
    