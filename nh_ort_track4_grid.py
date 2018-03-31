#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 00:05:08 2018

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

from nh_ort_track3_plot_trajectory import nh_ort_track3_plot_trajectory
from nh_ort_track3_read            import nh_ort_track3_read
from nh_ort_track3_read            import stretch_hbt, stretch_hbt_invert  # Some stretch routines good for trajectories
from nh_ort_track3_read            import plot_flattened_grids_table


# =============================================================================
# Define a class to handle all of the details of the Track-4 grids.
# This class does not create the grids, but it does read, write, plot, and fly thru.
# The ultimate goal of this class is to create the time-dependent particle densities that
# get given to Doug Mehoke for Track 5.
#     
# =============================================================================
    
class nh_ort_track4_grid:
    
    def __init__(self, grid_4d):
    
        """
        Save the grid to the class, and do any initialization.
        """
        
        self.grid_4d = grid_4d
        
        self.num_grids = hbt.sizex(grid_4d)
        
        self.name_trajectory = 'primary'
        
        self.name_test       = 'ort2-ring' 
        
        self.axis_sum = 0   # Which axis do we sum along for images? Should be as visible from Sun and/or SC.
        
# =============================================================================
# Set all ring parameters based on passed-in values
# =============================================================================
    
    def set_parameters(self, 
                             albedo=None, speed=None, q=None, rho=None,   # These are the physical params,
                                                                          # and they go into the filename and DPH sims.
                             b=None, beta=None, s=None,                   # These are the particle sizes we sum over.
                                                                          # They are used for plotting, but not else. 
                             name_trajectory=None, 
                             name_test=None):
    
        """
        Save the parameters corresponding to a run.

        parameters:
            albedo, speed, q, rho    -- These are for the runs. One per 4D grid.
            b, beta, s    -- These are for particle size, one per grid.
            
            *** Note that for Python logic to work, the  arrays b, beta, s must be passed as lists, 
                 and *not* NumPy arrays. Testing against None works on lists, but chokes on NumPy arrays.
        """
                                       
        if albedo:
            self.albedo = albedo

        if speed:
            self.speed = speed
                                                                                      
        if q:
            self.q = q

        if speed:
            self.rho = rho
            
        # Now process the keywords related to particle size. These should each have the same length as .num_grids.
        
        if b:
            self.b = b

        if beta:
            self.beta = beta
                                                                                      
        if s:
            self.s = s

        if name_test:
            self.beta = beta
                                                                                      
        if name_test:
            self.s = s

# =============================================================================
# Make plots of flattened grids, one plot for each particle size
# =============================================================================
                                                         
    def plot(self, stretch_percent=98):
        
        """
        Make a plot of all of the sub-grids in a grid. Each one is flattened, and they are all ganged up.
        """
        
        stretch = astropy.visualization.PercentileInterval(stretch_percent)
    
        hbt.figsize((15,15))
        for i in range(self.num_grids):
            b_i     = np.amax(self.b) - i     # b goes in opposite order from size, so plot it backwards.
            beta_i  = self.beta[i]              # Calc beta so we can get size. MRS slide 5.5
            s_i     = self.s[i]                 # Particle size. MRS slide 6.4. Millimeters.
   
            plt.subplot(1, self.num_grids, i+1)
            plt.imshow(stretch(np.sum(self.grid_4d[i], axis=self.axis_sum)))
            plt.title(r's={:.2f} mm, $\beta$={:.4f}'.format(s_i, beta_i))
            plt.tight_layout()
        plt.show()
        print(f'  albedo={self.albedo}, q={self.q}, rho={self.rho}, speed={self.speed}')


# =============================================================================
# Create the output filename automatically, based on the parameter values
# =============================================================================

    def create_filename(self):
               
        str_traj = self.name_trajectory
        
        str_test = 'ort2-ring'
        
        str_speed = 'v2.2'
        
        str_q = 'q{:3.1f}'.format(self.q)
        
        str_albedo = 'pv{:4.2f}'.format(self.albedo)
        
        str_rho = 'rho{:4.2f}'.format(self.rho)
        
#        str_inc = 'inc{:4.2f}'.format(self.inclination)
        
        file_out = f"{str_test}_{str_traj}_{str_speed}_{str_q}_{str_albedo}_{str_rho}.grid4d"
                        
        return file_out

# =============================================================================
# Write out the 4D grid to disk
# =============================================================================
    
    def write(self, file=None, dir=None):
        """
        Write the 4D grid itself to a file. The run parameters (albedo, rho, q, speed) are encoded into the filename.
        The remaining parameters (b, beta, s) are written.
        
        This file is only used for HBT's convenience. It is not given to Doug Mehoke nor is it the output of any 
        step.
        
        Format of the file is a pickle tuple, with (grid_4d, albedo, rho, q, speed, b, beta, s).
        
        ** Warning: These are really big files! 400+ MB by default. Will 'compress' down to 7 MB.
        
        So, in general, we should not plan on using these.
        
        Optional Parameters:
            file: Filename. By default, it is auto-generated from the parameters.
            
            dir:  Directory
        """
        
        if not dir:
            dir = os.path.expanduser('~/Data/ORT2/throop/track4/')
        
        if not file:
            file = self.create_filename()
            
        lun = open(os.path.join(dir, file), 'wb')
        
        pickle.dump((self.grid_4d, self.albedo, self.rho, self.q, self.speed, self.b, self.beta, self.s), lun)
        lun.close()
        print("Wrote: " + os.path.join(dir,file))
        
        return

# =============================================================================
# Read a 4D grid file from disk
# =============================================================================
    
    def read(self, file, dir=None):
        """
        Read a 4D grid file from disk

        Format of the file is a pickle tuple, with (grid_4d, albedo, rho, q, v, b, beta, s).
        
        Optional Parameters:
            dir: Directory
            
        """
        
        if not dir:
            dir = '~/Data/ORT2/throop/track4/'

        print("Reading: " + file)
        lun = open(os.path.join(dir, file), 'rb')
        (self.stack_field, self.img_field, self.stack_haz, self.img_haz, self.wcs_haz) = pickle.load(lun)
        lun.close()

        return

# =============================================================================
# Fly a trajectory through the grids
# =============================================================================
        
    def fly_trajectory(self):
        
        pass 
    
        return