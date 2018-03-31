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

#import cPickle
import gzip

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
    
    def __init__(self, arg):
    
        """
        Save the grid to the class, and do any initialization.
        
        Takes one input. It can be either of:
            
            density: A 4D grid: (num_particle_sizes, x, y, z). Typically (7, 200, 200, 200).
            
        or
          
            file: A filename, in `.grid4d.gz` format.
        
        """

        
        if (type(arg) == str):
            self.read(arg)

        else:
            self.density = arg
            
        self.num_grids = hbt.sizex(self.density)
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
                             name_test=None,
                             resolution_km=None):
    
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

        if rho:
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
            
        if resolution_km:
            self.resolution_km = resolution_km

# =============================================================================
# Make plots of flattened grids, one plot for each particle size
# =============================================================================
                                                         
    def plot(self, stretch_percent=98):
        
        """
        Make a plot of all of the sub-grids in a grid. Each one is flattened, and they are all ganged up.
        """

        (_, nx, ny, nz) = np.shape(self.density)
        
        xpos = hbt.frange(-nx/2, nx/2)*self.resolution_km[0]
        ypos = hbt.frange(-ny/2, ny/2)*self.resolution_km[1]
        zpos = hbt.frange(-nz/2, nz/2)*self.resolution_km[2]

        extent = (ypos[0], ypos[-1], zpos[0], zpos[-1])
        
        stretch = astropy.visualization.PercentileInterval(stretch_percent)
    
#        hbt.figsize((15,15))
        for i in range(self.num_grids):
            b_i     = np.amax(self.b) - i     # b goes in opposite order from size, so plot it backwards.
            beta_i  = self.beta[i]              # Calc beta so we can get size. MRS slide 5.5
            s_i     = self.s[i]                 # Particle size. MRS slide 6.4. Millimeters.
   
            plt.subplot(1, self.num_grids, i+1)
            plt.imshow(stretch(np.sum(self.density[i], axis=self.axis_sum)), extent=extent)
            plt.title(r's={:.2f} mm, $\beta$={:.4f}'.format(s_i, beta_i))
            plt.tight_layout()
            if (self.axis_sum == 0):
                plt.xlabel('Z [km]')
                plt.ylabel('Y [km]')
        plt.show()
        print(f'  albedo={self.albedo}, q={self.q}, rho={self.rho}, speed={self.speed}')


# =============================================================================
# Create the output filename automatically, based on the parameter values
# =============================================================================

    def create_filename(self):
               
        str_traj = self.name_trajectory
        
        str_test = 'ort2-ring'
        
        str_speed = 'v{:3.1f}'.format(self.speed)
        
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
        
        Format of the file is a pickle tuple, with fields:
              (density, albedo, rho, q, speed, b, beta, s, (resolution_km)).
              
        'density' is the large 4D array, with format (num_grids, x, y, z) -- typically (7, 200, 200, 200).      
        
        ** Warning: These are really big files! 400+ MB by default. Will 'compress' down to 7 MB.
        
        https://stackoverflow.com/questions/18474791/decreasing-the-size-of-cpickle-objects
        
        So, in general, we should not plan on using these.
        
        Optional Parameters:
            file: Filename. By default, it is auto-generated from the parameters.
            
            dir:  Directory
        """
        
        if not dir:
            dir = os.path.expanduser('~/Data/ORT2/throop/track4/')
        
        if not file:
            file = self.create_filename()
                
        do_compress = True
        
        if do_compress:
            
            file = file + '.gz'
            
            protocol=-1    
            with gzip.open(os.path.join(dir,file), 'wb') as f:
                pickle.dump((self.density, 
                             self.albedo, self.rho, self.q, self.speed, 
                             self.b, self.beta, self.s, 
                             self.resolution_km), 
                             f, 
                             protocol)

        else:
            lun = open(os.path.join(dir, file), 'wb')        
            pickle.dump((self.density, 
                             self.albedo, self.rho, self.q, self.speed, 
                             self.b, self.beta, self.s, 
                             self.resolution_km), lun)
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
            dir = os.path.expanduser('~/Data/ORT2/throop/track4/')

        print("Reading: " + file)

        with gzip.open(os.path.join(dir, file), 'rb') as f:
             (self.density, 
                             self.albedo, self.rho, self.q, self.speed, 
                             self.b, self.beta, self.s, 
                             self.resolution_km) = pickle.load(f)
        f.close()

        return        
    
# =============================================================================
# Fly a trajectory through the grids and sample it
# =============================================================================
        
    def fly_trajectory(self, name_observer, et_start, et_end, dt):

        """
        Sample the ring along a flight path.
        
        Parameters
        ----
        
        et_start:
            
        et_end:
        
        dt:
            
        """
    
        # Save the passed-in parameters in case we need them 
        
        self.et_start = et_start
        self.et_end   = et_end
        self.dt       = dt
        
        self.name_observer = name_observer
        
        # Set up the output time array
        
        num = math.ceil( (et_end - et_start) / dt.to('s').value ) + 1
        
        et_t = hbt.frange(int(et_start), int(et_end), num)
        
        # Calc offset in DT from C/A time
            
        delta_et_t = et_t - np.mean(et_t)
        
        # Loop over et
        
        radius_t = []
        x_t      = []
        y_t      = []
        z_t      = []
        v_t      = [] # Velocity
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
            
            v_t.append( sp.vnorm(st[3:6]) )
                
        # Convert into bin values. This is vectorized.
        # ** If values exceed the largest bin, then the index returned will be too large for the density lookup!
        
        bin_x_t = np.digitize(x_t, self.x_1d)
        bin_y_t = np.digitize(y_t, self.y_1d)
        bin_z_t = np.digitize(z_t, self.z_1d)
        
        v_t = np.array(v_t) * u.km/u.s
        
        # If indices are too large, drop them by one. This just handles the edge cases so the code doesn't crash.
        
        bin_x_t[bin_x_t >= len(self.x_1d)] = len(self.x_1d)-1
        bin_y_t[bin_y_t >= len(self.y_1d)] = len(self.y_1d)-1
        bin_z_t[bin_z_t >= len(self.z_1d)] = len(self.z_1d)-1
        
        # Now that we have all the XYZ bins that the s/c travels in, get the density for each one.
        # We should be able to do this vectorized -- but can't figure it out, so using a loop.
        
        number_t = 0. * et_t
        
        for i in range(len(et_t)):
            number_t[i] = self.number_arr[bin_x_t[i], bin_y_t[i], bin_z_t[i]]  # Number per km3

        # ** This gives us number density at the smallest binsize. We then need to apply n(r), to get 
        # density at all other bin sizes.
        # Make a cumulative sum of the density. We don't really need this, but just for fun.
        
        number_cum_t       = np.cumsum(number_t)
    
        # Now for fun, calculate the number of grains that intercept a s/c of a given area.
        # We do this calc very crudely, just by taking the area of s/c vs area of a bin edge. We are ignoring
        # the slant path of the s/c, so we could underestimate the flux by up to sqrt(3).
        
        # Calc the fraction of the bin volume that the s/c sweeps up during its passage. 
        
        # Get the binvolume, in km3, as an integer
        
        binvol = (self.binsize_3d[0] * self.binsize_3d[1] * self.binsize_3d[2])
    
        # Calc the fractional ratio between volume of a bin, and volume that s/c sweeps up in its path thru the bin.
        
        fracvol_t = ( (self.area_sc * v_t * self.dt) / (binvol) ).to(1).value
    
        number_sc_t = number_t * fracvol_t
        
        number_sc_cum_t = np.cumsum(number_sc_t)
            
        # Save all variables so they can be retrieved
        
        self.radius_t      = radius_t
        self.et_t          = et_t
        self.delta_et_t    = delta_et_t
        self.bin_x_t       = bin_x_t
        self.bin_y_t       = bin_y_t
        self.bin_z_t       = bin_z_t
        self.number_t      = number_t
        self.number_sc_t   = number_sc_t
        self.number_sc_cum_t=number_sc_cum_t
        self.number_cum_t  = number_cum_t
        self.lon_t         = lon_t
        self.lat_t         = lat_t
        self.x_t           = x_t
        self.y_t           = y_t
        self.z_t           = z_t
        self.v_t           = v_t
        
        self.radius_t = radius_t  # This is distance from the center, in 3D coords. Not 'ring radius.'
        
#        return()

        
        pass 
    
        return
    
# =============================================================================
# Output the trajectory and particle intercept info to a file
# =============================================================================
        
    def output_trajectory(self, suffix=None, do_positions=True):

        """
        Write an output table. The filename is auto-generated based on the parameters already supplied.
        
        The output value is typically in particles per km3.
          ** That is my memory from MRS. But, his email 5-Dec-2017 says 
             "Once the numbers are calibrated to particles per size bin per cubic meter, we tabulate those numbers 
             vs. time along each trajectory and hand them off to Doug M"
          ** Wiki says "Units are TBD"
          ** MRS's code nhdust.py (?) uses #/km3 at one point.
             
             ** For now I will assume # per km3, but this can be trivially changed later.
          
        Parameters
        -----
        do_positions:
            If True, output three columns which describe the XYZ position of the spacecraft as function of time.
        str:
            A string to be inserted into the filename, as part of the filename.
        """

        # Choose the method of table writing. Astropy is easier to code, but I can't use my own header to match MRS's.
        
        do_write_astropy = False
        
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
        
        # Get the binvolume, in km3, as an integer
        
        binvol_km3 = (self.binsize_3d[0] * self.binsize_3d[1] * self.binsize_3d[2]).to('km^3').value
        
        # For each bin, output the # of particles, divided by the bin volume, to get particles/km3.
        
        for i in range(len(self.n_dust)):
                t['n_{}, {:.3f}, # per km3'.format(i, self.r_dust[i])] = \
                                       np.array(self.number_t * self.n_dust[i] / binvol_km3).astype(int)
        
        # Create the output filename
        
        dir_out = '/Users/throop/Data/ORT2/throop/track4/'
        
        file_out = self.create_filename_track4()
        
#        file_out = 'ort1_q{}_i{}_a{}.txt'.format(self.q_dust, self.inclination, self.albedo)
#        if do_positions:
#            file_out = file_out.replace('.txt', '_pos.txt')
#       
        path_out = os.path.join(dir_out, file_out)
        
        if do_write_astropy:

            # Write the file using Astropy
            
            t.write(path_out, format = 'ascii', overwrite=True)
            print("Wrote: {} using astropy writer".format(path_out))

        else:

            # Write the file using explicit manual formatting
            # There are seven size radial bins, as per MRS e-mail 2-Feb-2018.
            
            lun = open(path_out, "w")
            lun.write("#    First line is '0' and then size bins, in mm\n")
            lun.write("#    Remaining are delta_ET, n(r_1), n(r_2), n(r_3), n(r_4), n(r_5), n(r_6), n(r_7)\n")
            lun.write("#    n(r) are number per km3\n")
            lun.write("#    Henry Throop {}\n".format(str(datetime.now())))          
            lun.write("{} {} {} {} {} {} {} {}\n".format(0, 
                                                    self.r_dust[0].to('mm').value,
                                                    self.r_dust[1].to('mm').value,
                                                    self.r_dust[2].to('mm').value,
                                                    self.r_dust[3].to('mm').value,
                                                    self.r_dust[4].to('mm').value,
                                                    self.r_dust[5].to('mm').value,
                                                    self.r_dust[6].to('mm').value))
            for i in range(len(t)):                                                
                lun.write("{} {} {} {} {} {} {} {}\n".format(
                        t[i][1],
                        t[i][2],
                        t[i][3],
                        t[i][4],
                        t[i][5],
                        t[i][6],
                        t[i][7],
                        t[i][8]))
            lun.close()
            print("Wrote: {} using manual writer".format(path_out))
           
         
        return
    