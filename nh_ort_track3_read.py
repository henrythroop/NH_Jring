#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:38:48 2018

This reads the NH MU69 ORT 'Track 3' data, which is the output of N-Body simulations done by
Doug Hamilton + David Kaufmann.

Incorporates code fragment from DK 
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

class nh_ort_track3_read:  # Q: Is this required to be same name as the file?

# =============================================================================
# Initialize the method.
# Read the specified track3 file into memory, from disk
# =============================================================================
    
    def __init__(self, name, dir_base='/Users/throop/data/ORT2/hamilton/deliveries', binfmt = 2):
        
        """
        Initialize the object. Read in data from a specified file, and do initial processing on it.
        
        Parameters
        ------
            
        name:
            The name of the run. This is a slash-separated string, such as 
            `'ort1-0001/1000010/speed1/beta_2.6e-04/subset00'`. These slashes are not the best way 
            to do this, and it's not very compact. But the names are inconsistent, and this is the 
            most reliable way to deal with it for now.
        
        Optional parameters
        ------
        
        binfmt:
            Type of file to read. 1 = old file (pure grid). 2 = new file (sparse matrix). For ORT1 and beyond, 
            the default is `binfmt=2`.
            
        dir_base:
            The top-level directory that I download files from ixion into. Typically this ends in `deliveries`.
        """    
        
        do_kaufmann  = ('kaufmann' in dir_base)
        do_hamilton  = ('hamilton' in dir_base) 
        
        self.name = name
    
        parts = name.split('/')

        dir = name
        
        # Construct the filename for the 'header', which is a short text file. It is basically the same
        # as the entire directory tree, but with '/' → '_'.
        # There is one difference -- probably a bug -- in that 'speed' term is off by one, as per DK email. OBOB.
      
        # Also, construct the filename for the data array itself. It is trivial.
        # As per the wiki, the data and header files should be identically named. They are not!
        # https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Integrations-to-Density
        
        if do_kaufmann:
            file_header = 'header.txt'
            file_data   = f'model.array{binfmt}'

        if do_hamilton:
            file_header = 'header.txt'
            file_data   = f'grid.array{binfmt}'
                
        file_header_full = os.path.join(dir_base, dir, file_header)
        file_data_full   = os.path.join(dir_base, dir, file_data)

        self.file_header_full = file_header_full
        self.file_data_full   = file_data_full
        
        # Read the datafiles. Using code fragment from DK email 13-Feb-2018.

        # Read header file
        
        f = open(file_header_full, 'r')
        text = f.read()
        f.close()
        lines = text.split('\n')
        
        self.header = lines
        
        nx = int(lines[0].split()[0])
        ny = int(lines[1].split()[0])
        nz = int(lines[2].split()[0])
        km_per_cell_x = float(lines[3].split()[0])
        km_per_cell_y = float(lines[4].split()[0])
        km_per_cell_z = float(lines[5].split()[0])
        
        # Read and process binary file
        
        density = np.zeros((nx, ny, nz), np.float32)
        f = open(file_data_full, 'rb')
        data = f.read()
        f.close()
                
        if binfmt == 1:
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        start = 4*(i*ny*nz + j*nz + k)
                        (density[k,j,i],) = struct.unpack('f', data[start:start+4])
        else: # binfmt == 2
            statinfo = os.stat(file_data_full)
            nentries = int(statinfo.st_size / 10) # record size = 10 bytes (3 2-byte ints + 4-byte float)
                                                  # int / int = float, in Python -- not sure why
                   
            for i in range(nentries):
                start = 10*i
                (ix,iy,iz) = struct.unpack('hhh', data[start:start+6])
                (density[iz-1,iy-1,ix-1],) = struct.unpack('f', data[start+6:start+10])
    
        self.density = density
        
        # Read a bunch of values from the header, and save them into the object
        
        self.km_per_cell_x = km_per_cell_x
        self.km_per_cell_y = km_per_cell_y
        self.km_per_cell_z = km_per_cell_z
        
        self.body_id   = self.get_header_val('body_id')
        self.duration  = self.get_header_val('duration')*u.year
        self.time_step = self.get_header_val('time_step')*u.s
        self.beta      = self.get_header_val('beta')
        self.speed     = self.get_header_val('speed')    # -3 or -2.2
        self.subset    = self.get_header_val('subset')
        self.grains    = self.get_header_val('grains')
        self.state_file= self.get_header_val('state_file')
        self.obj_file  = self.get_header_val('obj_file')
        
        # Do a calculation for particle radius
        
#        Beta = 5.7e-4 * Q_PR / (rho * s/(10mm))  ← s = radius in cm [sic]. rho = g/cm3.
#        Q_PR = 1 + 1.36 * p_v # ← p_v = albedo     From Showalter 3.6. A simple relationship. 

        
        # Below we list the possible values for p_v (albedo) and q_pr
        # These two below are linked and 1:1!
        
        p_v = np.array([0.7, 0.3, 0.1, 0.05])
        q_pr = np.array([1.95, 1.41, 1.14, 1.07])

        # This is a list of all the distinct allowable values of beta
        # This list is just taken by examination of DK's files.
        # Beta has 16 different values, 5.7e-8 .. 5.7e-3. Factor of 2.1 between them. Total range of 1e5.
        
        beta = [5.7e-8, 5.7e-7, 5.7e-6, 5.7e-5, 5.7e-4, 5.7e-3, 
                         2.6e-7, 2.6e-6, 2.6e-5, 2.5e-4, 2.6e-3, 
                         1.2e-7, 1.2e-6, 1.2e-5, 1.2e-4, 1.2e-3]

        
#        plt.plot(sorted(beta), marker = 'o')
#        plt.yscale('log')
#        plt.show()
        
        rho_dust = np.array([1, 0.46, 0.21, 0.1])*u.g/(u.cm)**3
        
        f = 5.7e-5 * u.g / (u.cm)**3
        
        r_dust = (5.7e-4 * q_pr) / (self.beta * rho_dust) * u.mm
        
        return  # Return the object so we can chain calls.

# =============================================================================
# Search through the header, and find the value corresponding to a given keyword
# =============================================================================

    def get_header_val(self, keyword):
        
        """ 
        Get a header value from a list of header lines.
        
        As per https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Integrations-to-Density,
        the header is supposed to be a fixed format, with one specific value per line. However, I am building
        this function more robustly, in case the specs change in the future.
        
        Parameters
        -----
        
        keyword:
            A string that matches the keyword for a specified entry in the header. `keyword` should be chosen so that
            it will not match multiple lines in the header.
        
        """
       
                                                 
        out = ''
        for line in self.header:                    # Loop over every line
            if keyword in line:
                out = line.split('#')[0].strip()
                if hbt.is_number(out):                  # Convert from string to float, if possible
                    out = float(out)
        return out
       
# =============================================================================
# Plot the distributions
# =============================================================================
    
    def plot(self):
        
#        stretch_percent = 90    
#        stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
 
        (nx, ny, nz) = np.shape(self.density)
        xpos = hbt.frange(-nx/2, nx/2)*self.km_per_cell_x
        ypos = hbt.frange(-ny/2, ny/2)*self.km_per_cell_y
        zpos = hbt.frange(-nz/2, nz/2)*self.km_per_cell_z

        extent = (ypos[0], ypos[-1], zpos[0], zpos[-1])
        
        for axis in (0,1,2):
            plt.subplot(2,2,axis+1)
            plt.imshow(stretch(np.sum(self.density, axis=axis)), extent=extent)
            plt.ylabel('Pos [km]')
            plt.xlabel('Pos [km]')
            if (axis == 2):
                plt.colorbar()
            plt.title(self.name)
        
        plt.subplot(2,2,4)
        plt.imshow((np.sum(self.density, axis=0)), extent=extent)
        plt.xlim((-10000,10000))
        plt.ylim((-10000,10000))
        
        plt.show()
        
        return self

# =============================================================================
# Print some info about the run
# =============================================================================

    def print_info(self):

        print(f'Name:       {self.name}')
        print(f'Beta:       {self.beta}')
        print(f'Duration:   {self.duration}')
        print(f'Speed:      {self.speed}')
        print(f'Timestep:   {self.time_step}')
        print(f'Body ID:    {self.body_id}')
        print(f'Subset:     {self.subset}')
        print(f'# grains:   {self.grains}')
        print(f'State file: {self.state_file}')
        print(f'Obj   file: {self.obj_file}')
        
        return

    def find_matching_indices(self, beta = None, speed = None, body_id = None, subset = None):
        
        mask_start = np.ones(self.num_)
        mask_beta  = (beta  == self.beta)
        mask_speed = (speed == self.speed)
# =============================================================================
# Make a plot showing all of the runs in the grid. 
# =============================================================================

def plot_flattened_grids(arr):
    
    """
    Plot all of the grids in an array of grids. This is just to give a quick 'family portrait' overview.
    They are plotted in sequential order, but no attempt is made to list by beta, etc.

    parameters
    -----
    
    arr: 
        3D array, with each plane to be plotted.

    """
    
#    stretch_percent = 98  # Since bg is mostly black, this requires a slightly different stretch than normal!
#    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
 
    num = hbt.sizex(arr) # Number of runs in the array
    num_grid_x = math.ceil(math.sqrt(num))
    num_grid_y = num_grid_x
    dx_pix_grid = hbt.sizey(arr)
    dy_pix_grid = dx_pix_grid
    
    img = np.zeros((num_grid_x*dx_pix_grid, num_grid_y*dy_pix_grid))
    i = 0  # Row
    j = 0  # Column
    for k in range(len(runs_full)):
        
        # For each one, copy the grid, normalize it, and put it into out
        
        img[i*dx_pix_grid:(i+1)*dx_pix_grid, j*dy_pix_grid:(j+1)*dy_pix_grid] = \
            hbt.normalize( density_flattened_arr[k,:,:] )                         
        i += 1
        
        if (i > (num_grid_x-1)):
            i = 0
            j += 1
            
    plt.imshow(stretch(img))
    plt.show()    
    
# =============================================================================
# Make a plot showing all of the runs in the grid. 
# *** This is the same as the function above -- just uses input in table, rather than array! ***
# =============================================================================

def plot_flattened_grids_table(t, stretch_percent=96):
    
    """
    Plot all of the grids in an array of grids. This is just to give a quick 'family portrait' overview.
    They are plotted in sequential order, but no attempt is made to list by beta, etc.
    
    parameters
    -----
    
    t: 
        Astropy table, with column 'img_2d' to be plotted.
    """
    
 
    num = hbt.sizex(t['img_2d']) # Number of runs in the array
    num_grid_x = math.ceil(math.sqrt(num))
    num_grid_y = num_grid_x
    dx_pix_grid = hbt.sizey(t['img_2d'][0])
    dy_pix_grid = dx_pix_grid
    
    img = np.zeros((num_grid_x*dx_pix_grid, num_grid_y*dy_pix_grid))
    i = 0  # Row
    j = 0  # Column
    for k in range(num):
        
        # For each one, copy the grid, normalize it, and put it into out
        
        img[i*dx_pix_grid:(i+1)*dx_pix_grid, j*dy_pix_grid:(j+1)*dy_pix_grid] = \
            hbt.normalize( t['img_2d'][k] )                         
        i += 1
        
        if (i > (num_grid_x-1)):
            i = 0
            j += 1
            
    plt.imshow(stretch(img))
    plt.show()    

# =============================================================================
# Define a good stretch to use. Same idea as the Astropy Visualization stretch, just better for these data.
# =============================================================================

def stretch(arr):
    return np.log(10 + arr)
        
# =============================================================================
# Run the file
# =============================================================================
    
if __name__ == '__main__':

    plt.set_cmap('plasma')
#    stretch_percent = 96
#    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
#    stretch = astropy.visualization.LinearStretch()

    pi = math.pi
    
# For Hamilton, sample is ~/Data/ORT2/hamilton/deliveries/sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/
#                          ort2-0003/y2.2/beta2.2e-01/subset04/grid.array2'
    
    dir_base='/Users/throop/data/ORT2/hamilton/deliveries'    
    runs_full = glob.glob(os.path.join(dir_base, '*', '*', '*', '*', '*', '*')) # Hamilton - ORT actual

    # Set the axis to sum along, based on who create the file. There is an inconsistency from DPH vs. DK in ORT2.
    # This will probably be fixed in ORT3.
    # If everyone uses Sunflower frame, then we should sum in the Y dir, which should be axis=1.
    
    if ('hamilton' in dir_base):
        axis_sum = 0
        
    do_short = False
    
    if do_short:
        runs_full = runs_full[7:9]
        
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
    
    run = runs_full[0].replace(dir_base, '')[1:]  # Populate one, just for cut & paste ease

    rings = []                                    # Make a list with all of the ring objects in it.
    
    hbt.figsize((15,15))
    
    # Loop over all of the files, and read each one in.
    
    for i,run_full in enumerate(runs_full):

        run  = run_full.replace(dir_base, '')[1:]  # Remove the base pathname from this, and initial '/'
        ring = nh_ort_track3_read(run)             # Read the data array itself
        ring.print_info()

        halfwidth_km = ring.km_per_cell_x * hbt.sizex(ring.density) / 2
        extent = [-halfwidth_km, halfwidth_km, -halfwidth_km, halfwidth_km]  # Make calibrated labels for X and Y axes
        
        beta_arr[i]                  = ring.beta  # Take the params from individual runs, and stuff into an array   
        time_step_arr[i]             = ring.time_step.value
        duration_arr[i]              = ring.duration.value
        speed_arr[i]                 = ring.speed
        body_id_arr[i]               = ring.body_id
        subset_arr[i]                = ring.subset
        grains_arr[i]                = ring.grains

        # Take the 3D arrays and flatten them, just so we can visualize more easily.
        # To get a full picture, we sum along all three axes (X, Y, Z) and show each.

        density_flattened_arr[i,:,:] = np.sum(ring.density, axis=axis_sum)
        
        do_plot_xyz_views = False
        
        if do_plot_xyz_views:
            plt.subplot(1,3,1)
            plt.imshow(stretch(np.sum(ring.density, axis=0)), extent=extent)
            plt.title('Summed along X')
            plt.subplot(1,3,2)
            plt.imshow(stretch(np.sum(ring.density, axis=1)), extent=extent)
            plt.title('Summed along Y')
            plt.subplot(1,3,3)
            plt.imshow(stretch(np.sum(ring.density, axis=2)), extent=extent)
            plt.title('Summed along Z')
            plt.show()
        
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
    #  Albedo:     No units
    #  ds:         Bin width. mm
    #  R:          Moon radius. km.     [Slide 6.4 inferred from E0.] [Also, 4.3 says that area is delivered in km2.]
    
    albedo = [0.05, 0.1, 0.3, 0.7]   # Q: Are albedo and beta independent? A: Yes, completely indep -- see eq @ 3.7. 
    q      = [-2, -3.5]              # Should be negative, since exponent in eq @ 6.4 is positive
    rho    = [1, 0.46, 0.21, 0.10]
    speed  = [-2.2, -3]
    orb_sol= [1]                     # Just one case here for 14-Mar-2018 ORT2 case, but might have more in future.
        
    b      = hbt.frange(-12,0)  # Exponent for beta, and radius, and bin width. 
                                #  b=[-12,0] → beta=[1e-4,1], which is the range that DPH uses.
                                # We iterate over 'b', rather than over beta or size.
                                
#    s      = 10.**b/3            # Particle size, some units?
    
    
    epsilon = 1e-5               # A small value for testing equality of floats.    
    
    # We create an astropy table for the output.
    
    t = Table(names=['Albedo', 'q', 'rho', 'speed', 'img_max', 'img_2d', 'E_0', 'profile'],
              dtype = [float, float, float, float,    float,   'object',     float, 'object'])
    
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
                        # They have already been flattened -- we are not summing in XYZ space again.
                        
                        img_i = (sarea * albedo_i * pi * s_i**2 * ds_i * s_i**q_i *
                                 np.sum(density_flattened_arr[is_good], axis=0) /
                                 (ring.km_per_cell_x * ring.km_per_cell_y) * 1e-12) / num_subsets
                        
                        # Add this image (with one beta) to the existing image (with a dfft beta)
                        
                        img += img_i
                        
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
                    
                    # Take a radial profile
                    
                    (radius_pix, profile) = get_radial_profile_circular(img)
                    plt.imshow(scale(ring))
                    
                    # We are now finished with summing this array over particle size (ie, beta). Save it!
                    
                    t.add_row((albedo_i, q_i, rho_i, speed_i, np.amax(img), img, 0, profile))
                
    # Now that all combinations have been done, make some plots

    hbt.figsize((8,6))
    for i in range(len(t)):        
        plt.plot(radius_pix * ring.km_per_cell_x, t['profile'][i])
    plt.yscale('log')    
    plt.xlabel('Radius [km]')
    plt.ylabel('I/F (non-normalized)')
    plt.ylim((1e-12, 1e-6))
    plt.show()
    
    
    columns = ['Albedo', 'q', 'rho', 'speed']
    for column in columns:
        plt.plot(t[column], t['img_max'], marker = 'o', markersize=3, ls = 'none')
        plt.yscale('log')
        plt.xlabel(column)
        plt.ylabel('Max I/F')
        plt.show()
    
    # And make a plot showing all of the flattened disks
    
    hbt.figsize((25,25))
    
    plot_flattened_grids_table(t)
    
    # Make radial profiles

    binwidth = 1

    radius = hbt.frange(1,100)
    for 
    profile 
    
    
    # - *Average* all of the subsets, 0 .. 7
    # - 
    # Iterate over beta       {-12 .. 1}
    # Iterate over two speeds {-3, -2.2}
    # 
    # The output for each combination will a value of E0.
        
    
    # Finally, do a weighted sum of all the planes, and plot it.
    # Based on the prelim data, we should 
    
#    plot_flattened_grids(np.sum(density_flattened_arr * weighting, axis=0))
#    
#    plt.imshow(stretch(np.sum(density_flattened_arr * weighting, axis=0)))
#    plt.imshow((np.sum(density_flattened_arr * weighting, axis=0)))
#    plt.show()
    
# =============================================================================
# Now do a one-off test to make some plots to validate the XYZ orientation
# =============================================================================
    
def tester():

    dir =  '/Users/throop/data/ORT2/hamilton/deliveries/'
    runs_full = [ dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003/y3.0/beta1.0e-04/subset07',
                  dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003/y3.0/beta1.0e-03/subset02',
                  dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003/y3.0/beta2.2e-02/subset02',
                  dir + 'sbp_ort2_ba2pro4_v1_DPH/sbp_ort2_ba2pro4_v1/ort2-0003/y3.0/beta2.2e-02/subset05']
    
    num_axis = {}
    num_axis = {'X' : 0, 'Y' : 1, 'Z' : 2}
    axes     = ['X',     'Y',     'Z']
    
    for run_full in runs_full:
        i = 1
        ring.print_info()
    
        for axis in axes:
            run  = run_full.replace(dir_base, '')[1:]  # Remove the base pathname from this, and initial '/'           
            ring = nh_ort_track3_read(run)
            plt.subplot(1,3,i)
            plt.imshow(stretch(np.sum(ring.density, axis=num_axis[axis])), extent=extent)
            plt.title(f'Summed along {axis}')
            i+=1
        plt.show()
    
