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

        #### START OF DK CODE FRAGMENT ####
        
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
                # Loop over list of non-empty cells. For each one, read xyz position (6 bytes), and density (4 bytes)
                start = 10*i
                (ix,iy,iz) = struct.unpack('hhh', data[start:start+6])   # h = short int
#                (density[iz-1,iy-1,ix-1],) = struct.unpack('f', data[start+6:start+10])  # Original order ZYX
                (density[ix-1,iy-1,iz-1],) = struct.unpack('f', data[start+6:start+10])   # HBT modified order XYZ, 
                                                                                          # to match wiki API.
                
    
        #### END OF DK CODE FRAGMENT ####
        
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
            plt.imshow(stretch_hbt(np.sum(self.density, axis=axis)), extent=extent)
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
            
    plt.imshow(stretch_hbt(img))
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
 
    stretch = astropy.visualization.PercentileInterval(stretch_percent)
    
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

def stretch_hbt(arr):
    return np.log(10 + arr)

def stretch_hbt_invert(arr):
    return np.exp(arr) - 10
        

