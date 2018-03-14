#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:38:48 2018

This reads the NH MU69 ORT 'Track 3' data, which is the output of N-Body simulations done by
Doug Hamilton + David Kaufmann.

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

# HBT imports

import hbt

class nh_ort_track3_read:  # Q: Is this required to be same name as the file?

# =============================================================================
# Initialize the method.
# Read the specified track3 file into memory, from disk
# =============================================================================
    
    def __init__(self, name, dir_base='/Users/throop/data/ORT1/kaufmann/deliveries', binfmt = 2):
        
        """
        Initialize the method. Read in data from a specified file, and do initial processing on it.
        
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
            
        self.name = name
    
        parts = name.split('/')

        dir = name
        
        # Construct the filename for the 'header', which is a short text file. It is basically the same
        # as the entire directory tree, but with '/' → '_'.
        # There is one difference -- probably a bug -- in that 'speed' term is off by one, as per DK email. OBOB.
        
        file_header =(parts[0].replace('-', '_') + '_' + parts[1] + '_' + 
                      parts[2].replace('1', '0').replace('2', '1') + '_' +
                      parts[3] + '_' + parts[4] +
                      '.header')

        # Construct the filename for the data array itself. It is trivial.
        # As per the wiki, the data and header files should be identically named. They are not!
        # https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Integrations-to-Density
        
        file_data = f'model.array{binfmt}'
        
        file_header_full = os.path.join(dir_base, dir, file_header)
        file_data_full   = os.path.join(dir_base, dir, file_data)

        self.file_header_full = file_header_full
        self.file_data_full   = file_data_full
        
        # Read the datafiles. Using code fragment from DK email 13-Feb-2018

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
        
        # Read a bunch of values from the header, and save them into the method
        
        self.km_per_cell_x = km_per_cell_x
        self.km_per_cell_y = km_per_cell_y
        self.km_per_cell_z = km_per_cell_z
        
        self.body_id   = self.get_header_val('body_id')
        self.duration  = self.get_header_val('duration')*u.year
        self.time_step = self.get_header_val('time_step')*u.s
        self.beta      = self.get_header_val('beta')
        self.weight    = self.get_header_val('weight')   # Statistical weight of this solution
        self.speed     = self.get_header_val('speed')    # Should be 1 or 2 as per wiki. But is really 0 or 1.
        self.ejected   = self.get_header_val('ejected')
        self.gm        = self.get_header_val('gm')
        
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
        
        stretch_percent = 90    
        stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
 
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

        print(f'Name:     {self.name}')
        print(f'Beta:     {self.beta}')
        print(f'Duration: {self.duration}')
        print(f'Speed:    {self.speed}')
        print(f'Timestep: {self.time_step}')
        print(f'Body ID:  {self.body_id}')
        print(f'GM:       {self.gm}')
        print(f'ejected:  {self.ejected}')
        
        return
        
# =============================================================================
# Run the file
# =============================================================================
    
if __name__ == '__main__':

    do_short = False

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    dir_base='/Users/throop/data/ORT1/kaufmann/deliveries'
    
    runs_full = glob.glob(os.path.join(dir_base, '*', '*', '*', '*', '*'))
    if do_short:
        runs_full = runs_full[0:10]
        
    run = runs_full[0].replace(dir_base, '')[1:]  # Populate one, just for cut & paste ease
    
    num_files             = len(runs_full)
    density_flattened_arr = np.zeros((num_files, 200, 200))
    beta_arr              = np.zeros(num_files)
    time_step_arr         = np.zeros(num_files)
    duration_arr          = np.zeros(num_files)
    speed_arr             = np.zeros(num_files)
    body_id_arr           = np.zeros(num_files)
    ejected_arr           = np.zeros(num_files)
    gm_arr                = np.zeros(num_files) 
    
    weighting = np.zeros((len(runs_full), 1, 1)) # Weight each individual run. We make this into a funny shape
                                                 # so we can broadcast it with the full array easily.
    
    for i,run_full in enumerate(runs_full):
        run = run_full.replace(dir_base, '')[1:]  # Remove the base pathname from this, and initial '/'
        ring = nh_ort_track3_read(run)            # Read the data array itself
#        ring.plot()
        ring.print_info()
        density_flattened_arr[i,:,:] = np.sum(ring.density, axis=0)
        beta_arr[i] = ring.beta
        time_step_arr[i] = ring.time_step.value
        duration_arr[i] = ring.duration.value
        speed_arr[i] = ring.speed
        body_id_arr[i] = ring.body_id
        gm_arr[i] = ring.gm
        ejected_arr[i] = ring.ejected
        
        print('-----')
    
    print(f"Read {num_files} files.")
    # Assign some (very arbitrary) weights to each plane
    
    weighting[:,0,0] = beta_arr[:]
    weighting[:,0,0] = (body_id_arr == 4)
    
    # Finally, do a weighted sum of all the planes, and plot it.
    # Based on the prelim data, we should 
    
    plt.imshow(stretch(np.sum(density_flattened_arr * weighting, axis=0)))
    plt.show()
    