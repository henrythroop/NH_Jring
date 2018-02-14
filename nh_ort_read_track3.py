#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:38:48 2018

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

class nh_ort_read_track3:  # Q: Is this required to be same name as the file?

# =============================================================================
# Initialize the method.
# Read the specified track3 file into memory, from disk
# =============================================================================
    
    def __init__(self, name, dir_base='/Users/throop/data/ORT1/kaufmann/deliveries', binfmt = 2):
        
        """
        Initialize the method.
        
        Paramateters
        ------
            
        name:
            The name of the run. This is a slash-separated string, such as 
            `'ort1-0001/1000010/speed1/beta_2.6e-04/subset00'`
             
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
        self.km_per_cell_x = km_per_cell_x
        self.km_per_cell_y = km_per_cell_y
        self.km_per_cell_z = km_per_cell_z
        
        self.body_id   = self.get_header_val('body_id')
        self.duration  = self.get_header_val('duration')*u.year
        self.time_step = self.get_header_val('time_step')*u.s
        self.beta      = self.get_header_val('beta')
        self.weight    = self.get_header_val('weight')   # Statistical weight of this solution
        self.speed     = self.get_header_val('speed')    # Should be 1 or 2 as per wiki. But is really 0 or 1. 
        
        # We should read the full header in to the class as well. But this will require more parsing.
        
        return

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
        
        return

# =============================================================================
# Run the file
# =============================================================================
    
if __name__ == '__main__':
    
    dir_base='/Users/throop/data/ORT1/kaufmann/deliveries'
    
    runs_full = glob.glob(os.path.join(dir_base, '*', '*', '*', '*', '*'))
    run = runs_full[0].replace(dir_base, '')[1:]
    
    for run_full in runs_full:
        run = run_full.replace(dir_base, '')[1:]  # Remove the base pathname from this, and initial '/'
        ring = nh_ort_read_track3(run)
        ring.plot()
    
     