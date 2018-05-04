#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 13:28:06 2016

@author: throop
"""

# Read LORRI outbound Pluto rings mosaics and calculate I/F limit from them.
# This is not a rigorous reduction -- just want to get an answer so I can write an MT
# about using LORRI vs. MVIC for unresolved KEM observations. 
# HBT 7-Dec-2016

# General python imports

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
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
from   scipy.optimize import curve_fit
#from   pylab import *  # So I can change plot size.
                       # Pylab defines the 'plot' command
import spiceypy as sp # was cspice
import skimage
from   astropy.wcs import WCS
from   astroquery.vo_conesearch import conesearch
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils import daofind
import time
import sys  # For stdout.write, without newline
from scipy.interpolate import griddata

from mpl_toolkits.axes_grid1 import host_subplot # For adding a second axis to a plot
import mpl_toolkits.axisartist as AA             # For adding a second axis to a plot

import imreg_dft as ird
import re # Regexp
import pickle # For load/save
import spiceypy as sp

# Imports for Tk

from astropy import units as u
from astropy.coordinates import SkyCoord

# HBT imports
import hbt

dir = '/Users/throop/Data/NH_LORRI_Ring/'
stretch_percent = 95

files = glob.glob(dir + '*/*_header*')

hdulist = []
a       = [] # List of image arrays, one per file

for file in files:
    hdulist.append(fits.open(file))
    a.append(hdulist[-1]['PRIMARY'].data)

hbt.set_plot_defaults
plt.set_cmap('Greys_r')
 
stretch = astropy.visualization.PercentileInterval(stretch_percent)

for i in range(np.size(files)): 
    arr = np.array(a[i])

    hbt.figsize((6,25))
    
    plt.imshow(stretch(arr))
    plt.title(hdulist[i]['PRIMARY'].header['VISITNAM'])
    plt.show()

    binning = 100
    
    hbt.figsize((5,5))    
    plt.plot(arr[800,:], label = 'Raw')
    plt.plot(hbt.smooth_boxcar(arr[800,:],binning), label='Binning = {}'.format(binning), color='red')
    plt.title(hdulist[i]['PRIMARY'].header['VISITNAM'])
    plt.ylabel('DN')
    plt.xlabel('Column')
    plt.legend()

    plt.show()
    

#==============================================================================
# A one-off function to add valid headers to Tod's headerless mosaic files
#==============================================================================

def merge_mosaic_header():
   
    files_mos = glob.glob(dir + '*/*mos*fit*') # Read Tod's mosaics
    files_raw = glob.glob(dir + '*/lor*fit*')  # Read the initial image in each of the two sequences 
                                               # O_RING_DEP_LORRI_202, _305
    
    hdulist = []
    for file in files:
        hdulist.append(fits.open(file))
    
    # Now merge the headers and the mosaics
    
    files_out = []
    for file in files_mos:
        files_out.append(file.replace('_mos_', '_mos_header_'))
    
    hbt.merge_fits(files_mos[0], files_raw[0], files_out[0])
    hbt.merge_fits(files_mos[1], files_raw[1], files_out[1])
    
#==============================================================================
# A one-off function to take the files with headers, and wcs, and add backplanes to them.
# ring_lor202_mos_v1.fits               -- from Tod
# ring_lor202_mos_v1_header.fits        -- with NH header from individual LORRI added 
# ring_lor202_mos_v1_header_wcs.fits    -- with WCS info added
# ring_lor202_mos_v1_header_wcs_pl.fits -- with backplanes added
#==============================================================================
    
def add_backplanes():
  
    dir = '/Users/throop/data/NH_LORRI_Ring/'
    
    files = glob.glob(dir + '/*_header_wcs.fits')
    
    file_tm           = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel
 
    # Make sure SPICE is running
    
    sp.furnsh(file_tm)
    
    for i,file in enumerate(files):

        file_out = file.replace('.fits', '_pl.fits')
        
        print("{}/{}: Generating backplane for {}".format(i, np.size(files), file))
    
        # Create the backplanes
    
        print("Generating backplane...")
                
        (planes,descs) = hbt.compute_backplanes(
                             file, frame = 'IAU_PLUTO', name_target='Pluto', name_observer='New Horizons')
            
    # Create a new FITS file, made of an existing file plus these new planes
    
        hdulist = fits.open(file)  # Read in the file we've just created, which has full fixed header
    
        type_out = 'float32'        # Write backplane in single precision, to save space.
        hdu_out = fits.HDUList()    # Create a brand new FITS file
        hdu_out.append(hdulist['PRIMARY'])
        hdu_out.append(fits.ImageHDU(planes['RA'].astype(type_out), name = 'RA'))
        hdu_out.append(fits.ImageHDU(planes['Dec'].astype(type_out), name = 'Dec'))
        hdu_out.append(fits.ImageHDU(planes['Phase'].astype(type_out), name = 'Phase'))
        hdu_out.append(fits.ImageHDU(planes['Longitude_eq'].astype(type_out), name = 'Longitude_eq'))
        hdu_out.append(fits.ImageHDU(planes['Radius_eq'].astype(type_out), name = 'Radius_eq'))
        
        hdu_out.writeto(file_out, clobber=True)
        print("Wrote file: " + file_out)