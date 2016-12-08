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
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
#from   photutils import datasets
from   astropy.stats import sigma_clipped_stats
from   scipy.stats import mode
from   scipy.stats import linregress
from   photutils import daofind
import wcsaxes
import time
import sys  # For stdout.write, without newline
from scipy.interpolate import griddata

from mpl_toolkits.axes_grid1 import host_subplot # For adding a second axis to a plot
import mpl_toolkits.axisartist as AA             # For adding a second axis to a plot

import imreg_dft as ird
import re # Regexp
import pickle # For load/save

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
#==============================================================================
    
def add_backplanes():
    
    dir_out = '/Users/throop/data/NH_LORRI_Ring/'
    
    files_merged = glob.glob(dir + '*/*mos_header*fit*')
    
    for file in files_merged:
        plane = hbt.create_backplane(file)
        file_out = file_out.replace('.fit', '_planes.pkl')
        
for i,file in enumerate(files):

    file_short = file.split('/')[-1]
    file_out = dir_out + file_short    
    file_out = file_out.replace('.fit', '_planes.pkl')
    print("{}/{}: Generating backplane for {}".format(i, np.size(files), file_short))

    plane = hbt.create_backplane(file)

    
    # Create the backplanes. 

    print "Generating backplanes..."
    
    planes = hbt.create_backplane(file_out, frame = 'IAU_PLUTO', name_target='Pluto', name_observer='New Horizons')

    file_header_pl = file_header.replace('.fits', '_pl.fits')     # File for backplanes

# Create backplaned FITS file.

# Create a new FITS file, made of an existing file plus these new planes

    hdulist = fits.open(file_header_out)  # Read in the file we've just created, which has full fixed header

    type_out = 'float32'  # Write backplane in single precision, to save space.
    hdu_out = fits.HDUList() # Create a brand new FITS file
    hdu_out.append(hdulist['PRIMARY'])
    hdu_out.append(fits.ImageHDU(planes['RA'].astype(type_out), name = 'RA'))
    hdu_out.append(fits.ImageHDU(planes['Dec'].astype(type_out), name = 'Dec'))
    hdu_out.append(fits.ImageHDU(planes['Phase'].astype(type_out), name = 'Phase'))
    hdu_out.append(fits.ImageHDU(planes['Longitude_eq'].astype(type_out), name = 'Longitude_eq'))
    hdu_out.append(fits.ImageHDU(planes['Radius_eq'].astype(type_out), name = 'Radius_eq'))
    
    hdu_out.writeto(file_header_pl, clobber=True)
    print "Wrote file: " + file_header_pl   
    
