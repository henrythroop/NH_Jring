#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

File to create the backplanes for MU69 encounter images.
Written for hazard search ORTs, but should work more generally.

Created on Tue Jan  9 13:19:09 2018

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
from   image_stack import image_stack
from   create_backplane import create_backplane

def nh_create_backplanes_fits(file_in = None, 
                              file_out = None, 
                              name_target = 'MU69', 
                              name_observer = 'New Horizons', 
                              clobber = False,
                              type = 'Sunflower'):
    
    """
    Function to create backplanes and add them to a FITS file.
    
    Idea is that this should be runnable as a commandline script.
    
    Parameters
    ----
    
    file_in:
        Input FITS file. Should be fully navigated and have good header info (correct WCS, etc)
    file_out:
        Output FITS file to be written. If None, then no file is written.
    name_target:
        String. Name of central body of the backplanes. For now, must be 'MU69'.
    name_observer:
        String. Name of observer. Usually 'New Horizons'
    clobber:
        Boolean. The original file will be overwritten, iff file_out == file_in and file==True.
    type:
        Type of orbit to assume for the backplanes. Can be 'Sunflower'. Other orbits might be added later, like ones 
        that face north pole, or are tilted relative to ecliptic, rather than face the Sun.
    
    Output
    ----    
    
    tuple:
        A tuple of(radius, azimuth, d_ra, d_dec, etc). This has all of the backplanes in it.
    
    Return status
    ----
        0 or 1. 0 = no problems.
        
    """

#- Open the file
#- Read the time, mode (1x1 vs 4x4), image scale, etc.
#- Compute the backplanes
#- Copy the FITS fileconso
#- Write the new FITS file
#- Return output to user

# Intialize graphics
    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

    plt.set_cmap('Greys_r')

    frame = 'J2000'
    name_target = 'MU69'
    name_observer = 'New Horizons'        
    abcorr = 'LT+S'
    
# Start up SPICE
    
    file_kernel = 'kernels_kem.tm'
    sp.furnsh(file_kernel)
    
# Load the image
        
    hdu = fits.open(file_in) 
    data = hdu['PRIMARY'].data
    header = hdu['PRIMARY'].header

# Grab some critical info from header

    et      = header['SPCSCET']
    sformat = header['SFORMAT']
    dx_pix  = header['NAXIS1']
    dy_pix  = header['NAXIS2']
    mode    = header['SFORMAT']
    exptime = header['EXPTIME']
    

# Do some processing on the header

    utc = sp.et2utc(et, 'C', 0)
    w = WCS(header)
    
# Look up position of MU69 in pixels.
    
    (st,lt) = sp.spkezr(name_target, et, frame, abcorr, name_observer)
    (_, ra, dec) = sp.recrad(st[0:3])
    (pos_pix_x, pos_pix_y) = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)
    
# Display the image, and MU69.
    
    plt.imshow(stretch(data))
    plt.plot(pos_pix_x, pos_pix_y, color='red', marker = 'o', ms=5, alpha=0.5)
    plt.show()

# Call a routine to actually create the backplanes, and return as a tuple.

    planes = create_backplane(file_in,
                              type = 'Sunflower',
                              name_target = name_target,
                              name_obserer = name_observer
                              )
    
# Return to user
    
    return(0)

# 

if (__name__ == '__main__'):
    
    file_in = os.path.join(os.path.expanduser('~'), 'Data', 'NH_KEM_Hazard', 'ORT1_Jan18', 
                               'lor_0406731132_0x633_sci_HAZARD_test1.fit')

    planes = nh_create_backplanes_fits(file_in, None)