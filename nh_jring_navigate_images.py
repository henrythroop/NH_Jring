#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 16:15:18 2017

@author: throop
"""

"""
Library routine to navigate an image based on stars.
  o Use DAOPhot to locate stars in the image
  o Use a star catalog to look up stars that should be in the image
  o Calculate the offset between these.
  o Return the result, in terms of a revised WCS for the frame, and a pixel offset.
  o Does not rewrite the FITS header.
"""

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
from   astropy.utils import data

from   scipy.optimize import curve_fit
                       # Pylab defines the 'plot' command
import spiceypy as sp
from   astropy.wcs import WCS
from   astropy.vo.client import conesearch # Virtual Observatory, ie star catalogs
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
from   scipy.stats import mode
from   scipy.stats import linregress
import wcsaxes
import time
from   scipy.interpolate import griddata
#import cv2

import re # Regexp
import pickle # For load/save

import cProfile # For profiling

from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure
import warnings
from   importlib import reload
from   time import gmtime, strftime

# HBT imports

import hbt

#def nh_jring_navigate_images():

"""
Navigate all the NH Jring images
"""

#    file = dir + 'lor_0034765323_0x630_sci_1.fit' # This one is faint -- hard to see much. But algo works -both.
##    file_in = dir + 'lor_0034602123_0x630_sci_1.fit'  # Algo works. Both g_t_i and fft.
##    file_in = dir + 'lor_0034613523_0x630_sci_1.fit'  # Fails fft. Works g_t_i.
##    file_in = dir + 'lor_0034676528_0x630_sci_1.fit' # Good test case, very simple, good stars. g_t_i works great.
##    file_in = dir + 'lor_0034676528_0x630_sci_1_opnav.fit' # Post-navigation
##    file_in = dir + 'lor_0034676528_0x630_sci_1_opnav_opnav.fit' # Post-navigation


dir =  '/Users/throop/Dropbox/Data/NH_Jring/data/jupiter/level2/lor/all/' 

files = glob.glob(dir + '*[123456].fit')  # Exclude any '*_opnav.fit'

plt.set_cmap('Greys_r')            
hbt.figsize((10,10))
do_plot = True

for i,file in enumerate(files):
    
    file_out = file.replace('.fit', '_opnav.fit')
    
    im = hbt.read_lorri(file) # Read the image, and process it a bit I think
    
    hdulist = fits.open(file) 
    header = hdulist['PRIMARY'].header
    mode = header['SFORMAT']     
    hdulist.close()           
    
# Want to only do this if it is LORRI 1x1, not 4x4.
# Also, skip if this file has already been done?

    if ((mode == '1X1') and (not(os.path.isfile(file_out)))):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")    
            w_orig  = WCS(file)           # Look up the WCS coordinates for this frame
            w       = WCS(file)           # Get a copy of it, which we'll change
        
        crval_orig = w_orig.wcs.crval
    
    # Do the navigation call
    
        method_opnav = 'fft'
        (w, (dy_pix, dx_pix)) = hbt.navigate_image_stellar(im, w, method = method_opnav,
                                 title = file.split('/')[-1], do_plot=do_plot)
    
        crval = w.wcs.crval
    
        print("Center was at RA {}, Dec {}".format(crval_orig[0], crval_orig[1]))
        print("Center is now at RA {}, Dec {}".format(crval[0], crval[1]))
        print("Determined OpNav offset: dx = {} pix, dy = {} pix".format(dx_pix, dy_pix))
    
    # Read in the existing FITS file to memory, so we can modify it and write it back out
    # Q: Is it just necessary to update these fields, or do I need to do more?
    
    # A: Looks like just updating CRVAL1/2 in the text header itself should be enough.
    #    See http://www.stsci.edu/hst/HST_overview/documents/multidrizzle/ch44.html
    
        header_wcs = w.to_header()   # NB: w.to_header, w.to_header(), and w.wcs.to_header() all do different things!
        
        hdulist = fits.open(file)
        header  = hdulist[0].header
        header['CRVAL1'] = crval[0] # RA, deg
        header['CRVAL2'] = crval[1] # Dec, deg
        header['OPNAVDX'] = (dx_pix, '[pix] Offset determined from stellar navigation, X')
        header['OPNAVDY'] = (dy_pix, '[pix] Offset determined from stellar navigation, Y')
        header['OPNAVMTH'] = (method_opnav, 'Method for OpNav stellar navigation')
        header['comment'] = 'OpNav CRVAL via navigate_image_stellar, ' + \
                             strftime("%Y-%m-%d %H:%M:%S UTC", gmtime())
    
    # Write out a revised FITS file
    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")            # Ignore warning about a long comment.
            hdulist.writeto(file_out, overwrite=True)  # Write to a new file (open, write, and close, in one command)
            
        hdulist.close()                            # Close the original file
        
        print("Wrote file {}/{}: ".format(i,len(files), file_out))