#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 22:00:03 2017

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

import re # Regexp
import pickle # For load/save


from   matplotlib.figure import Figure

# HBT imports

import hbt

file_tm = "/Users/throop/gv/dev/gv_kernels_new_horizons.txt"  # SPICE metakernel

files = glob.glob(dir + '/*/*.fits')

stretch_percent = 90    
stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

plt.set_cmap('Greys_r')

# Start up SPICE

#file_tm = 'kernels_nh_jupiter.tm'  # SPICE metakernel
sp.furnsh(file_tm) # Commented out for testing while SPICE was recopying kernel pool.

mode = []
exptime = []
filename_short = []
exptime = []
visitnam = []
sapname = []
sapdesc = []
reqid   = []
et      = []
utc     = []
target = []

t = Table()

t = Table(  [[],              [],          [],         [],        [],       [],       [],      [],   []],
    names=('filename_short', 'exptime', 'visitname', 'sapname', 'sapdesc', 'target', 'reqid', 'et',  'utc'),
    dtype = ('S30',           'float64', 'U30',      'U30',     'U30',     'U30',    'U30',   'float64', 'U30'))


for i,file in enumerate(files):
    hdulist        = fits.open(file)
    arr            = hdulist[0].data
    err            = hdulist[1].data
    quality        = hdulist[2].data
    
    filename_short = file.split('/')[-1]
    exptime = hdulist[0].header['EXPTIME']
    visitnam = hdulist[0].header['VISITNAM']
    sapname = hdulist[0].header['SAPNAME']
    sapdesc = hdulist[0].header['SAPDESC']
    target = hdulist[0].header['TARGET']
    reqid = hdulist[0].header['REQID']
    et = hdulist[0].header['SPCSCET']
    utc = sp.et2utc(et, 'C', 1)

    print("Read {}/{} {}".format(i, len(files), filename_short))
    hdulist.close()
    
# Put the FITS parameters into an astropy table row

    t.add_row([filename_short, exptime, visitnam, sapname, sapdesc, target, reqid, et, utc])

# Plot one image. Plot MU69 position on top of it to be sure that WCS is correct.

for i in range(120,150):
    
    hdulist        = fits.open(files[i])
    arr            = hdulist[0].data
    err            = hdulist[1].data
    quality        = hdulist[2].data
    
    w = WCS(files[i])
    
    vec,lt = sp.spkezr('2014 MU69', et, 'J2000', 'LT', 'New Horizons')
    vec_sc_targ = vec[0:3]
    (junk,ra,dec) = sp.recrad(vec_sc_targ) # Get the RA / Dec of the object
    x_pix, y_pix    = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0) # Convert to pixels
        
    plt.imshow(stretch(arr))
    plt.title("{}, exptime {}".format(t['filename_short'][i], t['exptime'][i]))
    plt.plot(x_pix, y_pix, marker = '+', color='red')
    
    plt.ylabel('Y pixels')
    plt.xlabel('X pixels')
    
    plt.show()
    
    
    hdulist.close()
