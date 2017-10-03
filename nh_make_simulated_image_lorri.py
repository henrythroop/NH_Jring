#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 10:12:26 2017

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

# =============================================================================
# Make a simulated image
# =============================================================================

def nh_make_simulated_image_lorri(do_ring = False, do_mu69 = False, 
                                  mode = '1X1', exptime = 10, 
                                  dist_solar = 40*u.au, dist_target = 0.01*u.au,
                                  a_ring = (1000*u.km, 3000*u.km), do_psf = True,
                                  pos = (None, None)):

    
    if (mode == '1X1'):
        naxis = 1024     # Number of pixels in the output array. Assumed to be square.
    else:
        naxis = 256      # 4X4 mode
    
    # Create the output array
    
    arr = np.zeros((naxis1, naxis1))    
  
    # Compute the 
    
    if pos[0] == None:
        pos = (naxis/2, naxis/2)
    
    
    