#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 08:52:54 2018

@author: throop
"""


import pdb
import glob
                  # but apparently it still must be imported.
from   subprocess import call
import warnings
import os.path
import os

import astropy
from   astropy.io import fits
import astropy.table   # I need the unique() function here. Why is in in table and not Table??
from   matplotlib.figure import Figure

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

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular

from   get_radial_profile_circular import get_radial_profile_circular
from   plot_img_wcs import plot_img_wcs
from   wcs_translate_pix import wcs_translate_pix

# HBT imports

import hbt

def nh_ort_directorize_by_reqid(dir=None):
    """
    This is just a one-off program to take the ORT images, and put them into folders sorted by reqid.
    Simon organizes them this way by default. Marc does not.
    """

    dir = '/Users/throop/Data/ORT3/throop/backplaned'
    
    files = glob.glob(os.path.join(dir, '*.fit'))
    
    for file in files:
        hdulist        = fits.open(file)
        reqid          = hdulist['PRIMARY'].header['REQID']
        
        # Put OPNAVs in their own directory. However, I see that Simon didn't actually do this, so we can comment out.
        if ('OPNAV' in reqid):
#            dir_out = os.path.join(dir, 'K1LR_OPNAV', reqid)
            dir_out = os.path.join(dir, reqid)
        else:
            dir_out = os.path.join(dir, reqid)
        if not os.path.isdir(dir_out):
            print(f'mkdir {dir_out}')
            
        print(f'mv {file} {dir_out}')    
    
    