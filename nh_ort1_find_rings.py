#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 02:47:12 2018

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

from astropy.coordinates import SkyCoord

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes

def nh_ort1_find_rings():
    
    plt.set_cmap('Greys_r')

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

    dir = '/Users/throop/Data/ORT1/throop/backplaned'
    files = glob.glob(os.path.join(dir, '*', '*fits'))
    
    hbt.figsize((15,8))
    
    # Set up output arrays
    
    ra_arr    = []
    dec_arr   = []
    reqid_arr = []
    exptime_arr= []
    et_arr    = []
    utc_arr   = [] 

    # Start up SPICE
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem.tm')
    
    for ii,file in enumerate(files):
        
        hdu = fits.open(file)
        print('Reading {}/{}: {}'.format(ii, len(files), os.path.basename(file)))
        img = hdu[0].data
        
        header = hdu[0].header
        
        ra_arr.append(header['CRVAL1'])
        dec_arr.append(header['CRVAL2'])
        exptime_arr.append(header['EXPTIME'])
        reqid_arr.append(header['REQID'])
        et_arr.append(header['SPCSCET'])
        utc_arr.append(sp.et2utc(header['SPCSCET'], 'C', 0))
        
        radius_eq = hdu['RADIUS_EQ'].data
        
        dradius = 1000
        num_bins_radius = 100
        
        bins_radius = hbt.frange(0, np.amax(radius), num_bins_radius)
        dn_median_arr = np.zeros(num_bins_radius)
        dn_mean_arr   = np.zeros(num_bins_radius)
        
        for i in range(num_bins_radius-1):
            is_good = np.logical_and(radius_eq > bins_radius[i], radius_eq < bins_radius[i+1])
            dn_median_arr[i] = np.nanmedian(img[is_good])
            dn_mean_arr[i]   = np.nanmean(img[is_good])

        do_plot = False

        if do_plot:
            
            plt.subplot(1,2,1)
            plt.plot(bins_radius, dn_median_arr, label = 'median')
            plt.plot(bins_radius, dn_mean_arr,   label = 'mean')
            plt.legend(loc = 'upper right')
            plt.title("{}/{}  {}".format(ii,len(files), os.path.basename(file)))
           
            
            plt.subplot(1,2,2)
            plt.imshow(stretch(img))
            plt.show()
        
        hdu.close()
        
# =============================================================================
# Read the values into NumPy arrays
# =============================================================================

    ra   = np.array(ra_arr)
    dec  = np.array(dec_arr)
    reqid = np.array(reqid_arr)
    et = np.array(et_arr)
    exptime = np.array(exptime_arr)
    utc  = np.array(utc_arr)
    
    plt.plot(ra, dec, ls='none', marker = 'o', ms=2)
    
    # Put them all into a table
    
    t = Table(          [ra, dec, et, utc, exptime, reqid], 
              names = ('RA', 'Dec', 'ET', 'UTC', 'EXPTIME', 'ReqID'))
    
    t = Table([a, b, c], names=('a', 'b', 'c'), meta={'name': 'first table'})
    
    w_haz0 = (t['ReqID'] == 'K1LR_HAZ00')
    w_haz1 = (t['ReqID'] == 'K1LR_HAZ01')
    w_haz2 = (t['ReqID'] == 'K1LR_HAZ02')
    w_haz3 = (t['ReqID'] == 'K1LR_HAZ03')
    w_haz4 = (t['ReqID'] == 'K1LR_HAZ04')
    
    plt.plot(ra[w_haz0], dec[w_haz0], marker='o', ls='none')
    plt.plot(ra[w_haz1], dec[w_haz1], marker='o', ls='none')
    plt.plot(ra[w_haz2], dec[w_haz2], marker='o', ls='none')
    plt.plot(ra[w_haz3], dec[w_haz3], marker='o', ls='none')
    plt.plot(ra[w_haz4], dec[w_haz4], marker='o', ls='none')
    plt.show()
    
    plt.plot(et[w_haz0], marker='o', ls='none')
    plt.plot(et[w_haz1], marker='o', ls='none')
    plt.plot(et[w_haz2], marker='o', ls='none')
    plt.plot(et[w_haz3], marker='o', ls='none')
    plt.plot(et[w_haz4], marker='o', ls='none')

# =============================================================================
# Run the function
# =============================================================================

if (__name__ == '__main__'):
    nh_ort1_find_rings()
    