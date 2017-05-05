#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:59:29 2017

@author: throop
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
from   astropy.table import Table, vstack
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
from astropy.utils import data

# HBT imports

import hbt

def nh_jring_create_objectlist(file_in, do_stars=True, bodies = [], num_stars_max = 100):
    
    ''' 
    Creates a text file which lists all the stars in a file, with lines like
    
        'star', <xpos>, <ypos>, mag [optional]

    <xpos> and <ypos> are the x and y coordinate centers, in pixels.
    The output is sorted by magnitude, if it is available. 
    
    <bodies> is a list, like ['Adrastea', 'Amalthea'], etc.
    
    It is OK to have xpos and ypos be outside the screen. We might want to know that sometimes, for satellites.
    '''

#    file_in    = 'lor_0034962025_0x630_sci_1_opnav.fit'
    dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all/'
    dir_out    = '/Users/throop/data/NH_Jring/out/'
    
    file_in_base = file_in.split('/')[-1]   # Strip the pathname
    file_out_base = file_in_base.replace('.fit', '_objects.txt')
    file_out   = dir_out + file_out_base

# If we were passed a

    file = dir_images + file_in_base
   
    header = hbt.get_image_header(dir_images + file_in_base)

    dx_pix = header['NAXIS1']
    dy_pix = header['NAXIS2']
    
    radius_search = 0.2 * u.deg
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = WCS(file)
        
    name_cat = u'Guide Star Catalog v2 1'# Works on gobi only (no tomato)
    url_cat = 'http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=GSC23&' # Works always

    with data.conf.set_temp('remote_timeout', 30): # This is the very strange syntax to set a timeout delay.
                                                   # The default is 3 seconds, and that times out often.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stars = conesearch.conesearch(w.wcs.crval, radius_search, cache=True, catalog_db = url_cat)

    # NB: the returned value is a Table, but I have to access via .array[] -- not sure why.
    
    ra_stars  = np.array(stars.array['ra'])*hbt.d2r # Convert to radians
    dec_stars = np.array(stars.array['dec'])*hbt.d2r # Convert to radians

    mag_stars = np.array(stars.array['Mag'])
    
    print("Stars downloaded: {}; mag = {} .. {}".format(np.size(mag_stars), np.nanmin(mag_stars), np.nanmax(mag_stars)))
    print("RA = {} .. {}".format(np.nanmin(ra_stars)*hbt.r2d, np.nanmax(ra_stars)*hbt.r2d))
    
    # Now sort by magnitude, and keep the 100 brightest

    order = np.argsort(mag_stars)
    order = np.array(order)[0:num_stars_max]

    ra_stars        = ra_stars[order]
    dec_stars       = dec_stars[order]
    mag_stars       = mag_stars[order]
    
    radec_stars        = np.transpose(np.array((ra_stars,dec_stars)))
    x_stars, y_stars   = w.wcs_world2pix(radec_stars[:,0]*hbt.r2d, radec_stars[:,1]*hbt.r2d, 0)
  
    is_good = np.logical_and( np.logical_and(x_stars >=0, x_stars <= dx_pix),
                              np.logical_and(y_stars >=0, y_stars <= dy_pix) )
    
# Now make a table
    
    t_stars          = Table()
    t_stars['name']  = np.zeros(np.shape(mag_stars[is_good]), dtype='U30')
    t_stars['name'][:] = u'star'    

    t_stars['x_pix'] = x_stars[is_good]
    t_stars['y_pix'] = y_stars[is_good]
    t_stars['mag']   = mag_stars[is_good]


#==============================================================================
# Now find the satellite locations
#==============================================================================

    if np.size(bodies) > 0:
        
# Look up satellite positions

        et = header['SPCSCET']
        utc = sp.et2utc(et, 'C', 0)
        x_bodies, y_bodies = hbt.get_pos_bodies(et, bodies, units='pixels', wcs=w)
        t_sats = Table()
        t_sats['x_pix'] = x_bodies
        t_sats['y_pix'] = y_bodies
        t_sats['name']  = np.array(bodies).astype('U30')

#==============================================================================
# Merge the stars and sats into one table
#==============================================================================
                 
        t_merged = vstack([t_stars, t_sats])

    else:
        t_merged = t_stars
              
#==============================================================================
# And write the table to disk
#==============================================================================

    t_merged.write(file_out, format = 'csv', overwrite = True)
    print("Wrote: {} ({} objects)".format(file_out, np.shape(t_merged)[0]))
#    
    return t_merged

#==============================================================================
# A few lines for diagnostic testing
#==============================================================================

def test():
    
    file_in    = 'lor_0034962025_0x630_sci_1_opnav.fit'
    
    sats = ['Adrastea', 'Thebe', 'Metis', 'Io', 'Europa', 'Amalthea']
    t = nh_jring_create_objectlist(file_in, bodies=sats)
