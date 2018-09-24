#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 16:15:18 2017

@author: throop
"""

"""
User-level program to interactively navigate all LORRI rings images.
Allows user to list all images, navigate one or a series automatically, etc.
New WCS coords are written out to .fits file. 
This is essentially a user interface to hbt.navigate_image_stellar().
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
from   astropy import units as u           # Units library
from   astropy.coordinates import SkyCoord # To define coordinates to use in star search
import time
import sys
#import cv2

import re # Regexp
import pickle # For load/save

import cProfile # For profiling

from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from   matplotlib.figure import Figure
import warnings
from   importlib import reload
from   time import gmtime, strftime
import shutil

# HBT imports

import hbt 

#nh_jring_navigate_images():

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


###

dir_images = '/Users/throop/data/NH_Jring/data/jupiter/level2/lor/all'
dir_out = '/Users/throop/data/NH_Jring/out/'

plt.set_cmap('Greys_r')            
hbt.figsize((15,15))
do_plot           = True
DO_SKIP_NAVIGATED = False
DO_SKIP_4X4       = False
DO_INTERACTIVE    = True
method_opnav      = 'fft'
is_success        = False
index_group       = 6
index_image       = 8

# Load the main observation table

filename_save = 'nh_jring_read_params_571.pkl' # Filename to save parameters in
print("Loading file: {}".format(filename_save))    
lun = open(dir_out + filename_save, 'rb')
t = pickle.load(lun)
lun.close()

# Extract the one group we are looking at 
# NB: This creates an entirely new table -- *not* a view into the original table. 
# If we modify values here, we need to explicitly write the changes back to the original.

groups = astropy.table.unique(t, keys=(['Desc']))['Desc']

groupmask = t['Desc'] == groups[index_group]
t_group = t[groupmask]

num_images_group = np.size(t_group)

###
        
i = 0   # i must be the current image number
        # ii is index within the list
        # k is the keyboard string
ii = 0

# Now start a keyboard loop and run with input from the user

while True:
    
    file_short = t_group[i]['Shortname']
    
    if DO_INTERACTIVE:
        k = input("File {}-{} ({} = {}): ".format(0, np.size(t_group)-1, i, file_short))
    else:
        k = repr(list_batch[ii])         # Load the next element, as a string
        if (ii == np.size(list_batch)-1):  # If we've hit the last element
            DO_INTERACTIVE = True    # And go back to interactive mode

    if (k in ['x', 'q']):            # QUIT
        sys.exit(0)

    if ('g' in k):               # Enter a group number
        k = input(f'Group number: ({index_group}): ')
        if k:
            try:
                index_group = int(k)
                groupmask = t['Desc'] == groups[index_group]
                t_group = t[groupmask]
            except:
                pass
            print(f'Group number = {index_group}')
        k = ''  # Reset it so that loop will not run
            
    if ('-' in k):
        (i1, i2) = np.array(k.split('-')).astype('int')
        list_batch = hbt.frange(int(i1), int(i2)).astype('int')
        ii = 0                        # Current element number, within our list
        k = repr(list_batch[ii])      # Extract one element from it
        print ("Running from {} to {}...".format(i1, i2))
        DO_INTERACTIVE = False

    if ('*' in k):                   # Wildcard search
        searchstr = k.replace('*', '')
        for ii,file in enumerate(files):
            if (searchstr in file):
                print("{}. {}".format(ii, file.split('/')[-1]))
                
    if (k == 'l'):                  # List all files
        print(t_group)
        for i in range(len(t_group)):
            print(f"{i} {t_group[i]['Shortname']}")
            
            # print("{}. {}".format(ii, file.split('/')[-1]))
    
    if (k == '?'):                  # Get help
        print(" <#> = navigate, <#-#> = navigate range, l = list, n = next, x = exit\n" + 
              f' sn = toggle Skip_Navigated ({DO_SKIP_NAVIGATED})\n' +
              f' s4 = toggle Skip_4x4 ({DO_SKIP_4X4})\n' + 
              f' g  = group number ({index_group})')
    
    if (k == 'sn'):
        DO_SKIP_NAVIGATED = not(DO_SKIP_NAVIGATED)
        print("DO_SKIP_NAVIGATED = {}".format(DO_SKIP_NAVIGATED))

    if (k == 's4'):
        DO_SKIP_4X4 = not(DO_SKIP_4X4)
        print("DO_SKIP_4X4 = {}".format(DO_SKIP_4X4))
        
    if (k == 'n') :                 # Next
        k = repr(i)

    # if (k == ''):                   # Next -- do this automatically if <cr> hit
    #     k = repr(i)
    
    if hbt.is_number(k):            # If a number was entered () 

        
        # file = t_group[i]['Filename']
        # file_short = file.split('/')[-1]
        
        file_out   = file.replace('.fit', '_opnav.fit')
        
        file = t_group[int(k)]['Filename']
        file_short = t_group[int(k)]['Shortname']
        im = hbt.read_lorri(file) # Read the image, and process it a bit I think
        print(f'Reading image {index_group}/{index_image}: {file_short}')
        hdulist = fits.open(file) 
        header  = hdulist['PRIMARY'].header
        mode    = header['SFORMAT']
        exptime = header['EXPTIME']
        hdulist.close()           
        
        is_navigated = os.path.isfile(file_out)  # Test for the presence of the _opnav.fits file

# If it's a 4x4 file, it's probably saturated and lots of things don't work. So doing navigation 
# it hopeless. But we don't want to lose track of the file (for the numbering scheme), 
# so we tag it and copy anyhow.
        
        if (mode == '4X4') and (DO_SKIP_4X4):
            print("{}/{}: Skipping OpNav due to 4x4".format(i, np.size(files)))
            print("Copying to {}".format(file_out))
            shutil.copyfile(file, file_out)

        elif (is_navigated and DO_SKIP_NAVIGATED):
            print ("{}/{}: Skipping since already navigated".format(i, np.size(files)))
    
        else:
    
            with warnings.catch_warnings():   # Block the warnings
                warnings.simplefilter("ignore")    
                w_orig  = WCS(file)           # Look up the WCS coordinates for this frame
                w       = WCS(file)           # Get a copy of it, which we'll change
            
            crval_orig = w_orig.wcs.crval
        
    # Do the navigation call
                    
            (w, (dy_pix, dx_pix)) = hbt.navigate_image_stellar(im, w, method = method_opnav,
                                     title = "{}, {} s".format(file.split('/')[-1], repr(exptime)), 
                                     do_plot=do_plot)
        
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
            
            print(f'Wrote file {index_group}/{index_image}: {file_out}')
            
            is_success = True
            
    # Now, if requested, do a check on the new navigation by loading that brand-new FITS file, and navigating *it*

            do_opnav_verify = True
            
            if do_opnav_verify:

                file    = file_out

                with warnings.catch_warnings():   # Block the warnings
                    warnings.simplefilter("ignore")    
                    w       = WCS(file)           # Get a copy of it, which we'll change
                    w_orig       = WCS(file)           # Get a copy of it, which we'll change

                (w, (dy_pix, dx_pix)) = hbt.navigate_image_stellar(im, w, method = method_opnav,
                                         title = "{}, {} s".format(file.split('/')[-1], repr(exptime)), 
                                         do_plot=do_plot)
                
                crval_orig = w_orig.wcs.crval
                print('Checking navigated frame! Answers below should be (0,0)')
                print("Center was at RA {}, Dec {}".format(crval_orig[0], crval_orig[1]))
                print("Center is now at RA {}, Dec {}".format(crval[0], crval[1]))
                print("Determined OpNav offset: dx = {} pix, dy = {} pix".format(dx_pix, dy_pix))


        
            
    else:
        # print("Error!")
        is_success = False
        
    if (DO_INTERACTIVE and is_success):
        i += 1        # Go to next image in master list
    else:
        ii += 1       # Go to next image in the sublist
