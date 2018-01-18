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
import imreg_dft as ird                    # Image translation

import re # Regexp
import pickle # For load/save

import scipy

from   matplotlib.figure import Figure

# HBT imports

import hbt

class image_stack:

    """ 
    This class stacks images. It takes a list of files, and it does various processing on them.
    
    Written for NH Hazard imaging.
    
    """
   
    def __init__(self, dir, do_force=False, do_verbose=False, nmax=None) :   

        """
        Init method: load the index and all files.
        
        This does not align, register, stack, or anything like that. It just loads the images.
        
        If a saved .pkl file is available, then it will be used instead of reading individual images.
        
        Parameters
        ----
        dir:
            The directory with all the files
            
        Optional keyword parameters
        ----
        do_force:
            Force loading stack from original files, rather than from a .pkl file.
        do_verbose:
            List each file explicitly when reading
            
        """
                
        files1 = glob.glob(os.path.join(dir,      '*.fit*'))     # Look in dir
        files2 = glob.glob(os.path.join(dir, '*', '*.fit*'))     # Look in subdirs
        
        files = files1 + files2
        
        # Truncate the list, if requested
        
        if (nmax):
            files = files[0:nmax]
            
        num_files = len(files)

        self.file_save = os.path.join(dir, 'image_stack_n{}.pkl'.format(num_files))
        
        # Initialize the center of this image. The shfits of each image are taken to be relative to this.
        # It could be that we just keep this at zero. Time will tell.
        
        self.shift_x_pix_center = 0
        self.shift_y_pix_center = 0
        
        # Set the internal zoom level, which will be used when flattening
        
        self.zoom = 1
        
        # If a save file exists, then load it, and immediately return
        
        if (os.path.isfile(self.file_save)) and not(do_force):
            self.load()
            return
        
        mode     = []
        exptime  = []
        filename_short = []
        exptime  = []
        visitnam = []
        sapname  = []
        sapdesc  = []
        reqid    = []
        et       = []
        utc      = []
        target   = []
        
        # Set up a dictionary for the data itself
        
        self.data     = {}
        
        # Set up the table. 
        # The standard method of representing Python 3 strings in numpy is via the unicode 'U' dtype."
        # Disadvantage of this is memory space, but not an issue for me.
        
        # Fields in the table are:
        #   filename_short
        #   exptime
        #   visitname
        #   sapname
        #   sapdesc
        #   target
        #   reqid
        #   et
        #   utc
        #   shift_x_pix -- the shift of this image, relative to the zero point (tbd)
        #   shift_y_pix -- the shift of this image, relative to the zero point (tbd)
        #   ra    -- 
        #   dec
        #   angle
        #   dx_pix -- x dimension
        #   dy_pix -- y dimension
        
        
        self.t = Table(  [[],              [],          [],         [],        [],       [],       [],      [],   [], 
                                                    [], [], [], [], [], [] ],
            names=('filename_short', 'exptime', 'visitname', 'sapname', 'sapdesc', 'target', 'reqid', 'et',  'utc', 
                              'shift_x_pix', 'shift_y_pix', 'ra_center', 'dec_center', 'angle', 'dx_pix'),
            dtype = ('U50',           'float64', 'U50',      'U50',     'U50',     'U50',    'U50',   'float64', 'U50', 
                                                    'float64', 'float64', 'float64', 'float64', 'float64', 'float64'  ))
        
        print("Reading {} files from {}".format(len(files), dir))
        
        for i,file in enumerate(files):
        
            # Read the data
        
            hdulist        = fits.open(file)
            arr            = hdulist[0].data
            err            = hdulist[1].data
            quality        = hdulist[2].data
            
            dx_pix         = hbt.sizex(arr)   # Usually 256 or 1024
            dy_pix         = hbt.sizey(arr)
            
            
            filename_short = os.path.basename(file).replace('.fits', '').replace('lor_', '')\
                     .replace('_0x633_pwcs','')\
                     .replace('_0x630_pwcs','')
            exptime = hdulist[0].header['EXPTIME']
            visitnam= hdulist[0].header['VISITNAM']
            sapname = hdulist[0].header['SAPNAME']
            sapdesc = hdulist[0].header['SAPDESC']
            target  = hdulist[0].header['TARGET']
            reqid   = hdulist[0].header['REQID']
            et      = hdulist[0].header['SPCSCET']
            angle   = hdulist[0].header['SPCEMEN']*hbt.d2r # Boresight roll angle
            utc     = sp.et2utc(et, 'C', 1)
        
            if do_verbose:
                print("Read {}/{} {}".format(i, len(files), filename_short))
            else:
                print(".", end="")

            hdulist.close()
        
            # Read the WCS coords
            
            w = WCS(file)
            
            # Get the RA/Dec location of the central pixel, in radians
            
            ra_center  = hdulist[0].header['CRVAL1']*hbt.d2r
            dec_center = hdulist[0].header['CRVAL2']*hbt.d2r
            
            pixscale_x = hdulist[0].header['CD1_1']*hbt.d2r  # radians per pixel
            pixscale_y = hdulist[0].header['CD2_2']*hbt.d2r 
            
            # Initialize the shift amount for this image, in pixels
            
            shift_x_pix = 0
            shift_y_pix = 0
            
            
        # Load the values for this image into a row of the astropy table
        
            self.t.add_row(
                      [filename_short, exptime, visitnam, sapname, sapdesc, target, reqid, 
                       et, utc, 
                       shift_x_pix, shift_y_pix,
                       ra_center,   
                       dec_center,  
                       angle,
                       dx_pix])     
        
        
        # Save the image from this file into a dictionary. We can access it later using the filename as a key.
            
            self.data[filename_short] = arr
        
        # End of loop over files
        
        # Sort by ET.
            
        self.t.sort('et')
        
        # Save the pixel scale, from most recent image. We assume that the pixel scale of all frames is identical.
        
        self.pixscale_x = pixscale_x
        self.pixscale_y = pixscale_y
        
        # Save the image size, from most recent image
        
        self.dx_pix     = dx_pix
        self.dy_pix     = dy_pix
        
        # Save the image stack size so it can be easily retrieved
        
        self.size       = (len(files), dx_pix, dy_pix) 
        
    # Finally, remove a few columns that we don't need, or that are wrong.
        
        self.t.remove_column('sapdesc')
        self.t.remove_column('sapname')
        self.t.remove_column('target')
        self.t.remove_column('visitname')

    # Initialize the 'indices' vector, which indicates which planes we use for flattening
    
        self.indices = np.ones(len(self.t), dtype=bool)
        
    # Return. Looks like an init method should not return anything.
    
        print()


# =============================================================================
# Calculate alignments between all frames
# =============================================================================

    def align(self, method = 'radec', center = None):
        """
        Take a loaded stack, and set the shift amounts for each image.
        
        Does not actually perform the shifts, but updates internal shift settings
        for each image. The shifting only happens when image is flattened.
        
        Updates the internal shift_x_pix, shift_y_pix.
        
        Optional arguments
        ----
        method:
            'radec' : Center based on RA / Dec positions in WCS
            
            'body'  : Center based on a given body name
        center:
            Tuple or string specifying the center to use
        """
        
        if (method.upper() == 'RADEC'):
            
            # If 'radec' is passed, then set the shifts s.t. the right amount is applied to put 
            # specified RA/Dec in the center.
            
            # Extract the ra and dec in radians. These specify the desired RA and dec.
            
            (ra, dec) = center
            
            t = self.t
            
            shift_x_pix = self.t['shift_x_pix']
            shift_y_pix = self.t['shift_y_pix']
            
            ra_center   = self.t['ra_center']   # RA of each individual image
            dec_center  = self.t['dec_center']
            
            num_images = len(shift_x_pix)
            
            for i in range(num_images):
                shift_x_pix[i] = (ra_center[i]  -  ra)/self.pixscale_x
                shift_y_pix[i] = (dec_center[i] - dec)/self.pixscale_y
            
            # Save these values back to the table, where they live
            
            self.t['shift_pix_x'] = shift_x_pix
            self.t['shift_pix_y'] = shift_y_pix
            
# =============================================================================
# Load a pickle file from disk. 
# Running this is just faster than reading and processing the files explicitly -- it duplicates no functionality
# =============================================================================

    def load(self):
        
        """
        Load a saved .pkl file from disk.
        Running this is just faster than reading and processing the files explicitly -- it duplicates no functionality.
        """
        print("Reading: " + self.file_save)

        lun = open(self.file_save, 'rb')
#        self = pickle.load(lun)
        (self.t, self.data, self.indices, self.pixscale_x, self.pixscale_y, self.dx_pix, self.dy_pix, self.size) \
          = pickle.load(lun)
        lun.close() 

# =============================================================================
# Save current state into a pickle file.
# =============================================================================

    def save(self):
        
        """
        Save current stack into a .pkl file, which includes contents of all images.
        The filename is held internally and user does not need it.
        """

        lun = open(self.file_save, 'wb')
#        pickle.dump(self, lun)
        pickle.dump((self.t, self.data, self.indices, 
                     self.pixscale_x, self.pixscale_y, self.dx_pix, self.dy_pix, self.size), lun)
        lun.close()
        print("Wrote: " + self.file_save)     
        
# =============================================================================
# Print the table
# =============================================================================

    def print(self):
        """
        Print the table. This lists all of the tabulated information about each
        image in the stack.
        """
        
        self.t.pprint(max_lines = -1, max_width = -1)
        
# =============================================================================
# Set the indices to extract
# =============================================================================

    def set_indices(self, indices):
        self.indices = indices
        
# =============================================================================
# Flatten a stack of images as per the currently-set indices
# =============================================================================

    def flatten(self, method='mean', zoom=1, do_subpixel = False):
        
        """
        Flatten a stack of images as per the currently-set indices.
        
        Optional parameters
        ----
        method: 
            'mean':   Faster but doesn't remove cosmic rays well.
            
            'median': Slower but better
            
        zoom:
            Value to scale output by -- e.g., zoom=4 → each pixel turns into 4x4=16 pixels.     
        """    
        
        self.num_planes = np.sum(self.indices)

        # Create a 3D output array for all of the images to go into. 
        # This will temporarily store the shifted+zoomed frames.
        
        self.arr = np.zeros((self.num_planes, self.dx_pix*zoom, self.dy_pix*zoom))
        
        w = np.where(self.indices)[0]  # List of all the indices
        print(f"Looping over {w} planes")

        for j,w_i in enumerate(w):  # Loop over all the indices. 'j' counts from 0 to N. w_i is each individual index.
            im = self.data[self.t['filename_short'][w_i]]
            im_expand = scipy.ndimage.zoom(im,zoom) # Expand the image. This increases total flux; keeps DN/px fixed.

            dx = self.t['shift_x_pix'][w_i]
            dy = self.t['shift_y_pix'][w_i]
            
            shift = (dy*zoom, dx*zoom)
        
            # Apply the proper shift in X and Y. What I am calling 'x' and 'y'             
            # Sub-pixel shifting is in theory better. But in reality it makes a trivial difference over
            # integer shifting.
            
            # Apply the sub-pixel shifting. I tried for hours to get scipy.ndimage.shift() to work, but 
            # it kept being screwy. fourier_shift works fine, though slightly slower.
            
            if do_subpixel:
                im_expand = np.fft.ifft2(scipy.ndimage.fourier_shift(np.fft.fft2(im_expand.copy()),shift)).real
            
            else:
                im_expand = np.roll(np.roll(im_expand, int(round(shift[0])), axis=0), int(round(shift[1])), axis=1)
  
            # Copy the shifted, zoomed image into place
            
            self.arr[j,:,:] = im_expand              
        
        # Merge all the individual frames, using mean or median
        
        print("Flattening array with dimension {} using {}".format(np.shape(self.arr), method))
        
        if (method == 'mean'):
            self.img_flat   = np.nanmean(self.arr,0)    # Fast
            
        if (method == 'median'):
            self.img_flat = np.nanmedian(self.arr,0)  # Slow -- about 15x longer
            
        return self.img_flat

# =============================================================================
# End of method definition
# =============================================================================

#=============================================================================
# Extract a radial profile. This is the simplest possible case: circular
# ** This is NOT part of the method! It is a standalone function.        
#=============================================================================

def get_radial_profile_circular(arr, pos = (0,0), width=1, method='mean'):
  
    dx = hbt.sizex(arr) 
    dy = hbt.sizey(arr)
    
    xx, yy = np.mgrid[:dx, :dy]  # What is this syntax all about? That is weird.
                                 # A: mgrid is a generator. np.meshgrid is the normal function version.
    
    dist_pix_2d = np.sqrt( ((xx - pos[0]) ** 2) + ((yy - pos[1]) ** 2) )

    dist_pix_1d = hbt.frange(0, int(np.amax(dist_pix_2d)/width))*width
    
    profile_1d    = 0. * dist_pix_1d.copy()
    
#    profile_1d_median  = 0. * dist_pix_1d.copy()
    
    for i in range(len(dist_pix_1d)-2):

    # Identify the pixels which are at the right distance
    
        is_good = np.logical_and(dist_pix_2d >= dist_pix_1d[i],
                                 dist_pix_2d <= dist_pix_1d[i+1]) 
    
        if (method == 'mean'):
            profile_1d[i]   = np.mean(arr[is_good])
    #        profile_1d = profile_1d_mean
    
        if (method == 'median'):
            profile_1d[i] = np.median(arr[is_good])
    #        profile_1d = profile_1d_median
        
    return (dist_pix_1d, profile_1d)

# =============================================================================
# Run the function
# =============================================================================

if (__name__ == '__main__'):
    
#    file_orig ='/Users/throop/Dropbox/Data/ORT1/spencer/ort1/lor_0405175932_0x633_sci_HAZARD_ort1.fit'
#    file_wcs = '/Users/throop/Dropbox/Data/ORT1/porter/pwcs_ort1/K1LR_HAZ00/lor_0405175932_0x633_pwcs.fits'
#    file_bp = '/Users/throop/Dropbox/Data/ORT1/throop/backplaned/K1LR_HAZ00/lor_0405175932_0x633_pwcs_backplaned.fits'
#    
#    hdu_wcs  = fits.open(file_wcs)
#    hdu_orig = fits.open(file_orig)
#    hdu_bp   = fits.open(file_bp)
#    
#    header_wcs  = hdu_wcs[0].header
#    header_orig = hdu_orig[0].header
#    header_bp   = hdu_bp[0].header
#    
#    print(fits.FITSDiff(hdu_wcs, hdu_orig).report())
    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    plt.set_cmap('Greys_r')

    zoom = 4
    
    reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
    reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
    
    dir_data    = '/Users/throop/Data/ORT1/throop/backplaned/'

    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem.tm')
    
    hbt.figsize((12,12))
    
    # Load and stack the field images
    
    stack_field = image_stack(os.path.join(dir_data, reqid_field))
    stack_haz0  = image_stack(os.path.join(dir_data, reqids_haz[0]))
    stack_haz4  = image_stack(os.path.join(dir_data, reqids_haz[4]))
    
    # Set the rough position of MU69, from GV
    
    et = stack_haz0.t['et'][0] # Look up ET for first image in the Hazard stack
    
    (st, lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    vec = st[0:3]
    (_, ra, dec) = sp.recrad(vec)
    
#    radec_mu69 = (274.711*hbt.d2r, -20.867*hbt.d2r)
    radec_mu69 = (ra*hbt.d2r, dec*hbt.d2r)
    
    # Align the frames
    
    stack_field.align(method = 'radec', center = (radec_mu69))
    stack_haz0.align(method  = 'radec', center = (radec_mu69))
    stack_haz4.align(method  = 'radec', center = (radec_mu69))
    
    img_field = stack_field.flatten(do_subpixel=False, method='median',zoom=1)
    img_haz0  = stack_haz0.flatten(do_subpixel=False,  method='median',zoom=1)
    img_haz4  = stack_haz4.flatten(do_subpixel=False,  method='median',zoom=1)

#    img_field_s = stack_field.flatten(do_subpixel=True, method='median',zoom=zoom)
#    img_haz0_s  = stack_haz0.flatten(do_subpixel=True,  method='median',zoom=zoom)
#    img_haz4_s  = stack_haz4.flatten(do_subpixel=True,  method='median',zoom=zoom)

    plt.imshow(stretch(img_field))
    plt.imshow(stretch(img_haz0))
    plt.show()

    diff_s = img_haz4_s - img_field_s
    diff   = img_haz4   - img_field
    
    plt.imshow(stretch( img_haz0 - img_field ))
    plt.imshow(stretch( img_haz4_s - img_field_s ))
    plt.show()
    
    # Manually get the offset
    
    out = ird.translation(img_field, img_haz4)
    tvec = np.round(np.array(out['tvec'])).astype(int)
    print("Detected shift = {} pix".format(out['tvec']))
    
    plt.imshow(stretch(np.roll(np.roll(img_haz4,round(tvec[1]),1),round(tvec[0]),0) - img_field))

#    stack_field.set_indices(hbt.frange(0,10))
    
    # Make a plot of the offset in pixels between frames
    
#    plt.plot(stack_field.t['shift_x_pix'], stack_field.t['shift_y_pix'])
#    plt.xlabel('RA pixels')
#    plt.ylabel('Dec pixels')
#    plt.show()
