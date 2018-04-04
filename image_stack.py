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
from   get_radial_profile_circular import get_radial_profile_circular

from   get_radial_profile_circular import get_radial_profile_circular
from   plot_img_wcs import plot_img_wcs
from   wcs_translate_pix import wcs_translate_pix

# HBT imports

import hbt

class image_stack:

    """ 
    This class stacks images. It takes a list of files, and it does various processing on them.
    
    Written for NH Hazard imaging.
    
    """
   
    def __init__(self, dir, do_force=False, do_verbose=False, nmax=None, prefix='lor', do_save=False) :   

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
        prefix:
            A text prefix which each filename must match (e.g., 'lor' to match only LORRI files).
        do_save:
            If set, save the results of this stacking as a .pkl file
            
        """
        
        # If we are passed a pickel file, then restore the file from pickle, rather than by reading in explicity.
        # I haven't tested this -- not sure if it works.
        
        if 'pkl' in dir:
            self.file_save = dir
            self.load()
            return
        
        name_target = 'MU69'
            
        do_lorri_destripe = True  # I didn't use this at first, but it is a clear improvement.
        
        files1 = glob.glob(os.path.join(dir,      prefix + '*.fit*'))     # Look in dir
        files2 = glob.glob(os.path.join(dir, '*', prefix + '*.fit*'))     # Look in subdirs
        
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
                
        # Set up the table 't'. This is an astropy table within the stack, that has a list of all of the 
        # useful image parameters taken from the FITS header.
        
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
                                                    [], [], [], [], [], [], [], [], [], [], [], [] ],
            names=('filename_short', 'exptime', 'visitname', 'sapname', 'sapdesc', 'target', 'reqid', 'et',  'utc', 
                  'shift_x_pix', 'shift_y_pix', 'ra_center', 'dec_center', 'angle', 'dx_pix', 'dy_pix', 
                  'pixscale_x_km', 'pixscale_y_km', 'dist_target_km', 'wcs', 'data'),
            dtype = ('U50',           'float64', 'U50',      'U50',     'U50',     'U50',    'U50',   'float64', 'U50', 
                    'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 
                    'float64', 'float64', 'float64', 'float64', 'object', 'object'  ))
        
        if (len(files)):
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
            
                # Destripe if requested (aka remove jailbars)
                
                if do_lorri_destripe:
                    arr = hbt.lorri_destripe(arr)
                    
                # Read the WCS coords of this file
                
                w = WCS(file)
                
                # Get the RA/Dec location of the central pixel, in radians
                
                ra_center  = hdulist[0].header['CRVAL1']*hbt.d2r
                dec_center = hdulist[0].header['CRVAL2']*hbt.d2r
                
                pixscale_x = abs(hdulist[0].header['CD1_1']*hbt.d2r)  # radians per 4x4 pixel
                pixscale_y = abs(hdulist[0].header['CD2_2']*hbt.d2r)  # Sometimes this is negative?
                
                # Initialize the shift amount for this image, in pixels
                
                shift_x_pix = 0
                shift_y_pix = 0
                
                # Calc the distance to MU69, and the pixel scale (non-zoomed)
                
                (st, lt) = sp.spkezr(name_target, et, 'J2000', 'LT', 'New Horizons')
                dist_target_km = sp.vnorm(st[0:2])
                pixscale_x_km = dist_target_km * pixscale_x
                pixscale_y_km = dist_target_km * pixscale_y
                
                # Load the values for this image into a row of the astropy table
            
                self.t.add_row(
                          [filename_short, exptime, visitnam, sapname, sapdesc, target, reqid, 
                           et, utc, 
                           shift_x_pix, shift_y_pix,
                           ra_center,   
                           dec_center,  
                           angle,
                           dx_pix, dy_pix, # X and Y dimensions
                           pixscale_x_km, pixscale_y_km,
                           dist_target_km,
                           w,              # WCS object 
                           arr])           # Actual imaging data 
            
            # End of loop over files
        else:
            print(f"No files found in {dir}!")
            return
        
#        if do_verbose:
        
        print("\n") # Print CR at end of "...."
            
        # Sort by ET.
            
        self.t.sort('et')
        
        # Save the pixel scale, from most recent image. We assume that the pixel scale of all frames is identical.
        
        self.pixscale_x = pixscale_x
        self.pixscale_y = pixscale_y
        
        self.pixscale_x_km = pixscale_x_km
        self.pixscale_y_km = pixscale_y_km
        
        self.dist_target_km   = dist_target_km
        
        self.et               = et 
        
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
        
        # Initialize the num_planes vector to set the size
        
        self.num_planes = len(files)
    
        # If we generated the files manually (not by reloading), and flag is not set, then offer to save them
        
        if not(do_save):
            answer = input("Save to pickle file? ")
            if ('y' in answer):
                do_save = True
                
        if do_save:
            self.save()            
            
        # Return. Looks like an init method should not return anything.
    
        print()

## =============================================================================
## Return size of the stack
## =============================================================================
#        
#        
#    def size(self):
#        """
#        Return size of the stack.
#        
#        Q: Why do I need this? Aren't a method's variables exposed to the outside world?
#        I thought that was the point of not needing getter/setters in python.
#        """
#        
#        return self.size
    
# =============================================================================
# Calculate alignments between all frames
# =============================================================================

    def align(self, method = 'wcs', center = None):
        """
        Take a loaded stack, and set the shift amounts for each image.
        
        Does not actually perform the shifts, but updates internal shift settings
        for each image. The shifting only happens when image is flattened.
        
        Updates the internal shift_x_pix, shift_y_pix, accessible thru self.t[<number>]['shift_pix_x'].
        These pixel shifts are raw, and do not reflect any zoom.
        
        Optional arguments
        ----
        method:
            'wcs' : Center based on full WCS object. This should usually be used.

            'wcs_simple' : Center based on RA / Dec positions in WCS. Does not use full WCS -- just CRVAL.
            
            'body'  : Center based on a given body name. Not yet implemented.
            
        center:
            Tuple or string specifying the center to use. For instance, center = (ra_radians, dec_radians)
            
        """
        
        if (method.upper() == 'WCS_SIMPLE'):
            
            # If 'radec' is passed, then set the shifts s.t. the right amount is applied to put 
            # specified RA/Dec in the center.
            
            # Extract the ra and dec in radians. These specify the desired RA and dec.
            
            (ra, dec) = center
                        
            shift_x_pix = self.t['shift_x_pix']
            shift_y_pix = self.t['shift_y_pix']
            
            ra_center   = self.t['ra_center']   # RA of each individual image
            dec_center  = self.t['dec_center']
            
            num_images = self.size[0]
            
            for i in range(num_images):
                shift_x_pix[i] = (ra_center[i]  -  ra)/self.pixscale_x
                shift_y_pix[i] = (dec_center[i] - dec)/self.pixscale_y
            
            # Save these values back to the table, where they live
            
            self.t['shift_pix_x'] = shift_x_pix
            self.t['shift_pix_y'] = shift_y_pix
       
        if (method.upper() == 'WCS'):
            
            # If 'WCS' is passed, do the same as 'radec', but use the whole WCS object explicitly.
            # This should give very nearly the same results.
            
            (ra, dec) = center
            
            num_images = self.size[0]
            
            shift_x_pix = self.t['shift_x_pix']
            shift_y_pix = self.t['shift_y_pix']
            
            for i in range(num_images):
                wcs = self.t['wcs'][i]  # Get WCS for this image
                                
                (pos_pix_x, pos_pix_y) = wcs.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)
                
                shift_x_pix[i] = self.dx_pix/2 - pos_pix_x
                shift_y_pix[i] = self.dy_pix/2 - pos_pix_y
           
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
        (self.t, self.indices, 
                     self.et,
                     self.pixscale_x_km,
                     self.pixscale_y_km,
                     self.pixscale_x,
                     self.pixscale_y,
                     self.dist_target_km,
                     self.dx_pix, self.dy_pix, self.size) =\
          pickle.load(lun)
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
        pickle.dump((self.t, self.indices, 
                     self.et,
                     self.pixscale_x_km,
                     self.pixscale_y_km,
                     self.pixscale_x,
                     self.pixscale_y,
                     self.dist_target_km,
                     self.dx_pix, self.dy_pix, self.size), lun)
        lun.close()
        print("Wrote: " + self.file_save)     
        
# =============================================================================
# Print the table
# =============================================================================

    def print(self):
        """
        Print the table. This lists all of the tabulated information about each
        image in the stack. 
        
        This function is only moderately useful, because the wcs and data objects
        do not print well in a table. It would be more useful if these columns were omitted.
        """
        
        self.t.pprint(max_lines = -1, max_width = -1)
        
# =============================================================================
# Set the indices to extract
# =============================================================================

    def set_indices(self, indices):
        """
        Set the indices of planes to use when flattening.
        
        Parameters
        ----
        indices:
            Boolean array, one per plane. Indicates whether this plane will be used during flatten().
        """
        
        self.indices = indices

# =============================================================================
# Calculate the amount of padding required for an image
# =============================================================================
        
    def calc_padding(self):

        """
        Returns the padding required for an image, based on the shifts of individual frames.
    
        return: (max, ((mix_x, max_x), (min_y, max_y)))
        """
        
        shift_x_pix = self.t['shift_x_pix']
        shift_y_pix = self.t['shift_y_pix']
        
        # Calculate the total shift range, from negative to positive.
        # Use ceil() to assure that we have enough headroom for rounding.
        
        shift_x_pix_min = hbt.ceilout(np.amin(shift_x_pix))
        shift_x_pix_max = hbt.ceilout(np.amax(shift_x_pix))
        shift_y_pix_min = hbt.ceilout(np.amin(shift_y_pix))
        shift_y_pix_max = hbt.ceilout(np.amax(shift_y_pix))
        
        # Calculate the padding on each edge to add
        pad_xy = np.amax(np.abs([shift_x_pix_min, shift_x_pix_max,
                                 shift_y_pix_min, shift_y_pix_max]))

        return (pad_xy, 
                ((shift_x_pix_min, shift_x_pix_max), 
                 (shift_y_pix_min, shift_y_pix_max)))
        
# =============================================================================
# Flatten a stack of images as per the currently-set indices
# =============================================================================

    def flatten(self, method='median', zoom=1, do_subpixel=False, do_wrap=False, padding = 'Auto'):
        
        """
        Flatten a stack of images as per the currently-set shift values and indices.
        
        Sets the value of self.wcs to reflect the output image.
        
        As per Amanda Zangari's rules for the ORT, the parameters should be chosen to match the central
        value of the stack -- not the first or last elements.
        
        Return values
        ----
        
        img_flat:
            A 2D flattened image array. Not an object, just a simple 2D array. This image is padded.
        
        wcs:
            The WCS parameters for this image.
            
        Optional parameters
        ----
        method: 
            'median': Slower but better. Recommended.

            'mean':   Faster but doesn't remove cosmic rays well.
                        
        zoom:
            Value to scale images by. Zoom is applied before the shift is performed. The output image
            retains the zoomed values. Zoom=4 → each pixel turns into 4x4=16 pixels.
            
        do_subpixel:
            Boolean. For shifts, use an FFT method which simulates sub-pixel shifts. This is OK for experimentation,
            but leaves a lot of ringing artifacts. Better to scale the output using zoom parameter.
            
        do_wrap:
            Boolean. If True, then the flattened output array will be the same size as the input, with each plane rolled
            s.t. they wrap around the edge. If False, the output array is padded so no wrapping occurs.
            
        padding:
            Amount of padding to add to edge. If 'Auto', then it will be calculated automatically based on the actual 
            shifts in this stack. However, this may be larger than necessary, *and* the shifts might not match other
            stacks.
            
            If int, then a fixed number of bins is added to each edge (R, L, T, B). Padding is added *before* the zoom.
            (Padding=10, Zoom=4) will increase the size of the output array by 80 pixels in each direction.
            
        """    
        
        self.num_planes = np.sum(self.indices)

        # Create a 3D output array for all of the images to go into. 
        # This will temporarily store the shifted+zoomed frames.
        
        w = np.where(self.indices)[0]  # List of all the indices
        
        shift_x_pix = self.t['shift_x_pix']
        shift_y_pix = self.t['shift_y_pix']

        # Calculate the padding on each edge to add

        if type(padding) == str:
            if (padding == 'Auto'):
                pad_xy = (self.calc_padding())[0]

        else:        
            pad_xy = padding
            
        if not(do_wrap):
            arr = np.zeros( (self.num_planes, (self.dx_pix + pad_xy*2)*zoom, 
                                              (self.dy_pix + pad_xy*2)*zoom) )
        else:    
            arr = np.zeros((self.num_planes, self.dx_pix*zoom, self.dy_pix*zoom))
        
        for j,w_i in enumerate(w):  # Loop over all the indices. 'j' counts from 0 to N. w_i is each individual index.
            im = self.t['data'][w_i]
            im_expand = scipy.ndimage.zoom(im,zoom) # Expand the image. This increases total flux; keeps DN/px fixed.

            if not(do_wrap):
                pad = np.array( ((pad_xy, pad_xy),(pad_xy, pad_xy)) ) * zoom
                im_expand = np.pad(im_expand, pad, mode = 'constant')   # This adds a const
                
            # Look up the shift amount, and zoom it as necessary
            
            dx = shift_x_pix[w_i]
            dy = shift_y_pix[w_i]
            shift = (dy*zoom, dx*zoom) # For roll(), axis=0 is y, and axis=1 is x! Opp. from normal, but easy to verify.
        
            # Apply the proper shift in X and Y. What I am calling 'x' and 'y'             
            # Sub-pixel shifting is in theory better. But in reality it makes a trivial difference over
            # integer shifting, and leaves a lot of ringing artifacts.
            
            # Apply the sub-pixel shifting. I tried for hours to get scipy.ndimage.shift() to work, but 
            # it kept being screwy. fourier_shift works fine, though slightly slower.
            
            if do_subpixel:
                im_expand = np.fft.ifft2(scipy.ndimage.fourier_shift(np.fft.fft2(im_expand.copy()),shift)).real
            
            else:
                im_expand = np.roll(np.roll(im_expand, int(round(shift[0])), axis=0), int(round(shift[1])), axis=1)
#                print(f'Adding image {j} with shifts dx {int(round(shift[0]))}, dy {int(round(shift[1]))}')
  
            # Copy the shifted, zoomed image into place
            
            arr[j,:,:] = im_expand              
        
        # Merge all the individual frames, using mean or median
        
        print(f'Flattening array with dimension {np.shape(arr)} using {method}; zoom={zoom}')
        print(f'  Pad_xy = {pad_xy}')
        
        if (method == 'mean'):
            arr_flat   = np.nanmean(arr,0)    # Fast
            
        if (method == 'median'):
            arr_flat = np.nanmedian(arr,0)  # Slow -- about 15x longer
        
        # Copy the WCS over. Use the one from the 0th image -- which is arbitrary, but as good as any.
        
        wcs = self.t[0]['wcs']
        
        # Plot Plane 0 to verify accuracy of its WCS
        
        plot_img_wcs(self.t['data'][0], wcs, title = 'Original, plane 0')
        
        print(f'WCS for above: {wcs}')
        
        # Adjust the WCS 
        
        # Change the center position, in RA/Dec, by applying the known offset in pixels
        # Offset is in raw pixels (ie, not zoomed)
        # Take shift amount from 0th image.
        
        dx_wcs_pix =  self.t['shift_pix_x'][0]  # Experiment with these to get them right.
        dy_wcs_pix =  self.t['shift_pix_y'][0]
        
        print(f'Calling wcs_translate_pix({-dy_wcs_pix} + {pad_xy}; {-dx_wcs_pix} + {pad_xy})')
        
        wcs_translate_pix(wcs, dy_wcs_pix + pad_xy, dx_wcs_pix + pad_xy)

        print(f'WCS after above: {wcs}')
        
        # Change the image size, in the WCS, to reflect the larger image size. This is just the CRPIX value
        
        wcs.wcs.crpix = [(hbt.sizex(im_expand)-1)/2, (hbt.sizey(im_expand)-1)/2]
        
        # Change the pixel scale to reflect the zoom
        
        wcs.wcs.pc = wcs.wcs.pc / zoom
        
        # And return the WCS along with the image

        plot_img_wcs(arr_flat, wcs, title = f'After zoom x{zoom} + adjust')

        print(f'WCS after above: {wcs}')
        
        return (arr_flat, wcs)

# =============================================================================
# Align and retrieve a single image from the stack.
# =============================================================================

    def image_single(self, num, zoom=1, do_wrap=False, padding = 'Auto'):

        """
        Return a single plane of an image stack, according to a given zoom and wrap.
        """
        
        if type(padding) == str:
            if (padding == 'Auto'):
                pad_xy = (self.calc_padding())[0]

        else:        
            pad_xy = padding

        if (do_wrap):
            raise(ValueError('do_wrap = True not supported'))
            
        # Get the image
            
        im = self.t['data'][num]
        
        # Zoom the image
        
        im_expand = scipy.ndimage.zoom(im,zoom) # Expand the image. This increases total flux; keeps DN/px fixed.

        # Pad the image 
        
        pad = np.array( ((pad_xy, pad_xy),(pad_xy, pad_xy)) ) * zoom
        im_expand = np.pad(im_expand, pad, mode = 'constant')

        # Shift the image

        dx = self.t['shift_x_pix'][num]
        dy = self.t['shift_y_pix'][num]
        shift = (dy*zoom, dx*zoom)
    
        im_expand = np.roll(np.roll(im_expand, int(round(shift[0])), axis=0), int(round(shift[1])), axis=1)

        return(im_expand)

    
# =============================================================================
# End of method definition
# =============================================================================

def dn2iof(val, mode):

    if (mode == '4X4'):
        # Convert DN values in array, to I/F values
            
        profile = val
    
def wcs_translate_radec(wcs, dra, ddec):
    pass

# =============================================================================
# Define a test function which is called when we run this file. This is just an example of using the class.
# =============================================================================
    
if (__name__ == '__main__'):

    # Restore a bunch of stacks from a pickle file. This is not the usual way to to this.
    
#    file = '/Users/throop/Data/ORT2/throop/backplaned/stacks_blink_ORT2_n5_z4.pkl'
#    print("Reading: " + file)           
#    lun = open(file, 'rb')
#    (stack_field, img_field, stack_haz, img_haz, wcs_haz) = pickle.load(lun)
#    lun.close()
    
    # Load a stack
    
    dir = '/Users/throop/Data/ORT2/throop/backplaned/K1LR_HAZ03/'
    stack = image_stack(dir)
    
    # Plot a set of images from it
    
    i = 1
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.imshow(stretch(stack.image_single(i)))
        plt.title(f"i={i}, et = {stack.t['et'][i]}")
    plt.show()

    # Align the stack
    
    radec_mu69 = (4.794979838984583, -0.3641418801015417)    
    stack.align(method = 'WCS', center = radec_mu69)

    # Flatten the stack
    
    stack.padding = 55
    zoom = 2
    (arr_flat, wcs) = stack.flatten(zoom=zoom)
    
    
    # Flatten it
    
    stack_flat = stack.flatten()
    
#    dir = /Users/throop/Data/ORT2/throop/backplaned/'
    
#            lun = open(self.file_save, 'rb')
#        self = pickle.load(lun)
#        (self.t, self.indices, 
#                     self.et,
#                     self.pixscale_x_km,
#                     self.pixscale_y_km,
#                     self.pixscale_x,
#                     self.pixscale_y,
#                     self.dist_target_km,
#                     self.dx_pix, self.dy_pix, self.size) =\
#          pickle.load(lun)
#        lun.close() 



#    s = image_stack(dir)
#        /Users/throop/Data/ORT2/throop/backplaned/stacks_blink_ORT2_n5_z4.pkl            
#
#    stacks_blink_ORT2_n5_z4.pkl
#    dir = /Users/throop/Data/ORT2/throop/backplaned/stacks_blink_ORT2_n5_z4.pkl
