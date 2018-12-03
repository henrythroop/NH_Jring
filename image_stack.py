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
from   astroquery.vo_conesearch import conesearch
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
from   wcs_translate_pix import wcs_translate_pix, wcs_zoom

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
        
        # If we are passed a pickle file, then restore the file from pickle, rather than by reading in explicity.
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

        self.file_save      = os.path.join(dir, 'image_stack_n{}.pkl'.format(num_files))
        
        # Initialize the center of this image. The shfits of each image are taken to be relative to this.
        # It could be that we just keep this at zero. Time will tell.
        
        self.shift_x_pix_center = 0
        self.shift_y_pix_center = 0
        
        # Set the internal zoom level, which will be used when flattening
        
        self.zoom = 1
        
        # Set a flag to indicate if flattened or not
        
        self.flattened = False
        
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
        
        
        self.t = Table(  [[],  [],            [],          [],         [],        [],       [],       [],      [],   [], 
                                                    [], [], [], [], [], [], [], [], [], [], [], [] ],
            names=('filename', 'filename_short', 'exptime', 'visitname', 'sapname', 'sapdesc', 'target', 'reqid', 'et',  'utc', 
                  'shift_x_pix', 'shift_y_pix', 'ra_center', 'dec_center', 'angle', 'dx_pix', 'dy_pix', 
                  'pixscale_x_km', 'pixscale_y_km', 'dist_target_km', 'wcs', 'data'),
            dtype = ('U150', 'U50',           'float64', 'U50',      'U50',     'U50',     'U50',    'U50',   'float64', 'U50', 
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
                backplane_radius = hdulist['RADIUS_EQ'].data
                
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
                    
                # Read the WCS coords of this file.
                # XXX Suppress warnings about WCS SIP coords which Buie's files get.
                # However, looks like AstroPy doesn't use proper warning mechanism??
                
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
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
                # *** Need to add backplane data here, as planes[] or backplanes[] ***
            
                self.t.add_row(
                          [file, filename_short, exptime, visitnam, sapname, sapdesc, target, reqid, 
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
            # print(f'File to save: {self.file_save}')
            answer = input(f'Save to pickle file {self.file_save.split("/")[-1]}? ')
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
        
        if (method.upper() == 'WCS_SIMPLE'):  # Do not use this!
            
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
       
        if (method.upper() == 'WCS'):  # This is the normal method
            
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
# Load a pickle file from disk. 
# Running this is just faster than flattening the array explicitly -- it duplicates no functionality
# =============================================================================

    def load_flat(self):
        
        """
        Load a saved .pkl file from disk.
        Running this is just faster than flattening explicity -- it duplicates no functionality.
        """
        print("Reading: " + self.file_flat_save)

        lun = open(self.file_flat_save, 'rb')
#        self = pickle.load(lun)
        (            self.flattened,
                     self.arr_flat,
                     self.wcs_flat,
                     self.num_planes ) = pickle.load(lun)
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
# Save flattened stack into a pickle file.
# =============================================================================

    def save_flat(self):
        
        """
        Save flattened stack into a .pkl file.
        The filename is held internally and user does not need it.
        """

        lun = open(self.file_flat_save, 'wb')

        pickle.dump((self.flattened, 
                     self.arr_flat, 
                     self.wcs_flat,
                     self.num_planes), lun)
        lun.close()
        print("Wrote: " + self.file_flat_save)     
        
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
        
        # XXX bug: hbt.ceilout(0.0) = -1. This is causing a problem!
        
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
# Flatten a stack of images as per the currently set indices
# =============================================================================

    def flatten(self, method='median', zoom=1, do_subpixel=False, do_wrap=False, padding = 'Auto', do_plot=True,
                do_force=False, do_save=False):
        
        """
        Flatten a stack of images as per the currently-set shift values and indices.
        
        Sets the value of self.wcs to reflect the output image.
        
        As per Amanda Zangari's rules for the ORT, the parameters should be chosen to match the central
        value of the stack -- not the first or last elements. **What parameters?? This is just an image
        stack, which is a single image array, and a WCS array.**
        
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
            
        do_force:
            Force explicit flattening from memory, rather than from a .pkl file.

        do_save:
            If set, save the results of this flattening as a .pkl file
            
        """    
        
        
        if self.flattened and (not do_force):
            print('Already flattened!')
            return (self.arr_flat, self.wcs_flat)       # If we try to flatten an already flattened stack, just return
                                                        # the already computed results... unless do_force is set,
                                                        # in which case we manually do it.
        
        # Create the filename to save to, in case we need it
        
        self.file_flat_save      = self.file_save.replace('.pkl', f'_flat_z{zoom}.pkl')
        # os.path.join(dir, f'image_stack_n{num_files}_flat_z{zoom}.pkl')
        
        if (os.path.isfile(self.file_flat_save)) and not(do_force):
            self.load_flat()
            return (self.arr_flat, self.wcs_flat)
        
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
                print(f'Calculated padding: {pad_xy} pixels required')

        else:        
            pad_xy = padding
            
        if not(do_wrap):
            arr = np.zeros( (self.num_planes, (self.dx_pix + pad_xy*2)*zoom, 
                                              (self.dy_pix + pad_xy*2)*zoom) )
        else:    
            arr = np.zeros((self.num_planes, self.dx_pix*zoom, self.dy_pix*zoom))
        
# =============================================================================
#         Now loop and do the actual shifting and zooming of each image
# =============================================================================
        
        for j,w_i in enumerate(w):  # Loop over all the indices. 'j' counts from 0 to N. w_i is each individual index.
            im = self.t['data'][w_i]
            im_expand = scipy.ndimage.zoom(im,zoom) # Expand the image. This increases total flux; keeps DN/px fixed.

            if not(do_wrap):
                pad = np.array( ((pad_xy, pad_xy),(pad_xy, pad_xy)) ) * zoom
                im_expand = np.pad(im_expand, pad, mode = 'constant')   # This pads each axis by a constant (eg, 5 pix)
                
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
  
            # Copy the shifted, zoomed image into place in the output array
            
            arr[j,:,:] = im_expand              
        
# =============================================================================
#       Flatten all the individual frames down to one, using mean or median
# =============================================================================
        
        print(f'Flattening array with dimension {np.shape(arr)} using {method}; zoom={zoom}')
        print(f'  Pad_xy = {pad_xy}')
        
        if (method == 'mean'):
            arr_flat   = np.nanmean(arr,0)    # Fast
            
        if (method == 'median'):
            arr_flat = np.nanmedian(arr,0)  # Slow -- about 15x longer

# =============================================================================
#      Modify the WCS structure to fit the new image that has been created
#      We want to create one WCS which matches the final single output image --
#      *not* a WCS for each plane.            
# =============================================================================
     
        # Copy the WCS over. Start with the one from the 0th image (which is arbitrary),
        # and we modify it to be correct.
        
        wcs = self.t[0]['wcs'].deepcopy()
        
        # Plot Plane 0 to verify accuracy of its WCS
    
        if do_plot:
            plot_img_wcs(self.t['data'][0], wcs, title = 'Original, plane 0',
                         name_target='MU69', et = self.t['et'][0], name_observer='New Horizons') 
               # XXX I want to plot the stack name here!
        
#        print(f'WCS for above: {wcs}')
        
        # Adjust the WCS to reflect that of the flattened image.
        
        # Change the center position, in RA/Dec, by applying the known offset in pixels
        # Offset is in raw pixels (ie, not zoomed)
        # Take shift amount from 0th image, because that is where I grabbed the WCS from.
        
        dx_wcs_pix =  self.t['shift_pix_x'][0]  # Experiment with these to get them right. -- probably wrong as is. XXX
        dy_wcs_pix =  self.t['shift_pix_y'][0]
        
        print('WCS before shifting')
        print(f'Shift_pix_x = {dx_wcs_pix} pixels horizontal')
        print(f'Shift_pix_y = {dy_wcs_pix} pixels vertical')
        print(f'Pad_xy      = {pad_xy}     pixels vertical and pixels horizontal')
        print()
        print(f'WCS original: {wcs}')  
         
        # print(f'Calling wcs_translate_pix({pad_xy:.1f}, {pad_xy:.1f}) for pad only')
        # wcs_pad = wcs.deepcopy()

        wcs_translate_pix(wcs, pad_xy, pad_xy)
        
        # plot_img_wcs(arr_flat, wcs_pad, title = 'pad only') 
        # print(f'WCS_PAD: {wcs_pad}')  

        print(f'Calling wcs_translate_pix({dx_wcs_pix:.1f}, {dy_wcs_pix:.1f})')
        wcs_translate_pix(wcs, dx_wcs_pix, dy_wcs_pix) 
        
        # Zoom the WCS to match the already zoomed image
        
        wcs_zoom(wcs, zoom, np.array(np.shape(im)))
        plot_img_wcs(arr_flat, wcs, title = 'pad + trans + zoom',
                     name_target='MU69', et = self.t['et'][0], name_observer='New Horizons') 

        # And return the WCS along with the image
        
        if do_plot:
            plot_img_wcs(arr_flat, wcs, title = f'After zoom x{zoom} + adjust',
                         name_target='MU69', et = self.t['et'][0], name_observer='New Horizons')  # XXX data OK, but position of red
                                                               # point on plot is sometimes incorrect.

#        print(f'WCS after above: {wcs}')
        
        self.flattened = True
        
        self.arr_flat = arr_flat
        self.wcs_flat = wcs

        # Prompt to save the flattened image, if requested
        
        if not(do_save):
            print(f'File to save: {self.file_flat_save}')
            answer = input(f'Save to pickle file {self.file_flat_save}? ')
            if ('y' in answer):
                do_save = True
                
        if do_save:
            self.save_flat()   
            
        return (arr_flat, wcs)

# =============================================================================
# Return the central index of an image stack (that is, the index of the image closest to the median ET)
# =============================================================================

    @property
    def index_central(self):
        """
        Return the central index of a image stack.
        """
        
        ets = sorted(self.t['et'])
        
        med = np.median(ets)
        
        index = np.where(ets >= med)[0]   # Median value is not always one of the values itself -- med([9,10]) = 9.5
                                          # So, return a close one instead of an exact match. 
        return index[0]
    
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
        
        print(f'image_single: retrieved pad_xy = {pad_xy}')
        
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

    # Demo function to try out the stack functionality

    plt.set_cmap('Greys_r')
    
    hbt.figsize((10,10))
    
    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
        
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) 
    
    # Load a stack
    
    do_force = False
    
    # dir = '/Users/throop/Data/ORT2/throop/backplaned/K1LR_HAZ03/'
    # dir = '/Users/throop/Data/ORT2/throop/backplaned/K1LR_HAZ03/'
#    stack = image_stack(dir, do_force=True, do_save=True)
    
    dir = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018301/'
    stack = image_stack(dir, do_force=do_force)
    
    # Plot a set of images from it
    
    i = 1
    num_images = 4
    for i in range(num_images):
        plt.subplot(2,2,i+1)
        plt.imshow(stretch(stack.image_single(i)), origin='lower')
        plt.title(f"i={i}, et = {stack.t['et'][i]}")
    plt.show()

    # Align the stack.
    # The values for the shifts (in pixels) are calculated based 
    # on putting MU69's ra/dec in the center of each frame.
    # This assumes that MU69 does not move. And to a very good approximation, that is true.
    # Total MU69 motion 15-Aug .. 15-Dec < 1 pixel.
    
    ra_mu69  = 274.73344   # Value for 2-Nov-2018. Matches plot_img_wcs.py
    dec_mu69 = -20.86170
    
    radec_mu69 = (ra_mu69*hbt.d2r, dec_mu69*hbt.d2r)
    
    stack.align(method = 'WCS', center = radec_mu69)

    # Flatten the stack
    
    stack.padding = 65
    zoom = 2
    (arr_flat, wcs_flat) = stack.flatten(zoom=zoom, do_force=True, do_plot=True, do_save=True)
    
    plot_img_wcs(arr_flat, wcs_flat, title = f'Zoom {zoom}', width=50,
                 name_observer='New Horizons', name_target='MU69', et = stack.t['et'][0])
        
    # Make a plot of ET, showing it for each image and the whole stack
    # Highlight the central time in red.
    
    do_plot_et = False
    
    if do_plot_et:
        t0 = np.amin(stack.t['et'])
        plt.plot(stack.t['et'] - t0, marker = 'o')
        plt.xlabel('Image number')
        plt.ylabel('dET [sec]')
        plt.plot(stack.index_central, stack.t['et'][stack.index_central] - t0, marker = 'o', color = 'red')
        plt.show()    
    
# =============================================================================
#     Now try out stack alignment
# =============================================================================

    dir1 = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018284/'
    stack1 = image_stack(dir1, nmax=3)
    
    dir3 = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018301/'
    stack3 = image_stack(dir3)

    radec_mu69 = (4.794979838984583, -0.3641418801015417)    
    stack1.align(method = 'WCS', center = radec_mu69)  # This sets values stack1.t[0]['shift_pix_x']
    stack3.align(method = 'WCS', center = radec_mu69)
    
    print("Stack 1: padding = {}".format(stack1.calc_padding()[0]))
    print("Stack 1: shift_pix_xy = {}, {}".format(stack1.t[0]['shift_pix_x'], stack1.t[0]['shift_pix_y']))
    
    zoom = 2
    
    hbt.figsize((12,12))
    
    # Plot individual frames. These work.
    
    plot_img_wcs(stretch(stack1.t[0]['data']), stack1.t[0]['wcs'])
    plot_img_wcs(stretch(stack1.t[1]['data']), stack1.t[1]['wcs'])
    
    # Flatten the frames
    
    (arr_flat_1, wcs1) = stack1.flatten(zoom = 1, do_force=True, do_save=True)  
    (arr_flat_3, wcs3) = stack3.flatten(zoom = 1, do_force=True, do_save=True)
    (arr_flat_3, wcs3) = stack3.flatten(zoom = 1, do_force=True, do_save=True)
    
    # Now verify that WCS remains OK
    
    plot_img_wcs(arr_flat_1, wcs1, title = f'After zoom x{zoom} + adjust', width=50)
    plot_img_wcs(arr_flat_3, wcs3, title = f'After zoom x{zoom} + adjust', width=200)
    
# =============================================================================
# Now for debugging the small offsets, process some actual MU69 OpNav data.
# =============================================================================

    ra_mu69  = 274.73344   # Value for 2-Nov-2018. Matches plot_img_wcs.py
    dec_mu69 = -20.86170
    
    radec_mu69 = (ra_mu69*hbt.d2r, dec_mu69*hbt.d2r)
    
    dir1 = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018311/'
    stack1 = image_stack(dir1)  
                
    # Plot each individual frame. Check against each individual WCS.
    # Concl: These are correct. MU69 is in exactly the right location at each one.
    # It is just in the UL corner of the blob. 
    # Also, MU69 should not be at the center. I haven't done the shifts for that yet.

    hbt.figsize((12,12))
    
    for i in range(stack1.size[0]):
        plot_img_wcs(stretch(stack1.t[i]['data']), stack1.t[i]['wcs'], width=100)
    
    # Calculate how to align the frames

    stack1.align(method = 'WCS', center = radec_mu69)  # This sets values stack1.t[0]['shift_pix_x']

    # print("Stack 1: padding = {}".format(stack1.calc_padding()[0]))
    # print("Stack 1: shift_pix_xy = {}, {}".format(stack1.t[0]['shift_pix_x'], stack1.t[0]['shift_pix_y']))

    # Flatten

    zoom = 5
     
    (arr_flat_1, wcs1) = stack1.flatten(zoom = zoom, do_force=True, do_save=True)  

    # Plot the result after align and flatten
    # Concl: Yes, this looks definitely wrong. Off by ~0 pixels at zoom1, but a few at zoom4.
    # The *array* does look properly centered on MU69 (ie, MU69 is in the center of the plot). 
    # But the WCS red dot is shifted from where it should be.
    
    plot_img_wcs(stretch(arr_flat_1), wcs1,            title=f'Stack 1, zoom={zoom}')
    plot_img_wcs(stretch(arr_flat_1), wcs1, width=100, title=f'Stack 1, zoom={zoom}')
    plot_img_wcs(stretch(arr_flat_1), wcs1, width=30,  title=f'Stack 1, zoom={zoom}')
    
    # Flatten the frames
    
    (arr_flat_3, wcs3) = stack3.flatten(zoom = 1, do_force=True, do_save=True)
    (arr_flat_3, wcs3) = stack3.flatten(zoom = 1, do_force=True, do_save=True)
    
    # Now verify that WCS remains OK
    
    plot_img_wcs(arr_flat_1, wcs1, title = f'After zoom x{zoom} + adjust', width=50)
    plot_img_wcs(arr_flat_3, wcs3, title = f'After zoom x{zoom} + adjust', width=200)
