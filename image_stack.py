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
                
            # Read the WCS coords
            
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
        
        Updates the internal shift_x_pix, shift_y_pix.
        
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
        (self.t, self.indices, self.pixscale_x, self.pixscale_y, self.dx_pix, self.dy_pix, self.size) \
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
        pickle.dump((self.t, self.indices, 
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
                im_expand = np.pad(im_expand, pad, mode = 'constant')
                
            # Look up the shift amount, and zoom it as necessary
            
            dx = shift_x_pix[w_i]
            dy = shift_y_pix[w_i]
            shift = (dy*zoom, dx*zoom)
        
            # Apply the proper shift in X and Y. What I am calling 'x' and 'y'             
            # Sub-pixel shifting is in theory better. But in reality it makes a trivial difference over
            # integer shifting, and leaves a lot of ringing artifacts.
            
            # Apply the sub-pixel shifting. I tried for hours to get scipy.ndimage.shift() to work, but 
            # it kept being screwy. fourier_shift works fine, though slightly slower.
            
            if do_subpixel:
                im_expand = np.fft.ifft2(scipy.ndimage.fourier_shift(np.fft.fft2(im_expand.copy()),shift)).real
            
            else:
                im_expand = np.roll(np.roll(im_expand, int(round(shift[0])), axis=0), int(round(shift[1])), axis=1)
  
            # Copy the shifted, zoomed image into place
            
            arr[j,:,:] = im_expand              
        
        # Merge all the individual frames, using mean or median
        
        print("Flattening array with dimension {} using {}".format(np.shape(arr), method))
        
        if (method == 'mean'):
            arr_flat   = np.nanmean(arr,0)    # Fast
            
        if (method == 'median'):
            arr_flat = np.nanmedian(arr,0)  # Slow -- about 15x longer
            
        return arr_flat

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

def dn2iof(val, mode):

    if (mode == '4X4'):
        # Convert DN values in array, to I/F values
            
        profile = val
        

        
# =============================================================================
# End of method definition
# =============================================================================

# =============================================================================
# Run the function. This just goes thru the motion of making some stacks.
# It is not necessary to run this, but useful for diagnostics.
# =============================================================================

if (__name__ == '__main__'):
    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    plt.set_cmap('Greys_r')

    zoom = 1      # How much to magnify images by before shifting. 4 (ie, 1x1 expands to 4x4) is typical
                  # 1 is faster; 4 is slower but better.
    
    name_ort = 'ORT1'
    
    if (name_ort == 'ORT1'):
        dir_data    = '/Users/throop/Data/ORT1/throop/backplaned/'
        reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'

    if (name_ort == 'ORT2'):
        dir_data    = '/Users/throop/Data/ORT2/throop/backplaned/'
        reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02']
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
    
    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
    
    hbt.figsize((12,12))
    
    # Load and stack the field images
    
    do_force = True   # If set, reload the stacks from individual frames, rather than restoring from a pickle file.
    
    stack_field = image_stack(os.path.join(dir_data, reqid_field),   do_force=do_force, do_save=do_force)

    stack_haz = {}
    
    for reqid_i in reqids_haz:
        stack_haz[reqid_i] = image_stack(os.path.join(dir_data, reqid_i), do_force=do_force, do_save=do_force)
            
    # Set the rough position of MU69
    
    et = stack_haz[reqids_haz[0]].t['et'][0] # Look up ET for first image in the Hazard stack
    (st, lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    vec = st[0:3]
    (_, ra, dec) = sp.recrad(vec)
    radec_mu69 = (ra, dec) # Keep in radians
    
    # Align the frames
    
    stack_field.align(method = 'wcs', center = (radec_mu69))
    for reqid_i in reqids_haz:
        stack_haz[reqid_i].align(method  = 'wcs', center = (radec_mu69))
    
    # Calc the padding required. This can only be done after the images are loaded and aligned.

    pad_field = stack_field.calc_padding()[0]
    pad_haz = []
    for reqid_i in reqids_haz:
        pad_haz.append(stack_haz[reqid_i].calc_padding()[0])
    pad_haz.append(pad_field)
        
    pad = max(pad_haz)
       
    # Flatten the stacks into single output images
    # If we get an error here, it is probably due to a too-small 'pad' value.
    
    img_field = stack_field.flatten(do_subpixel=False, method='median',zoom=zoom, padding=pad)
    
    img_haz = {}
    img_haz_diff = {}
    for reqid_i in reqids_haz:
        img_haz[reqid_i]  = stack_haz[reqid_i].flatten(do_subpixel=False,  method='median',zoom=zoom, padding=pad)
        img_haz_diff[reqid_i] = img_haz[reqid_i] - img_field
        
    # Plot the trimmed, flattened images. This is just for show. They are not scaled very well for ring search.
    
    plt.imshow(stretch(img_field))
    for reqid_i in reqids_haz:
        plt.imshow(stretch(hbt.trim_image(img_haz_diff[reqid_i])))
        plt.title(reqid_i)
        plt.show()
        
        # Save as FITS
    
        file_out = '/Users/throop/Desktop/test_zoom{}.fits'.format(zoom)
        hdu = fits.PrimaryHDU(stretch(diff_trim))
        hdu.writeto(file_out, overwrite=True)
        print(f'Wrote: {file_out}')
    
    # Calculate and display radial profiles
    
    hbt.figsize((10,8))
    hbt.set_fontsize(20)
    pos =  np.array(np.shape(img_field))/2  # MU69 will be at the center of this array
    for reqid_i in reqids_haz:
        (radius, profile_dn) = get_radial_profile_circular(img_haz[reqid_i] - img_field, pos=pos, width=2)
        radius_km = radius * stack_haz[reqid_i].pixscale_x_km
        plt.plot(radius_km, profile_dn, label=reqid_i)
    plt.xlim(0,50_000)
    plt.xlabel('Expanded Pixels')
    plt.ylabel('DN Median')
    plt.legend()
    plt.show()
    
    # Now, convert from DN to I/F
    
    for reqid_i in reqids_haz:
 
        (radius, profile_dn) = get_radial_profile_circular(img_haz[reqid_i] - img_field, pos=pos, width=2)
        radius_km = radius * stack_haz[reqid_i].pixscale_x_km
        
        RSOLAR_LORRI_1X1 = 221999.98  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
        RSOLAR_LORRI_4X4 = 3800640.0  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
        
        C = profile_dn # Get the DN values of the ring. Typical value is 1 DN.
        
        # Define the solar flux, from Hal's paper.
        
        FSOLAR_LORRI  = 176.	     	    # We want to be sure to use LORRI value, not MVIC value!
        F_solar = FSOLAR_LORRI # Flux from Hal's paper
        
        RSOLAR = RSOLAR_LORRI_4X4
        
        # Calculate the MU69-Sun distance, in AU (or look it up). 
        
        km2au = 1 / (u.au/u.km).to('1')
        
        et = stack_haz[reqid_i].t['et'][0]
        (st,lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
        r_nh_mu69 = sp.vnorm(st[0:3]) * km2au # NH distance, in AU
        
        (st,lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'Sun')
        r_sun_mu69 = sp.vnorm(st[0:3]) * km2au # NH distance, in AU
        
        pixscale_km =  (r_nh_mu69/km2au * (0.3*hbt.d2r / 256)) / zoom # km per pix (assuming 4x4)
        
        TEXP = stack_haz[reqid_i].t['exptime'][0]
        
        I = C / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. All v similar, except for spectrum assumed.
        
        # Apply Hal's conversion formula from p. 7, to compute I/F and print it.
        
        profile_IoF = math.pi * I * r_sun_mu69**2 / F_solar # Equation from Hal's paper
        plt.plot(radius_km, profile_IoF, label=reqid_i)
        plt.xlim(0,50_000)
    plt.xlabel('Dist [km], 4x4 zoomed??')
    plt.ylabel('DN Median')
    plt.ylim((-2e-7,5e-7))
    plt.legend()
    plt.title(name_ort)
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
    plt.show()

    # Now see if I can take all fields, zoom them appropriately, and sum them.
    # We want to scale them all to the scale of the final frame.
    
    pixscale_km_out = stack_haz[reqids_haz[-1]].pixscale_x_km
#    out = 0*img_haz['K1LR_HAZ02'].copy()
    size_out = np.shape(img_haz[reqids_haz[-1]])
    
    img_rescale_3d = np.zeros((len(reqids_haz),size_out[0],size_out[1]))
    img_haz_rescale = {}
    for i,reqid_i in enumerate(reqids_haz):
        
        zoomfac = stack_haz[reqid_i].pixscale_x_km / pixscale_km_out
        arr = scipy.ndimage.zoom(img_haz[reqid_i] - img_field, zoomfac)
        
        size_in = np.shape(arr)
        edge_left = int( (size_in[0]-size_out[0])/2)
        edge_top = int( (size_in[1]-size_out[1])/2)
        
        arr_out = arr[ edge_left : edge_left+size_out[0], 
                       edge_top  : edge_left+size_out[1] ] 
        img_haz_rescale[reqid_i] = arr_out
        img_rescale_3d[i,:,:] = arr_out
    
    img_rescale_mean     = np.mean(img_rescale_3d, axis=0)
    img_rescale_median = np.median(img_rescale_3d, axis=0)
    
    plt.imshow( stretch(img_rescale))
    plt.title(f'restretched and summed, {name_ort}')
    plt.show()
    
    file_out = f'/Users/throop/Desktop/img_rescale_median_{name_ort}.fits'
    hdu = fits.PrimaryHDU(img_rescale_median)
    hdu.writeto(file_out, overwrite=True)
    print(f'Wrote: {file_out}')
        
    print(f'zooming by {zoomfac}, to size {np.shape(img_haz_rescale[reqid_i])}')
    
    (radius_med, profile_dn_median) = get_radial_profile_circular(img_rescale_median, pos=pos, width=2)
    (radius_med, profile_dn_mean)     = get_radial_profile_circular(img_rescale_mean,    pos=pos, width=2)
    
    hbt.figsize((8,6))
    plt.plot(radius, profile_dn_median,label = 'All curves, Median')
    plt.plot(radius, profile_dn_mean, label = 'All curves, Mean')
    plt.xlim((0,80))
    plt.legend()
    plt.ylabel('DN')
    plt.title(name_ort)
    plt.xlabel('pixels of some sort')
    plt.show()
    
#    stack_field.set_indices(hbt.frange(0,10))

# Make a plot of the offset in pixels between frames

#    plt.plot(stack_field.t['shift_x_pix'], stack_field.t['shift_y_pix'])
#    plt.xlabel('RA pixels')
#    plt.ylabel('Dec pixels')
#    plt.show()

 # Load the first HAZ4 frame. Plot MU69's position on it.
 
    file = '/Users/throop/Dropbox/Data/ORT1/throop/backplaned/K1LR_HAZ04/lor_0406990332_0x633_pwcs_backplaned.fits'
    hdu = fits.open(file)
    w = WCS(file)
    plt.imshow(stretch(hdu[0].data))
     
    et = hdu[0].header['SPCSCET']
    (st, lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    vec = st[0:3]
    (_, ra, dec) = sp.recrad(vec) 
    (pos_pix_x, pos_pix_y) = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)
     
    plt.plot(pos_pix_x, pos_pix_y, marker = 'o', ms=5, color='red')
    plt.show()
     
    hdu.close()

