#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 11:43:11 2018

THIS IS MY MAIN CODE TO MAKE SUPERSTACKS, IMPLANT RINGS, TAKE RADIAL PROFILES, ETC,
FOR THE MU69 FLYBY.

IT IS CALLED 'ORT' BUT THAT IS JUST HISTORICAL.

@author: throop
"""

import glob
import math
import os.path
import os

import astropy
from   astropy.io import fits
import astropy.table
import matplotlib
import matplotlib.pyplot as plt # pyplot
from   matplotlib.figure import Figure
import numpy as np
import astropy.modeling
import spiceypy as sp
from   astropy import units as u           # Units library
import pickle # For load/save

import scipy
import copy
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel

# HBT imports

import hbt

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular
from   get_radial_profile_backplane import get_radial_profile_backplane
from   get_radial_profile_backplane import get_radial_profile_backplane_quadrant
from   plot_img_wcs import plot_img_wcs
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes
from   scipy.optimize import curve_fit
from   wcs_translate_pix import wcs_translate_pix, wcs_zoom

# Define a gaussian function, for the fit
    
def gaus(x,a,x0,sigma):
    """
    Just a simple Gaussian function, to be used for a fit function.
    """
    
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

#%%%
    
def nh_ort_make_superstack(stack, 
                           img, 
                           img_field,
                           exptime,
                           exptime_field,
                           name_stack_base, 
                           do_save_all=False, 
                           do_backplanes=True, 
                           dir='', 
                           str_name='', 
                           method = 'wcs',
                           do_center=True,
                           frame='2014_MU69_SUNFLOWER_ROT',
                           wcs = '',
                           format_out = 'small',  # 'small' or 'full',
                           zoom = 1,
                           do_fast_backplanes = False, **kwargs):
    
    """
    This makes a superstack image. The output image has resolution matching the lowest resolution 
    of all the input images (that is, the closest to MU69).
    
    This is a standalone function -- not part of any class.
    
    If requested, the superstacks *and* the rescaled individual stacks, are all written to disk.
    
    Parameters
    -----
    
    stack:
        A dictionary of `image_stack` objects.
        
    img:
        A dictonary of 2D images, created by flattening the stacks above. The dictionary 
        keys of `img` and `stack` must be the same.
        
    img_field:
        A dictionary of 2D images, for the field frames. These images are aligned exactly with the 
        'img' imags. 
        If img_field = None, then the fields are ignnored entirely.
    
    name_stack_base:
        String, which identifies an item in the `stack`.' 
        Which one of the input stacks do we use for the 'base' for the superstack? Typically this will be 
        the highest-resolution one. This stack is used to set the output resolution, and the ET (ie, center position).
    
    format_out:
        string to describe the format of the output. 
        - default 'small' = crop to the size of the final frame (ie, where all the data overlap). 
        - 'full' = use all the data, showing all frames, even where some are not overlapping.
        
    Optional parameters
    -----
    
    zoom:
        The zoom value to use when stacking.
        
    index_reference:
        Index to which of the final images to use as the output scale.
        (Alternatively, we might allow the dictionary to be passed)
        
    do_save_all:
        Boolean. If set, write all of the rescaled stacks, and superstacks, to disk.
        
    str_name:
        String. A prefix to use in writing the filenames for the stacks.
        
    do_center:
        Boolean. If set, code will do a centroiding on the central region of each image and center each 
        stack based on that. This is Q&D, but works if the brightest object is moving, and near the center,
        such as is the case for MU69 ORT.
    
    method:
        Method used to align each individual stack. Can be 'brightest' or 'wcs'.
        I used to use 'brightest' (where center was set based on seeing MU69), but now should use 'wcs'.
        
    do_backplanes:
        Boolean. If set, will return a set of backplanes for the superstacks.
        NB: THIS REQUIRES DO_SAVE_ALL = True.
    
    do_fast_backplanes:
        Boolean. If set, and generating backplanes, then will generate only an abbreviated set.
        [NOT CURRENTLY IMPLEMENTED.]
        
    wcs:
        A pre-computed WCS which will be put into the output image. Should be already zoomed, centered, etc. 
        as needed. This is desired to have, because the routine will output FITS files, which should have full WCS.
        
    frame: 
        A SPICE frame to use for creating the backplane.
        
    """
        
    keys = list(stack.keys())
    
    # Set up a bunch of dictionaries, one per plane
    
    pixscale_km = {}   # Pixel scale. Note that this is the 
    magfac      = {}   # Magnification factor
    shape_mag   = {}   # Shape when magnified
            
    # Search for the highest-res stack, and make it so that the superstack matches the resolution of that stack.
    # Regardless of size, we always take the highest resolution one.    

    for key in keys:
        pixscale_km[key] = stack[key].pixscale_x_km
    pixscale_km_out = min(pixscale_km.values())
    for key in keys:
        magfac[key]      = stack[key].pixscale_x_km / pixscale_km_out
        shape_mag[key]   = (np.ceil(np.array(np.shape(img[key])) * magfac[key])).astype(int)

    # Look up the name of the 'base' plane -- that is, the highest resolution (smallest pixscale)
    
    name_stack_base = keys[np.where(pixscale_km_out == np.array(list(pixscale_km.values())))[0][0]]
        
#    pixscale = stack[name_stack_base].pixscale_x_km  # This is the pixel scale of the final stack.
                                                                 # Zoom has not been applied yet.
    if format_out == 'small':

        # Take the shape of the highest res frame, as is
    
        size_out = np.shape(img[name_stack_base])
    
    else:  # Take the shape of the largest frame, when mag'd. This is more complicated to determine than it should be!
        
        size_out_x = np.max(np.array(list(shape_mag.values()))[:,0]).astype(int)
        size_out_y = np.max(np.array(list(shape_mag.values()))[:,1]).astype(int)
        
        size_out = (size_out_x, size_out_y)
        
    print(f'Final output size = {size_out}')
    
    # Create an array (not a list) for the output images to go into
    
    img_rescale_3d = np.zeros((len(keys),size_out[0],size_out[1]))
    img_rescale = {}
    img_zoom    = {}
    
    # Loop over every individual stack, and zoom them.
    # Assuming we use WCS centering, each one is centered.
    
    for i,key in enumerate(keys):
        
        img_denan = img[key].copy()
        img_denan[np.isnan(img_denan)] = 0.  # Convert any NAN to 0, since ndimage.zoom propagates NAN
        img_zoom[key] = scipy.ndimage.zoom(img_denan, magfac[key])
        
        arr = img_zoom[key]
        
        print(f'Zoomed image {key} by {magfac[key]}')
        
        size_in   = np.shape(arr)
        
        # Calculate the position in the output array, based on centering the image
        
        if (format_out == 'small'):   # If we crop this image to fit it exactly to the output frame
            
            edge_left = int( (size_in[0]-size_out[0])/2)
            edge_top  = int( (size_in[1]-size_out[1])/2)
            
            arr_out = arr[ edge_left : edge_left+size_out[0], 
                           edge_top  : edge_left+size_out[1] ]
            
            print(f'Zoomed image cropped to shape {np.shape(arr_out)}')

        if (format_out == 'full'):  # If we take this entire frame, pad it, and then put that into the output
            
            edge_left = int( (size_out[0]-size_in[0])/2)
            edge_top  = int( (size_out[1]-size_in[1])/2)
            
            arr_out = np.zeros(size_out)
            
            arr_out[ edge_left : edge_left+size_in[0], edge_top:edge_top+size_in[1] ] = arr
            
            print(f'Zoomed image padded to shape {np.shape(arr_out)}')
            
        # if (method == 'wcs'):
        #     # Align based on WCS (ie, alignment has already been done). 
        #     # This is a much better method. If we have done things right, then MU69 should
        #     # be already centered, so we just copy the image right in.
        #     # The only possible issue is that for MU69, I am currently assuming no target motion. Might have to 
        #     # adjust my stack-creation routines later if this is an issue. And, this would not work properly for 
        #     # (e.g.) HST images.
            
        #     shift_x = 0
        #     shift_y = 0
        #     arr_out = arr_out
            
        # # Do a one-off to rotation 365A / 365B into proper orientation. This will invalidate their WCS info, but oh well.
        
        # if ('2018365B' in reqid_i) or ('2018365B' in reqid_i):
        #     arr_out = np.rot90(arr_out, 1)
 
        plt.imshow(stretch(arr_out))
        plt.title(f'zoomed frame {key}, x{magfac[key]}')
        plt.show()
                    
        img_rescale[reqid_i]  = arr_out
        img_rescale_3d[i,:,:] = arr_out

        # Write the scaled stacks to disk, if requested
        
        if do_save_all:
            file_out = f'stack_{key}_{str_name}_rescaled_hbt.fits'
            path_out = os.path.join(dir_out, file_out)
            hdu = fits.PrimaryHDU(arr_out) 
            hdu.writeto(path_out, overwrite=True)
            print(f'Wrote: {path_out}')
    
    # Do the actual image flattening
        
    img_rescale_median = np.median(img_rescale_3d, axis=0)  
    img_rescale_mean   = np.mean(  img_rescale_3d, axis=0)
    
    plt.imshow(stretch(img_rescale_median))
    
    # plt.plot(img_rescale_3d[1,800,:])
    # plt.plot(img_rescale_3d[2,800,:])
    # plt.plot(img_rescale_3d[3,800,:])

    # plt.ylim(0,50)
    # plt.xlim(500,1000)
    
    plt.show()
      
    # Put the WCS info in a header, which we will write to FITS file
    
    header = wcs.to_header()
    
    # Also, put an ET into the FITS header. Get this from the stack.
    
    et = stack_haz[name_stack_base].t['et'][0]
    header['SPCSCET'] = f'{et}'
    header['SPCUTCID'] = sp.et2utc(et, 'C', 0)
    
    # Write the superstacks (with WCS) to disk, if requested
    
    if do_save_all:
        
        # Write full-size images
        
        file_out = f'superstack_n{len(keys)}_{str_name}_median_wcs_hbt.fits'  # Write the median stack file
        hdu = fits.PrimaryHDU(img_rescale_median, header=header)
        path_out = os.path.join(dir_out, file_out)
        hdu.writeto(path_out, overwrite=True)
        print(f'Wrote median stack: {path_out}')
      
        path_out_main = path_out  # Strangely, we don't need to .copy() a string. Ahh - that is just for np objects!
        
        file_out = f'superstack_n{len(keys)}_{str_name}_mean_wcs_hbt.fits'                 # Write the mean stack file
        hdu = fits.PrimaryHDU(img_rescale_mean, header=header)
        path_out = os.path.join(dir_out, file_out)
        hdu.writeto(path_out, overwrite=True)
        print(f'Wrote mean   stack: {path_out}')
        
        if do_backplanes:
            # Now that we have written a FITS file to disk, compute backplanes for it
            # ** for compute_backplanes to work, use_sky_plane = False must be set in that function.
            
            # frame = '2014_MU69_SUNFLOWER_ROT'
            
            print(f'Computing backplanes using frame={frame}: {path_out_main}')
            planes = compute_backplanes(path_out_main, 'MU69', frame, 'New Horizons', angle3=30*hbt.d2r)

    if not do_backplanes:
        planes = None           # Make sure to set planes to None if none requested
        
    return((img_rescale_mean, img_rescale_median, planes, pixscale_km_out))


### Make a ring annulus mask
    
def mask_annulus(arr, a, da, pixscale=1, tau=1, shape='tophat', width_gauss = None):
    
    """
    Create a fake ring with a specified radius, width, and optical depth.
    
    Parameters
    -----
    
    arr:
        Input array. Used just as a size reference. Must be square.
        
    a:
        Ring radius, in km or pixels.
        
    da:
        Ring width, in km or pixels.
        
    pixscale:
        Pixels per km. Or, 1 to indicate that all measurements are in pixels.
        
    tau:
        Optical depth.
    
    shape:
        Ring profile. 'tophat', 'triangle', etc.
        
    width_gauss:
        Width of the gaussian convolution kernel, in pixels.
        
    """
    
    arr_out = np.zeros(np.shape(arr))
    
    dist    = hbt.dist_center(hbt.sizex(arr)) * pixscale
    
    # if (shape == 'tophat'):
    
    mask = ((dist > (a - da/2)) & (dist < (a + da/2))).astype(int) * tau
    
    if width_gauss:
        kernel = astropy.convolution.Gaussian2DKernel(width_gauss / pixscale)
        mask_convolved = astropy.convolution.convolve(mask, kernel)
        
        mask = mask_convolved
        
    return mask
    
#%%%

# =============================================================================
# End of function definitions
# =============================================================================
    
# =============================================================================
# Below is the one-off code (not a function) that does all the image stacking for the MU69 ORT's.
#   - Load stacks
#   - Align them
#   - Make superstacks
#   - Make radial profile     
# =============================================================================

if (__name__ == '__main__'):

#%%%    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    cmap_superstack = 'plasma'
    cmap_stack      = 'Greys_r'
    
    plt.set_cmap(cmap_stack)


    ########## SET PARAMETERS HERE #################
    
    zoom = 4     # How much to magnify images by before shifting. 4 (ie, 1x1 expands to 4x4) is typical
                  # 1 is faster; 4 is slower but better.

    width = 1  # Bin width for radial profiles
    
    do_tunacan = False

    # DO_FORCE_STACKS: If set, reload the stacks from individual frames, rather than restoring from a pickle file.
    # NB: When adding a new set of OpNavs to an existing run, there sometimes is not enough padding
    # to allow for room for the latest OpNav to be added, if it has more jitter than previous.
    # So, typically do do_force=True when adding a new OpNav visit.
  
    # DO_FORCE_FLATTEN: Same thing
      
    do_force_stacks_haz   = True      
    do_force_flatten_haz  = True    

    do_force_stacks_field = True  # Keep as False. Except if using new parameters  [well, not sure. Maybe better True]
    do_force_flatten_field= True  # Keep as False
     
    ########## END PARAMETERS HERE #################
    
    name_ort      = 'MU69_Approach'  # This is the actual encounter -- not ORT!
    initials_user = 'HBT'
    dir_data      = '/Users/throop/Data'
    
    if (name_ort == 'MU69_Approach'):
        dir_images    = os.path.join(dir_data, name_ort, 'throop', 'backplaned')
        dir_out       = os.path.join(dir_data, name_ort, 'throop', 'stacks')

        reqids_haz = ['KELR_MU69_APDEEP_L4_2018365A',
                       'KELR_MU69_APDEEP_L4_2018365B',
                       'KELR_MU69_APDEEP_L4_2018365E',
                       'KELR_MU69_APDEEP_L4_2018365F',
                      ]
        
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
   
    # Check for a tunacan profile
    
    if do_tunacan:
        frame = '2014_MU69_TUNACAN_ROT'
        frametype = 'Tunacan'
    else:
        frame = '2014_MU69_SUNFLOWER_ROT'
        frametype = 'Sunflower'

    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
        
    # Make sure output dir exists
    
    os.makedirs(dir_out, exist_ok=True)
    
    # Load and stack the field images

    # We create a field stack for each reqid. The contents are the same, but the shift values will be different
    
    stack_haz    = {}
    stack_field  = {}
    img_haz      = {}
    img_field    = {}
    wcs_haz      = {}
    wcs_field    = {}
    et_haz       = {}
    img_haz_diff = {}
    radec_mu69_haz={}    

    pad_haz      = []  # Array, not a dictionary, since all we do here is take a max
    pad_field    = []

#%%%
    
# =============================================================================
# Load each stack from disk, from the raw FITS files
# =============================================================================

    is_stack_ready = {}
    
    for i,reqid_i in enumerate(reqids_haz):    # Use the FIELD reqid, but store under the HAZ reqid. 
                                               # This is so that we will have a version of each FIELD 
                                               # shifted appropriately for each HAZ REQID.

        # Generate a filename which we'll use for saving files from this stack.
        # This is not the path of the original FITS files
           
        # Load the hazard stack
        
        num_frames = len(glob.glob(os.path.join(dir_images, reqid_i, '*_backplaned.fits')))

        file_stack_base = os.path.join(dir_out, 
                                     f'stack_{reqid_i}_{name_ort}_n{num_frames}_z{zoom}')

        file_stack_pkl = file_stack_base + '.pkl'

        # Check if we have a flattened, zoomed stack already on disk.
        # If we do, then reload it. This lets us skip this loop (to load), *and* the next loop (to flatten).
        
        if os.path.isfile(file_stack_pkl):
            print(f'Restoring flattened stack: {file_stack_pkl}')

            lun = open(file_stack_pkl, 'rb')
            
            (stack_haz[reqid_i], 
             # stack_field[reqid_i], 
             img_haz[reqid_i],   
             # img_field[reqid_i],
             wcs_haz[reqid_i],   
             # wcs_field[reqid_i],
             et_haz[reqid_i],
             # img_haz_diff[reqid_i],
             radec_mu69_haz[reqid_i],
             )                                              = pickle.load(lun)
            
            lun.close()             
            
            is_stack_ready[reqid_i] = True               # Set a flag to indicate flattened, zoomed stack is ready.

        else: # Reload raw files from disk

            stack_haz[reqid_i] = image_stack(os.path.join(dir_images, reqid_i), do_force=do_force_stacks_haz, 
                                  do_save=False)
    
            # Load the corresponding field stack
            
            # if i == 0:
            #     stack_field[reqid_i] = image_stack(os.path.join(dir_images, reqid_field), 
            #                           do_force = do_force_stacks_field, 
            #                           do_save = False)
            # else:
            #     stack_field[reqid_i] = copy.deepcopy(stack_field[reqids_haz[0]])
            
            stack_field = None
            
            # Look up ET and MU69 position in each stack
            
            et_haz[reqid_i] = stack_haz[reqid_i].t['et'][0] # Look up ET for each stack Hazard stack (0th frame in each)
            (st, lt) = sp.spkezr('MU69', et_haz[reqid_i], 'J2000', 'LT', 'New Horizons')
            (_, ra, dec) = sp.recrad(st[0:3])
            radec_mu69_haz[reqid_i] = (ra, dec) # Keep in radians
        
        # Align the field frames and hazard frames to MU69 position
        
            stack_haz[reqid_i].align(  method  = 'wcs', center = (radec_mu69_haz[reqid_i]))  # align all images in stack
            # stack_field[reqid_i].align(method  = 'wcs', center = (radec_mu69_haz[reqid_i]))
        
        # Calc the padding required. This can only be done after the images are loaded and aligned.
        
            pad_haz.append(  stack_haz  [reqid_i].calc_padding()[0])
            # pad_field.append(stack_field[reqid_i].calc_padding()[0])
            
            is_stack_ready[reqid_i] = False

#%%%
            
# Done with loading stacks. Now calculate the padding, or take a large fixed values
# To save memory we can calculate an optimum padding. But it is more flexible if we just take a large one.
            
    do_calc_padding = False

    if do_calc_padding:        
        pad = np.amax([pad_haz, pad_field])
    else:
        pad = 100 # For a mosaic, what is the appropriate padding value to use? Try 100.

#%%%
        
# =============================================================================
#     Loop over the reqids, and zoom and flatten each stack.
# =============================================================================

    # NB: If we get an error here, it is probably due to a too-small 'pad' value, often caused by adding new
    # files and not increasing pad() to accomodate them all.
    
# =============================================================================
#     Flatten the stacks and zoom them
# =============================================================================

    do_plot_individual_stacks = True
    
    for reqid_i in reqids_haz:
        
        # First check if we already loaded the stacks from above. If so, then the flag will 
        # be set, and we're good to go. The zooming is the time-sink here, so no need to repeat it.
        
        if is_stack_ready[reqid_i]:
            print(f'Stack {reqid_i} ready and zoomed')
            
        else:

            num_planes = stack_haz[reqid_i].num_planes
            
            exptime_haz_i   = np.median(stack_haz[reqid_i].t['exptime'])
            
            # exptime_field_i = np.median(stack_field[reqid_i].t['exptime'])
            
            # Flatten the hazard frames
            
            (img_haz[reqid_i], wcs_haz[reqid_i])  =\
                  stack_haz[reqid_i].flatten(do_subpixel=False,  method='median0',zoom=zoom, padding=pad, 
                           do_force=do_force_flatten_haz, do_save=False)
            
            # Flatten the field frames
            
            # (img_field[reqid_i], wcs_field[reqid_i])  =\
            #        stack_field[reqid_i].flatten(do_subpixel=False,  method='median0',zoom=zoom, padding=pad, 
            #                 do_force=do_force_flatten_field, do_save=False)
            
            # Zero out the field. We do not want to use it here, period.
            
            # img_field[reqid_i] *= 0
            
            # Do the image subtraction
            
            # img_haz_diff[reqid_i] = img_haz[reqid_i] - img_field[reqid_i]
                                  
            # Save stack as FITS
            # Put three planes in this: Haz, Field, Diff
            
            file_out_fits = os.path.join(dir_out, 
                                         f'stack_{reqid_i}_{name_ort}_n{num_planes}_z{zoom}.fits')
            
            hdu1 = fits.PrimaryHDU(img_haz[reqid_i])
            # hdu2 = fits.ImageHDU(img_field[reqid_i])
            # hdu3 = fits.ImageHDU(img_haz_diff[reqid_i])
            
            hdulist = fits.HDUList([hdu1])
            hdulist.writeto(file_out_fits, overwrite=True)

            # hdu.writeto(file_out_fits, overwrite=True)
            
            print(f'Wrote single-plane FITS: {file_out_fits}')
            
            # Plot an extracted image, and save as PNG
            
            file_out_png = file_out_fits.replace('fits', 'png')
            width_extract = 800
            
            hbt.figsize((15,15))
            plt.subplot(1,2,1)   # Plot image
            plot_img_wcs(img_haz[reqid_i], wcs_haz[reqid_i], title = reqid_i, 
                 name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[reqid_i], 
                 width=width_extract,
                 do_show=False,
                 cmap = 'plasma')

            # plt.subplot(1,2,2)  # Plot image - field
            # plot_img_wcs(img_haz_diff[reqid_i], wcs_haz[reqid_i], title = reqid_i + ' diff', 
            #      name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[reqid_i], 
            #      width=width_extract,
            #      do_show=False,
            #      cmap = 'plasma')
            
            plt.savefig(file_out_png)
            plt.show()
            print(f'Wrote: {file_out_png}')
            
            # Save stack as pickle. This is a large pickle file, that includes
            # the hazard stack *and its accompanying field stack*. They are in matched pairs 
            # with same shift values so it really makes sense to save them together.
            # Also save the WCS, and the ET, and the radec, all in one big pickle.
            
            file_stack_base = os.path.join(dir_out, 
                                 f'stack_{reqid_i}_{name_ort}_n{num_planes}_z{zoom}')
            
            file_stack_pkl = file_stack_base + '.pkl'
            lun = open(file_stack_pkl, 'wb')
            pickle.dump((stack_haz[reqid_i], 
                         # stack_field[reqid_i], 
                         img_haz[reqid_i],   
                         # img_field[reqid_i],
                         wcs_haz[reqid_i],   
                         # wcs_field[reqid_i],
                         et_haz[reqid_i],
                         # img_haz_diff[reqid_i],
                         radec_mu69_haz[reqid_i]), lun)
    
            pickle.dump(stack_haz[reqid_i], lun)
    
            lun.close()
            print("Wrote: " + file_stack_pkl)    
            
            # Finally, set a flat indicating that the stack is ready
            
            is_stack_ready[reqid_i] = True

# =============================================================================
# # Make a useful string for titles
# =============================================================================

    num_frames       = 0
    num_visits_opnav = 0
    num_visits_haz   = 0
    for reqid_i in reqids_haz:
        num_frames += stack_haz[reqid_i].num_planes
        if 'HAZ' in reqid_i.upper():
            num_visits_haz +=1
        if 'OPNAV' in reqid_i.upper():
            num_visits_opnav +=1
            
    if (num_visits_opnav + num_visits_haz) == 1:  # If only one reqid, say 'Stack of 10, 2018315'
        str_stack = 'Stack'
    else:
        str_stack = 'Superstack'
        # str_reqid = reqids_haz[0].split('_')[-1]
        # str_stack = f'Stack of {stack_haz[reqids_haz[0]].size[0]}'
    # else:    
    
    if len(reqids_haz) == 1:    
        str_reqid = reqids_haz[0].split('_')[-1]
    else:
        str_reqid = reqids_haz[0].split('_')[-1] + ' .. ' + reqids_haz[-1].split('_')[-1]
        
    str_reqid = str_reqid.replace('2018', '')  # Shorten REQID from '2018345' → '345'
    
    # str_stack = 'Superstack'
    
    if num_visits_haz >= 1:
        str_stack += f', {num_visits_haz} HAZ'
    if num_visits_opnav >= 1:
        str_stack += f', {num_visits_opnav} OPNAV'
        
    str_stack += f', {num_frames} frames'
            
#%%%
# =============================================================================
#     Make a 'superstack' -- ie stack of stacks, put on the same spatial scale and stacked
# =============================================================================
    
    keys = list(stack_haz.keys()) 
    
    pixscale_km = {}
    for key in keys:
        pixscale_km[key] = stack_haz[key].pixscale_x_km / zoom
        
    # Search for the highest-res stack, and make it so that the superstack matches the resolution of that stack.
    
    pixscale_km_out = min(pixscale_km.values())
    
    # Look up what the index is for this image
    
    name_stack_base = keys[np.where(pixscale_km_out == np.array(list(pixscale_km.values())))[0][0]]
    
    str_name = f'{name_ort}_z{zoom}'

    # Grab the WCS from the highest-res frame, and apply that to this frame.
    
    wcs_superstack = wcs_haz[name_stack_base]
    
    # Get the exposure times
    
    exptime_haz   = np.median(stack_haz[name_stack_base].t['exptime'])
    # exptime_field = np.median(stack_field[name_stack_base].t['exptime'])

    # Make the superstack, and get the backplanes.
    
    print('Making superstack...')
    
    do_backplanes = False  # IF this is false, then 'None' is returned for backplanes
            
    (img_superstack_mean, img_superstack_median, backplanes, pixscale_km_superstack) = \
                                                                          nh_ort_make_superstack(stack_haz, 
                                                                          img_haz, None, 
                                                                          exptime_haz, None,
                                                                          name_stack_base, 
                                                                          do_save_all=True, dir=dir_out,
                                                                          str_name=str_name, do_center=True,
                                                                          do_backplanes=do_backplanes,
                                                                          frame=frame,
                                                                          format_out = 'full',
                                                                          wcs = wcs_haz[name_stack_base],
                                                                          )
    
    
    # plot_img_wcs(img_superstack_median, wcs_haz[name_stack_base], width=1800)

    # Set the pixel scale of the superstack properly. Strangely, stacks themself list the original
    # resolution, not the mag'd resolution, so nh_ort_make_superstack has no way of knowing.
    
    pixscale_km_superstack /= zoom
    
    # Center the superstack on MU69
    # NB: We need to call this routine 2x: once with big box, and once with narrower, to converge in on it.

    (img_superstack_mean, shift)   = hbt.center_array_by_mass(img_superstack_mean,   boxwidth=300*zoom)
    (img_superstack_mean, shift)   = hbt.center_array_by_mass(img_superstack_mean,   boxwidth=100*zoom)
    (img_superstack_mean, shift)   = hbt.center_array_by_mass(img_superstack_mean,   boxwidth= 20*zoom)
    print(f'Superstack Shift mean = {shift}')
    
    (img_superstack_median, shift) = hbt.center_array_by_mass(img_superstack_median, boxwidth=300*zoom)
    (img_superstack_median, shift) = hbt.center_array_by_mass(img_superstack_median, boxwidth=100*zoom)
    (img_superstack_median, shift) = hbt.center_array_by_mass(img_superstack_median, boxwidth= 20*zoom)
    print(f'Superstack Shift median = {shift}')

#%%%    
    # Center each of the stacks on MU69, using simple centroid

    keys = img_haz.keys()
    
    for key in keys:
        for boxwidth_i in [20]:
            (img_shift,shift) = hbt.center_array_by_mass(img_haz[key], boxwidth=boxwidth_i*zoom)
            img_haz[key] = img_shift
            plt.imshow(stretch(crop(img_shift, (400,400))))
            plt.title(f'{key}, just shifted {shift}, boxwidth {boxwidth_i}')
            plt.show()      
    
    # Test centering algorithms
    
    
    # Check the centering, using WCS
    
    for key in keys:
            plot_img_wcs(img_haz[key], wcs_haz[key], cmap=cmap_superstack,
                  width=200,
                  scale_km = stack_haz[key].pixscale_x_km / zoom,
                  label_axis='km',
                  do_show=True, title=key)
    
    for key in keys:
            plot_img_wcs(img_haz[key], wcs_haz[key], cmap=cmap_superstack,
                  width=200,
                  scale_km = stack_haz[key].pixscale_x_km / zoom,
                  label_axis='km',
                  do_show=True, title=key)
            
    # Make a plot of all the stacks
    
    keys = list(keys)
    
    for key in keys:
        plot_img_wcs(img_haz[key], wcs_haz[key], cmap=cmap_superstack,
                  name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[key],
                  width=200,
                  scale_km = stack_haz[key].pixscale_x_km / zoom,
                  label_axis='km',
                  do_show=False, title=key)
        plt.show()

    # Construct a radius backplane. Make this from scratch, *not* from SPICE.
    
    plane_radius_superstack = hbt.dist_center(np.shape(img_superstack_median)[0]) * pixscale_km_superstack
         
#%%%

    # plane_radius_superstack = backplanes[0]['Radius_eq']
    # plane_longitude_superstack = backplanes[0]['Longitude_eq']

# =============================================================================
# Convert the stacks and superstacks from DN to I/F
# =============================================================================
    
    # Apply Hal's conversion formula from p. 7, to compute I/F and print it.

    RSOLAR_LORRI_1X1 = 221999.98  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
    RSOLAR_LORRI_4X4 = 3800640.0  # Diffuse sensitivity, LORRI 1X1. Units are (DN/s/pixel)/(erg/cm^2/s/A/sr)
    
    # Define the solar flux, from Hal's paper.
    FSOLAR_LORRI  = 176.	     	    # We want to be sure to use LORRI value, not MVIC value!
    F_solar       = FSOLAR_LORRI # Flux from Hal's paper
    RSOLAR        = RSOLAR_LORRI_4X4

    km2au = 1 / (u.au/u.km).to('1')
        
    # Calculate the MU69-Sun distance, in AU (or look it up).         
    
    et        = stack_haz[reqids_haz[-1]].t['et'][0] # ET for the final image stack
    (st,lt)   = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    r_nh_mu69 = sp.vnorm(st[0:3]) * km2au # NH distance, in AU
    (st,lt)   = sp.spkezr('MU69', et, 'J2000', 'LT', 'Sun')
    r_sun_mu69= sp.vnorm(st[0:3]) * km2au # NH distance, in AU
    pixscale_km =  (r_nh_mu69/km2au * (0.3*hbt.d2r / 256)) / zoom # km per pix (assuming LORRI 4x4)
    TEXP        = stack_haz[reqid_i].t['exptime'][0]  # Exposure time of the field frames. All are 29.967 sec.

    # Convert stacks
    
    img_haz_iof = {}
    for key in keys:
        img_haz_iof[key] = math.pi * img_haz[key]  * r_sun_mu69**2 / F_solar / stack_haz[key].t['exptime'][0] / RSOLAR

    # Convert superstacks

    I_median = img_superstack_median / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    I_mean   = img_superstack_mean   / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    
    img_superstack_median_iof = math.pi * I_median * r_sun_mu69**2 / F_solar # Equation from Hal's paper
    img_superstack_mean_iof   = math.pi * I_mean   * r_sun_mu69**2 / F_solar # Equation from Hal's paper


### Take radial profiles of each of the stacks
### Take radial profiles of the superstack
    
#%%%
# =============================================================================
#     Make a pair of final plots of the superstack, with and without aimpoints + sunflower rings
# =============================================================================
    
    width_pix_plot_arr = [1000, 400]
    
    # Define the ring sizes.
    # For the SUNFLOWER ring, we just take a radial profile outward, and plot it.
    # For the TUNACAN   ring, we need to assume a vertical thickness. I make two rings, assuming two thicknesses.
    
    a_ring_km = [250, 500]      
    da_ring_km = 10   # Width of the lines in the ring to plot, in km    
        
    hbt.figsize((20,8))
    hbt.fontsize(14)

    plt.subplot(1,2,1)
    
    # Set a min and max range for the plot.
    # We might change this, but for comparing images, it is good to keep this fixed.
    # NB: Simon's plots go -1e-6 .. 3e-6
    
    vmin = -1e-6
    vmax =  3e-5
        
    plot_img_wcs(img_superstack_median_iof, wcs_superstack, cmap=cmap_superstack, 
#                 title = f'{str_stack}, {str_reqid}, ring at {a_ring_km} km', 
                 title = f'{str_stack}, {str_reqid}', 
                 # width=width_pix_plot,
                 do_show=False,
                 scale_km=pixscale_km_superstack,
                 label_axis = 'KM',
                 name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                 do_colorbar=True, do_stretch=True, width=6000)
        
    # Construct the image of ring masks. This can be either from superstack, or constructed ad hoc.

    radius_km = plane_radius_superstack

    # Make the mask
    
    mask_ring0 = (radius_km > (a_ring_km[0] - da_ring_km/2)) & (radius_km < ((a_ring_km[0] + da_ring_km/2) ))
    mask_ring1 = (radius_km > (a_ring_km[1] - da_ring_km/2)) & (radius_km < ((a_ring_km[1] + da_ring_km/2) ))

    mask_ring = mask_ring0 | mask_ring1

    img_superstack_median_iof_masked = img_superstack_median_iof * (mask_ring == False).astype(int)

    for w_i in width_pix_plot_arr:     
        plt.subplot(1,2,1)
    
        vmin = 0e-5
        vmax = 2e-5
        plot_img_wcs(img_superstack_median_iof_masked, wcs_superstack, cmap=cmap_superstack, 
                     title = f'{str_stack}, {str_reqid}', 
                     width=w_i,
                     do_show=False,
                     name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                     do_colorbar=True, vmin=vmin, vmax=vmax, do_stretch=False,
                     scale_km=pixscale_km_superstack,
                     label_axis = '[km]')
    
        plt.subplot(1,2,2)
    
        plot_img_wcs(img_superstack_median_iof, wcs_superstack, cmap=cmap_superstack, 
        # plot_img_wcs(diff, wcs_superstack, cmap=cmap_superstack, 
                     title = f'{str_stack}, {str_reqid}', 
                     width=w_i,
                     do_show=False,
                     name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                     do_colorbar=True, vmin=vmin, vmax=vmax, do_stretch=False,
                     scale_km=pixscale_km_superstack,
                     label_axis = '[km]')
        
        plt.show()
    
#%%%    
# =============================================================================
# Make a set of plots showing the image stack and the radial profile, for each stack.
# =============================================================================
    
    hbt.figsize((12,6))
    hbt.fontsize(12)
    
    plot_img_wcs((img_superstack_median_iof), wcs_superstack, cmap=cmap_superstack, title='Superstack', width=200)
            
    for key in keys:
        keystr = key.split('_')[-1]
  
        pixscale_i = stack_haz[key].pixscale_x_km / zoom
        
        plt.subplot(1,2,1)
        plot_img_wcs(stretch(img_haz[key]), wcs_haz[key], do_stretch=True, cmap=cmap_superstack,
                     scale_km = pixscale_i, label_axis = '[km]', title = key, do_show=False)
        
        plt.subplot(1,2,2)
        width_profile_pix = 1
        (radius,profile) = get_radial_profile_circular(img_haz_iof[key], width=width_profile_pix)
        plt.xlabel('km')
        plt.ylabel('I/F')
        plt.ylim((0, 2e-5))
        bin_radius_max = np.where([np.isnan(profile)][0])[0][0]
        radius_max = radius[bin_radius_max] * pixscale_i
        iof_min = np.amin(profile[0:bin_radius_max])

        plt.plot(radius * pixscale_i, profile)
        plt.xlim((0, radius_max))
        plt.ylim((iof_min, 2e-5))
        plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
        plt.title(f'Profile, {keystr}, dr={width_profile_pix * pixscale_i:.1f} km')
        
        plt.tight_layout()
        plt.show()
        
    hbt.figsize() 
    hbt.fontsize()

#%%%
# =============================================================================
# Make plots showing the image stack and the radial profile, for superstack
# =============================================================================
    
    hbt.figsize((12,6))
    hbt.fontsize(12)
    
    keystr = 'superstack'
  
    pixscale_i = pixscale_km_superstack
    
    plt.subplot(1,2,1)
    plot_img_wcs((img_superstack_median_iof), wcs_superstack, cmap=cmap_superstack, title='Superstack', 
                 do_show=False)
    
    # plot_img_wcs(stretch(img_haz[key]), wcs_haz[key], do_stretch=True, cmap=cmap_superstack,
    #              scale_km = pixscale_i, label_axis = '[km]', title = key, do_show=False)
    
    plt.subplot(1,2,2)
    width_profile_pix = zoom
    
    # Get the profile. This is really slow for zoom=4! And, ideally we would do it for width_profile_pix = 1.
    
    (radius,profile) = get_radial_profile_circular(img_superstack_median_iof, width=width_profile_pix)

    # Make a plot
    
    plt.xlabel('km')
    plt.ylabel('I/F')
    plt.ylim((0, 2e-5))
    # radius_max = radius[bin_radius_max] * pixscale_i
    radius_max = 1400
    bin_radius_max = np.digitize(radius_max, radius)
    
    iof_min = np.median(profile[0:bin_radius_max-2]) # Doesn't work well.

    plt.plot(radius * pixscale_i, profile - iof_min)
    plt.xlim((0, radius_max))
    plt.ylim((0, 5e-6))
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
    plt.title(f'Profile, {keystr}, dr={width_profile_pix * pixscale_i:.1f} km')
    plt.xlabel('km')
    plt.ylabel('I/F')
    
    
    plt.tight_layout()
    plt.show()
        
    hbt.figsize() 
    hbt.fontsize()

# =============================================================================
# Now look up the pole position, for plot labels
# =============================================================================

    # Look up the pole position of MU69, based on the loaded frame kernel
    # This is just for labeling plots -- nothing else.
       
    name_target = 'MU69'
    
    if do_tunacan:
        frame = '2014_MU69_TUNACAN_ROT'
    else:
        frame = '2014_MU69_SUNFLOWER_ROT'
        
    name_observer = 'New Horizons'
    
    vec = [0, -1, 0]                                # -Y vector is the rotational pole
    mx_frame_j2k =  sp.pxform(frame, 'J2000', et_haz[name_stack_base])
    vec_j2k = sp.mxv(mx_frame_j2k, vec)
    (_, ra_pole, dec_pole) = sp.recrad(vec_j2k)

#%%%

# =============================================================================
# Now make some implants, and add these to the dictionary for analysis
# =============================================================================
    
# Create some implant rings
    
    tau_implant = 1e-6
    a_ring      = 300
    da_ring     = 5

    tau_implant = 2e-6
    a_ring      = 1000
    da_ring     = 10

    width_gauss_implant = zoom*2
    
    keys = img_haz_iof.keys()

    # First get rid of the old implants, if they are there
    # Iterate over all the keys in the dictionary. But, we can't change a dictionary in place while iterating,
    # so we make a *list* of the keys, and then iterate over that.
    
    for key in list(keys):
        if '_implant' in key:
            del img_haz_iof[key]
            del stack_haz[key]
            del wcs_haz[key]
            print(f'Deleted {key} from dictionary')

    # Now create and add the new implants. Just add each one to the dictionary.
            
    for key in list(keys):
        key_implant = key + '_implant'
        
        if key_implant not in keys and 'implant' not in key:
            pixscale_km = stack_haz[key].pixscale_x_km / zoom
            arr = img_haz_iof[key]

            arr_implant = (arr + 
               mask_annulus(arr, a_ring, da_ring, pixscale=pixscale_km, width_gauss=width_gauss_implant, tau=tau_implant) )
            
            img_haz[key_implant]     = arr_implant
            img_haz_iof[key_implant] = arr_implant
            stack_haz[key_implant]   = stack_haz[key]
            wcs_haz[key_implant]     = wcs_haz[key]
            
            print(f'Implanted {key} → {key_implant}')
        else:
            print(f'Already exists: {key_implant}')

    keys = img_haz_iof.keys()
    
#%%%    
# =============================================================================
# Now do some special processing, like subtracting Field A - Field, B, etc.
# Differential.    
# =============================================================================
    
    do_implants = True
    
    key1 = ['KELR_MU69_APDEEP_L4_2018365A', 'KELR_MU69_APDEEP_L4_2018365E']
    key2 = ['KELR_MU69_APDEEP_L4_2018365B', 'KELR_MU69_APDEEP_L4_2018365F']

    if do_implants:
        key1 = ['KELR_MU69_APDEEP_L4_2018365A_implant', 'KELR_MU69_APDEEP_L4_2018365E_implant']
        key2 = ['KELR_MU69_APDEEP_L4_2018365B_implant', 'KELR_MU69_APDEEP_L4_2018365F_implant']
    
    pairs = ['AB', 'EF'] 
    
    dxy   = {}
    profile_diff = {}
    radius_km_diff = {}
    
    dxy['1'] = [[9,0],   [50,   0]]  # zoom level 1. These shift values are btwn A and B, etc. Set these by hand.
    dxy['2'] = [[18, 0], [99,   0]]  # zoom level 2
    dxy['4'] = [[36, 0], [200,  1]]  # zoom level 4
    
    hbt.fontsize(15)
    width_profile_pix = 1
    
    for i in range(len(key1)):

        pair = pairs[i]
        
        pixscale_km = stack_haz[key1[i]].pixscale_x_km / zoom

        if 'implant' in key1[i]:
            keystr = f'Stack difference, {key1[i][-12:]} - {key2[i][-12:]}'
        else:    
            keystr = f'Stack difference, {key1[i][-4:]} - {key2[i][-4:]}'
 
        hbt.figsize((15,15)) 

        diff = (img_haz_iof[key1[i]] - np.roll(img_haz_iof[key2[i]],dxy[f'{zoom}'][i], axis=[0,1]))
        plt.subplot(1,2,1)
        
        plot_img_wcs(stretch(diff), wcs_haz[key1[i]], scale_km=pixscale_km, label_axis='km',title=keystr, 
                     do_show=False, cmap='Greys_r', do_colorbar=False)
        
        # Use this line below for testing shifts manually. remove the stretch() initially.

        # plt.imshow(stretch( (img_haz_iof[key1[i]] - np.roll(img_haz_iof[key2[i]],[9,0], axis=[0,1])) ))
        
        # End commented line!
        
        plt.subplot(1,2,2)
        
        plot_img_wcs(diff, wcs_haz[key1[i]], scale_km=pixscale_km, label_axis='km', width=200*zoom, 
                     do_stretch=True, title=keystr, do_show=False)
 
        # plot_img_wcs(img_haz_iof[key1[i]], wcs_haz[key1[i]], scale_km=pixscale_km, label_axis='km', width=200*zoom, 
        #              do_stretch=True, title=keystr, do_show=False)
        
        # plt.imshow(stretch(diff))
        # plt.title(keystr)
        plt.show()

        (radius,profile) = get_radial_profile_circular(diff, width=width_profile_pix)
        
        profile_diff[pair] = profile
        radius_km_diff[pair]  = radius * pixscale_km
        
        # plt.xlabel('km')
        # plt.ylabel('I/F')
        # plt.ylim((0, 2e-5))
        # radius_max = radius[bin_radius_max] * pixscale_i
        
        bin_radius_max = np.where(np.isnan(profile))[0][0]
        
        iof_min = np.nanmin(profile[0:bin_radius_max-2]) # Doesn't work well.
    
        hbt.figsize((15,5)) 
        plt.plot(radius*pixscale_km, profile - iof_min)
        plt.xlim((0, radius[bin_radius_max]*pixscale_km))
        plt.xlabel('km')
        plt.ylim((0, 5e-6))
        plt.ylabel('I/F')
        plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
        plt.title(f'Circular profile, {keystr}, dr={width_profile_pix * pixscale_km:.1f} km, zoom={zoom}')
        plt.show()

# Plot all the radial profiles together on one plot

    radius_km_max = 2000        
    for pair in pairs:
        
        bin_radius_max = int(np.digitize(1000, radius_km_diff[pair]))
        bias = np.amin(profile_diff[pair][0:bin_radius_max])
        
        plt.plot(radius_km_diff[pair], profile_diff[pair]-bias, label=pair)
        
    plt.legend()
    plt.title(f'Circular profile, {keystr}, dr={width_profile_pix * pixscale_km:.1f} km, zoom={zoom}')
    plt.ylim((-3e-7,5e-6))
    plt.xlabel('km')
    plt.ylabel('I/F')
    plt.xlim(0,radius_km_max)
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
    if do_implants:
        plt.axvline(a_ring, alpha=0.2, color='red', linestyle='--')
    plt.show()

#%%%
    
# Make an implated ring at a specified location, in the differential image
    
    diff_masked = diff * (1 - mask_annulus(diff, 250, 20, pixscale=pixscale_km, width_gauss=2))

    plot_img_wcs(diff_masked, wcs_superstack, cmap=cmap_superstack, 
                     title = f'{str_stack}, {str_reqid}', 
                     width=w_i,
                     do_show=False,
                     name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                     do_colorbar=True, vmin=-vmax, vmax=vmax, do_stretch=False,
                     scale_km=pixscale_km_superstack,
                     label_axis = '[km]')



    # if (shape == 'triangle'):
    #     mask_triangle = ()
    
# =============================================================================
#     Calculate the radial profile, in I/F units 
# =============================================================================
    
    # Take the radial profile of the superstack
    
    num_pts = 1000  # Number of radial points to use.
    
    profile_iof = {}  # Create a dictionary to save the radial profiles in
    
    # Take a sunflower radial profile  
                                             
    # (radius,  profile_iof_median)   = get_radial_profile_backplane(img_superstack_median_iof,
    #                                      plane_radius_superstack, method = 'median', num_pts = num_pts)

    # (radius_quadrant,  profile_iof_quadrant)   = get_radial_profile_backplane_quadrant(img_superstack_mean_iof,
    #                                      plane_radius_superstack, plane_azimuth_superstack, method = 'median', 
    #                                      num_pts = num_pts/2)
        
    # profile_iof[0] = profile_iof_median  # Just pick one -- they both are similar

    # # Calculate the bias level, crudely

    # radius_max_km = 1200
    # bin_radial_end = np.digitize(radius_max_km*0.9, radius)
    # bias       = np.nanmin(profile_iof[0][10:bin_radial_end])
    # bias = 8e-6
    
# =============================================================================
# Fit a gaussian to the radial profile, if requested
# =============================================================================

    # do_fit_profile = False
    
    # if do_fit_profile:
    #     r_0         = 4000      # Inner radius to ignore in gauss fit, in km
    #     radius_ring = 9000      # Starting point for gaussfit for ring position, in km
    #     hw_ring     = 1000      # Starting point for ring halfwidth, in km
        
        
    #     bin_0 = hbt.x2bin(r_0, radius)
    #     x = radius[bin_0:]
    #     y = profile_iof[bin_0:]            
        
    #     popt,pcov = curve_fit(gaus,x,y,p0=[radius_ring,0,hw_ring])
    
# =============================================================================
# Plot the radial profile
# =============================================================================
    
    # hbt.figsize((16,6))
    # hbt.fontsize(12)
    # ylim_profile = (-2e-6, 2e-5)
    # xlim_profile = (0, radius_max_km)

    # plt.subplot(1,2,1)
    # plt.plot(radius, profile_iof[0] - bias,
    #      label = f'Median, pole = ({ra_pole*hbt.r2d:.0f}°, {dec_pole*hbt.r2d:.0f}°), ' + 
    #              f'dr={round(np.diff(radius)[0]):.0f} km')

    # do_plot_errorbars = False
    
    # if do_plot_errorbars:
    #     plt.errorbar(radius, profile_iof - bias, yerr = profile_iof_std,
    #          label = f'Median, pole = ({ra_pole*hbt.r2d:.0f}, {dec_pole*hbt.r2d:.0f}) deg')
        
    
    # do_plot_mean=False
    # if do_plot_mean:            # Can plot the mean. But the median is much more useful
    #     plt.plot(radius, profile_iof_mean - bias, label = 'mean')
    
    # if do_fit_profile:
    #     plt.plot(x,gaus(x,*popt),'ro:', marker = None, ls = '--', lw=1.5, 
    #              label = f'Fit, radius={popt[1]:.0f} km, FWHM={2.35 * popt[2]:.0f} km')
    
    # FWHM = 2.35 * sigma: https://ned.ipac.caltech.edu/level5/Leo/Stats2_3.html
        
#     plt.xlim(xlim_profile)
#     plt.ylim(ylim_profile)
#     plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
#     plt.xlabel('Radius [km]')
#     plt.ylabel('I/F')
#     plt.legend(loc = 'upper right')
#     plt.title(f'Radial profile, {frametype}, {str_stack} {str_reqid}')

# # Plot the radial profile, for quadrants

#     plt.subplot(1,2,2)
#     plt.plot(radius_quadrant, profile_iof_quadrant[0] - bias,
#          label = f'Median, pole = ({ra_pole*hbt.r2d:.0f}°, {dec_pole*hbt.r2d:.0f}°), ' + 
#                  f'dr={round(np.diff(radius_quadrant)[0]):.0f} km', color='white')
    
#     for i in range(4):
#         plt.plot(radius_quadrant, profile_iof_quadrant[i] - bias, label=f'Quadrant {i}', alpha=0.7)

#     plt.xlim((0, radius_max_km))
#     plt.ylim(ylim_profile)
#     plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
#     plt.xlabel('Radius [km]')
#     plt.ylabel('I/F')
#     plt.legend(loc = 'upper right')
#     plt.title(f'Radial profile, {frametype}, {str_stack} {str_reqid}')
#     plt.show()

#%%%
                         
# =============================================================================
#      Write the ring profile to a text file as per wiki
#      https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Detection-to-Ring  
# =============================================================================
    
    version = 1
    
    file_out = os.path.join(dir_out, f'{initials_user.lower()}_{name_ort.lower()}_v{version}.ring')
    lun = open(file_out,"w")
    
    # This needs to be generalized to print as many profiles as we have. 
    # I don't know how to write a print format line with a variable number of fields!
    
    if do_tunacan:
        for i in range(len(profile_iof[i])):
                lun.write('{:10.3f} {:11.3e} {:11.3e}\n'.format(radius[i], profile_iof[0][i], profile_iof[1][i]))
    else:
        for i in range(len(profile_iof)):
            lun.write('{:10.3f} {:11.3e}\n'.format(radius[i], profile_iof[0][i]-bias))
    lun.close()
    
    do_ring_out = True

    if do_ring_out:
        print(f'Wrote: {file_out}')
        print(f' scp {file_out} ixion:\~/MU69_Approach/astrometry' )  # We don't copy it, but put up string for user.

        
######
        

# 