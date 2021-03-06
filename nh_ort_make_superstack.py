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
                           do_fast_backplanes = False, **kwargs):
    
    """
    This makes a superstack image. The output image has resolution matching the lowest resolution 
    of all the input images (that is, the closest to MU69).
    
    This is a standalone function -- not part of any class.
    al
    If requested, the superstacks *and* the rescaled individual stacks, are all written to disk.
    
    Parameters
    -----
    
    stack:
        A dictionary of `image_stack` objects.
        
    img:
        A dictonary of 2D images, created by flattening the stacks above. The dictionary 
        keys of `img` and `stack` must be the same.
    
    name_stack_base:
        String, which identifies an item in the `stack`.' 
        Which one of the input stacks do we use for the 'base' for the superstack? Typically this will be 
        the highest-resolution one. This stack is used to set the output resolution, and the ET (ie, center position).
        
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
    
    pixscale_km = {}
    for key in keys:
        pixscale_km[key] = stack[key].pixscale_x_km
        
    # Search for the highest-res stack, and make it so that the superstack matches the resolution of that stack.
    
    pixscale_km_out = min(pixscale_km.values())
    
    # Look up what the index is for this image
    
    name_stack_base = keys[np.where(pixscale_km_out == np.array(list(pixscale_km.values())))[0][0]]
    
#    pixscale = stack[name_stack_base].pixscale_x_km  # This is the pixel scale of the final stack.
                                                                 # Zoom has not been applied yet.
    size_out = np.shape(img[name_stack_base])
    
    img_rescale_3d = np.zeros((len(keys),size_out[0],size_out[1]))
    img_rescale = {}
    
    # Loop over every individual stack
    
    for i,key in enumerate(keys):
        
        magfac    = stack[key].pixscale_x_km / pixscale_km_out
        arr       = scipy.ndimage.zoom(img[key] - img_field[key], magfac)
        
        size_in   = np.shape(arr)
        edge_left = int( (size_in[0]-size_out[0])/2)
        edge_top  = int( (size_in[1]-size_out[1])/2)
        
        arr_out = arr[ edge_left : edge_left+size_out[0], 
                       edge_top  : edge_left+size_out[1] ]
        
        # If requested, do a centroid on the central region of the image. 
        # Then shift the entire image based on this, to put brightest object at center of frame.
        # This is the old style of aligning. It does *not* work for actual MU69 OpNav approach images, because
        # these are taken far enough out that MU69 is definitely not the brightest thing in the images.
    
        if (method == 'brightest'):
            extract = arr_out[int(size_out[0]/2)-50:int(size_out[0]/2+50),    # 100x100 pixel region. Hardcoded, ugh.
                              int(size_out[0]/2)-50:int(size_out[0]/2+50)]
            shift_x = 50 - hbt.wheremax(np.sum(extract, axis=0))
            shift_y = 50 - hbt.wheremax(np.sum(extract, axis=1))  # vertical axis on screen. 
                                                                  #  Axis=0. - means value is too large. 
            arr_out = np.roll(np.roll(arr_out, shift_x, axis=1), shift_y, axis=0)
            
            # print(f'Rolling by {shift_x}, {shift_y} to align {key} into superstack')
        
        if (method == 'wcs'):
            # Align based on WCS (ie, alignment has already been done). 
            # This is a much better method. If we have done things right, then MU69 should
            # be already centered, so we just copy the image right in.
            # The only possible issue is that for MU69, I am currently assuming no target motion. Might have to 
            # adjust my stack-creation routines later if this is an issue. And, this would not work properly for 
            # (e.g.) HST images.
            
            shift_x = 0
            shift_y = 0
            arr_out = arr_out
            
            # print(f'Aligned by WCS. Rolling by {shift_x}, {shift_y} to align {key} into superstack')
                       
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

    if do_backplanes:        
        return((img_rescale_mean, img_rescale_median, planes))
    else:
        return((img_rescale_mean, img_rescale_median))

# =============================================================================
# Implant a ring into an image
# =============================================================================

def ring_implant(img, file_ring_implant, resolution, iof_max_ring, backplanes, width_plot_pix=450, 
                 do_plot=True, **kwargs):
    """
    Implant a ring into an image.

    Parameters
    -----
     
    img: 
        Image array to implant in. Assumed units of I/F
     
    file_ring_implant:
        Filename of a 'grid file' defining the ring. This file is output by NH_ORT_TRACK4_FLYBY.PY

    resolution:
        Resolution of the grid image, km/pix. This is specified ultimately by the header.txt file.
         
    iof:
        The I/F of the peak of the smeared image
         
    do_plot:
        Boolean. Make plots of the implanted ring?
    
    Optional parameters
    -----    
    
    width_plot_pix:
        Size of the plot, in zoomed LORRI pixels.
        
    """
         
    lun = open(file_ring_implant, 'rb')
    print(f'Reading ring implant pickle file {file_ring_implant}')
    (density) = pickle.load(lun)
    lun.close()

    if 'sunflower' in file_ring_implant:
        axis_sum = 1

    if 'tunacan' in file_ring_implant:
        axis_sum = 1
            
    # Flatten, rotate, and scale the image appropriately, from 3D → 2D

    img_ring     = np.rot90(np.sum(density, axis_sum),1)

    img_ring_iof = iof_max_ring * img_ring / np.amax(img_ring) # Scale the I/F. We adjust this later, after convolve.
    
    # Get the km per pix of the superstack, and scale ring to match it
    
    scale_superstack = np.abs( (backplanes[0]['dRA_km'][0,1]) - (backplanes[0]['dRA_km'][0,0]) )

    scale_ring = resolution  # resolution of grid image, km/pix
    # if 'sunflower10k' in file_ring_pickle:
        
    magfac = scale_ring / scale_superstack
    
    img_ring_iof_zoom = scipy.ndimage.zoom(img_ring_iof, magfac)
    
    # Convolve a PSF with the ring
    
    psf = img_superstack_median_iof.copy()
    wheremax = np.where(psf == np.amax(psf))  # Returns 2D position
    dxy_psf = 12  # Approx fullwidth of MU69, in pixels, to extract the PSF. Even number.
    psf_extract = psf[int(wheremax[0]-dxy_psf/2):int(wheremax[0]+dxy_psf/2), 
                      int(wheremax[1]-dxy_psf/2):int(wheremax[1]+dxy_psf/2)] 
    psf_extract = psf_extract / np.sum(psf_extract)
    
    shape_ss = np.shape(img_superstack_median_iof)
    shape_ring = np.shape(img_ring_iof_zoom)

    # Scale the ring s.t. I/F is set to the proper value. We do this after convolve, rather than before, 
    # since that is what we did in NH_ORT_TRACK4_CALIBRATE.PY as well
    
    img_ring_iof_zoom_c = scipy.signal.convolve2d(img_ring_iof_zoom, psf_extract)    
    
    img_ring_iof_zoom_c *= (iof_max_ring / np.amax(img_ring_iof_zoom_c))
    
    # Tweak the scaled image so it's even
    
    if shape_ring[0] %2:
        img_ring_iof_zoom = img_ring_iof_zoom[:,0:-1]

    if shape_ring[1] %2:
        img_ring_iof_zoom = img_ring_iof_zoom[0:-1,:]
        
    shape_ring = np.shape(img_ring_iof_zoom_c)

    img_merged_iof = img_superstack_median_iof.copy()
    
    img_merged_iof[ int(shape_ss[0]/2-shape_ring[0]/2) : int(shape_ss[0]/2+shape_ring[0]/2),
                int(shape_ss[1]/2-shape_ring[1]/2) : int(shape_ss[1]/2+shape_ring[1]/2) ] += img_ring_iof_zoom_c
    
    if do_plot:               
        hbt.figsize((15,15))
        plt.subplot(2,2,1)
                
        plot_img_wcs(img_superstack_median_iof, wcs_superstack, width=width_plot_pix, do_show=False, title='Stack',
                     name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                     cmap='plasma', **kwargs)
        
        plt.subplot(2,2,2)
        plot_img_wcs(img_ring_iof_zoom, wcs_superstack, width=width_plot_pix, do_show=False, 
                     title=f'{file_ring_short}, I/F={iof_max_ring:.0e}', do_stretch=False, do_inhibit_axes=True, 
                     name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                     cmap='plasma')

        plt.subplot(2,2,3)
        plot_img_wcs(img_ring_iof_zoom_c, wcs_superstack, width=width_plot_pix, do_show=False, 
                     title=f'Convolved', do_stretch=False, do_inhibit_axes=True, 
                     name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                     cmap='plasma')

        plt.subplot(2,2,4)        
        plot_img_wcs(img_merged_iof, wcs_superstack, width=width_plot_pix, do_show=False, 
                     name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                     title=f'Stack + Convolved',
                     cmap='plasma', **kwargs)
        
        plt.show()
        hbt.figsize()
        
    return img_merged_iof
    
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
    
    zoom = 1     # How much to magnify images by before shifting. 4 (ie, 1x1 expands to 4x4) is typical
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
        reqids_haz  = [
                        'KALR_MU69_OpNav_L4_2018228', 
                        'KALR_MU69_OpNav_L4_2018258', 'KALR_MU69_OpNav_L4_2018264',
                        'KALR_MU69_OpNav_L4_2018267', 
                        'KALR_MU69_OpNav_L4_2018284', 'KALR_MU69_OpNav_L4_2018287', 'KALR_MU69_OpNav_L4_2018298',
                        'KALR_MU69_OpNav_L4_2018301', 'KALR_MU69_OpNav_L4_2018304',
                        'KALR_MU69_OpNav_L4_2018306', 'KALR_MU69_OpNav_L4_2018311',
                        'KALR_MU69_OpNav_L4_2018314', 
                        'KALR_MU69_OpNav_L4_2018315',
                        'KALR_MU69_OpNav_L4_2018316',
                        'KALR_MU69_OpNav_L4_2018317',
                        'KALR_MU69_OpNav_L4_2018319',
                        'KALR_MU69_OpNav_L4_2018325',
                        'KALR_MU69_Hazard_L4_2018325',  # 110 frames
                        'KALR_MU69_OpNav_L4_2018326',

                 # A pretty good stack: 330 .. 341
                 
                        'KALR_MU69_OpNav_L4_2018330',  # 10 frames
                        'KALR_MU69_OpNav_L4_2018331',  # 10 frames
                        'KALR_MU69_OpNav_L4_2018332',  # 10 frames
                        
                        'KALR_MU69_OpNav_L4_2018334',  # 10 frames
                        'KALR_MU69_OpNav_L4_2018335',  # 10 frames
                        'KALR_MU69_OpNav_L4_2018337',  # 10 frames
                        'KALR_MU69_Hazard_L4_2018334',  # 96 frames

                        'KALR_MU69_OpNav_L4_2018338',
                        'KALR_MU69_OpNav_L4_2018339',
                        'KALR_MU69_Hazard_L4_2018340',  # 61 frames
                        'KALR_MU69_OpNav_L4_2018340',
                        'KALR_MU69_OpNav_L4_2018341',
                        
                        'KALR_MU69_OpNav_L4_2018342', # 6 frames
                        'KALR_MU69_OpNav_L4_2018343',  # 6 frames 
                        'KALR_MU69_OpNav_L4_2018344', # 6 frames
                        'KALR_MU69_OpNav_L4_2018345', # 6 frames
                        'KALR_MU69_Hazard_L4_2018340',  # 61 frames
                        'KALR_MU69_Hazard_L4_2018344',  # 52 frames
                        'KALR_MU69_Hazard_L4_2018347',
                       ]

        # reqids_haz  = [          # For final MU69 Approach image
        #                             'KALR_MU69_Hazard_L4_2018334',
        #                             'KALR_MU69_Hazard_L4_2018340',
        #                             'KALR_MU69_Hazard_L4_2018344',
        #                             'KALR_MU69_Hazard_L4_2018347',
        # ]

   
        # reqids_haz  = [        

        #                 'KALR_MU69_OpNav_L4_2018342', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018343',  # 6 frames 
        #                 'KALR_MU69_OpNav_L4_2018344',  # 6 frames 
        #                 'KALR_MU69_OpNav_L4_2018345', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018347', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018348', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018349', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018350', # 6 frames

        # ]
        
        # reqids_haz += [        

        #                 'KALR_MU69_OpNav_L4_2018353', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018354',  # 6 frames 
        #                 'KALR_MU69_OpNav_L4_2018355', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018357', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018358', # 6 frames
        #                 'KALR_MU69_OpNav_L4_2018359', # 6 frames

        # ]
        
        reqids_haz = ['KELR_MU69_APDEEP_L4_2018365A',
                      # 'KELR_MU69_APDEEP_L4_2018365B',
                      # 'KELR_MU69_APDEEP_L4_2018365E',
                      # 'KELR_MU69_APDEEP_L4_2018365F',
                      ]

  
#####
                        
        # reqids_haz  = ['KALR_MU69_OpNav_L4_2018298','KALR_MU69_OpNav_L4_2018301']
        # reqids_haz  = ['KALR_MU69_OpNav_L4_2018301']
        # reqids_haz  = ['KALR_MU69_OpNav_L4_2018287']
        # reqids_haz  = ['KALR_MU69_OpNav_L4_2018284']
        
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
        
    if (name_ort == 'ORT1'):
        dir_images    = os.path.join(dir_data, name_ort, 'backplaned')
        dir_out       = os.path.join(dir_data, name_ort)
        reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'

    if (name_ort == 'ORT2'):
        dir_images    = os.path.join(dir_data, name_ort, 'throop', 'backplaned')
        dir_out       = os.path.join(dir_data, name_ort)
        reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
        reqids_haz  = ['K1LR_HAZ02', 'K1LR_HAZ03']
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
        
    if (name_ort == 'ORT2_OPNAV'):
        dir_images    = '/Users/throop/Data/ORT2/throop/backplaned/'
        dir_out       = os.path.join(dir_data, name_ort.split('_')[-1])
        dirs = glob.glob(dir_data + '/*LR_OPNAV*')         # Manually construct a list of all the OPNAV dirs
        reqids_haz = []
        for dir_i in dirs:
            reqids_haz.append(os.path.basename(dir_i))
        reqids_haz = sorted(reqids_haz)    
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
        
    if (name_ort == 'ORT3'):
        dir_images    = os.path.join(dir_data, name_ort, 'throop', 'backplaned')
        dir_out       = os.path.join(dir_data, name_ort)
        reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
#        reqids_haz  = ['K1LR_HAZ02', 'K1LR_HAZ03']
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
        
    if (name_ort == 'ORT4'):
        dir_images    = os.path.join(dir_data, name_ort, 'throop', 'backplaned')
        dir_out       = os.path.join(dir_data, name_ort, 'throop', 'stacks')
        reqids_haz  = ['K1LR_HAZ00', 'K1LR_HAZ01', 'K1LR_HAZ02', 'K1LR_HAZ03', 'K1LR_HAZ04']
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
        a_xy = (1, math.cos(hbt.d2r * 30))

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
            
            (stack_haz[reqid_i], stack_field[reqid_i], 
             img_haz[reqid_i],   img_field[reqid_i],
             wcs_haz[reqid_i],   wcs_field[reqid_i],
             et_haz[reqid_i],
             img_haz_diff[reqid_i],
             radec_mu69_haz[reqid_i],
             )                                              = pickle.load(lun)
            
            lun.close()             
            
            is_stack_ready[reqid_i] = True               # Set a flag to indicate flattened, zoomed stack is ready.

        else: # Reload raw files from disk

            stack_haz[reqid_i] = image_stack(os.path.join(dir_images, reqid_i), do_force=do_force_stacks_haz, 
                                  do_save=False)
    
            # Load the corresponding field stack
            
            if i == 0:
                stack_field[reqid_i] = image_stack(os.path.join(dir_images, reqid_field), 
                                      do_force = do_force_stacks_field, 
                                      do_save = False)
            else:
                stack_field[reqid_i] = copy.deepcopy(stack_field[reqids_haz[0]])
                    
            # Look up ET and MU69 position in each stack
            
            et_haz[reqid_i] = stack_haz[reqid_i].t['et'][0] # Look up ET for each stack Hazard stack (0th frame in each)
            (st, lt) = sp.spkezr('MU69', et_haz[reqid_i], 'J2000', 'LT', 'New Horizons')
            (_, ra, dec) = sp.recrad(st[0:3])
            radec_mu69_haz[reqid_i] = (ra, dec) # Keep in radians
        
        # Align the field frames and hazard frames to MU69 position
        
            stack_haz[reqid_i].align(  method  = 'wcs', center = (radec_mu69_haz[reqid_i]))  # align all images in stack
            stack_field[reqid_i].align(method  = 'wcs', center = (radec_mu69_haz[reqid_i]))
        
        # Calc the padding required. This can only be done after the images are loaded and aligned.
        
            pad_haz.append(  stack_haz  [reqid_i].calc_padding()[0])
            pad_field.append(stack_field[reqid_i].calc_padding()[0])
            
            is_stack_ready[reqid_i] = False

#%%%
            
# Done with loading stacks. Now calculate the padding, or take a large fixed values
# To save memory we can calculate an optimum padding. But it is more flexible if we just take a large one.
            
    do_calc_padding = False

    if do_calc_padding:        
        pad = np.amax([pad_haz, pad_field])
    else:
        pad = 80

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
            exptime_field_i = np.median(stack_field[reqid_i].t['exptime'])
            
            # Flatten the hazard frames
            
            (img_haz[reqid_i], wcs_haz[reqid_i])  =\
                  stack_haz[reqid_i].flatten(do_subpixel=False,  method='median0',zoom=zoom, padding=pad, 
                           do_force=do_force_flatten_haz, do_save=False)
            
            # Flatten the field frames
            
            (img_field[reqid_i], wcs_field[reqid_i])  =\
                  stack_field[reqid_i].flatten(do_subpixel=False,  method='median0',zoom=zoom, padding=pad, 
                           do_force=do_force_flatten_field, do_save=False)
            
            # Scale the field frames s.t. exptimes match the main exposure
            # We do not adjust the EXPTIME in the field header, but the DN values will be scaled.
            
            img_field[reqid_i] *= (exptime_haz_i / exptime_field_i)
            
            # Do the image subtraction
            
            img_haz_diff[reqid_i] = img_haz[reqid_i] - img_field[reqid_i]
                                  
            # Save stack as FITS
            # Put three planes in this: Haz, Field, Diff
            
            file_out_fits = os.path.join(dir_out, 
                                         f'stack_{reqid_i}_{name_ort}_n{num_planes}_z{zoom}.fits')
            
            hdu1 = fits.PrimaryHDU(img_haz[reqid_i])
            hdu2 = fits.ImageHDU(img_field[reqid_i])
            hdu3 = fits.ImageHDU(img_haz_diff[reqid_i])
            
            hdulist = fits.HDUList([hdu1, hdu2, hdu3])
            hdulist.writeto(file_out_fits, overwrite=True)

            # hdu.writeto(file_out_fits, overwrite=True)
            
            print(f'Wrote multi-plane FITS: {file_out_fits}')
            
            # Plot an extracted image, and save as PNG
            
            file_out_png = file_out_fits.replace('fits', 'png')
            width_extract = 100
            
            hbt.figsize((15,15))
            plt.subplot(1,2,1)   # Plot image
            plot_img_wcs(img_haz[reqid_i], wcs_haz[reqid_i], title = reqid_i, 
                 name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[reqid_i], 
                 width=width_extract,
                 do_show=False,
                 cmap = 'plasma')

            plt.subplot(1,2,2)  # Plot image - field
            plot_img_wcs(img_haz_diff[reqid_i], wcs_haz[reqid_i], title = reqid_i + ' diff', 
                 name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[reqid_i], 
                 width=width_extract,
                 do_show=False,
                 cmap = 'plasma')
            
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
            pickle.dump((stack_haz[reqid_i], stack_field[reqid_i], 
                         img_haz[reqid_i],   img_field[reqid_i],
                         wcs_haz[reqid_i],   wcs_field[reqid_i],
                         et_haz[reqid_i],
                         img_haz_diff[reqid_i],
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
        pixscale_km[key] = stack_haz[key].pixscale_x_km
        
    # Search for the highest-res stack, and make it so that the superstack matches the resolution of that stack.
    
    pixscale_km_out = min(pixscale_km.values())
    
    # Look up what the index is for this image
    
    name_stack_base = keys[np.where(pixscale_km_out == np.array(list(pixscale_km.values())))[0][0]]
    
    str_name = f'{name_ort}_z{zoom}'

    # Grab the WCS from the highest-res frame, and apply that to this frame.
    
    wcs_superstack = wcs_haz[name_stack_base]
    
    # Get the exposure times
    
    exptime_haz   = np.median(stack_haz[name_stack_base].t['exptime'])
    exptime_field = np.median(stack_field[name_stack_base].t['exptime'])

    # Make the superstack, and get the backplanes.
    
    print('Making superstack...')
    
    (img_superstack_mean, img_superstack_median, backplanes) = nh_ort_make_superstack(stack_haz, 
                                                                          img_haz, img_field, 
                                                                          exptime_haz,
                                                                          exptime_field,
                                                                          name_stack_base, 
                                                                          do_save_all=True, dir=dir_out,
                                                                          str_name=str_name, do_center=True,
                                                                          do_backplanes=True,
                                                                          frame=frame,
                                                                          wcs = wcs_haz[name_stack_base],
                                                                          )
    
    for key in keys:
        plot_img_wcs(img_haz[key], wcs_haz[key], cmap=cmap_superstack,
                  name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[key], width=100,
                  do_show=False, title=key)
        plt.show()
    
#%%%
# =============================================================================
#     Adjust the WCS and backplanes, if needed
# =============================================================================    
    
    # OK, now do a q&d adjust. The radius backplane is off by a few pixels, and I don't know why. 
    # Just go ahead and shift it and force it to be centered (that is, at center of array, which is where
    # MU69 usually should be). 
    # This is not the best solution, but it'll work for now.
    
    # Get the centroid of the radius 

    plane_radius_superstack = backplanes[0]['Radius_eq']
    r = plane_radius_superstack.copy()

    plane_longitude_superstack = backplanes[0]['Longitude_eq']

    centroid = scipy.ndimage.measurements.center_of_mass(r < 10000)
#    print(f'Centroid is {centroid}')

    dxy = np.array(np.shape(r))/2 - np.array(scipy.ndimage.measurements.center_of_mass(r < 5000))  
    r_roll = np.roll(r, np.round(dxy.astype(int)), axis=(0, 1))

    centroid = scipy.ndimage.measurements.center_of_mass(r_roll < 10000)
#    print(f'Centroid is {centroid}')
    
    # Save the shifted 'radius' backplane back to the array
    
    backplanes[0]['Radius_eq'] = r_roll

    plane_radius_superstack = backplanes[0]['Radius_eq']
    
    do_wcs_tweaks = False
        
    if do_wcs_tweaks:
        
        dx_wcs_tweak = 0.0   # Positive values shift to the right
        dy_wcs_tweak = 0.0  #   was 1.5, 1. Now 3.5, 6. Now (4, 3). Now (4,4). Now (0,0)
        
        wcs_superstack_tweak = wcs_superstack.deepcopy()
        wcs_translate_pix(wcs_superstack_tweak, dx_wcs_tweak, dy_wcs_tweak)
        
        # Tweak the backplane slightly as well. Not sure quite why this is necessary, since backplane should 
        # be exactly aligned w/ WCS.
        
        dx_backplane_tweak = 0
        dy_backplane_tweak = 0
        plane_radius_superstack_tweak = plane_radius_superstack.copy()
        plane_radius_superstack_tweak = np.roll(np.roll(plane_radius_superstack_tweak, dx_backplane_tweak, axis=1), 
                                          dy_backplane_tweak, axis=0)
        
        # Copy the tweaked versions over
        
        wcs_superstack          = wcs_superstack_tweak
        plane_radius_superstack = plane_radius_superstack_tweak # Do not copy this back to the backplane itself! Breaks.
    
# =============================================================================
# Convert the stacks from DN to I/F
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

    I_median = img_superstack_median / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    I_mean   = img_superstack_mean   / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    
    img_superstack_median_iof = math.pi * I_median * r_sun_mu69**2 / F_solar # Equation from Hal's paper
    img_superstack_mean_iof   = math.pi * I_mean   * r_sun_mu69**2 / F_solar # Equation from Hal's paper

#%%%    
# =============================================================================
# Implants: read a model ring image and place it over the data image, if requested
# =============================================================================
        
    do_implant = False
    
    if do_implant:
        dir_ring_img =    '/Users/throop/data/ORT5/throop/deliveries/'
    
        # Set the ring to read. Note that when changing this, need to also change the resolution 
        # in the function itself. Either 250 km/pix or 500 km/pix. That is the resolution used for the 
        # simulations by DPH/DK. If ring position looks off by 2x, then change this value in the function.
        
        # (file_ring_implant, resolution) = (dir_ring_img + 
        #      'dph-sunflower3.5k/ort5_None_y3.0_q2.5_pv0.05_rho0.46.dust_img.pkl', 250)
            
          
        # (file_ring_implant, resolution) = (dir_ring_img + \
        #     'dph-tunacan3.5kinc55/ort5_None_y2.2_q2.5_pv0.05_rho0.46.dust_img.pkl', 250)
            
        (file_ring_implant, resolution) = (dir_ring_img + \
            'dph-tunacan10kinc70/None_None_y3.0_q2.5_pv0.05_rho1.00.dust_img.pkl', 500)

        (file_ring_implant, resolution) = (dir_ring_img + \
            'dph-sunflower10k/ort5_None_y2.2_q2.5_pv0.05_rho0.46.dust_img.pkl', 250)
        
        file_ring_short = file_ring_implant.split('/')[7]

        # Set the I/F of the ring
        
        iof_max_ring = 5e-7

        # And implant it!
        
        img_merged_iof  = ring_implant(img_superstack_median_iof, file_ring_implant, resolution, iof_max_ring,
                                            backplanes,
                                            vmin=-1e-6, vmax=2e-6,
                                            width_plot_pix=700,
                                            do_plot=True)
            
#%%%
# =============================================================================
#     Make a pair of final plots of the superstack, with and without aimpoints + sunflower rings
# =============================================================================
    
    width_pix_plot = 400*zoom
    
    # Define the ring sizes.
    # For the SUNFLOWER ring, we just take a radial profile outward, and plot it.
    # For the TUNACAN   ring, we need to assume a vertical thickness. I make two rings, assuming two thicknesses.
    
    a_ring_km = [3500, 10000]         # For the sunflower *and* the tunacan, we use these as the central ring radii.
                                      # These values are used to *plot* both TUNA and SUNFLOWER. But they are not
                                      # used in computation of either one.
                                  
    thick_tuna_km = [1000, 2000]  # Midplane distance to limit to, in km. For TUNACAN only. This is used in 
                                     # computing the ring profile. 
    
    plt.show()  # for some reason the figsize() isn't working, so maybe flush things out??
    
    hbt.figsize((20,8))
    hbt.fontsize(14)

    plt.subplot(1,2,1)
    
    # Set a min and max range for the plot.
    # We might change this, but for comparing images, it is good to keep this fixed.
    # NB: Simon's plots go -1e-6 .. 3e-6
    
    vmin = -1e-6
    vmax =  1e-5
    
    plot_img_wcs(img_superstack_median_iof, wcs_superstack, cmap=cmap_superstack, 
#                 title = f'{str_stack}, {str_reqid}, ring at {a_ring_km} km', 
                 title = f'{str_stack}, {str_reqid}', 
                 width=width_pix_plot,do_show=False,
                 name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                 do_colorbar=True, do_stretch=False, vmin=vmin, vmax=vmax)
        
    # Construct the image of ring masks
    # We plot these the same in TUNACAN or SUNFLOWER -- they are just marked on directly
    
    mask_tuna = {}
    
    if do_tunacan:

        da_ring_km = 300    # Define the edge sizes of the tunacan ring, for plotting only 

        rad  = plane_radius_superstack
        vert = backplanes[0]['Altitude_eq']  # Extract the vertical height backplane
        
        for i in range(len(thick_tuna_km)):
        
            mask_tuna_topbot = ( ((np.abs(vert) > thick_tuna_km[i]) & (np.abs(vert) < (thick_tuna_km[i] + da_ring_km)))
                                  & (np.abs(rad) < a_ring_km[i]) )
            
            mask_tuna_rl     = ( ((np.abs(rad)  > a_ring_km[i]) & (np.abs(rad) < (a_ring_km[i] + da_ring_km)))
                                  & (np.abs(vert) < thick_tuna_km[i]) )
    
            mask_tuna[i] = mask_tuna_topbot | mask_tuna_rl
       
        mask_ring = mask_tuna[0] | mask_tuna[1]  # Combine inner and outer tunacans, to make one pretty plot with both.

    else: # Make the sunflower ring masks
        
        da_ring_km = 300   # Width of the lines in the ring to plot, in km    
        rad = plane_radius_superstack    
        mask_ring0 = (rad > a_ring_km[0]) & (rad < (a_ring_km[0] + da_ring_km) )
        mask_ring1 = (rad > a_ring_km[1]) & (rad < (a_ring_km[1] + da_ring_km) )
    
        mask_ring = mask_ring0 | mask_ring1
    
    plt.imshow(np.logical_not(mask_ring), alpha=0.12, origin='lower')
    
    # Mark the aimpoints on the plot
    
    ut_ca = '2019 1 Jan 05:33'
    et_ca = sp.utc2et(ut_ca)
    frame = 'J2000'
    abcorr = 'None' # Tried LT and NONE

    trajectories = ['alternate', 'prime']  # Make 'prime' second, so I don't screw things up for someone else.
    
    for i,trajectory in enumerate(trajectories):
        
        # Dump all kernels and load the proper kernel for this trajectory
        
        hbt.unload_kernels_all()
        file_tm = f'kernels_kem_{trajectory}.tm'
        sp.furnsh(file_tm)
        
        # Vector from MU69 → NH at C/A
        (st, lt) = sp.spkezr('New Horizons', et_ca, frame, abcorr, 'MU69')
        vec_mu69_nh_ca = st[0:3]
        print(f'Trajectory = {trajectory}, MU69-NH dist = {sp.vnorm(vec_mu69_nh_ca):7.6} km')
        
        # Vector from NH to MU69, at time of image
        (st, lt) = sp.spkezr('MU69', et_haz[name_stack_base], frame, abcorr, 'New Horizons')
        vec_nh_mu69_haz = st[0:3]
        
        # Add these two vectors, to get the vec from NH to aimpoint, at time of image
        vec_nh_mu69_haz_aimpoint = vec_nh_mu69_haz + vec_mu69_nh_ca
    
        # And convert into an xy position
        (_, ra_aimpoint, dec_aimpoint) = sp.recrad(vec_nh_mu69_haz_aimpoint)
        (x, y) = wcs_superstack.wcs_world2pix(ra_aimpoint*hbt.r2d, dec_aimpoint*hbt.r2d, 0)
        print(f'Trajectory = {trajectory}, CA RA/Dec = {ra_aimpoint*hbt.r2d, (dec_aimpoint*hbt.r2d)}')
        
        plt.plot(x, y, marker = 'X', markersize=10, color='lightgreen')
     
    plt.subplot(1,2,2)

    plot_img_wcs(img_superstack_median_iof, wcs_superstack, cmap=cmap_superstack, 
                 title = f'{str_stack}, {str_reqid}', 
                 width=width_pix_plot,do_show=False,
                 name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[name_stack_base],
                 do_colorbar=True, vmin=vmin, vmax=vmax, do_stretch=False)

    plt.show()
    
#%%%    
# =============================================================================
# Make a set of plots showing the field, the superstack, and all the individual stacks
#   This is a grid 3 x N.    
# =============================================================================
    
    hbt.figsize((8,8))
    hbt.fontsize(12)
    
    # plot_img_wcs(stretch(img_field), wcs_superstack, cmap=cmap_superstack, width=150, title='Field')
    
    stretchz=astropy.visualization.ManualInterval(vmin=2, vmax=40)
    
    width_postage = 120
    
    for key in keys:
        keystr = key.split('_')[-1]
        plt.subplot(1,3,1)

        plot_img_wcs(stretch(img_field[key]), wcs_field[key], cmap=cmap_superstack, width=width_postage, 
                     title = 'Field',
                     do_show=False, name_observer = 'New Horizons', name_target = 'Sun', et = et_haz[key],
                     )

        plt.subplot(1,3,2)
        plot_img_wcs(stretch(img_haz[key]), wcs_field[key], cmap=cmap_superstack, width=width_postage, 
                     title = f'Haz, {keystr}',
                     do_show=False, name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[key],
                     do_inhibit_axes=True,
                     )
        
        plt.subplot(1,3,3)
        plot_img_wcs(stretch(img_haz[key] - img_field[key]), 
                     wcs_field[key], cmap=cmap_superstack, width=width_postage, title = f'Haz-field, {keystr}',
                     do_show=False, name_observer = 'New Horizons', name_target = 'MU69', et = et_haz[key],
                     do_inhibit_axes=True,
                     )
        
        plt.show()
     
    hbt.figsize() 
    hbt.fontsize()
    
#%%%

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

# =============================================================================
#     Calculate the radial profile, in I/F units 
# =============================================================================
    
    # Take the radial profile of the superstack
    
    num_pts = 1000  # Number of radial points to use.
    
    profile_iof = {}  # Create a dictionary to save the radial profiles in
    
    if do_tunacan:
        if zoom == 2:
            num_pts = 300  # Number of radial points to use. Fewer for tunacan, since fewer pixels
        else:
            num_pts = 200
        
        # For tunacan, take two different radial profiles. One for an outer ring, and one for an inner.
        # My code here is so awful! This clearly should be a loop, and it should be made to be the same for both
        # tunacan and sunflower. Ugh!
    
        for i in range(len(thick_tuna_km)): # Iterate over tuna thickness:   
            is_good = (np.abs(backplanes[0]['Altitude_eq']) < thick_tuna_km[i])

            # Mask only the ring pixels in the 'mean' stack
            
            img_superstack_median_iof_masked = img_superstack_median_iof * is_good
            img_superstack_median_iof_masked[is_good == 0] = np.nan    # Set all pixel > vertical altitude to NaN

            # Mask only the ring the pixels in the 'median' stack
            
            img_superstack_mean_iof_masked = img_superstack_mean_iof * is_good
            img_superstack_mean_iof_masked[is_good == 0] = np.nan    # Set all pixel > vertical altitude to NaN

            (radius,  profile_iof_mean) = \
                                        get_radial_profile_backplane(img_superstack_mean_iof_masked,
                                                 plane_radius_superstack, method = 'mean', num_pts = num_pts, 
                                                 do_std=False)
                                        
            (radius,  profile_iof_median) = \
                                        get_radial_profile_backplane(img_superstack_median_iof_masked,
                                                 plane_radius_superstack, method = 'mean', num_pts = num_pts, 
                                                 do_std=False)
            # Save this profile to a dictionary
            
            profile_iof[i] = profile_iof_median
          
    else:  # Take a sunflower radial profile
        
                                             
        (radius,  profile_iof_median)   = get_radial_profile_backplane(img_superstack_median_iof,
                                             plane_radius_superstack, method = 'median', num_pts = num_pts)

        (radius_quadrant,  profile_iof_quadrant)   = get_radial_profile_backplane_quadrant(img_superstack_mean_iof,
                                             plane_radius_superstack, plane_longitude_superstack, method = 'mean', 
                                             num_pts = num_pts/4)
        if do_implant:
            (radius,  profile_merged_iof)   = get_radial_profile_backplane(img_merged_iof,
                                             plane_radius_superstack, method = 'median', num_pts = num_pts)

            
            
        profile_iof[0] = profile_iof_median  # Just pick one -- they both are similar

    # Calculate the bias level, crudely

    radius_max_km = 5000
        
    bin_radial_end = np.digitize(radius_max_km, radius)
    
    bias       = np.amin(profile_iof[0][0:bin_radial_end])

    bias_merged = bias 
# =============================================================================
# Fit a gaussian to the radial profile, if requested
# =============================================================================

    do_fit_profile = False
    
    if do_fit_profile:
        r_0         = 4000      # Inner radius to ignore in gauss fit, in km
        radius_ring = 9000      # Starting point for gaussfit for ring position, in km
        hw_ring     = 1000      # Starting point for ring halfwidth, in km
        
        
        bin_0 = hbt.x2bin(r_0, radius)
        x = radius[bin_0:]
        y = profile_iof[bin_0:]            
        
        popt,pcov = curve_fit(gaus,x,y,p0=[radius_ring,0,hw_ring])
    
# =============================================================================
# Plot the radial profile
# =============================================================================
    
    hbt.figsize((8,6))
    hbt.fontsize(12)
    
    if do_tunacan:
        for i,thick in enumerate(thick_tuna_km):
            plt.plot(radius, profile_iof[i] - np.nanmin(profile_iof[i][0:bin_radial_end]),
                 label = f'Median, pole=({ra_pole*hbt.r2d:.0f}°, {dec_pole*hbt.r2d:.0f}°), ' + 
                         f'dr={round(np.diff(radius)[0]):.0f} km, z/2={thick_tuna_km[i]} km')
        
    else:    
        plt.plot(radius, profile_iof[0] - bias,
             label = f'Median, pole = ({ra_pole*hbt.r2d:.0f}°, {dec_pole*hbt.r2d:.0f}°), ' + 
                     f'dr={round(np.diff(radius)[0]):.0f} km')
    
    do_plot_errorbars = False
    
    if do_plot_errorbars:
        plt.errorbar(radius, profile_iof - bias, yerr = profile_iof_std,
             label = f'Median, pole = ({ra_pole*hbt.r2d:.0f}, {dec_pole*hbt.r2d:.0f}) deg')
    
    plt.xlim(0,radius_max_km)
    if do_tunacan:
        plt.ylim((0,1e-5))
    else:
        plt.ylim((0,2e-6))
    
    do_plot_mean=False
    if do_plot_mean:            # Can plot the mean. But the median is much more useful
        plt.plot(radius, profile_iof_mean - bias, label = 'mean')
    
    if do_fit_profile:
        plt.plot(x,gaus(x,*popt),'ro:', marker = None, ls = '--', lw=1.5, 
                 label = f'Fit, radius={popt[1]:.0f} km, FWHM={2.35 * popt[2]:.0f} km')
    
    # FWHM = 2.35 * sigma: https://ned.ipac.caltech.edu/level5/Leo/Stats2_3.html
        
    plt.xlim((0, radius_max_km))
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
    plt.xlabel('Radius [km]')
    plt.ylabel('I/F')
    plt.legend(loc = 'upper right')
    plt.title(f'Radial profile, {frametype}, {str_stack} {str_reqid}')
    plt.show()


#%%%
# =============================================================================
# Plot the radial profile, for quadrants
# =============================================================================

    do_profile_quadrant = True
  
    if do_profile_quadrant:
      
        (radius_quadrant,  profile_iof_quadrant)   = get_radial_profile_backplane_quadrant(img_superstack_mean_iof,
                                             plane_radius_superstack, plane_longitude_superstack, method = 'median', 
                                             num_pts = num_pts/8)

        plt.plot(radius_quadrant, profile_iof_quadrant[0] - bias,
             label = f'Median, pole = ({ra_pole*hbt.r2d:.0f}°, {dec_pole*hbt.r2d:.0f}°), ' + 
                     f'dr={round(np.diff(radius)[0]):.0f} km', color='white')

        
        for i in range(4):
            plt.plot(radius_quadrant, profile_iof_quadrant[i] - bias, label=f'Quadrant {i}')

        plt.ylim((-1e-7, 2e-6))
        plt.xlim((0,10000))

        plt.xlim((0, radius_max_km))
        plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
        plt.xlabel('Radius [km]')
        plt.ylabel('I/F')
        plt.legend(loc = 'upper right')
        plt.title(f'Radial profile, {frametype}, {str_stack} {str_reqid}')
        plt.show()

#%%%
# =============================================================================
# Plot the radial profile, with implant
# =============================================================================
    
    if do_implant:
        hbt.figsize((8,6))
        hbt.fontsize(12)
        
        if do_tunacan:
            for i,thick in enumerate(thick_tuna_km):
                plt.plot(radius, profile_merged_iof[i] - np.nanmin(profile_iof[i][0:bin_radial_end]),
                      label = f'Median, pole=({ra_pole*hbt.r2d:.0f}°, {dec_pole*hbt.r2d:.0f}°), ' + 
                              f'dr={round(np.diff(radius)[0]):.0f} km, z/2={thick_tuna_km[i]} km')
            
        else:    
            plt.plot(radius, profile_merged_iof - bias_merged,
                 label = f'Median, pole = ({ra_pole*hbt.r2d:.0f}°, {dec_pole*hbt.r2d:.0f}°), ' + 
                         f'dr={round(np.diff(radius)[0]):.0f} km')
        
        do_plot_errorbars = False
        
        if do_plot_errorbars:
            plt.errorbar(radius, profile_iof - bias, yerr = profile_iof_std,
                 label = f'Median, pole = ({ra_pole*hbt.r2d:.0f}, {dec_pole*hbt.r2d:.0f}) deg')
            
        plt.xlim(0,50000)
        if do_tunacan:
            plt.ylim((0,1e-5))
        else:
            plt.ylim((0,2e-6))
        
        do_plot_mean=False
        if do_plot_mean:            # Can plot the mean. But the median is much more useful
            plt.plot(radius, profile_iof_mean - bias, label = 'mean')
        
        if do_fit_profile:
            plt.plot(x,gaus(x,*popt),'ro:', marker = None, ls = '--', lw=1.5, 
                     label = f'Fit, radius={popt[1]:.0f} km, FWHM={2.35 * popt[2]:.0f} km')
        
        # FWHM = 2.35 * sigma: https://ned.ipac.caltech.edu/level5/Leo/Stats2_3.html
            
        plt.xlim((0, radius_max_km))
        plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
        plt.xlabel('Radius [km]')
        plt.ylabel('I/F')
        plt.legend(loc = 'upper right')
        plt.title(f'With ring implant {file_ring_short}, I/F={iof_max_ring:.0e}')
        plt.show()
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
        

