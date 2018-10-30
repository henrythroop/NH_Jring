#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 11:43:11 2018

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

# HBT imports

import hbt

from   matplotlib.figure import Figure
from   get_radial_profile_circular import get_radial_profile_circular
from   get_radial_profile_backplane import get_radial_profile_backplane
from   plot_img_wcs import plot_img_wcs
from   image_stack import image_stack
from   compute_backplanes import compute_backplanes
from   scipy.optimize import curve_fit

# Define a gaussian function, for the fit
    
def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
    

def nh_ort_make_superstack(stack, img, img_field, 
                           name_stack_base, 
                           do_save_all=True, do_backplanes=True, dir='', str_name='', 
                           do_center=True, 
                           do_wcs = False):
    
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
        
    do_wcs:
        Boolean. If set, will include a WCS in the outut tuple.
        
    do_backplanes:
        Boolean. If set, will return a full set of backplanes for the superstacks.
        
    """
    
    cmap = 'plasma'
    
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
        arr       = scipy.ndimage.zoom(img[key] - img_field, magfac)
        
        size_in   = np.shape(arr)
        edge_left = int( (size_in[0]-size_out[0])/2)
        edge_top  = int( (size_in[1]-size_out[1])/2)
        
        arr_out = arr[ edge_left : edge_left+size_out[0], 
                       edge_top  : edge_left+size_out[1] ]
        
        # If requested, do a centroid on the central region of the image. 
        # Then shift the entire image based on this, to put brightest object at center of frame.
    
        if do_center:
            extract = arr_out[int(size_out[0]/2)-50:int(size_out[0]/2+50),    # 100x100 pixel region. Hardcoded, ugh.
                              int(size_out[0]/2)-50:int(size_out[0]/2+50)]
            shift_x = 50 - hbt.wheremax(np.sum(extract, axis=0))
            shift_y = 50 - hbt.wheremax(np.sum(extract, axis=1))  # vertical axis on screen. 
                                                                  #  Axis=0. - means value is too large. 
            arr_out = np.roll(np.roll(arr_out, shift_x, axis=1), shift_y, axis=0)
            print(f'Rolling by {shift_x}, {shift_y} to align {key} into superstack')
            
        img_rescale[reqid_i]  = arr_out
        img_rescale_3d[i,:,:] = arr_out

        # Write the scaled stacks to disk, if requested
        
        if do_save_all:
            file_out = f'stack_{key}_{str_name}_rescaled_hbt.fits'
            path_out = os.path.join(dir_out, file_out)
            hdu = fits.PrimaryHDU(arr_out) 
            hdu.writeto(path_out, overwrite=True)
            print(f'Wrote: {path_out}')
            
    img_rescale_median = np.median(img_rescale_3d, axis=0)  
    img_rescale_mean   = np.mean(img_rescale_3d, axis=0)
    
    # Get a WCS system, if requested. Just grab the one from the stack that we used.
    # However, there may be some small mis-alignment issues. 
    
    wcs = stack[name_stack_base].t['wcs'][0]
    
    # Now adjust the WCS so the center of the image is is at dx/2, dy/2, *and* RA/Dec is position of MU69.
    
    # Look up RA, Dec of MU69 for that date
    
    et = stack_haz[name_stack_base].t['et'][0]    # name_stack_base indicates which of the stacks is used for 
    (st,lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    (range_, ra, dec) = sp.recrad(st[0:3])
    
    wcs.wcs.crval = (ra*hbt.r2d, dec*hbt.r2d)
    wcs.wcs.crpix = (hbt.sizex(img_rescale_mean)/2, hbt.sizey(img_rescale_mean)/2)
    
    # Put this WCS info in a header, which we will write to FITS file
    
    header = wcs.to_header()
    
    # Also, put the correct ET into the FITS header.
    
    et = stack_haz[name_stack_base].t['et'][0]
    header['SPCSCET'] = f'{et}'
    header['SPCUTCID'] = sp.et2utc(et, 'C', 0)
    
    # Make a 'small' version of the extracted image. It will be just much easier to deal with. Put WCS info it it too.
    
    x0 = int(hbt.sizex(img_rescale_mean) * 0.375) # This reduces size from 1600 pixels to 400 pixels
    x1 = int(hbt.sizex(img_rescale_mean) * 0.625)
    y0 = int(hbt.sizey(img_rescale_mean) * 0.375)
    y1 = int(hbt.sizey(img_rescale_mean) * 0.625)
    
    img_rescale_mean_sm = img_rescale_mean[x0:x1, y0:y1]
    img_rescale_median_sm = img_rescale_median[x0:x1, y0:y1]
    
    wcs_sm = wcs.copy()
    wcs_sm.wcs.crpix = (hbt.sizex(img_rescale_mean_sm)/2, hbt.sizey(img_rescale_mean_sm)/2)

    header_sm = wcs_sm.to_header()
    header_sm['SPCSCET'] = f'{et}'
    header_sm['SPCUTCID'] = sp.et2utc(et, 'C', 0)

    # Write the superstacks (with WCS) to disk, if requested
    
    if do_save_all:
        
        # Write full-size images
        
        file_out = f'superstack_n{len(keys)}_{str_name}_median_wcs_hbt.fits'
        hdu = fits.PrimaryHDU(img_rescale_median, header=header)
        path_out = os.path.join(dir_out, file_out)
        hdu.writeto(path_out, overwrite=True)
        print(f'Wrote: {path_out}')
      
        file_out = f'superstack_{str_name}_mean_wcs_hbt.fits'
        hdu = fits.PrimaryHDU(img_rescale_mean, header=header)
        path_out = os.path.join(dir_out, file_out)
        hdu.writeto(path_out, overwrite=True)
        print(f'Wrote: {path_out}')

        # Write small images
        
        file_out_superstack_median_sm = f'superstack_{str_name}_median_wcs_sm_hbt.fits'
        hdu = fits.PrimaryHDU(img_rescale_median_sm, header=header_sm)
        path_out = os.path.join(dir_out, file_out_superstack_median_sm)
        hdu.writeto(path_out, overwrite=True)
        print(f'Wrote: {path_out}')
      
        file_out_superstack_mean_sm = f'superstack_{str_name}_mean_wcs_sm_hbt.fits'
        hdu = fits.PrimaryHDU(img_rescale_mean_sm, header=header_sm)
        path_out = os.path.join(dir_out, file_out_superstack_mean_sm)
        hdu.writeto(path_out, overwrite=True)
        print(f'Wrote: {path_out}')
        
        if do_backplanes:
            # Now that we have written a FITS file to disk, compute backplanes for it
            # ** for compute_backplanes to work, use_sky_plane = False must be set in that function.
            
            frame = '2014_MU69_SUNFLOWER_ROT'
            
            print(f'Computing backplanes for file {path_out}')
            planes = compute_backplanes(path_out, 'MU69', frame, 'New Horizons', angle3=30*hbt.d2r)
        
        # And plot the Radius backplane
        
#        plt.imshow(planes[0]['Radius_eq'])

    if do_backplanes:        
        return((img_rescale_mean, img_rescale_median, planes))
    else:
        return((img_rescale_mean, img_rescale_median))
 
# =============================================================================
# End of function definition
# =============================================================================
    
# =============================================================================
# One-off code (not a function) that does all the image stacking for the MU69 ORT's.
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

    zoom = 4     # How much to magnify images by before shifting. 4 (ie, 1x1 expands to 4x4) is typical
                  # 1 is faster; 4 is slower but better.

    width = 1  # Bin width for radial profiles
    
#    name_ort = 'ORT1'
#    name_ort = 'ORT2_OPNAV'
#    name_ort = 'ORT3'
#    name_ort = 'ORT4'
    
    name_ort = 'MU69_Approach'  # This is the actual encounter -- not ORT!
    initials_user = 'HBT'
    dir_data = '/Users/throop/Data'
    
    a_xy = (1,1)   # Projected ellipticity of ring

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

    if (name_ort == 'MU69_Approach'):
        dir_images    = os.path.join(dir_data, name_ort, 'throop', 'backplaned')
        dir_out       = os.path.join(dir_data, name_ort, 'throop', 'stacks')
        reqids_haz  = ['KALR_MU69_OpNav_L4_2018284', 'KALR_MU69_OpNav_L4_2018287', 'KALR_MU69_OpNav_L4_2018_298']
        reqids_haz  = ['KALR_MU69_OpNav_L4_2018298','KALR_MU69_OpNav_L4_2018301']
        # reqids_haz  = ['KALR_MU69_OpNav_L4_2018301']
        # reqids_haz  = ['KALR_MU69_OpNav_L4_2018287']
        # reqids_haz  = ['KALR_MU69_OpNav_L4_2018284']
        reqid_field = 'K1LR_MU69ApprField_115d_L2_2017264'
        a_xy = (1, math.cos(hbt.d2r * 30))
    
    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
    
    hbt.figsize((12,12))
    
    # Make sure output dir exists
    
    os.makedirs(dir_out, exist_ok=True)
    
    # Load and stack the field images
    
    do_force = False   # If set, reload the stacks from individual frames, rather than restoring from a pickle file.
    
    stack_field = image_stack(os.path.join(dir_images, reqid_field),   do_force=do_force, do_save=do_force)

    stack_haz = {}
    
    # Load and stack the data images
    
    for reqid_i in reqids_haz:
        stack_haz[reqid_i] = image_stack(os.path.join(dir_images, reqid_i), do_force=do_force, do_save=do_force)
            
    # Set the rough position of MU69
    
    et = stack_haz[reqids_haz[0]].t['et'][0] # Look up ET for first image in the Hazard stack
    (st, lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    vec = st[0:3]
    (_, ra, dec) = sp.recrad(vec)
    radec_mu69 = (ra, dec) # Keep in radians
    
    # Align the frames
    
    stack_field.align(method = 'wcs', center = (radec_mu69))
    for reqid_i in reqids_haz:
        stack_haz[reqid_i].align(method  = 'wcs', center = (radec_mu69))  # In each individual stack, align all images
    
    # Calc the padding required. This can only be done after the images are loaded and aligned.

    pad_field = stack_field.calc_padding()[0]
    pad_haz = []
    for reqid_i in reqids_haz:
        pad_haz.append(stack_haz[reqid_i].calc_padding()[0])
    pad_haz.append(pad_field)
        
    pad = max(pad_haz)
       
    # Flatten the stacks into single output images
    # If we get an error here, it is probably due to a too-small 'pad' value.
    # This step can be very slow. Right now results are not saved. In theory they could be.
    
    # Flatten the field stack
    
    (img_field, wcs_field) = stack_field.flatten(do_subpixel=False, method='median',zoom=zoom, padding=pad)

    # Verify that the WCS has been properly preserved after flattening here.
    
    hbt.figsize((12,12)) # This is ignored, for some reason.
    
    plot_img_wcs(img_field, wcs_field, title = 'Field Stack')
    
#%%%
    
    # Loop over and flatten the main hazard stacks (HAZ00, HAZ01, etc)
    # When we do this, the output image is shifted around within the padding amount, so that *all* 
    # ouput images have the same size. So, maybe flatten should return img *and* wcs?
    # And the expectation would be that all of these individual WCS would match. That would be the point.
    
# =============================================================================
#     Flatten the main image stacks
# =============================================================================
    
    img_haz = {}
    wcs_haz = {}
    img_haz_diff = {}

    do_plot_individual_stacks = True
    
    for reqid_i in reqids_haz:
        (img_haz[reqid_i], wcs_haz[reqid_i])  =\
              stack_haz[reqid_i].flatten(do_subpixel=False,  method='median',zoom=zoom, padding=pad)
        img_haz_diff[reqid_i] = img_haz[reqid_i] - img_field
        
        if do_plot_individual_stacks:
            plot_img_wcs(img_haz[reqid_i], wcs_haz[reqid_i], title = reqid_i)
        
# =============================================================================
#     Plot the stacks, differenced with the field image.
# =============================================================================
    
    hbt.figsize((12,12))
    hbt.set_fontsize(12)
    
    plt.imshow(stretch(img_field), origin='lower')
    for i,reqid_i in enumerate(reqids_haz):        
        diff_trim = hbt.trim_image(img_haz_diff[reqid_i])
        plt.imshow(stretch(diff_trim), origin='lower')
        plt.title(f'Difference stack {i}/{len(reqids_haz)}: {reqid_i}')
        plt.show()
        
        # Save as FITS
    
        file_out = os.path.join(dir_out, 'stack_{}_{}_z{}_hbt.fits'.format(reqid_i, name_ort, zoom))
        hdu = fits.PrimaryHDU(stretch(diff_trim))
        hdu.writeto(file_out, overwrite=True)
        print(f'Wrote: {file_out}')

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
    
    # Create and save the median superstack

    str_name = f'{name_ort}_z{zoom}'
    (img_superstack_mean, img_superstack_median) = nh_ort_make_superstack(stack_haz, img_haz, img_field, 
                                                                          name_stack_base, 
                                                                          do_save_all=True, dir=dir_out,
                                                                          str_name = str_name, do_center=True,
                                                                          do_backplanes=False)

    # Make the superstack, and return backplanes.
    
    (img_superstack_mean, img_superstack_median,backplanes) = nh_ort_make_superstack(stack_haz, img_haz, img_field, 
                                                                          name_stack_base, 
                                                                          do_save_all=True, dir=dir_out,
                                                                          str_name = str_name, do_center=True,
                                                                          do_backplanes=True)

    
    # Plot the superstack to screen
    
    hbt.figsize(14,14)
    plt.subplot(1,2,1)
    plt.imshow( stretch(img_superstack_median), origin='lower', cmap=cmap_superstack)
    plt.title(f'restretched and medianed, {name_ort}, n={len(keys)}')

    plt.subplot(1,2,2)
    plt.imshow( stretch(img_superstack_mean), origin='lower', cmap=cmap_superstack)
    plt.title(f'restretched and mean, {name_ort}, n={len(keys)}')

    plt.show()
    hbt.figsize()
 
# =============================================================================
# Now start to calculate radial profiles
# =============================================================================

# Create a backplane for this superstack
# Look up the Radius backplane, in the chosen stack. Extracting the RADIUS_EQ backplane.
    
    file_backplanes = stack_haz[name_stack_base].t['filename'][0]
    f = fits.open(file_backplanes)
    plane_radius_orig = f['RADIUS_EQ'].data  # Read the RADIUS backplane value
    f.close()
    
    # Do a really brute-force thing: center this backplane, by just centroiding it. 
    # I should be able to center it better, using shift coordinates from above.
    # However, I don't actually ever compute these. Even in assembling the superstacks, I just center on the 
    # brightest thing in the frame, rather than follow MU69's motion. I should fix that, but for now this works.
    
    plane_radius  = scipy.ndimage.zoom(plane_radius_orig, zoom)
    wheremin_x    = hbt.wheremin(np.amin(plane_radius,axis=1))
    wheremin_y    = hbt.wheremin(np.amin(plane_radius,axis=0))
    plane_radius  = np.roll(plane_radius, 512-wheremin_y, axis=1)
    plane_radius  = np.roll(plane_radius, 512-wheremin_x, axis=0)

    # Pad this so it is the same size as the superstack
    
    dpad = int( (hbt.sizex(img_superstack_mean) - hbt.sizex(plane_radius) )/2)
    
    plane_radius = np.pad(plane_radius, dpad, mode='maximum')    
        
    name_target = 'MU69'
    frame = '2014_MU69_SUNFLOWER_ROT'
    name_observer = 'New Horizons'

    file_out_superstack_median_sm = f'superstack_{str_name}_median_wcs_sm_hbt.fits'

    # Look up the pole position of MU69, based on the loaded frame kernel
    # That is, the RA + Dec of the pole, in degrees, as derived from NH_ORT_FIND_RING_POLE.PY.
    # I put this value into the .tf frame file for MU69, and grab it here.
    # This is just for labeling plots -- nothing else.
    
    vec = [0, -1, 0]                                # -Y vector is the rotational pole
    mx_frame_j2k =  sp.pxform(frame, 'J2000', et)
    vec_j2k = sp.mxv(mx_frame_j2k, vec)
    (_, ra_pole, dec_pole) = sp.recrad(vec_j2k)

# =============================================================================
#     Make a radial profile, in I/F units 
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
    TEXP        = stack_haz[reqid_i].t['exptime'][0]

    I_median = img_superstack_median / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    I_mean   = img_superstack_mean   / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    
    img_superstack_median_iof = math.pi * I_median * r_sun_mu69**2 / F_solar # Equation from Hal's paper
    img_superstack_mean_iof   = math.pi * I_mean   * r_sun_mu69**2 / F_solar # Equation from Hal's paper
    
    # Take the radial profile of the superstack
    
    num_pts = 800  # Number of radial points to use.
    
    (radius,  profile_iof_mean)   = get_radial_profile_backplane(img_superstack_mean_iof,
                                         plane_radius, method = 'mean', num_pts = num_pts)
                                         
    (radius,  profile_iof_median)   = get_radial_profile_backplane(img_superstack_median_iof,
                                         plane_radius, method = 'median', num_pts = num_pts)

    profile_iof = profile_iof_median  # Just pick one -- they both are similar
    
  # Fit a gaussian to the radial profile. Set the initial guesses.
    
    r_0         = 4000      # Inner radius to ignore in gauss fit, in km
    radius_ring = 9000      # Starting point for gassfit for ring position, in km
    hw_ring     = 1000      # Starting point for ring halfwidth, in km
    
    bin_0 = hbt.x2bin(r_0, radius)
    x = radius[bin_0:]
    y = profile_iof[bin_0:]            
    
    popt,pcov = curve_fit(gaus,x,y,p0=[radius_ring,0,hw_ring])

    hbt.figsize((8,6))
 
    bias = 0.000004 +  7.5e-7
    plt.plot(radius, profile_iof - bias, 
             label = f'OpNav, backplaned, pole = ({ra_pole*hbt.r2d:.0f}, {dec_pole*hbt.r2d:.0f}) deg')
    
    do_plot_fit = False
    
    if do_plot_fit:
        plt.plot(x,gaus(x,*popt),'ro:', marker = None, ls = '--', lw=1.5, 
                 label = f'Fit, radius={popt[1]:.0f} km, FWHM={2.35 * popt[2]:.0f} km')
    
    # FWHM = 2.35 * sigma: https://ned.ipac.caltech.edu/level5/Leo/Stats2_3.html
    
    if (len(reqids_haz)) == 1:
        str_reqid = reqids_haz[0]
    else:    
        str_reqid = reqids_haz[0] + f' + {len(reqids_haz)-1} more'

    plt.ylim((0, 2e-6))
    plt.xlim((0, 20000))
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))
    plt.xlabel('Radius [km]')
    plt.ylabel('I/F')
    plt.legend(loc = 'upper right')
    plt.title(f'Radial profile, superstack, backplane, {str_reqid}')
    plt.show()
    
# =============================================================================
#      Write the ring profile to a text file as per wiki
#      https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Detection-to-Ring  
# =============================================================================
    
    version = 1
    
    file_out = os.path.join(dir_out, f'{initials_user.lower()}_{name_ort.lower()}_v{version}.ring')
    lun = open(file_out,"w")
    for i in range(len(profile_iof)):
        lun.write('{:10.3f} {:11.3e}\n'.format(radius[i], profile_iof[i]))
    lun.close()
    
    do_ring_out = False

    if do_ring_out:
        print(f'Wrote: {file_out}')
        print(f' scp {file_out} ixion:\~/astrometry' )  # We don't copy it, but we put up the string so user can.


# =============================================================================
# Make an image overlaying a radius mask with the ring. This is just to visualize if the geometry is close.
# =============================================================================
    
    do_ring_mask_plot = False
    
    if do_ring_mask_plot:
        sigma_fit = popt[2]
        radius_fit = popt[1]
        
        radius_fit = 10000
        sigma_fit   = 2000
        
        
        hbt.figsize(12,12)
        alpha = 0.3
    
        plt.subplot(1,2,1)    
        dpad = (hbt.sizex(img_superstack_mean) - hbt.sizex(plane_radius))/2
        plane_radius = np.pad(plane_radius, int(dpad), 'constant')
        plt.imshow(stretch( (plane_radius < 120000) * img_superstack_mean)[600:1000, 600:1000], origin='lower') 
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.imshow( (plane_radius < 0)[600:1000, 600:1000], cmap = 'Reds', alpha=alpha,
                     origin='lower')
        plt.title('Superstack')
        
        plt.subplot(1,2,2)    
        dpad = (hbt.sizex(img_superstack_mean) - hbt.sizex(plane_radius))/2
        plane_radius = np.pad(plane_radius, int(dpad), 'constant')
        plt.imshow(stretch( (plane_radius < 120000) * img_superstack_mean)[600:1000, 600:1000], origin='lower') 
        
        mask = ((plane_radius < (radius_fit+sigma_fit/2)) & (plane_radius > (radius_fit-sigma_fit/2)))
        
            
        plt.imshow( ((plane_radius < (radius_fit+sigma_fit/2)) & 
                     (plane_radius > (radius_fit-sigma_fit/2)))[600:1000, 600:1000], cmap = 'Reds', alpha=alpha,
                     origin='lower')
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.title('Superstack + derived ring image')
        plt.show()
    
# =============================================================================
# Now, as a one-off for ORT3, try to measure the I/F of this edge-on boxy looking ring
# =============================================================================

    if (name_ort == 'ORT3'):
#        ring_extract_large      = img_superstack_mean_iof[325:440,310:500]
        
        hbt.figsize((8,8))
        hbt.fontsize(12)
        center_x = int(hbt.sizex(img_superstack_mean)/2)
        center_y = int(hbt.sizex(img_superstack_mean)/2)
        
        pos_mu69 = ( center_x, center_y )   # x, y 
        dy = 200
        dx = 200
        img_extract_large      = img_superstack_median_iof[pos_mu69[1]-int(dy/2):pos_mu69[1]+int(dy/2),
                                                         pos_mu69[0]-int(dx/2):pos_mu69[0]+int(dx/2)]  
        
                                                                # MU69 is centered on 100, 100

        mask = hbt.dist_center(dx) > 0  # Mask anything this far from MU69 as True
        m = img_extract_large.copy()
        m[mask == False] = np.nan
        plt.imshow(stretch(m), origin='lower')
        
        img_extract_large = m.copy()     # Copy the masked image back to the full extracted image
        
        plt.imshow(stretch(img_extract_large), origin='lower', cmap='plasma')
        plt.show()
        
        plt.imshow(stretch(img_superstack_mean_iof), origin='lower')

        width_ring = 100
        height_ring = 30
        padding_x = 20         # Padding on each edge (half-padding)
        padding_y = 30
        
        img_extract_yprofile   = img_extract_large[
                          int(dy/2 - height_ring/2-padding_y):int(dy/2+height_ring/2+padding_y),
                          int(dx/2 - width_ring/2):int(dx/2+width_ring/2)]
        plt.imshow(stretch(img_extract_yprofile), origin='lower') 
        plt.show()
                          
        img_extract_xprofile   = img_extract_large[
                          int(dy/2 - height_ring/2):int(dy/2+height_ring/2),
                          int(dx/2 - width_ring/2-padding_x):int(dx/2+width_ring/2+padding_x)]
        plt.imshow(stretch(img_extract_xprofile), origin='lower')
        plt.show()
                
        # Make profile in X dir
        
        profile_x = np.nanmedian(img_extract_xprofile, axis=0)
    
        x_km = np.arange(len(profile_x)) * pixscale_km
        x_km -= np.mean(x_km)
        
        plt.subplot(2,1,1)
        plt.plot(x_km, profile_x, color='pink', lw=5)
        plt.xlabel('km')
        plt.title('Profile, X dir')
        plt.ylabel('I/F')    
        plt.ylim((-4e-7, 6e-7))
            
        # Make profile in Y dir
        
        profile_y = np.nanmedian(img_extract_yprofile, axis=1)
        
        y_km = np.arange(len(profile_y)) * pixscale_km
        y_km -= np.mean(y_km)

        plt.subplot(2,1,2)        
        plt.plot(y_km, profile_y, color='pink', lw=5)

        plt.xlabel('km')
        plt.title('Profile, Y dir')
        plt.ylabel('I/F')
    
        plt.ylim((-4e-7, 6e-7))
        plt.tight_layout()
        plt.show()
        

### Testing for ORT3
        
        path_out = '/Users/throop/Data/ORT4/superstack_ORT4_z1_mean_wcs_hbt.fits'
        planes = compute_backplanes(path_out, 'MU69', frame, 'New Horizons')
        plt.imshow(planes[0]['Radius_eq'])
        
       