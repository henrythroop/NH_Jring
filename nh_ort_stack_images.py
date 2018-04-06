#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 11:43:11 2018

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

# HBT imports

import hbt
from image_stack import image_stack

def make_img_superstack(stack, img, img_field):
    
    """
    This makes a superstack image. The output image has resolution matching the lowest resolution 
    of all the input images (that is, the closest to MU69).
    
    This is a standalone function -- not part of any class.
    
    Parameters
    -----
    
    stack:
        A dictionary of `image_stack` objects.
        
    img:
        A dictonary of 2D images, created by flattening the stacks above. The dictionary 
        keys of `img` and `stack` must be the same.
    
    Optional parameters
    -----
    
    zoom:
        The zoom value to use when stacking.
        
    index_reference:
        Index to which of the final images to use as the output scale.
        (Alternatively, we might allow the dictionary to be passed)
        
    """
    
    keys = list(stack.keys())
    
    pixscale_km = {}
    for key in keys:
        pixscale_km[key] = stack[key].pixscale_x_km
        
    # Search for the lowest-res image, and output to that.
    
    pixscale_km_out = min(pixscale_km.values())
    
    # Look up what the index is for this image
    
    key_out = keys[np.where(pixscale_km_out == np.array(list(pixscale_km.values())))[0][0]]
    
    pixscale = stack[key_out].pixscale_x_km  # This is the pixel scale of the final stack.
                                                                 # Zoom has not been applied yet.
    size_out = np.shape(img[key_out])
    
    img_rescale_3d = np.zeros((len(keys),size_out[0],size_out[1]))
    img_rescale = {}
    for i,key in enumerate(keys):
        
        magfac = stack[key].pixscale_x_km / pixscale_km_out
        arr = scipy.ndimage.zoom(img[key] - img_field, magfac)
        
        size_in = np.shape(arr)
        edge_left = int( (size_in[0]-size_out[0])/2)
        edge_top  = int( (size_in[1]-size_out[1])/2)
        
        arr_out = arr[ edge_left : edge_left+size_out[0], 
                       edge_top  : edge_left+size_out[1] ] 
        img_rescale[reqid_i] = arr_out
        img_rescale_3d[i,:,:] = arr_out

        img_rescale_median = np.median(img_rescale_3d, axis=0)  
        img_rescale_mean   = np.mean(img_rescale_3d, axis=0)  
    
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

if (True):
    
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.
    
    plt.set_cmap('Greys_r')

    zoom = 2      # How much to magnify images by before shifting. 4 (ie, 1x1 expands to 4x4) is typical
                  # 1 is faster; 4 is slower but better.
    
#    name_ort = 'ORT1'
#    name_ort = 'ORT2_OPNAV'
    name_ort = 'ORT3'
    initials_user = 'HBT'
    dir_data = '/Users/throop/Data'

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
 
    if (name_ort == 'ORT3'):
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
    
    # Start up SPICE if needed
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
    
    hbt.figsize((12,12))
    
    # Load and stack the field images
    
    do_force = False   # If set, reload the stacks from individual frames, rather than restoring from a pickle file.
    
    stack_field = image_stack(os.path.join(dir_images, reqid_field),   do_force=do_force, do_save=do_force)

    stack_haz = {}
    
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
    
    # Flatten the main hazard stacks
    # When we do this, the output image is shifted around within the padding amount, so that *all* 
    # ouput images have the same size. So, maybe flatten should return img *and* wcs?
    # And the expectation would be that all of these individual WCS would match. That would be the point.
    
    img_haz = {}
    wcs_haz = {}
    img_haz_diff = {}
    
    for reqid_i in reqids_haz:
        (img_haz[reqid_i], wcs_haz[reqid_i])  =\
              stack_haz[reqid_i].flatten(do_subpixel=False,  method='median',zoom=zoom, padding=pad)
        img_haz_diff[reqid_i] = img_haz[reqid_i] - img_field
        
        plot_img_wcs(img_haz[reqid_i], wcs_haz[reqid_i], title = reqid_i)
        
    # Plot the trimmed, flattened images. This is just for show. They are not scaled very well for ring search.
    
    plt.imshow(stretch(img_field))
    for reqid_i in reqids_haz:        
        diff_trim = hbt.trim_image(img_haz_diff[reqid_i])
        plt.imshow(stretch(diff_trim))
        plt.title(reqid_i)
        plt.show()
        
        # Save as FITS
    
        file_out = os.path.join(dir_out, '{}_z{}.fits'.format(reqid_i, zoom))
        hdu = fits.PrimaryHDU(stretch(diff_trim))
        hdu.writeto(file_out, overwrite=True)
        print(f'Wrote: {file_out}')
    
# =============================================================================
#      Calculate and display radial profiles
# =============================================================================
    
    hbt.figsize((10,8))
    hbt.set_fontsize(12)
    pos =  np.array(np.shape(img_field))/2  # MU69 will be at the center of this array
    for reqid_i in reqids_haz:
        (radius, profile_dn) = get_radial_profile_circular(img_haz[reqid_i] - img_field, pos=pos, width=2)
        radius_km = radius * stack_haz[reqid_i].pixscale_x_km
        plt.plot(radius_km, profile_dn, label=reqid_i)
    plt.xlim(0,50_000)
    plt.xlabel('Expanded Pixels')
    plt.ylabel('DN Median')
    plt.legend(loc = 'upper right')
    plt.show()
    
    # Convert individual radial profiles from DN to I/F
    
    width = 2  # Bin width for radial profiles
    
    for reqid_i in reqids_haz:
        
        # Take the radial profile
        
        (radius, profile_dn) = get_radial_profile_circular(img_haz[reqid_i] - img_field, pos=pos, width=width)
        radius_km = radius * stack_haz[reqid_i].pixscale_x_km / zoom
        
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
    plt.ylabel('I/F Median')
    plt.ylim((-2e-7,5e-7))
    plt.legend()
    plt.title(f'{name_ort}, binwidth={width} pix, zoom={zoom}')
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
    plt.show()

# =============================================================================
#     Make a 'superstack' -- ie stack of stacks, put on the same spatial scale and stacked
# =============================================================================
    
    (img_superstack_mean, img_superstack_median) = make_img_superstack(stack_haz, img_haz, img_field)
  
    # Create, display, and save the median superstack
    
    plt.imshow( stretch(img_superstack_median))
    plt.title(f'restretched and medianed, {name_ort}')
    plt.show()
    file_out = os.path.join(dir_out, f'img_superstack_median_{name_ort}.fits')
    hdu = fits.PrimaryHDU(img_superstack_median)
    hdu.writeto(file_out, overwrite=True)
    print(f'Wrote: {file_out}')
 
    # Create, display, and save the mean superstack

    plt.imshow( stretch(img_superstack_mean))
    plt.title(f'restretched and mean, {name_ort}')
    plt.show()    
    file_out = os.path.join(dir_out, f'img_superstack_mean_{name_ort}.fits')
    hdu = fits.PrimaryHDU(img_superstack_mean)
    hdu.writeto(file_out, overwrite=True)
    print(f'Wrote: {file_out}')
       
#    print(f'zooming by {magfac}, to size {np.shape(img_haz_rescale[reqid_i])}')
    
# =============================================================================
#     Make a radial profile, in I/F units 
# =============================================================================
    
    # Apply Hal's conversion formula from p. 7, to compute I/F and print it.

    # Calculate the MU69-Sun distance, in AU (or look it up).         
    et = stack_haz[reqids_haz[-1]].t['et'][0] # ET for the final image stack
    (st,lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    r_nh_mu69 = sp.vnorm(st[0:3]) * km2au # NH distance, in AU
    (st,lt) = sp.spkezr('MU69', et, 'J2000', 'LT', 'Sun')
    r_sun_mu69 = sp.vnorm(st[0:3]) * km2au # NH distance, in AU
    pixscale_km =  (r_nh_mu69/km2au * (0.3*hbt.d2r / 256)) / zoom # km per pix (assuming LORRI 4x4)
    TEXP = stack_haz[reqid_i].t['exptime'][0]

    I_median = img_superstack_median / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    I_mean   = img_superstack_mean   / TEXP / RSOLAR   # Could use RSOLAR, RJUPITER, or RPLUTO. Dfft spectra.
    
    img_superstack_median_iof = math.pi * I_median * r_sun_mu69**2 / F_solar # Equation from Hal's paper
    img_superstack_mean_iof   = math.pi * I_mean   * r_sun_mu69**2 / F_solar # Equation from Hal's paper
    
    binwidth = 1
    (radius_median,profile_iof_median) = get_radial_profile_circular(img_superstack_median_iof, pos=pos, width=binwidth)
    (radius_mean,  profile_iof_mean)   = get_radial_profile_circular(img_superstack_mean_iof,   pos=pos, width=binwidth)
    
    hbt.figsize((8,6))
    radius_profile_km = radius_median*pixscale_km
    plt.plot(radius_profile_km, profile_iof_median,label = 'All curves, Median')
    plt.plot(radius_profile_km, profile_iof_mean,  label = 'All curves, Mean')
    plt.xlim((0,20000))
    plt.ylim((-0.5e-7,3e-7))
    plt.legend(loc = 'upper right')
    plt.ylabel('I/F')
    plt.title(f'{name_ort}, binwidth={binwidth} pix, zoom={zoom}')
    plt.xlabel('Radius [km]')
    plt.show()
    
# =============================================================================
#      Write the ring profile to a text file as per wiki
#      https://www.spaceops.swri.org/nh/wiki/index.php/KBO/Hazards/Pipeline/Detection-to-Ring  
# =============================================================================
    version = 1
    
    file_out = os.path.join(dir_out, f'{initials_user.lower()}_{name_ort.lower()}_v{version}.ring')
    lun = open(file_out,"w")
    for i in range(len(profile_iof_median)):
        lun.write('{:10.3f} {:11.3e}\n'.format(radius_profile_km[i], profile_iof_median[i]))
    lun.close()
    print(f'Wrote: {file_out}')
    print(f' scp {file_out} ixion:\~/astrometry' )  # We don't copy it, but we put up the string so user can.

    