#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 16:31:13 2018

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
from   astropy.wcs import WCS

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
from   wcs_translate_pix import wcs_translate_pix, wcs_zoom

# This is a one-off code to load some stacks, and check them to make sure that MU69 is centered properly.
# Some stacks seem to have a centering error, so I will investigate it here.

# =============================================================================
# Define a test function which is called when we run this file. This is just an example of using the class.
# =============================================================================
    
if (__name__ == '__main__'):

    #%%%

    # Demo function to try out the stack functionality

    plt.set_cmap('Greys_r')
    
    hbt.figsize((10,10))
    
    # Start up SPICE if needed
    
    hbt.unload_kernels_all()
    
    if (sp.ktotal('ALL') == 0):
        sp.furnsh('kernels_kem_prime.tm')
        
    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) 
    
    # Load a stack
    
    do_force = False
    do_save = True
    
    dirs = [
       # '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018331', # Now working!
        # '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018332', # Now working!
        # '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018334', # Now working!
        '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018335', # Now working!
        # '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018258',
       
       # '/Users/throop/Data/MU69_Approach/throop/backplaned/K1LR_MU69ApprField_115d_L2_2017264'
       ]
 
    #%%%
    
    for dir in dirs:
        
        str_reqid = dir.split('/')[-1]
        
        stack = image_stack(dir, do_force=do_force)
        
        # Plot a set of images from it
        
        i = 1
        num_images_plot = 4       # Max number of images to plot
        num_images = len(stack.t) # Number of images in this sequence 
        
        hbt.figsize((15,15))
        
        for i in range(num_images_plot):
            plt.subplot(2,int(np.ceil(num_images_plot/2)),i+1)
            filename_short = os.path.basename(stack.t['filename'][i])
            plot_img_wcs(stack.image_single(i), stack.t['wcs'][i], 
                           name_target='MU69', name_observer='New Horizons', et = stack.t['et'][i],
                           title = f'{str_reqid}: {i}/{num_images}',
                           width=50, do_show=False)
            print(f'{i} {stack.t["filename"][i]}')
            # plt.imshow(stretch(stack.image_single(i)), origin='lower')
            # plt.title(f"i={i}, et = {stack.t['et'][i]}")
        plt.show()
#%%%
        
    # Align the stack.
    # The values for the shifts (in pixels) are calculated based 
    # on putting MU69's ra/dec in the center of each frame.
    # This assumes that MU69 does not move within the stack itself.
    
    # ra_mu69  = 274.73344   # Value for 2-Nov-2018. Matches plot_img_wcs.py
    # dec_mu69 = -20.86170
    
    et = stack.t['et'][0]
    (st, _) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    (_, ra, dec) = sp.recrad(st[0:3])
    radec_mu69 = (ra, dec)
    
    stack.align(method = 'WCS', center = radec_mu69)  # Pass position in radians

    
#%%%    
    # Flatten the stack
    
    stack.padding = 65
    zoom = 4
    
    self = stack
    
    # Do the flattening
    
    (arr_flat, wcs_flat) = stack.flatten(zoom=zoom, do_force=True, do_plot=False, do_save=True)
    
    # Plot the newly flattened stack
    
    plot_img_wcs(arr_flat, wcs_flat, title = f'Zoom {zoom}, {str_reqid}', width=100,
                 name_target='MU69', name_observer='New Horizons', et = stack.t['et'][i],
                 )
        
    
# =============================================================================
# Now for debugging the small offsets, process some actual MU69 OpNav data.
# =============================================================================

#%%%
    
    # Do the most simple WCS possible
    # Concl: This one works!!
    
    do_simple = False
    
    if do_simple:
        file = '/Users/throop/Data/MU69_Approach/throop/backplaned/' + \
               'KALR_MU69_OpNav_L4_2018332/lor_0405688048_0x633_pwcs2_backplaned.fits'
        
        hdu = fits.open(file)
        img = hdu['PRIMARY'].data
        w   = WCS(file)
        et  = hdu['PRIMARY'].header['SPCSCET']
        
        plt.imshow(stretch(img), origin='lower')
        
        (st, lt)     = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
        (_, ra, dec) = sp.recrad(st[0:3])
        
        ra  *= hbt.r2d
        dec *= hbt.r2d
        
        (x,y) = w.wcs_world2pix(ra, dec, 0)
        
        hdu.close()
        
        plt.plot(x,y, marker = 'o', ms=10, alpha=0.5)
        plt.show
        