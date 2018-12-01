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
    
    dir = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018332'
    # dir = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018328'
    
    # dir = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018316'
    # dir = '/Users/throop/Data/MU69_Approach/throop/backplaned/KALR_MU69_OpNav_L4_2018301'
    
    str_reqid = dir.split('/')[-1]
    
    stack = image_stack(dir, do_force=do_force)
    
    # Plot a set of images from it
    
    i = 1
    num_images = 4
    
    hbt.figsize((15,15))
    
    for i in range(num_images):
        plt.subplot(2,int(np.ceil(num_images/2)),i+1)
        filename_short = os.path.basename(stack.t['filename'][i])
        plot_img_wcs(stack.image_single(i), stack.t['wcs'][i], 
                       name_target='MU69', name_observer='New Horizons', et = stack.t['et'][i],
                       title = f'{i}',
                       width=70, do_show=False)
        
        # plt.imshow(stretch(stack.image_single(i)), origin='lower')
        # plt.title(f"i={i}, et = {stack.t['et'][i]}")
    plt.show()

    # Align the stack.
    # The values for the shifts (in pixels) are calculated based 
    # on putting MU69's ra/dec in the center of each frame.
    # This assumes that MU69 does not move. And to a very good approximation, that is true.
    # Total MU69 motion 15-Aug .. 15-Dec < 1 pixel.
    
    # ra_mu69  = 274.73344   # Value for 2-Nov-2018. Matches plot_img_wcs.py
    # dec_mu69 = -20.86170
    
    et = stack.t['et'][0]
    (st, _) = sp.spkezr('MU69', et, 'J2000', 'LT', 'New Horizons')
    (_, ra, dec) = sp.recrad(st[0:3])
    radec_mu69 = (ra, dec)
    
    stack.align(method = 'WCS', center = radec_mu69)  # Pass position in radians
    
    # Flatten the stack
    
    stack.padding = 65
    zoom = 4
    (arr_flat, wcs_flat) = stack.flatten(zoom=zoom, do_force=True, do_plot=True, do_save=True)
    
    plot_img_wcs(arr_flat, wcs_flat, title = f'Zoom {zoom}, {str_reqid}', width=50,
                 name_target='MU69', name_observer='New Horizons', et = stack.t['et'][i],
                 )
        
    
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