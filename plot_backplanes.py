#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 22:16:22 2018

@author: throop
"""


import math      
import astropy
from   astropy.io import fits
import numpy as np
import spiceypy as sp
from   astropy.visualization import wcsaxes
import hbt
from   astropy.wcs import WCS
import os
import matplotlib.pyplot as plt

from plot_img_wcs import plot_img_wcs

def plot_backplanes(file,
                       name_target = None,
                       name_observer = None):
    
    """
    This is a simple function to take a FITS file, and make a plot to the screen of all the backplanes.
    
    If name_target is passed, then code will also superimpose the position of the target, from the observer.
    
    SPICE must be up and running.
    
    This is a general function and will work on many FITS files.

    Parameters
    -----
    file: 
        FITS file to plot
        
    Optional keyword parameters
    -----
    name_target:
        Name of SPICE target
    name_observer:
        name of observer
        
    If both of the optional parameters are passed, then the code will plot the position of target seen from observer.
        
    """
    # Plot all of the planes to the screen, for validation

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

    plt.set_cmap('Greys_r')
    
    hdu = fits.open(file)
    
    num = len(hdu)
    
    # Make a grid nxy x nxy images, one for each backplane
    
    nxy = math.ceil(math.sqrt(num))
  
    # Plot all the backplanes, one by one

    hbt.figsize((12,12))
    fig = plt.subplots()
    
    for i in range(num):
        plt.subplot(nxy,nxy,i+1)
        plt.imshow(stretch(hdu[i].data), origin='lower')
        # plt.gca().get_xaxis().set_visible(False)   # Inhibit axes labels, since they get squashed.
        # plt.gca().get_yaxis().set_visible(False)
        plt.title(hdu[i].name)
        i+=1

    plt.tight_layout()
    plt.show()

# Plot the image, using plot_img_wcs

    hbt.figsize((8,8)) 
    
    w = WCS(file)
    plot_img_wcs(hdu[0].data, w, do_show=False, name_target='MU69', name_observer = 'New Horizons',
         et = hdu[0].header['SPCSCET'], width=100)
    
    radius_ring = 10_000  # This needs to be adjusted for different distances.
    radius_arr = hdu['Radius_eq'].data
    radius_good = np.logical_and(radius_arr > radius_ring*0.95, radius_arr < radius_ring*1.05)
    plt.imshow(radius_good, alpha=0.3, origin='lower', cmap='plasma')
    plt.show()    
    
#    Plot the image itself
 
    # plt.imshow(stretch(hdu[0].data), origin='lower')

    # # Over the images, superimpose one of the planes

    # plt.imshow(stretch(hdu['Longitude_eq'].data), alpha=0.5, cmap=plt.cm.Reds_r, origin='lower')

    # # If requested, look up position of target.

    # if (name_target):
    #     et = hdu[0].header['SPCSCET']
    #     utc = sp.et2utc(et, 'C', 0)
    #     abcorr = 'LT'
    #     frame = 'J2000'  

        
        
        
        
        
        
        
        
        
    #     (st,lt) = sp.spkezr(name_target, et, frame, abcorr, name_observer)
    #     vec_obs_target = st[0:3]
    #     (_, ra, dec) = sp.recrad(vec_obs_target)
    #     (pos_pix_x, pos_pix_y) = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)

    #     # Plot a ring, if we are looking at MU69. Use the backplane to filter by radius
     
    #     if name_target == 'MU69': 
    #         radius_ring = 10_000  # This needs to be adjusted for different distances.
    #         radius_arr = hdu['Radius_eq'].data
    #         radius_good = np.logical_and(radius_arr > radius_ring*0.95, radius_arr < radius_ring*1.05)
    #         plt.imshow(radius_good, alpha=0.3, origin='lower')
        
    #     # Plot the center of the system, as determined from the radius backplane
    #     # This shows up just as a bright pixel. It is not otherwise marked. It is really just for testing.

    #     if name_target == 'MU69': 
    #         radius_arr = hdu['Radius_eq'].data
    #         radius_center_mask = (radius_arr == np.amin(radius_arr))
    #         # plt.imshow(radius_center_mask, alpha=0.3, origin='lower')
    #         plt.imshow(radius_arr < 3000, alpha=0.3, origin='lower')
            
    # MAKE A SET OF CROSSHAIRS. THIS IS COOL FOR DEBUGGING, BUT NOT NEEDED NOW
    
    #         ddec_arr = hdu['dDec_km'].data
    #         ddec_mask = (np.abs(ddec_arr) < 1000)
    #         dra_arr = hdu['dRA_km'].data
    #         dra_mask = (np.abs(dra_arr) < 1000)
    #         dradec_mask = np.logical_or(dra_mask,ddec_mask)
    #         plt.imshow(dradec_mask, alpha=0.3, origin='lower')
              
        # Plot target body
        
        # plt.plot(pos_pix_x, pos_pix_y, ms=2, marker = 'o', color='green')    
        # plt.title("{}, {}".format(os.path.basename(file), utc))

        # print()
    
    hdu.close()

# =============================================================================
# End of function definition
# =============================================================================

# =============================================================================
# Make a test call to the function
# =============================================================================

if (__name__ == '__main__'):
        
#    file = '/Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ00/lor_0405178272_0x633_pwcs_backplaned.fits'

    file  = '/Users/throop/Data/MU69_Approach/throop/test/lor_0405121318_0x633_pwcs2_backplane.fits'
    
    name_observer = 'New Horizons'
    name_target = 'MU69'
    plot_backplanes(file, name_target = name_target, name_observer = name_observer)
    
    