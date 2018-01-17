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

def plot_backplanes(file,
                       name_target = None,
                       name_observer = None):
    
    """
    This is a simple function to take a FITS file, and make a plot to the screen of all the backplanes.
    
    If name_target is passed, then code will also superimpose the position of the target, from the observer.
    
    SPICE must be up and running.
    
    This is a general function and will work on any 
    
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

    fig = plt.subplots()
    hbt.figsize((12,12))
    
    for i in range(num):
        plt.subplot(nxy,nxy,i+1)
        plt.imshow(stretch(hdu[i].data))
        plt.title(hdu[i].name)
        i+=1

    plt.show()
    
#    Plot the image itself
 
    hbt.figsize((10,10)) 
    plt.imshow(stretch(hdu[0].data))

    # Plot one of the planes

    plt.imshow(stretch(hdu['Longitude_eq'].data), alpha=0.5, cmap=plt.cm.Reds_r)

    # If requested, look up position of target.

    if (name_target):
        et = hdu[0].header['SPCSCET']
        utc = sp.et2utc(et, 'C', 0)
        abcorr = 'LT'
        frame = 'J2000'  

        w = WCS(file)
        
        (st,lt) = sp.spkezr(name_target, et, frame, abcorr, name_observer)
        vec_obs_target = st[0:3]
        (_, ra, dec) = sp.recrad(vec_obs_target)
        (pos_pix_x, pos_pix_y) = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)

        # Plot a ring, if we are looking at MU69. Use the backplane to filter by radius
     
        if name_target == 'MU69': 
            radius_ring = 100_000  # This needs to be adjusted for different distances.
            radius_arr = hdu['Radius_eq'].data
            radius_good = np.logical_and(radius_arr > radius_ring*0.95, radius_arr < radius_ring*1.05)
            plt.imshow(radius_good, alpha=0.3)
        
        # Plot target body
        
        plt.plot(pos_pix_x, pos_pix_y, ms=10, marker = 'o', color='green')    
        plt.title("{}, {}".format(os.path.basename(file), utc))
        plt.show() 
    
    hdu.close()

# =============================================================================
# End of function definition
# =============================================================================

# =============================================================================
# Make a test call to the function
# =============================================================================

if (__name__ == '__main__'):
    
    file = '/Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ00/lor_0405178272_0x633_pwcs_backplaned.fits'
    
    name_observer = 'New Horizons'
    name_target = 'MU69'
    plot_backplanes(file, name_target = name_target, name_observer = name_observer)
    
    