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

def nh_plot_backplanes(file):
    
    """
    This is a simple function to take a FITS file, and make a plot to the screen of all the backplanes.
    
    Will also plot location of MU69.
    
    Developed for NH KEM ORT-1.
    
    """
    # Plot all of the planes to the screen, for validation

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent) # PI(90) scales to 5th..95th %ile.

    plt.set_cmap('Greys_r')
    
    hdu = fits.open(file)
    
    num = len(hdu)

    fig = plt.subplots()

    hbt.figsize((12,12))
    
    for i in range(num):
        plt.subplot(4,4,i+1)
        plt.imshow(stretch(hdu[i].data))
        plt.title(hdu[i].name)
        i+=1

    plt.show()
        
    # Start up SPICE
    
    file_kernel = '/Users/throop/git/NH_rings/kernels_kem.tm'
    sp.furnsh(file_kernel)
        
    # Look up position of MU69 in pixels.

    et = hdu[0].header['SPCSCET']
    utc = sp.et2utc(et, 'C', 0)
    abcorr = 'LT'
    frame = 'J2000'
    name_target = 'MU69'
    name_observer = 'New Horizons'
    w = WCS(file)
    
    (st,lt) = sp.spkezr(name_target, et, frame, abcorr, name_observer)
    vec_obs_mu69 = st[0:3]
    (_, ra, dec) = sp.recrad(vec_obs_mu69)
    (pos_pix_x, pos_pix_y) = w.wcs_world2pix(ra*hbt.r2d, dec*hbt.r2d, 0)
    
#    Plot the image itself
 
    hbt.figsize((10,10)) 
    plt.imshow(stretch(hdu[0].data))

    # Plot one of the planes

    plt.imshow(stretch(hdu['Longitude_eq'].data), alpha=0.5, cmap=plt.cm.Reds_r)

    # Plot the ring
 
    radius_ring = 100_000  # This needs to be adjusted for different distances.
    radius_arr = hdu['Radius_eq'].data
    radius_good = np.logical_and(radius_arr > radius_ring*0.95, radius_arr < radius_ring*1.05)
    plt.imshow(radius_good, alpha=0.3)
    
    # Plot MU69
    
    plt.plot(pos_pix_x, pos_pix_y, ms=10, marker = 'o', color='green')    
    plt.title("{}, {}".format(os.path.basename(file), utc))
    plt.show()    

# =============================================================================
# End of function definition
# =============================================================================

if (__name__ == '__main__'):
    
    file = '/Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ00/lor_0405178272_0x633_pwcs_backplaned.fits'
    
    nh_plot_backplanes(file)
    
    