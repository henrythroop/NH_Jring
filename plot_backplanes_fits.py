#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 16:18:05 2018

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

def plot_backplanes_fits(file):

    """
    This file takes an image, and plots all of the backplanes for it.
    """

    # Start up SPICE
    
    file_kernel = '/Users/throop/git/NH_rings/kernels_kem.tm'
    sp.furnsh(file_kernel)

    hdulist = fits.open(file)

    # Loop over all of the planes, and plot each one
    
    i=1
    fig = plt.subplots()
    for hdu in hdulist:
        plt.subplot(3,4,i)
        plt.imshow(hdu.data)
        plt.title("{} / {}".format(i-1, hdu.name))
        i+=1

    plt.show()

    stretch_percent = 90    
    stretch = astropy.visualization.PercentileInterval(stretch_percent)
        
    # Look up position of MU69 in pixels.

    et = hdu[0].header['SPCSCET']
    utc = sp.et2utc(et, 'C', 0)
    abcorr = 'LT'
    frame = 'J2000'
    name_target = 'MU69'
    name_observer = 'New Horizons'
    w = WCS(file_new)
    
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
    plt.title("{}, {}".format(os.path.basename(file_new), utc))
    plt.show()  

    # Close the file
    
    hdu.close()

if (__name__ == '__main__'):
    
    file = '/Users/throop/Data/NH_KEM_Hazard/ORT1_Jan18/lor_0406731132_0x633_sci_HAZARD_test1-1_backplaned.fit'
    plot_backplanes_fits(file)
    